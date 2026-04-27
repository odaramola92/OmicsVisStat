"""
Two-Way ANOVA Implementation Matching R Approach

This implementation matches the R script approach exactly:
1. Extracts factors (Diet, Treatment) from column names
2. Creates 4-level Group factor (Diet_Treatment combinations)
3. Runs two-way ANOVA (Factor1 * Factor2)
4. Runs Tukey HSD on Group factor
5. Applies FDR correction once across metabolites for ANOVA p-values
6. Returns: Metabolite info + abundance columns + ANOVA stats + Tukey comparisons
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Sequence, Callable
from collections import Counter
import warnings
import re
from concurrent.futures import ThreadPoolExecutor

try:
    import statsmodels.api as _sm
    import statsmodels.formula.api as _smf
except Exception:
    _sm = None
    _smf = None

try:
    import statsmodels.stats.multicomp as _sm_mc
    from statsmodels.stats.multitest import multipletests
except Exception:
    _sm_mc = None
    multipletests = None


def _floor_pval(value: Optional[float]) -> Optional[float]:
    """Return p-value with only a machine-tiny floor.

    Final display flooring is applied per comparison column later using
    _apply_dynamic_pvalue_floor().
    """
    try:
        if value is None:
            return value
        v = float(value)
        if np.isnan(v):
            return v
        if v < 0:
            return np.nan
        # Keep raw p-value (including 0 from numeric underflow);
        # dynamic flooring is applied later, after analysis.
        return v
    except Exception:
        return np.nan


def _apply_dynamic_pvalue_floor(series: pd.Series, scale: float = 0.001) -> pd.Series:
    """Apply dynamic floor per p-value column.

    Floor rule requested:
      floor = smallest_positive_pvalue_in_column * scale
    where scale defaults to 0.01 (two orders of magnitude smaller).
    """
    s = pd.to_numeric(series, errors='coerce')
    finite_mask = np.isfinite(s)
    positive = s[finite_mask & (s > 0)]
    if positive.empty:
        return s

    dynamic_floor = float(positive.min()) * float(scale)
    if not np.isfinite(dynamic_floor) or dynamic_floor <= 0:
        return s

    # Replace non-positive finite values (e.g., 0 from underflow) and anything below floor.
    s = s.mask(finite_mask & (s <= 0), dynamic_floor)
    s = s.mask(np.isfinite(s) & (s < dynamic_floor), dynamic_floor)
    return s


def _safe_fc(m1: float, m2: float) -> Tuple[float, float]:
    """Compute fold-change safely. Returns (fc, log2fc)."""
    if m2 == 0 or np.isnan(m1) or np.isnan(m2):
        return np.nan, np.nan
    fc = m1 / m2
    if fc <= 0:
        return fc, np.nan
    return fc, np.log2(fc)


def perform_two_way_anova_new(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    factor_a_map: Dict[str, str],
    factor_b_map: Dict[str, str],
    *,
    group_map: Optional[Dict[str, str]] = None,
    group_order: Optional[List[str]] = None,
    factor_a_name: str = 'Diet',
    factor_b_name: str = 'Treatment',
    drop_zeros: bool = True,
    min_per_cell: int = 1,
    fdr: bool = True,
    typ: int = 2,
    n_jobs: int = 3,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Run two-way ANOVA matching R approach exactly.
    
    R Approach:
    1. Extract factors from column names (HFD/NC -> Diet, TBI/Sham -> Treatment)
    2. Create Group = paste(Diet, Treatment, sep = "_")
    3. Run: aov(value ~ Diet * Treatment)
    4. Run: TukeyHSD(aov(value ~ Group))
    5. FDR correction across metabolites for ANOVA p-values
    6. Output: Metabolite + abundance + ANOVA p/adj_p + Tukey comparisons
    
    Parameters
    ----------
    df : pd.DataFrame
        Input data with metabolite info columns + abundance columns
    sample_cols : Sequence[str]
        Column names for abundance data
    factor_a_map : dict
        Sample -> Factor A level (e.g., 'HFD' or 'NC')
    factor_b_map : dict
        Sample -> Factor B level (e.g., 'TBI' or 'Sham')
    group_map : dict, optional
        Sample -> Group name (e.g., 'HFD_1' -> 'HFD_TBI')
    group_order : list, optional
        Desired order of groups for output
    factor_a_name : str
        Name for Factor A (default 'Diet')
    factor_b_name : str
        Name for Factor B (default 'Treatment')
    
    Returns
    -------
    (complete_result_df, individual_pairwise_sheets_dict)
        - complete_result_df: Metabolite info + abundance + ANOVA stats + Tukey comparisons
        - individual_pairwise_sheets_dict: One sheet per pairwise comparison
    """
    
    if _sm is None or _smf is None or _sm_mc is None:
        raise ImportError("statsmodels is required for two-way ANOVA")

    # Validation
    missing_a = [c for c in sample_cols if c not in factor_a_map]
    missing_b = [c for c in sample_cols if c not in factor_b_map]
    if missing_a or missing_b:
        raise ValueError(f"Missing factor mappings for samples")

    # Find metabolite ID column
    metabolite_id_col: Optional[str] = None
    for candidate in ['Name', 'Protein', 'Metabolite', 'Molecule', 'Compound', 'LipidID', 'Lipid_ID', 'Lipid_Class']:
        if candidate in df.columns:
            metabolite_id_col = candidate
            break
    if metabolite_id_col is None and len(df.columns) > 0:
        metabolite_id_col = df.columns[0]

    # Get non-abundance columns (metabolite info columns)
    info_cols = [c for c in df.columns if c not in sample_cols]

    # Build Group factor from Factor A + Factor B
    # Group naming: If group_map provided, use it; else combine factors with underscore
    sample_to_group = {}
    for s in sample_cols:
        if group_map and s in group_map:
            sample_to_group[s] = group_map[s]
        else:
            # Default: FactorA_FactorB
            sample_to_group[s] = f"{factor_a_map[s]}_{factor_b_map[s]}"
    
    # Get unique groups in order
    if group_order:
        unique_groups = [g for g in group_order if g in sample_to_group.values()]
        # Add any groups not in group_order
        for g in set(sample_to_group.values()):
            if g not in unique_groups:
                unique_groups.append(g)
    else:
        unique_groups = sorted(set(sample_to_group.values()))
    
    # Generate all pairwise combinations
    pairwise_combos = []
    for i, g1 in enumerate(unique_groups):
        for g2 in unique_groups[i+1:]:
            pairwise_combos.append((g1, g2))

    print(f"\n🔬 Two-Way ANOVA Configuration:")
    print(f"   Factor A ({factor_a_name}): {sorted(set(factor_a_map.values()))}")
    print(f"   Factor B ({factor_b_name}): {sorted(set(factor_b_map.values()))}")
    print(f"   Groups: {unique_groups}")
    print(f"   Pairwise comparisons: {len(pairwise_combos)}")
    print(f"   Metabolites: {len(df)}")
    print(f"   Formula: value ~ {factor_a_name} * {factor_b_name}")
    print(f"   Tukey HSD: value ~ Group\n")

    def compute_entry(item):
        """Compute ANOVA + Tukey for one metabolite (matching R approach)"""
        idx, row = item
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        
        # Progress indicator
        if n_jobs == 1 and idx % 50 == 0:
            try:
                print(f"  Processing metabolite {idx}/{len(df)}...")
            except:
                pass
        
        # Initialize result with metabolite info + abundance values
        result = {}
        for col in info_cols:
            result[col] = row[col]
        for s in sample_cols:
            try:
                result[s] = float(pd.to_numeric(row[s], errors='coerce'))
            except:
                result[s] = np.nan
        
        # Initialize ANOVA columns
        result[f'p_{factor_a_name}'] = np.nan
        result[f'p_{factor_b_name}'] = np.nan
        result[f'p_interaction'] = np.nan
        result[f'adj_p_{factor_a_name}'] = np.nan
        result[f'adj_p_{factor_b_name}'] = np.nan
        result[f'adj_p_interaction'] = np.nan
        
        # Initialize Tukey columns for all pairwise comparisons
        for g1, g2 in pairwise_combos:
            comp_name = f"{g1}_vs_{g2}"
            result[f'{comp_name}_adj_p'] = np.nan
            result[f'{comp_name}_FC'] = np.nan
            result[f'{comp_name}_log2FC'] = np.nan
            result[f'{comp_name}_neg_log10_adj_p'] = np.nan
            result[f'{comp_name}_model_effect'] = np.nan
            result[f'{comp_name}_model_se'] = np.nan
            result[f'{comp_name}_ci_lower_95'] = np.nan
            result[f'{comp_name}_ci_upper_95'] = np.nan
        
        # Collect values and factor assignments
        values: List[float] = []
        factor_a_levels: List[str] = []
        factor_b_levels: List[str] = []
        group_levels: List[str] = []
        
        for s in sample_cols:
            try:
                v = float(pd.to_numeric(row[s], errors='coerce'))
            except:
                v = np.nan
            
            if np.isnan(v):
                continue
            if drop_zeros and v == 0:
                continue
            
            values.append(v)
            factor_a_levels.append(factor_a_map[s])
            factor_b_levels.append(factor_b_map[s])
            group_levels.append(sample_to_group[s])
        
        # Check minimum requirements
        if len(values) < 3:
            return result
        if len(set(factor_a_levels)) < 2 or len(set(factor_b_levels)) < 2:
            return result
        
        # Check min_per_cell if required
        if min_per_cell > 1:
            cell_counts = Counter(zip(factor_a_levels, factor_b_levels))
            if any(ct < min_per_cell for ct in cell_counts.values()):
                return result
        
        # ============================================================
        # TWO-WAY ANOVA: value ~ Factor_A * Factor_B
        # ============================================================
        long_df = pd.DataFrame({
            'value': values,
            'Factor_A': factor_a_levels,
            'Factor_B': factor_b_levels,
            'Group': group_levels
        })
        
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                # Run two-way ANOVA (matching R: aov(value ~ Diet * Treatment))
                model = _smf.ols('value ~ C(Factor_A) * C(Factor_B)', data=long_df).fit()
                anova_tbl = _sm.stats.anova_lm(model, typ=typ)
            
            # Extract p-values (matching R ANOVA table)
            result[f'p_{factor_a_name}'] = _floor_pval(anova_tbl.loc['C(Factor_A)', 'PR(>F)'])
            result[f'p_{factor_b_name}'] = _floor_pval(anova_tbl.loc['C(Factor_B)', 'PR(>F)'])
            result[f'p_interaction'] = _floor_pval(anova_tbl.loc['C(Factor_A):C(Factor_B)', 'PR(>F)'])
        except Exception as e:
            if n_jobs == 1:
                try:
                    print(f"  ⚠️ ANOVA failed for {metabolite_id}: {e}")
                except:
                    pass
        
        # ============================================================
        # TUKEY HSD: on Group factor (matching R: TukeyHSD(aov(value ~ Group)))
        # ============================================================
        try:
            # Build group means for fold-change calculation
            group_means = {}
            for g in unique_groups:
                g_vals = [v for v, gl in zip(values, group_levels) if gl == g]
                if g_vals:
                    group_means[g] = np.mean(g_vals)
            
            # Run Tukey HSD on Group factor
            if len(set(group_levels)) >= 2:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')
                    mc = _sm_mc.MultiComparison(values, group_levels)
                    tukey = mc.tukeyhsd(alpha=0.05)
                    
                    # Also fit OLS model on Group for model-based CI extraction
                    long_df.loc[:, 'Group'] = group_levels
                    model_group = _smf.ols('value ~ C(Group)', data=long_df).fit()
                
                tukey_df = pd.DataFrame(
                    data=tukey._results_table.data[1:],
                    columns=tukey._results_table.data[0]
                )
                
                # Extract Tukey results for each pairwise comparison
                for _, row_t in tukey_df.iterrows():
                    g1 = str(row_t['group1'])
                    g2 = str(row_t['group2'])
                    p_adj = _floor_pval(row_t.get('p-adj', np.nan))
                    
                    # Match comparison to our pairwise_combos order
                    comp_key = None
                    if (g1, g2) in pairwise_combos:
                        comp_key = (g1, g2)
                    elif (g2, g1) in pairwise_combos:
                        comp_key = (g2, g1)
                        g1, g2 = g2, g1  # swap to match our order
                    
                    if comp_key:
                        comp_name = f"{g1}_vs_{g2}"
                        result[f'{comp_name}_adj_p'] = p_adj
                        
                        # Calculate fold-change: FC = g2/g1 (second/first)
                        m1 = group_means.get(g1, np.nan)
                        m2 = group_means.get(g2, np.nan)
                        fc, log2fc = _safe_fc(m2, m1)
                        
                        result[f'{comp_name}_FC'] = fc
                        result[f'{comp_name}_log2FC'] = log2fc
                        result[f'{comp_name}_neg_log10_adj_p'] = -np.log10(p_adj) if p_adj > 0 else 15
                        
                        # Extract model-based effect and SE
                        try:
                            g1_col = f'C(Group)[T.{g1}]'
                            g2_col = f'C(Group)[T.{g2}]'
                            
                            # Handle both groups as dummies vs one as reference
                            if g1_col in model_group.params.index and g2_col in model_group.params.index:
                                # Both are dummies - contrast is g2 - g1
                                model_effect = float(model_group.params[g2_col] - model_group.params[g1_col])
                                se_g2 = float(model_group.bse[g2_col])
                                se_g1 = float(model_group.bse[g1_col])
                                model_se = np.sqrt(se_g2**2 + se_g1**2)
                            elif g2_col in model_group.params.index:
                                # g1 is reference, g2 is dummy
                                model_effect = float(model_group.params[g2_col])
                                model_se = float(model_group.bse[g2_col])
                            elif g1_col in model_group.params.index:
                                # g2 is reference, g1 is dummy (flip sign)
                                model_effect = -float(model_group.params[g1_col])
                                model_se = float(model_group.bse[g1_col])
                            else:
                                model_effect = np.nan
                                model_se = np.nan
                            
                            result[f'{comp_name}_model_effect'] = model_effect
                            result[f'{comp_name}_model_se'] = model_se
                            try:
                                # Use scipy/statsmodels-independent t critical via scipy if available through statsmodels internals.
                                from scipy import stats as _scipy_stats
                                t_crit = _scipy_stats.t.ppf(0.975, model_group.df_resid) if model_group.df_resid > 0 else np.nan
                            except Exception:
                                t_crit = np.nan
                            if np.isfinite(t_crit) and np.isfinite(model_effect) and np.isfinite(model_se):
                                result[f'{comp_name}_ci_lower_95'] = model_effect - t_crit * model_se
                                result[f'{comp_name}_ci_upper_95'] = model_effect + t_crit * model_se
                        except Exception:
                            result[f'{comp_name}_model_effect'] = np.nan
                            result[f'{comp_name}_model_se'] = np.nan
        
        except Exception as e:
            if n_jobs == 1:
                try:
                    print(f"  ⚠️ Tukey HSD failed for {metabolite_id}: {e}")
                except:
                    pass
        
        return result

    # Process all metabolites
    items = list(df.iterrows())
    total = len(items)
    
    print(f"⏳ Processing {total} metabolites...")
    
    if n_jobs and n_jobs > 1:
        with ThreadPoolExecutor(max_workers=max(1, int(n_jobs))) as ex:
            rows = []
            for i, result in enumerate(ex.map(compute_entry, items), 1):
                rows.append(result)
                if progress_callback and (i % max(1, total // 20) == 0 or i == total):
                    try:
                        metab_name = result.get(metabolite_id_col, 'Unknown')
                        progress_callback(i, total, str(metab_name))
                    except:
                        pass
    else:
        rows = []
        for i, it in enumerate(items, 1):
            result = compute_entry(it)
            rows.append(result)
            if progress_callback and (i % max(1, total // 20) == 0 or i == total):
                try:
                    metab_name = result.get(metabolite_id_col, 'Unknown')
                    progress_callback(i, total, str(metab_name))
                except:
                    pass

    complete_result = pd.DataFrame(rows)
    print(f"✅ Metabolite processing complete: {len(complete_result)} rows\n")

    # ============================================================
    # FDR CORRECTION (BH method) - Applied once across metabolites
    # Matching R: p.adjust(p_values, method = "BH")
    # ============================================================
    if fdr and not complete_result.empty and multipletests is not None:
        print(f"🔄 Applying FDR correction (Benjamini-Hochberg) across {len(complete_result)} metabolites...")
        
        # Correct ANOVA p-values (matching R approach)
        for p_col in [f'p_{factor_a_name}', f'p_{factor_b_name}', 'p_interaction']:
            adj_col = p_col.replace('p_', 'adj_p_')
            
            mask = complete_result[p_col].notna()
            if mask.any():
                try:
                    _, padj_vec, _, _ = multipletests(
                        complete_result.loc[mask, p_col],
                        method='fdr_bh'
                    )
                    complete_result.loc[mask, adj_col] = padj_vec
                    n_sig = (padj_vec < 0.05).sum()
                    print(f"   {p_col}: {n_sig}/{mask.sum()} significant after FDR (α=0.05)")
                except Exception as e:
                    print(f"   ⚠️ FDR correction failed for {p_col}: {e}")
        
        print(f"✅ FDR correction complete\n")
    
    # Note: Tukey HSD p-values are already adjusted (p-adj from TukeyHSD)
    # R script uses these directly without additional correction
    # Our implementation matches this: we don't apply additional FDR to Tukey p-values

    # ============================================================
    # DYNAMIC P-VALUE FLOORING (per comparison column)
    # floor = 0.01 * smallest positive p-value in that column
    # ============================================================
    pval_cols = []
    for col in complete_result.columns:
        if col.startswith('p_') or col.startswith('adj_p_') or col.endswith('_adj_p'):
            pval_cols.append(col)

    for col in pval_cols:
        complete_result[col] = _apply_dynamic_pvalue_floor(complete_result[col])

    # Recompute -log10 columns from floored adj p-values.
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        p_col = f'{comp_name}_adj_p'
        neg_log_col = f'{comp_name}_neg_log10_adj_p'
        if p_col in complete_result.columns and neg_log_col in complete_result.columns:
            complete_result[neg_log_col] = -np.log10(
                pd.to_numeric(complete_result[p_col], errors='coerce').replace(0, np.finfo(float).tiny)
            )

    # ============================================================
    # CREATE INDIVIDUAL PAIRWISE SHEETS
    # Each sheet contains: Metabolite info + 2 groups' samples + comparison stats
    # ============================================================
    print(f"📊 Creating individual pairwise comparison sheets...")
    individual_sheets = {}
    
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        
        # Get samples from each group
        samples_g1 = [s for s in sample_cols if sample_to_group[s] == g1]
        samples_g2 = [s for s in sample_cols if sample_to_group[s] == g2]
        all_samples = samples_g1 + samples_g2
        
        if not all_samples:
            continue
        
        # Build sheet with: info cols + samples + comparison stats
        cols_to_include = info_cols.copy()
        cols_to_include.extend([s for s in all_samples if s in complete_result.columns])
        
        # Add comparison-specific stats
        stat_cols = [
            f'{comp_name}_adj_p',
            f'{comp_name}_FC',
            f'{comp_name}_log2FC',
            f'{comp_name}_neg_log10_adj_p',
            f'{comp_name}_model_effect',
            f'{comp_name}_model_se',
            f'{comp_name}_ci_lower_95',
            f'{comp_name}_ci_upper_95'
        ]
        cols_to_include.extend([c for c in stat_cols if c in complete_result.columns])
        
        # Create sheet
        sheet_df = complete_result[cols_to_include].copy()
        
        # Keep sheet even if p-values are missing so users can inspect raw abundances
        # (previously dropped all-NaN rows, which resulted in zero pairwise sheets)
        individual_sheets[comp_name] = sheet_df
        print(f"   ✓ {comp_name}: {len(sheet_df)} metabolites")
    
    print(f"✅ Created {len(individual_sheets)} pairwise sheets\n")

    # ============================================================
    # COLUMN ORDERING (match R output structure)
    # R output: Metabolite info + abundance columns + ANOVA stats + Tukey comparisons
    # ============================================================
    # Order columns logically
    ordered_cols = []
    
    # 1. Metabolite info columns (first)
    ordered_cols.extend(info_cols)
    
    # 2. Abundance columns
    ordered_cols.extend([s for s in sample_cols if s in complete_result.columns])
    
    # 3. ANOVA columns (p-values then adjusted p-values)
    anova_cols = [
        f'p_{factor_a_name}',
        f'p_{factor_b_name}',
        f'p_interaction',
        f'adj_p_{factor_a_name}',
        f'adj_p_{factor_b_name}',
        f'adj_p_interaction'
    ]
    ordered_cols.extend([c for c in anova_cols if c in complete_result.columns])
    
    # 4. Tukey pairwise comparisons
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        tukey_cols = [
            f'{comp_name}_adj_p',
            f'{comp_name}_FC',
            f'{comp_name}_log2FC',
            f'{comp_name}_neg_log10_adj_p',
            f'{comp_name}_model_effect',
            f'{comp_name}_model_se',
            f'{comp_name}_ci_lower_95',
            f'{comp_name}_ci_upper_95'
        ]
        ordered_cols.extend([c for c in tukey_cols if c in complete_result.columns])
    
    # Add any remaining columns not already included
    for c in complete_result.columns:
        if c not in ordered_cols:
            ordered_cols.append(c)
    
    complete_result = complete_result[ordered_cols]

    print(f"📋 Output Summary:")
    print(f"   Rows: {len(complete_result)}")
    print(f"   Columns: {len(complete_result.columns)}")
    print(f"   - Info columns: {len(info_cols)}")
    print(f"   - Abundance columns: {len([c for c in sample_cols if c in complete_result.columns])}")
    print(f"   - ANOVA columns: {len([c for c in anova_cols if c in complete_result.columns])}")
    print(f"   - Tukey columns: {len([c for c in complete_result.columns if '_vs_' in c])}")
    print()

    return complete_result, individual_sheets
