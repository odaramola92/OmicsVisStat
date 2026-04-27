"""
One-Way ANOVA Implementation Following R Standard Approach

This implementation matches R's one-way ANOVA approach:
1. Run one-way ANOVA (F-test) or Kruskal-Wallis per metabolite
2. Run post-hoc test (Tukey HSD for ANOVA, Dunn for Kruskal) for pairwise comparisons
3. Apply FDR correction once across metabolites
4. Return: Complete result + individual pairwise sheets

Matches the clean pattern from two_way_anova_new_format.py
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Sequence, Callable
from collections import Counter
import warnings
from concurrent.futures import ThreadPoolExecutor

try:
    import statsmodels.api as _sm
    import statsmodels.formula.api as _smf
    import statsmodels.stats.multicomp as _sm_mc
    from statsmodels.stats.multitest import multipletests
    from scipy import stats
except Exception:
    _sm = None
    _smf = None
    _sm_mc = None
    multipletests = None
    stats = None


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


def _welch_effect_ci(v1: np.ndarray, v2: np.ndarray) -> Tuple[float, float, float, float]:
    """Return Welch-style mean-difference effect, SE, and 95% CI for (g2 - g1)."""
    try:
        n1 = len(v1)
        n2 = len(v2)
        if n1 < 2 or n2 < 2 or stats is None:
            return np.nan, np.nan, np.nan, np.nan

        m1 = float(np.mean(v1))
        m2 = float(np.mean(v2))
        s1 = float(np.var(v1, ddof=1))
        s2 = float(np.var(v2, ddof=1))

        effect = m2 - m1
        se2 = (s1 / n1) + (s2 / n2)
        if not np.isfinite(se2) or se2 <= 0:
            return effect, np.nan, np.nan, np.nan

        se = float(np.sqrt(se2))

        # Welch-Satterthwaite degrees of freedom.
        num = se2 ** 2
        den = ((s1 / n1) ** 2 / (n1 - 1)) + ((s2 / n2) ** 2 / (n2 - 1))
        if not np.isfinite(den) or den <= 0:
            return effect, se, np.nan, np.nan

        df = num / den
        if not np.isfinite(df) or df <= 0:
            return effect, se, np.nan, np.nan

        t_crit = float(stats.t.ppf(0.975, df))
        if not np.isfinite(t_crit):
            return effect, se, np.nan, np.nan

        ci_low = effect - t_crit * se
        ci_high = effect + t_crit * se
        return effect, se, float(ci_low), float(ci_high)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def _mannwhitney_effect_ci(v1: np.ndarray, v2: np.ndarray) -> Tuple[float, float, float, float]:
    """Return rank-biserial effect, asymptotic SE, and 95% CI for (g2 vs g1).

    Effect is rank-biserial correlation in [-1, 1], oriented so positive means
    higher tendency in g2 than g1.
    """
    try:
        n1 = len(v1)
        n2 = len(v2)
        if n1 < 1 or n2 < 1 or stats is None:
            return np.nan, np.nan, np.nan, np.nan

        u1, _ = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        total_pairs = float(n1 * n2)
        if total_pairs <= 0:
            return np.nan, np.nan, np.nan, np.nan

        # Convert U for v1 to orientation for g2-vs-g1.
        u2 = total_pairs - float(u1)
        effect = (2.0 * u2 / total_pairs) - 1.0

        # Large-sample approximation without tie correction.
        var_u = total_pairs * (n1 + n2 + 1.0) / 12.0
        if not np.isfinite(var_u) or var_u <= 0:
            return float(effect), np.nan, np.nan, np.nan

        se = 2.0 * np.sqrt(var_u) / total_pairs
        if not np.isfinite(se) or se <= 0:
            return float(effect), np.nan, np.nan, np.nan

        z_crit = 1.959963984540054
        ci_low = max(-1.0, effect - z_crit * se)
        ci_high = min(1.0, effect + z_crit * se)
        return float(effect), float(se), float(ci_low), float(ci_high)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def perform_one_way_anova(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Dict[str, str],
    *,
    group_order: Optional[List[str]] = None,
    overall_test: str = 'anova',
    pairwise_test: Optional[str] = None,
    drop_zeros: bool = True,
    min_group_size: int = 2,
    fdr: bool = True,
    typ: int = 2,
    n_jobs: int = 3,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Run one-way ANOVA/Kruskal-Wallis matching R approach.
    
    R Approach:
    1. Run: aov(value ~ group) or kruskal.test(value ~ group)
    2. Run post-hoc: TukeyHSD(fit) or dunn.test()
    3. FDR correction across metabolites for overall p-values
    4. Output: Metabolite + abundance + overall stats + pairwise stats
    
    Parameters
    ----------
    df : pd.DataFrame
        Input data with metabolite info columns + abundance columns
    sample_cols : Sequence[str]
        Column names for abundance data
    group_map : dict
        Sample -> Group name (e.g., 'Sample_1' -> 'Control')
    group_order : list, optional
        Desired order of groups for output
    overall_test : str
        'anova' for one-way ANOVA (parametric) or 'kruskal' for Kruskal-Wallis (non-parametric)
    pairwise_test : str, optional
        'tukey' for Tukey HSD (ANOVA only), 'dunn' for Dunn test (Kruskal only),
        'welch' for Welch's t-test, 'mannwhitney' for Mann-Whitney U test
        If None, defaults to 'tukey' for ANOVA and 'dunn' for Kruskal
    drop_zeros : bool
        Whether to drop zero values
    min_group_size : int
        Minimum samples required per group
    fdr : bool
        Apply FDR correction (Benjamini-Hochberg)
    n_jobs : int
        Number of parallel workers
    progress_callback : callable, optional
        Progress callback function(current, total, metabolite_name)
    
    Returns
    -------
    (complete_result_df, individual_pairwise_sheets_dict)
        - complete_result_df: Metabolite info + abundance + overall stats + pairwise comparisons
        - individual_pairwise_sheets_dict: One sheet per pairwise comparison
    """
    
    if stats is None or _sm_mc is None:
        raise ImportError("scipy and statsmodels are required for one-way ANOVA")

    # Validation
    missing = [c for c in sample_cols if c not in group_map]
    if missing:
        raise ValueError(f"Missing group mappings for samples: {missing[:5]}")

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

    # Get unique groups in order
    if group_order:
        unique_groups = [g for g in group_order if g in group_map.values()]
        # Add any groups not in group_order
        for g in set(group_map.values()):
            if g not in unique_groups:
                unique_groups.append(g)
    else:
        unique_groups = sorted(set(group_map.values()))
    
    # Filter groups by minimum size
    group_samples = {g: [s for s in sample_cols if group_map[s] == g] for g in unique_groups}
    valid_groups = [g for g in unique_groups if len(group_samples[g]) >= min_group_size]
    
    if len(valid_groups) < 2:
        raise ValueError(f"Need at least 2 groups with ≥{min_group_size} samples. Found {len(valid_groups)}")
    
    unique_groups = valid_groups
    
    # Generate all pairwise combinations
    pairwise_combos = []
    for i, g1 in enumerate(unique_groups):
        for g2 in unique_groups[i+1:]:
            pairwise_combos.append((g1, g2))

    # Determine pairwise test method
    if pairwise_test is None:
        pairwise_test = 'tukey' if overall_test == 'anova' else 'dunn'

    def compute_entry(item):
        """Compute one-way ANOVA + post-hoc for one metabolite"""
        idx, row = item
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        
        # Progress indicator removed for cleaner output
        
        # Initialize result with metabolite info + abundance values
        result = {}
        for col in info_cols:
            result[col] = row[col]
        for s in sample_cols:
            try:
                result[s] = float(pd.to_numeric(row[s], errors='coerce'))
            except:
                result[s] = np.nan
        
        # Initialize overall test columns
        result['overall_statistic'] = np.nan
        result['overall_p'] = np.nan
        result['overall_adj_p'] = np.nan
        
        # Initialize pairwise columns for all comparisons
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
        
        # Collect values per group
        values: List[float] = []
        group_levels: List[str] = []
        group_values_dict = {}
        
        for g in unique_groups:
            g_samples = group_samples[g]
            g_vals = []
            for s in g_samples:
                try:
                    v = float(pd.to_numeric(row[s], errors='coerce'))
                except:
                    v = np.nan
                
                if np.isnan(v):
                    continue
                if drop_zeros and v == 0:
                    continue
                
                g_vals.append(v)
                values.append(v)
                group_levels.append(g)
            
            group_values_dict[g] = np.array(g_vals)
        
        # Check minimum requirements
        if len(values) < 3:
            return result
        
        # Filter to groups that meet minimum size requirement
        valid_groups = [g for g in unique_groups if len(group_values_dict[g]) >= min_group_size]
        
        # Need at least 2 valid groups to perform comparison
        if len(valid_groups) < 2:
            return result
        
        # Update group lists to only include valid groups
        valid_group_values = [group_values_dict[g] for g in valid_groups]
        
        # ============================================================
        # ONE-WAY ANOVA or KRUSKAL-WALLIS (on valid groups only)
        # ============================================================
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                
                if overall_test == 'anova':
                    # One-way ANOVA (matching R: aov(value ~ group))
                    statistic, p_value = stats.f_oneway(*valid_group_values)
                else:  # kruskal
                    # Kruskal-Wallis (matching R: kruskal.test(value ~ group))
                    statistic, p_value = stats.kruskal(*valid_group_values)
                
                result['overall_statistic'] = float(statistic)
                result['overall_p'] = _floor_pval(float(p_value))
        except Exception as e:
            pass  # Silently skip failed tests
        
        # ============================================================
        # POST-HOC PAIRWISE TESTS (only between valid groups)
        # ============================================================
        # Filter pairwise_combos to only include comparisons between valid groups
        valid_pairwise_combos = [(g1, g2) for g1, g2 in pairwise_combos if g1 in valid_groups and g2 in valid_groups]
        
        try:
            # Build group means for fold-change calculation
            group_means = {g: np.mean(group_values_dict[g]) if len(group_values_dict[g]) > 0 else np.nan 
                          for g in unique_groups}
            
            if pairwise_test == 'tukey' and overall_test == 'anova' and len(set(group_levels)) >= 2:
                # Tukey HSD (matching R: TukeyHSD(aov(value ~ group)))
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')
                    mc = _sm_mc.MultiComparison(values, group_levels)
                    tukey = mc.tukeyhsd(alpha=0.05)
                    
                    # Also fit OLS model to extract model-based CI
                    long_df = pd.DataFrame({'value': values, 'group': group_levels})
                    model = _smf.ols('value ~ C(group)', data=long_df).fit()
                
                tukey_df = pd.DataFrame(
                    data=tukey._results_table.data[1:],
                    columns=tukey._results_table.data[0]
                )
                
                # Extract Tukey results for each pairwise comparison (only valid groups)
                for _, row_t in tukey_df.iterrows():
                    g1 = str(row_t['group1'])
                    g2 = str(row_t['group2'])
                    p_adj = _floor_pval(row_t.get('p-adj', np.nan))
                    
                    # Match comparison to our valid pairwise_combos order
                    comp_key = None
                    if (g1, g2) in valid_pairwise_combos:
                        comp_key = (g1, g2)
                    elif (g2, g1) in valid_pairwise_combos:
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
                            g1_col = f'C(group)[T.{g1}]'
                            g2_col = f'C(group)[T.{g2}]'
                            
                            # Handle both groups as dummies vs one as reference
                            if g1_col in model.params.index and g2_col in model.params.index:
                                # Both are dummies - contrast is g2 - g1
                                model_effect = float(model.params[g2_col] - model.params[g1_col])
                                se_g2 = float(model.bse[g2_col])
                                se_g1 = float(model.bse[g1_col])
                                model_se = np.sqrt(se_g2**2 + se_g1**2)
                            elif g2_col in model.params.index:
                                # g1 is reference, g2 is dummy
                                model_effect = float(model.params[g2_col])
                                model_se = float(model.bse[g2_col])
                            elif g1_col in model.params.index:
                                # g2 is reference, g1 is dummy (flip sign)
                                model_effect = -float(model.params[g1_col])
                                model_se = float(model.bse[g1_col])
                            else:
                                model_effect = np.nan
                                model_se = np.nan
                            
                            result[f'{comp_name}_model_effect'] = model_effect
                            result[f'{comp_name}_model_se'] = model_se
                            try:
                                t_crit = stats.t.ppf(0.975, model.df_resid) if model.df_resid > 0 else np.nan
                                if np.isfinite(t_crit) and np.isfinite(model_effect) and np.isfinite(model_se):
                                    result[f'{comp_name}_ci_lower_95'] = model_effect - t_crit * model_se
                                    result[f'{comp_name}_ci_upper_95'] = model_effect + t_crit * model_se
                            except Exception:
                                pass
                        except Exception:
                            result[f'{comp_name}_model_effect'] = np.nan
                            result[f'{comp_name}_model_se'] = np.nan
            
            elif pairwise_test == 'dunn' and overall_test == 'kruskal':
                # Dunn test for Kruskal-Wallis post-hoc
                # Note: Requires scikit-posthocs package
                try:
                    import scikit_posthocs as sp
                    dunn_df = sp.posthoc_dunn(values, group_levels, p_adjust='holm')
                    
                    for g1, g2 in valid_pairwise_combos:
                        if g1 in dunn_df.index and g2 in dunn_df.columns:
                            p_adj = _floor_pval(dunn_df.loc[g1, g2])
                            comp_name = f"{g1}_vs_{g2}"
                            result[f'{comp_name}_adj_p'] = p_adj
                            
                            # Calculate fold-change
                            m1 = group_means.get(g1, np.nan)
                            m2 = group_means.get(g2, np.nan)
                            fc, log2fc = _safe_fc(m2, m1)
                            
                            result[f'{comp_name}_FC'] = fc
                            result[f'{comp_name}_log2FC'] = log2fc
                            result[f'{comp_name}_neg_log10_adj_p'] = -np.log10(p_adj) if p_adj > 0 else 15

                            # Nonparametric effect summary for consistency in output columns.
                            v1 = group_values_dict.get(g1, np.array([]))
                            v2 = group_values_dict.get(g2, np.array([]))
                            eff, se_eff, ci_low, ci_high = _mannwhitney_effect_ci(v1, v2)
                            result[f'{comp_name}_model_effect'] = eff
                            result[f'{comp_name}_model_se'] = se_eff
                            result[f'{comp_name}_ci_lower_95'] = ci_low
                            result[f'{comp_name}_ci_upper_95'] = ci_high
                except ImportError:
                    # Fallback to pairwise Mann-Whitney with Holm correction
                    pairwise_ps = []
                    for g1, g2 in valid_pairwise_combos:
                        v1 = group_values_dict[g1]
                        v2 = group_values_dict[g2]
                        if len(v1) >= min_group_size and len(v2) >= min_group_size:
                            _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                            pairwise_ps.append(_floor_pval(p))
                        else:
                            pairwise_ps.append(np.nan)
                    
                    # Holm correction
                    valid_mask = np.array([not np.isnan(p) for p in pairwise_ps])
                    if valid_mask.any():
                        _, padj_vec, _, _ = multipletests([p for p, v in zip(pairwise_ps, valid_mask) if v], method='holm')
                        padj_idx = 0
                        for i, (g1, g2) in enumerate(valid_pairwise_combos):
                            if valid_mask[i]:
                                p_adj = padj_vec[padj_idx]
                                padj_idx += 1
                            else:
                                p_adj = np.nan
                            
                            comp_name = f"{g1}_vs_{g2}"
                            result[f'{comp_name}_adj_p'] = p_adj
                            
                            # Calculate fold-change
                            m1 = group_means.get(g1, np.nan)
                            m2 = group_means.get(g2, np.nan)
                            fc, log2fc = _safe_fc(m2, m1)
                            
                            result[f'{comp_name}_FC'] = fc
                            result[f'{comp_name}_log2FC'] = log2fc
                            result[f'{comp_name}_neg_log10_adj_p'] = -np.log10(p_adj) if p_adj > 0 else 15

                            # Nonparametric effect summary for consistency in output columns.
                            v1 = group_values_dict.get(g1, np.array([]))
                            v2 = group_values_dict.get(g2, np.array([]))
                            eff, se_eff, ci_low, ci_high = _mannwhitney_effect_ci(v1, v2)
                            result[f'{comp_name}_model_effect'] = eff
                            result[f'{comp_name}_model_se'] = se_eff
                            result[f'{comp_name}_ci_lower_95'] = ci_low
                            result[f'{comp_name}_ci_upper_95'] = ci_high
            
            else:
                # Pairwise Welch t-test or Mann-Whitney U test (only valid groups)
                for g1, g2 in valid_pairwise_combos:
                    v1 = group_values_dict[g1]
                    v2 = group_values_dict[g2]
                    
                    if len(v1) < min_group_size or len(v2) < min_group_size:
                        continue
                    
                    if pairwise_test == 'welch':
                        _, p = stats.ttest_ind(v1, v2, equal_var=False)
                        eff, se_eff, ci_low, ci_high = _welch_effect_ci(v1, v2)
                    else:  # mannwhitney
                        _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                        eff, se_eff, ci_low, ci_high = _mannwhitney_effect_ci(v1, v2)
                    
                    p_adj = _floor_pval(p)  # Will be corrected globally later if FDR enabled
                    comp_name = f"{g1}_vs_{g2}"
                    result[f'{comp_name}_adj_p'] = p_adj
                    
                    # Calculate fold-change
                    m1 = np.mean(v1)
                    m2 = np.mean(v2)
                    fc, log2fc = _safe_fc(m2, m1)
                    
                    result[f'{comp_name}_FC'] = fc
                    result[f'{comp_name}_log2FC'] = log2fc
                    result[f'{comp_name}_neg_log10_adj_p'] = -np.log10(p_adj) if p_adj > 0 else 15
                    result[f'{comp_name}_model_effect'] = eff
                    result[f'{comp_name}_model_se'] = se_eff
                    result[f'{comp_name}_ci_lower_95'] = ci_low
                    result[f'{comp_name}_ci_upper_95'] = ci_high
        
        except Exception as e:
            pass  # Silently skip failed tests
        
        return result

    # Process all metabolites
    items = list(df.iterrows())
    total = len(items)

    # Stability guard: Tukey/Dunn posthoc can become unresponsive under threaded
    # execution in some environments; force serial mode for these paths.
    effective_n_jobs = max(1, int(n_jobs)) if n_jobs else 1
    if pairwise_test in ('tukey', 'dunn'):
        effective_n_jobs = 1
    
    if effective_n_jobs > 1:
        with ThreadPoolExecutor(max_workers=effective_n_jobs) as ex:
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

    # ============================================================
    # FDR CORRECTION (BH method) - Applied once across metabolites
    # Matching R: p.adjust(p_values, method = "BH")
    # ============================================================
    if fdr and not complete_result.empty and multipletests is not None:
        # Correct overall test p-values
        mask = complete_result['overall_p'].notna()
        if mask.any():
            try:
                _, padj_vec, _, _ = multipletests(
                    complete_result.loc[mask, 'overall_p'],
                    method='fdr_bh'
                )
                complete_result.loc[mask, 'overall_adj_p'] = padj_vec
            except Exception as e:
                pass
        
        # For pairwise tests: Tukey/Dunn p-values are already adjusted
        # Welch/Mann-Whitney may need global correction
        if pairwise_test in ('welch', 'mannwhitney'):
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                p_col = f'{comp_name}_adj_p'
                mask = complete_result[p_col].notna()
                if mask.any():
                    try:
                        _, padj_vec, _, _ = multipletests(
                            complete_result.loc[mask, p_col],
                            method='fdr_bh'
                        )
                        complete_result.loc[mask, p_col] = padj_vec
                    except:
                        pass
    
    # Note: Tukey HSD and Dunn test p-values are already adjusted
    # Our implementation matches R: we don't apply additional FDR to these

    # ============================================================
    # DYNAMIC P-VALUE FLOORING (per comparison column)
    # floor = 0.01 * smallest positive p-value in that column
    # ============================================================
    pval_cols = []
    for col in complete_result.columns:
        if col in ('overall_p', 'overall_adj_p'):
            pval_cols.append(col)
        elif col.endswith('_adj_p'):
            pval_cols.append(col)

    for col in pval_cols:
        complete_result[col] = _apply_dynamic_pvalue_floor(complete_result[col])

    # Recompute -log10 columns from the floored adj p-values.
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
    # ============================================================
    individual_sheets = {}
    
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        
        # Get samples from each group
        samples_g1 = group_samples[g1]
        samples_g2 = group_samples[g2]
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
        
        # Drop rows where comparison has no result
        p_adj_col = f'{comp_name}_adj_p'
        if p_adj_col in sheet_df.columns:
            sheet_df = sheet_df.dropna(subset=[p_adj_col])
        
        if not sheet_df.empty:
            individual_sheets[comp_name] = sheet_df

    # ============================================================
    # COLUMN ORDERING
    # Order: Metabolite info + abundance + overall stats + pairwise comparisons
    # ============================================================
    ordered_cols = []
    
    # 1. Metabolite info columns
    ordered_cols.extend(info_cols)
    
    # 2. Abundance columns
    ordered_cols.extend([s for s in sample_cols if s in complete_result.columns])
    
    # 3. Overall test columns
    overall_cols = ['overall_statistic', 'overall_p', 'overall_adj_p']
    ordered_cols.extend([c for c in overall_cols if c in complete_result.columns])
    
    # 4. Pairwise comparisons
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        pairwise_cols = [
            f'{comp_name}_adj_p',
            f'{comp_name}_FC',
            f'{comp_name}_log2FC',
            f'{comp_name}_neg_log10_adj_p',
            f'{comp_name}_model_effect',
            f'{comp_name}_model_se',
            f'{comp_name}_ci_lower_95',
            f'{comp_name}_ci_upper_95'
        ]
        ordered_cols.extend([c for c in pairwise_cols if c in complete_result.columns])
    
    # Add any remaining columns not already included
    for c in complete_result.columns:
        if c not in ordered_cols:
            ordered_cols.append(c)
    
    complete_result = complete_result[ordered_cols]

    return complete_result, individual_sheets
