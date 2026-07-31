"""
Covariate Adjustment Module for MetaboGraph

This module provides functionality to adjust metabolite intensities for covariates
using linear models (OLS regression), similar to MetaboAnalyst's approach.

Features:
- Support for continuous and categorical covariates
- Automatic detection of covariates from main data or separate file
- Group comparison with covariate adjustment
- FDR correction (Benjamini-Hochberg)
- Export of model coefficients and adjusted p-values

Design:
- Can detect covariates embedded in the main metabolite matrix
- Can load covariates from a separate file
- Handles missing values appropriately
- Provides comprehensive output including model diagnostics
"""

from __future__ import annotations
import pandas as pd
import numpy as np
import os
from typing import Dict, List, Optional, Tuple, Any, Sequence
from dataclasses import dataclass
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', message='.*kurtosistest.*')

# Try to import statsmodels
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
    sm = None
    smf = None

# Try to import scipy for limma implementation (pure Python)
try:
    from scipy import stats as scipy_stats
    from scipy.special import digamma, polygamma
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    scipy_stats = None
    digamma = None
    polygamma = None


@dataclass
class CovariateAnalysisResult:
    """Container for covariate-adjusted analysis results"""
    metabolite_results: pd.DataFrame  # Main results with p-values per metabolite
    model_diagnostics: pd.DataFrame   # R², F-statistic, etc. per metabolite
    coefficient_table: pd.DataFrame   # All coefficients from all models
    adjusted_intensities: Optional[pd.DataFrame] = None  # Covariate-adjusted values
    summary_stats: Optional[Dict[str, Any]] = None  # Overall summary


def detect_covariate_columns(
    df: pd.DataFrame,
    sample_cols: List[str],
    exclude_cols: Optional[List[str]] = None
) -> Tuple[List[str], List[str]]:
    """
    Automatically detect potential covariate columns in the dataframe.
    
    Covariates are typically:
    - Age, Sex, Gender, BMI, Weight, Height
    - PMI (Post-mortem interval)
    - Batch, Run, Date
    - Any non-metabolite metadata columns
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    sample_cols : List[str]
        Known sample/intensity columns to exclude
    exclude_cols : Optional[List[str]]
        Additional columns to exclude from covariate detection
        
    Returns
    -------
    Tuple[List[str], List[str]]
        (continuous_covariates, categorical_covariates)
    """
    if exclude_cols is None:
        exclude_cols = []
    
    # Common covariate names
    common_continuous = ['age', 'pmi', 'bmi', 'weight', 'height', 'time', 'duration']
    common_categorical = ['sex', 'gender', 'batch', 'group', 'cohort', 'site', 'race', 'ethnicity']
    
    continuous_covars = []
    categorical_covars = []
    
    # Exclude sample columns and metabolite identifiers
    excluded = set(sample_cols + exclude_cols)
    metadata_candidates = [c for c in df.columns if c not in excluded]
    
    for col in metadata_candidates:
        col_lower = col.lower()
        
        # Skip if likely a metabolite identifier
        if col_lower in ['name', 'metabolite', 'compound', 'molecule', 'hmdb', 'hmdb_id', 
                         'kegg', 'pubchem', 'inchi', 'smiles', 'formula', 'mass']:
            continue
        
        # Check if it's a known covariate type
        is_continuous = any(keyword in col_lower for keyword in common_continuous)
        is_categorical = any(keyword in col_lower for keyword in common_categorical)
        
        # Analyze column data
        series = df[col].dropna()
        if len(series) == 0:
            continue
            
        # Check if numeric
        try:
            numeric_series = pd.to_numeric(series, errors='coerce')
            numeric_ratio = numeric_series.notna().sum() / len(series)
            
            if numeric_ratio > 0.8:  # Mostly numeric
                unique_count = numeric_series.nunique()
                if unique_count <= 10 and not is_continuous:
                    # Few unique values -> categorical (unless explicitly continuous)
                    categorical_covars.append(col)
                else:
                    # Many unique values -> continuous
                    continuous_covars.append(col)
            else:
                # Not numeric -> categorical
                unique_count = series.nunique()
                if unique_count <= 20:  # Reasonable for categorical
                    categorical_covars.append(col)
        except Exception:
            # If conversion fails, treat as categorical if unique count is low
            unique_count = series.nunique()
            if unique_count <= 20:
                categorical_covars.append(col)
    
    return continuous_covars, categorical_covars


def load_covariate_file(
    file_path: str,
    sample_id_col: Optional[str] = None
) -> pd.DataFrame:
    """
    Load covariate data from a separate file.
    
    Supports Excel (.xlsx, .xls), CSV (.csv), and TSV (.txt, .tsv) formats.
    
    Parameters
    ----------
    file_path : str
        Path to the covariate file
    sample_id_col : Optional[str]
        Column name containing sample IDs (default: first column)
        
    Returns
    -------
    pd.DataFrame
        Covariate dataframe with sample IDs as index
    """
    # Determine file type and read
    if file_path.endswith(('.xlsx', '.xls')):
        df_cov = pd.read_excel(file_path)
    elif file_path.endswith('.csv'):
        df_cov = pd.read_csv(file_path)
    elif file_path.endswith(('.txt', '.tsv')):
        df_cov = pd.read_csv(file_path, sep='\t')
    else:
        raise ValueError(f"Unsupported file format: {file_path}")
    
    # Set index to sample ID column
    if sample_id_col and sample_id_col in df_cov.columns:
        df_cov = df_cov.set_index(sample_id_col)
    elif sample_id_col is None and len(df_cov.columns) > 0:
        # Use first column as index
        df_cov = df_cov.set_index(df_cov.columns[0])
    
    return df_cov


def prepare_design_matrix(
    sample_cols: List[str],
    group_map: Dict[str, str],
    covariate_data: pd.DataFrame,
    covariate_cols: List[str],
    group_var_name: str = 'Group'
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """
    Create design matrix (X) for regression with group and covariates.
    
    Parameters
    ----------
    sample_cols : List[str]
        Sample column names from metabolite data
    group_map : Dict[str, str]
        Mapping of sample_col -> group label
    covariate_data : pd.DataFrame
        DataFrame with covariates (samples as rows or columns)
    covariate_cols : List[str]
        Names of covariate columns to include
    group_var_name : str
        Name for the group variable in design matrix
        
    Returns
    -------
    Tuple[pd.DataFrame, List[str], List[str]]
        (design_matrix, categorical_columns, included_covariates)
    """
    # Create DataFrame with sample IDs
    design_df = pd.DataFrame({
        'Sample': sample_cols,
        group_var_name: [group_map[s] for s in sample_cols]
    })
    
    # Add covariates
    included_covariates = []
    for cov in covariate_cols:
        if cov in covariate_data.columns:
            # Covariates are in columns (samples as rows)
            values = []
            for sample in sample_cols:
                if sample in covariate_data.index:
                    values.append(covariate_data.loc[sample, cov])
                else:
                    values.append(np.nan)
            design_df[cov] = values
            included_covariates.append(cov)
        elif cov in covariate_data.index:
            # Covariates are in rows (samples as columns)
            values = []
            for sample in sample_cols:
                if sample in covariate_data.columns:
                    values.append(covariate_data.loc[cov, sample])
                else:
                    values.append(np.nan)
            design_df[cov] = values
            included_covariates.append(cov)
    
    # Identify categorical columns
    categorical_cols = [group_var_name]
    for col in included_covariates:
        if col in design_df.columns:
            # Check if categorical
            try:
                pd.to_numeric(design_df[col], errors='raise')
            except (ValueError, TypeError):
                categorical_cols.append(col)
            else:
                # Numeric but few unique values -> treat as categorical
                if design_df[col].nunique() <= 5:
                    categorical_cols.append(col)
    
    return design_df, categorical_cols, included_covariates


def _apply_group_reference_coding(
    design_df: pd.DataFrame,
    group_var_name: str,
    reference_group: Optional[str]
) -> Optional[str]:
    """Apply explicit group baseline coding and return the effective reference."""
    if group_var_name not in design_df.columns:
        return None

    observed_groups = [g for g in pd.unique(design_df[group_var_name]) if pd.notna(g)]
    if not observed_groups:
        return None

    # Keep behavior deterministic: requested reference first, otherwise alphabetical first.
    sorted_groups = sorted(observed_groups, key=lambda x: str(x))
    if reference_group in observed_groups:
        categories = [reference_group] + [g for g in sorted_groups if g != reference_group]
    else:
        categories = sorted_groups

    design_df[group_var_name] = pd.Categorical(
        design_df[group_var_name],
        categories=categories,
        ordered=True
    )
    return categories[0]


def run_covariate_adjusted_analysis(
    df_intensities: pd.DataFrame,
    sample_cols: List[str],
    group_map: Dict[str, str],
    covariate_data: Optional[pd.DataFrame] = None,
    covariate_cols: Optional[List[str]] = None,
    *,
    group_var_name: str = 'Group',
    reference_group: Optional[str] = None,
    apply_fdr: bool = True,
    fdr_method: str = 'BH',
    alpha: float = 0.05,
    metabolite_id_col: Optional[str] = None,
    return_adjusted_intensities: bool = False,
    group_order: Optional[List[str]] = None,
    min_samples_per_group: int = 2,
    min_samples_type: str = 'absolute'
) -> CovariateAnalysisResult:
    """
    Run covariate-adjusted analysis for metabolite data.
    
    For each metabolite, fits an OLS model:
        metabolite_intensity ~ Group + Covariate1 + Covariate2 + ...
    
    Tests the significance of the Group effect while controlling for covariates.
    Supports multiple groups and generates all pairwise comparisons.
    
    Can also run without covariates (simple group comparison with OLS).
    
    Parameters
    ----------
    df_intensities : pd.DataFrame
        Metabolite intensity matrix (metabolites as rows, samples as columns)
    sample_cols : List[str]
        Sample column names to include
    group_map : Dict[str, str]
        Mapping of sample_col -> group label
    covariate_data : Optional[pd.DataFrame]
        Covariate values (samples as rows or columns). Can be None for no covariates.
    covariate_cols : Optional[List[str]]
        Names of covariates to include in model. Can be None or empty for no covariates.
    group_var_name : str
        Name for group variable (default: 'Group')
    reference_group : Optional[str]
        Reference group for dummy coding (if None, uses alphabetical first)
    apply_fdr : bool
        Apply FDR correction (if False, uses uncorrected p-values)
    fdr_method : str
        FDR correction method: 'BH' (Benjamini-Hochberg), 'Bonferroni', 'Holm', 'Hochberg', 'BY', 'None'
    alpha : float
        Significance threshold
    metabolite_id_col : Optional[str]
        Column name for metabolite identifiers
    return_adjusted_intensities : bool
        Whether to compute and return residual-adjusted intensities
    group_order : Optional[List[str]]
        Order of groups for pairwise comparisons (default: sorted)
    min_samples_per_group : int
        Minimum valid samples required per group (default: 2)
    min_samples_type : str
        Type of threshold: 'absolute' (count) or 'percentage' (default: 'absolute')
        
    Returns
    -------
    CovariateAnalysisResult
        Complete results including p-values, coefficients, diagnostics
    """
    if not STATSMODELS_AVAILABLE:
        raise ImportError(
            "statsmodels is required for covariate adjustment. "
            "Install it with: pip install statsmodels"
        )
    
    # Handle case with no covariates
    if covariate_cols is None or len(covariate_cols) == 0:
        covariate_cols = []
        covariate_data = pd.DataFrame(index=sample_cols)  # Empty dataframe
    
    # Determine metabolite ID column
    if metabolite_id_col is None:
        id_candidates = ['Name', 'Metabolite', 'Compound', 'Molecule']
        for candidate in id_candidates:
            if candidate in df_intensities.columns:
                metabolite_id_col = candidate
                break
        if metabolite_id_col is None:
            metabolite_id_col = df_intensities.columns[0]
    
    # Get unique groups in specified order
    unique_groups = list(dict.fromkeys(group_map.values()))  # Preserve order
    if group_order:
        # Use specified order, adding any missing groups at the end
        ordered = [g for g in group_order if g in unique_groups]
        ordered += [g for g in unique_groups if g not in ordered]
        unique_groups = ordered
    else:
        unique_groups = sorted(unique_groups)
    
    # Generate pairwise combinations
    pairwise_combos = []
    for i, g1 in enumerate(unique_groups):
        for g2 in unique_groups[i+1:]:
            pairwise_combos.append((g1, g2))
    
    # Prepare design matrix
    if len(covariate_cols) > 0:
        design_df, categorical_cols, included_covariates = prepare_design_matrix(
            sample_cols, group_map, covariate_data, covariate_cols, group_var_name
        )
    else:
        # No covariates - simple group design
        design_df = pd.DataFrame({
            'Sample': sample_cols,
            group_var_name: [group_map[s] for s in sample_cols]
        })
        categorical_cols = [group_var_name]
        included_covariates = []

    missing_covariates = [c for c in covariate_cols if c not in included_covariates]

    effective_reference_group = _apply_group_reference_coding(
        design_df,
        group_var_name,
        reference_group
    )
    
    # Convert categorical variables to dummy variables
    X = pd.get_dummies(
        design_df.drop(columns=['Sample']),
        columns=categorical_cols,
        drop_first=True,
        dtype=float
    )
    
    # Add constant term
    X = sm.add_constant(X)
    
    # Identify group coefficient columns
    group_cols = [col for col in X.columns if col.startswith(f'{group_var_name}_')]
    
    # Run model for each metabolite
    results_list = []
    diagnostics_list = []
    all_coefficients = []
    adjusted_intensities_dict = {}
    pairwise_results = {comp: [] for comp in pairwise_combos}
    
    n_metabolites = len(df_intensities)
    n_excluded_insufficient = 0
    n_excluded_invalid_values = 0  # OLS handles invalid values per-row, so this stays 0
    
    for idx, row in df_intensities.iterrows():
        metabolite_id = row[metabolite_id_col] if (metabolite_id_col and metabolite_id_col in df_intensities.columns) else f"Metabolite_{idx}"
        
        # Initialize result row with metabolite metadata and sample values
        # Only include metabolite_id_col as a key if it's not None
        result_row = {}
        if metabolite_id_col:
            result_row[metabolite_id_col] = metabolite_id
        else:
            result_row['metabolite_id'] = metabolite_id
        
        # Add all non-sample columns from original data
        for col in df_intensities.columns:
            if col not in sample_cols and col != metabolite_id_col:
                result_row[col] = row[col]
        
        # Add sample intensity values
        for sample in sample_cols:
            result_row[sample] = row[sample]
        
        # PRE-FILTER: Check if each pairwise comparison has sufficient valid data
        # If min_samples_type is None (imputation mode), skip filtering (already done in imputation step)
        skip_metabolite = False
        group_valid_counts = {}
        
        if min_samples_type is not None and min_samples_per_group is not None:
            # Standard filtering mode: apply min_samples thresholds
            for group in unique_groups:
                group_samples = [s for s in sample_cols if group_map.get(s) == group]
                # Count only valid (non-zero, non-NaN) values
                valid_values = [row[s] for s in group_samples 
                               if pd.notna(row[s]) and row[s] != 0]
                group_valid_counts[group] = len(valid_values)
            
            # Check each pairwise comparison using user-defined threshold
            for g1, g2 in pairwise_combos:
                # Calculate required minimum samples for each group
                g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
                g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
                
                if min_samples_type == 'percentage':
                    # Percentage-based threshold
                    g1_min_required = max(1, int(np.ceil(len(g1_samples) * min_samples_per_group / 100.0)))
                    g2_min_required = max(1, int(np.ceil(len(g2_samples) * min_samples_per_group / 100.0)))
                else:
                    # Absolute count threshold
                    g1_min_required = min_samples_per_group
                    g2_min_required = min_samples_per_group
                
                if group_valid_counts.get(g1, 0) < g1_min_required or group_valid_counts.get(g2, 0) < g2_min_required:
                    skip_metabolite = True
                    break
        else:
            # Imputation mode: filtering already applied, skip here
            for group in unique_groups:
                group_samples = [s for s in sample_cols if group_map.get(s) == group]
                valid_values = [row[s] for s in group_samples 
                               if pd.notna(row[s]) and row[s] != 0]
                group_valid_counts[group] = len(valid_values)
        
        # Store actual valid n per group
        for group in unique_groups:
            result_row[f'n_{group}'] = group_valid_counts.get(group, 0)
        
        if skip_metabolite:
            # Insufficient valid data - skip statistics but KEEP in complete results
            n_excluded_insufficient += 1
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': 'insufficient_valid_data'
            })
            # Add NaN for all pairwise stats
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
            # Add to results_list so it appears in complete sheet
            results_list.append(result_row)
            # Skip statistical computation
            continue
        
        # Extract intensity values for model
        y = row[sample_cols].values.astype(float)
        
        # Remove missing/zero values for regression (zeros are treated as missing)
        valid_mask = ~(np.isnan(y) | (y == 0) | np.isnan(X).any(axis=1))
        y_clean = y[valid_mask]
        X_clean = X[valid_mask]
        
        min_obs_required = X.shape[1] + 1  # one residual df beyond full rank model
        if len(y_clean) < min_obs_required:
            # Insufficient data for model
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': 'insufficient_data'
            })
            # Add NaN for all pairwise stats
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_ci_lower_95'] = np.nan
                result_row[f'{comp_name}_ci_upper_95'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
            results_list.append(result_row)
            continue
        
        try:
            # Fit OLS model
            model = sm.OLS(y_clean, X_clean).fit()
            
            # Store model quality metrics
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'r_squared': model.rsquared,
                'adj_r_squared': model.rsquared_adj,
                'f_statistic': model.fvalue,
                'f_pvalue': model.f_pvalue,
                'status': 'success'
            })
            
            # Calculate pairwise comparisons
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                
                # Get sample indices for each group
                g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
                g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
                
                # Get intensity values for each group - ONLY valid values (non-zero, non-NaN)
                g1_values = [row[s] for s in g1_samples 
                            if pd.notna(row[s]) and row[s] != 0]
                g2_values = [row[s] for s in g2_samples 
                            if pd.notna(row[s]) and row[s] != 0]
                
                # Calculate group means from valid values only
                g1_mean = np.mean(g1_values) if g1_values else np.nan
                g2_mean = np.mean(g2_values) if g2_values else np.nan
                result_row[f'{g1}_Mean'] = g1_mean
                result_row[f'{g2}_Mean'] = g2_mean
                
                # Calculate FC and log2FC
                if not np.isnan(g1_mean) and not np.isnan(g2_mean) and g1_mean > 0:
                    fc = g2_mean / g1_mean
                    result_row[f'{comp_name}_FC'] = fc
                    result_row[f'{comp_name}_log2FC'] = np.log2(fc) if fc > 0 else np.nan
                else:
                    result_row[f'{comp_name}_FC'] = np.nan
                    result_row[f'{comp_name}_log2FC'] = np.nan
                
                # Test pairwise comparison using t-test on coefficients.
                # Handle scalar/array statsmodels outputs robustly across versions.
                try:
                    from scipy import stats as scipy_stats
                    
                    g1_col = f'{group_var_name}_{g1}'
                    g2_col = f'{group_var_name}_{g2}'
                    conf_int = model.conf_int(alpha=0.05)  # 95% CI
                    
                    # Initialize model-based variables
                    model_effect = np.nan
                    model_se = np.nan
                    model_ci_lower = np.nan
                    model_ci_upper = np.nan

                    if g1_col in model.params.index and g2_col in model.params.index:
                        # Both are dummy variables - test the contrast (g2 - g1).
                        contrast = np.zeros(len(model.params))
                        contrast[model.params.index.get_loc(g2_col)] = 1
                        contrast[model.params.index.get_loc(g1_col)] = -1
                        t_test = model.t_test(contrast)

                        p_raw = np.asarray(t_test.pvalue)
                        pvalue = float(p_raw.reshape(-1)[0]) if p_raw.size else np.nan

                        ci_raw = np.asarray(t_test.conf_int(alpha=0.05))
                        if ci_raw.ndim == 2 and ci_raw.shape[1] >= 2:
                            ci_lower = float(ci_raw[0, 0])
                            ci_upper = float(ci_raw[0, 1])
                        elif ci_raw.ndim == 1 and ci_raw.size >= 2:
                            ci_lower = float(ci_raw[0])
                            ci_upper = float(ci_raw[1])
                        else:
                            ci_lower = np.nan
                            ci_upper = np.nan
                        
                        # Extract model-based effect and SE for contrast
                        model_effect = float(np.asarray(t_test.effect_size).reshape(-1)[0]) if hasattr(t_test, 'effect_size') else np.nan
                        if np.isnan(model_effect):
                            # Compute effect manually as g2 - g1
                            model_effect = float(model.params[g2_col] - model.params[g1_col])
                        
                        # SE for contrast: use t_test object which properly accounts for covariance
                        # The t_test.sd attribute gives the standard deviation of the contrast
                        if hasattr(t_test, 'sd'):
                            model_se = float(np.asarray(t_test.sd).reshape(-1)[0])
                        else:
                            # Fallback: approximate SE (ignores covariance, less accurate)
                            se_g2 = float(model.bse[g2_col])
                            se_g1 = float(model.bse[g1_col])
                            model_se = np.sqrt(se_g2**2 + se_g1**2)
                        
                        # Compute CI from effect ± t_crit * SE
                        df_resid = model.df_resid
                        t_crit = scipy_stats.t.ppf(0.975, df_resid)  # two-tailed 0.05
                        model_ci_lower = model_effect - t_crit * model_se
                        model_ci_upper = model_effect + t_crit * model_se
                        
                    elif g2_col in model.params.index:
                        # g1 is reference, g2 is dummy - test g2 coefficient.
                        pvalue = float(model.pvalues[g2_col])
                        ci_lower = conf_int.loc[g2_col, 0]
                        ci_upper = conf_int.loc[g2_col, 1]
                        
                        # Model-based: use coefficient and its SE/CI directly
                        model_effect = float(model.params[g2_col])
                        model_se = float(model.bse[g2_col])
                        model_ci_lower = conf_int.loc[g2_col, 0]
                        model_ci_upper = conf_int.loc[g2_col, 1]
                        
                    elif g1_col in model.params.index:
                        # g2 is reference, g1 is dummy - test g1 coefficient.
                        pvalue = float(model.pvalues[g1_col])
                        ci_lower = conf_int.loc[g1_col, 0]
                        ci_upper = conf_int.loc[g1_col, 1]
                        
                        # CRITICAL FIX: g1 coefficient represents (μ_g1 - μ_g2) since g2 is reference
                        # But g1_vs_g2 should represent (μ_g2 - μ_g1), so NEGATE the effect
                        model_effect = -float(model.params[g1_col])
                        model_se = float(model.bse[g1_col])
                        # Flip CI bounds due to negation of effect
                        model_ci_lower = -conf_int.loc[g1_col, 1]
                        model_ci_upper = -conf_int.loc[g1_col, 0]
                        
                    else:
                        # Both are reference (shouldn't happen)
                        pvalue = np.nan
                        ci_lower = np.nan
                        ci_upper = np.nan
                except Exception as e:
                    pvalue = np.nan
                    ci_lower = np.nan
                    ci_upper = np.nan
                    model_effect = np.nan
                    model_se = np.nan
                    model_ci_lower = np.nan
                    model_ci_upper = np.nan
                
                result_row[f'{comp_name}_pvalue'] = pvalue
                result_row[f'{comp_name}_ci_lower_95'] = ci_lower
                result_row[f'{comp_name}_ci_upper_95'] = ci_upper
                result_row[f'{comp_name}_adj_p'] = np.nan  # FDR filled later
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan  # Filled after FDR
                
                # Add model-based columns (coefficient + SE from model)
                result_row[f'{comp_name}_model_effect'] = model_effect
                result_row[f'{comp_name}_model_se'] = model_se
                
                # Store for pairwise sheet
                pairwise_row = result_row.copy()
                pairwise_row[f'{comp_name}_pvalue'] = pvalue
                pairwise_results[(g1, g2)].append(pairwise_row)
            
            results_list.append(result_row)
            
            # Store diagnostics
            diagnostics_list.append({
                metabolite_id_col: metabolite_id,
                'r_squared': model.rsquared,
                'adj_r_squared': model.rsquared_adj,
                'f_statistic': model.fvalue,
                'f_pvalue': model.f_pvalue,
                'aic': model.aic,
                'bic': model.bic,
                'log_likelihood': model.llf,
                'n_obs': int(model.nobs),
                'df_model': int(model.df_model),
                'df_resid': int(model.df_resid)
            })
            
            # Store coefficients with confidence intervals
            conf_int = model.conf_int(alpha=0.05)  # 95% confidence interval
            for coef_name, coef_value in model.params.items():
                # Get confidence interval bounds for this coefficient
                ci_lower = conf_int.loc[coef_name, 0] if coef_name in conf_int.index else np.nan
                ci_upper = conf_int.loc[coef_name, 1] if coef_name in conf_int.index else np.nan
                
                all_coefficients.append({
                    metabolite_id_col: metabolite_id,
                    'coefficient': coef_name,
                    'value': coef_value,
                    'std_err': model.bse[coef_name] if coef_name in model.bse else np.nan,
                    't_statistic': model.tvalues[coef_name] if coef_name in model.tvalues else np.nan,
                    'pvalue': model.pvalues[coef_name] if coef_name in model.pvalues else np.nan,
                    'ci_lower_95': ci_lower,
                    'ci_upper_95': ci_upper
                })
            
            # Compute adjusted intensities (residuals + group effect)
            if return_adjusted_intensities:
                # Remove covariate effects but keep group effect
                covar_cols = [c for c in X_clean.columns if c not in group_cols and c != 'const']
                residuals = y_clean - model.predict(X_clean)
                
                # CRITICAL FIX: Correct indexing for group effects
                # group_effect[valid_mask][i] creates a copy, so assignment doesn't update original
                group_effect = np.zeros(len(y_clean))
                valid_indices = np.where(valid_mask)[0]  # Get indices where valid_mask is True
                
                for j, i in enumerate(valid_indices):
                    # j: index in the valid subset (y_clean, residuals, group_effect)
                    # i: index in original sample_cols
                    sample = sample_cols[i]
                    group = group_map[sample]
                    # Find corresponding group dummy column
                    group_col = f'{group_var_name}_{group}'
                    if group_col in model.params.index:
                        # j indexes position in valid arrays (group_effect, residuals)
                        group_effect[j] = model.params[group_col]
                
                adjusted_vals = residuals + group_effect + model.params['const']

                # Re-expand to the original sample order so the output always matches sample_cols.
                adjusted_full = np.full(len(sample_cols), np.nan, dtype=float)
                adjusted_full[valid_mask] = adjusted_vals
                adjusted_intensities_dict[metabolite_id] = adjusted_full
                
        except Exception as e:
            # Model fitting failed
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': f'error: {str(e)}'
            })
            # Add NaN for all pairwise stats
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_ci_lower_95'] = np.nan
                result_row[f'{comp_name}_ci_upper_95'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{comp_name}_model_effect'] = np.nan
                result_row[f'{comp_name}_model_se'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
            results_list.append(result_row)
    
    # Convert to DataFrame
    df_results = pd.DataFrame(results_list)
    df_diagnostics = pd.DataFrame(diagnostics_list)
    df_coefficients = pd.DataFrame(all_coefficients)
    
    # Apply p-value correction for each pairwise comparison
    # Always create adj_p columns - if apply_fdr=False, just copy raw p-values
    if len(df_results) > 0:
        for g1, g2 in pairwise_combos:
            comp_name = f"{g1}_vs_{g2}"
            pval_col = f'{comp_name}_pvalue'
            adj_pval_col = f'{comp_name}_adj_p'
            neg_log_col = f'{comp_name}_neg_log10_adj_p'
            
            if pval_col in df_results.columns:
                valid_pvals = df_results[pval_col].dropna()
                if len(valid_pvals) > 0:
                    # Apply selected correction method
                    if not apply_fdr or fdr_method == 'None':
                        # No correction - use raw p-values
                        df_results.loc[valid_pvals.index, adj_pval_col] = valid_pvals.values
                    elif fdr_method == 'BH':
                        try:
                            from scipy.stats import false_discovery_control
                            adjusted = false_discovery_control(valid_pvals.values, method='bh')
                            df_results.loc[valid_pvals.index, adj_pval_col] = adjusted
                        except Exception:
                            # Fallback to manual BH correction
                            n = len(valid_pvals)
                            sorted_idx = valid_pvals.argsort()
                            sorted_pvals = valid_pvals.iloc[sorted_idx]
                            ranks = np.arange(1, n + 1)
                            adjusted = np.minimum.accumulate((sorted_pvals * n / ranks)[::-1])[::-1]
                            adjusted = np.clip(adjusted, 0, 1)
                            
                            # Map back to original order
                            adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                            df_results.loc[valid_pvals.index, adj_pval_col] = [
                                adj_dict[i] for i in valid_pvals.index
                            ]
                    elif fdr_method == 'Bonferroni':
                        adjusted = np.minimum(valid_pvals.values * len(valid_pvals), 1.0)
                        df_results.loc[valid_pvals.index, adj_pval_col] = adjusted
                    elif fdr_method == 'Holm':
                        n = len(valid_pvals)
                        sorted_idx = valid_pvals.argsort()
                        sorted_pvals = valid_pvals.iloc[sorted_idx]
                        ranks = np.arange(1, n + 1)
                        adjusted = np.maximum.accumulate(np.minimum(sorted_pvals.values * (n - ranks + 1), 1.0))
                        adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                        df_results.loc[valid_pvals.index, adj_pval_col] = [adj_dict[i] for i in valid_pvals.index]
                    elif fdr_method == 'Hochberg':
                        # Hochberg step-up: sort ascending, multiply by (n-i+1), then cummin from right
                        n = len(valid_pvals)
                        sorted_idx = valid_pvals.argsort()
                        sorted_pvals = valid_pvals.iloc[sorted_idx]
                        ranks = np.arange(1, n + 1)
                        # Multiply by (n - rank + 1) for step-up
                        adjusted = np.minimum(sorted_pvals.values * (n - ranks + 1), 1.0)
                        # Cumulative minimum from right to left (enforce monotonicity)
                        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
                        adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                        df_results.loc[valid_pvals.index, adj_pval_col] = [adj_dict[i] for i in valid_pvals.index]
                    elif fdr_method == 'BY':
                        # Benjamini-Yekutieli: like BH but with harmonic constant for dependent tests
                        n = len(valid_pvals)
                        sorted_idx = valid_pvals.argsort()
                        sorted_pvals = valid_pvals.iloc[sorted_idx]
                        ranks = np.arange(1, n + 1)
                        # Harmonic constant: sum(1/i) for i=1..n
                        cm = np.sum(1.0 / np.arange(1, n + 1))
                        # Multiply by (n * cm / rank)
                        adjusted = sorted_pvals.values * n * cm / ranks
                        # Cumulative minimum from right to left
                        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
                        adjusted = np.clip(adjusted, 0, 1)
                        adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                        df_results.loc[valid_pvals.index, adj_pval_col] = [adj_dict[i] for i in valid_pvals.index]
                    
                    # Calculate -log10(adj_p)
                    adj_pvals = df_results[adj_pval_col].replace(0, np.finfo(float).eps)
                    df_results[neg_log_col] = -np.log10(adj_pvals)
    
    # Mark significant metabolites for each comparison
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        adj_pval_col = f'{comp_name}_adj_p'
        sig_col = f'{comp_name}_significant'
        if adj_pval_col in df_results.columns:
            df_results[sig_col] = df_results[adj_pval_col] < alpha
    
    # Prepare adjusted intensities dataframe
    df_adjusted = None
    if return_adjusted_intensities and adjusted_intensities_dict:
        df_adjusted = pd.DataFrame(adjusted_intensities_dict).T
        df_adjusted.columns = sample_cols
    
    # Build pairwise comparison DataFrames - Use df_results (with FDR) instead of pairwise_results
    pairwise_dfs = {}
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        
        # Extract relevant columns from df_results (which has FDR-corrected values)
        if not df_results.empty:
            # Start with complete results
            pairwise_df = df_results.copy()
            
            # FILTER: Remove metabolites with NaN p-values or unsuccessful status
            pval_col = f'{comp_name}_pvalue'
            if pval_col in pairwise_df.columns:
                pairwise_df = pairwise_df[pairwise_df[pval_col].notna()]
            if 'status' in pairwise_df.columns:
                pairwise_df = pairwise_df[pairwise_df['status'] == 'success']
            
            if pairwise_df.empty:
                continue
            
            # Get samples for this comparison only
            g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
            g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
            comp_samples = g1_samples + g2_samples
            
            # Keep only relevant columns for this comparison
            metadata_cols = [c for c in df_results.columns 
                           if c not in sample_cols and c != metabolite_id_col]
            keep_cols = [metabolite_id_col] + comp_samples + metadata_cols
            
            # Filter to only columns that exist and are relevant to this comparison
            keep_cols = [c for c in keep_cols if c in pairwise_df.columns and 
                        (c == metabolite_id_col or c in comp_samples or 
                         comp_name in c or c.endswith('_Mean') or c.startswith('n_') or
                         c in ['r_squared', 'adj_r_squared', 'f_statistic', 'f_pvalue', 'status', 'n_samples'])]
            keep_cols = list(dict.fromkeys(keep_cols))  # Remove duplicates
            
            # Rename columns to remove comparison prefix for cleaner sheet
            pairwise_df = pairwise_df[keep_cols].copy()
            rename_map = {
                f'{comp_name}_pvalue': 'pvalue',
                f'{comp_name}_adj_p': 'adj_p',
                f'{comp_name}_FC': 'FC',
                f'{comp_name}_log2FC': 'log2FC',
                f'{comp_name}_neg_log10_adj_p': 'neg_log10_adj_p',
                f'{comp_name}_significant': 'significant'
            }
            pairwise_df = pairwise_df.rename(columns=rename_map)
            
            # Only add if there are metabolites left after filtering
            if not pairwise_df.empty:
                pairwise_dfs[comp_name] = pairwise_df
    
    # Summary statistics with actual sample counts per group
    sig_counts = {}
    tested_counts = {}
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        sig_col = f'{comp_name}_significant'
        if sig_col in df_results.columns:
            # Count only rows with valid statistics
            valid_rows = df_results[df_results['status'] == 'success']
            sig_counts[comp_name] = int(valid_rows[sig_col].sum())
            tested_counts[comp_name] = len(valid_rows)
    
    # Calculate actual sample sizes per group (from actual data)
    group_sample_counts = {}
    for group in unique_groups:
        group_samples = [s for s in sample_cols if group_map.get(s) == group]
        group_sample_counts[group] = len(group_samples)
    
    n_tested = len([r for r in results_list if r.get('status') == 'success'])
    
    summary = {
        'n_metabolites_total': n_metabolites,
        'n_metabolites_tested': n_tested,
        'n_metabolites_excluded': n_excluded_insufficient,
        'n_metabolites_excluded_invalid_values': n_excluded_invalid_values,
        'n_covariates_requested': len(covariate_cols),
        'n_covariates': len(included_covariates),
        'covariates': ', '.join(included_covariates) if included_covariates else 'None',
        'covariates_requested': ', '.join(covariate_cols) if covariate_cols else 'None',
        'covariates_missing': ', '.join(missing_covariates) if missing_covariates else 'None',
        'n_groups': len(unique_groups),
        'groups': ', '.join(unique_groups),
        'group_sample_counts': group_sample_counts,
        'pairwise_comparisons': pairwise_combos,
        'n_tested_per_comparison': tested_counts,
        'n_significant_per_comparison': sig_counts,
        'alpha': alpha,
        'fdr_method': 'Benjamini-Hochberg' if apply_fdr else 'None',
        'method': 'Linear Model (OLS)',
        'reference_group': effective_reference_group if effective_reference_group is not None else 'None',
        'design_matrix_columns': ', '.join(X.columns.astype(str).tolist())
    }
    
    # Store pairwise DataFrames in the result object
    result = CovariateAnalysisResult(
        metabolite_results=df_results,
        model_diagnostics=df_diagnostics,
        coefficient_table=df_coefficients,
        adjusted_intensities=df_adjusted,
        summary_stats=summary
    )
    
    # Add pairwise DataFrames as additional attribute
    result.pairwise_sheets = pairwise_dfs
    
    return result


def run_limma_covariate_analysis(
    df_intensities: pd.DataFrame,
    sample_cols: List[str],
    group_map: Dict[str, str],
    covariate_data: Optional[pd.DataFrame] = None,
    covariate_cols: Optional[List[str]] = None,
    *,
    group_var_name: str = 'Group',
    reference_group: Optional[str] = None,
    apply_fdr: bool = True,
    fdr_method: str = 'BH',
    alpha: float = 0.05,
    metabolite_id_col: Optional[str] = None,
    return_adjusted_intensities: bool = False,
    group_order: Optional[List[str]] = None,
    min_samples_per_group: int = 2,
    min_samples_type: str = 'absolute'
) -> CovariateAnalysisResult:
    """
    Run limma-style empirical Bayes covariate-adjusted analysis for metabolite data.
    
    This is a pure Python implementation of the limma approach (Smyth 2004) which uses 
    empirical Bayes to moderate (shrink) variance estimates across metabolites. This is 
    particularly useful when sample sizes are small, as it borrows information across 
    metabolites to produce more stable variance estimates.
    
    For each metabolite, fits a linear model:
        metabolite_intensity ~ Group + Covariate1 + Covariate2 + ...
    
    Then applies empirical Bayes moderation to:
    1. Estimate prior variance from all metabolites
    2. Shrink individual variance estimates toward the prior
    3. Compute moderated t-statistics with increased degrees of freedom
    
    Parameters
    ----------
    df_intensities : pd.DataFrame
        Metabolite intensity matrix (metabolites as rows, samples as columns)
    sample_cols : List[str]
        Sample column names to include
    group_map : Dict[str, str]
        Mapping of sample_col -> group label
    covariate_data : Optional[pd.DataFrame]
        Covariate values (samples as rows or columns). Can be None for no covariates.
    covariate_cols : Optional[List[str]]
        Names of covariates to include in model. Can be None or empty for no covariates.
    group_var_name : str
        Name for group variable (default: 'Group')
    reference_group : Optional[str]
        Reference group for contrasts (if None, uses alphabetical first)
    apply_fdr : bool
        Apply FDR correction (default: True)
    fdr_method : str
        FDR correction method: 'BH' (Benjamini-Hochberg), 'Bonferroni', etc.
    alpha : float
        Significance threshold
    metabolite_id_col : Optional[str]
        Column name for metabolite identifiers
    return_adjusted_intensities : bool
        Whether to compute and return residual-adjusted intensities
    group_order : Optional[List[str]]
        Order of groups for pairwise comparisons (default: sorted)
    min_samples_per_group : int
        Minimum valid samples required per group (default: 2)
    min_samples_type : str
        Type of threshold: 'absolute' (count) or 'percentage' (default: 'absolute')
        
    Returns
    -------
    CovariateAnalysisResult
        Complete results including moderated p-values, coefficients, diagnostics
        
    References
    ----------
    Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing 
    differential expression in microarray experiments. Statistical Applications 
    in Genetics and Molecular Biology, 3(1), Article 3.
    """
    if not STATSMODELS_AVAILABLE:
        raise ImportError(
            "statsmodels is required for limma analysis. "
            "Install it with: pip install statsmodels"
        )
    
    # Import scipy stats for t-distribution
    from scipy import stats as scipy_stats
    from scipy.special import digamma, polygamma
    
    # Handle case with no covariates
    if covariate_cols is None or len(covariate_cols) == 0:
        covariate_cols = []
        covariate_data = pd.DataFrame(index=sample_cols)
    
    # Determine metabolite ID column
    if metabolite_id_col is None:
        id_candidates = ['Name', 'Metabolite', 'Compound', 'Molecule']
        for candidate in id_candidates:
            if candidate in df_intensities.columns:
                metabolite_id_col = candidate
                break
        if metabolite_id_col is None:
            metabolite_id_col = df_intensities.columns[0]
    
    # Get unique groups in specified order
    unique_groups = list(dict.fromkeys(group_map.values()))
    if group_order:
        ordered = [g for g in group_order if g in unique_groups]
        ordered += [g for g in unique_groups if g not in ordered]
        unique_groups = ordered
    else:
        unique_groups = sorted(unique_groups)
    
    # Generate pairwise combinations
    pairwise_combos = []
    for i, g1 in enumerate(unique_groups):
        for g2 in unique_groups[i+1:]:
            pairwise_combos.append((g1, g2))
    
    # Prepare design matrix
    if len(covariate_cols) > 0:
        design_df, categorical_cols, included_covariates = prepare_design_matrix(
            sample_cols, group_map, covariate_data, covariate_cols, group_var_name
        )
    else:
        design_df = pd.DataFrame({
            'Sample': sample_cols,
            group_var_name: [group_map[s] for s in sample_cols]
        })
        categorical_cols = [group_var_name]
        included_covariates = []

    missing_covariates = [c for c in covariate_cols if c not in included_covariates]

    effective_reference_group = _apply_group_reference_coding(
        design_df,
        group_var_name,
        reference_group
    )
    
    # Convert categorical variables to dummy variables
    X = pd.get_dummies(
        design_df.drop(columns=['Sample']),
        columns=categorical_cols,
        drop_first=True,
        dtype=float
    )
    X = sm.add_constant(X)
    
    # Identify group coefficient columns
    group_cols = [col for col in X.columns if col.startswith(f'{group_var_name}_')]
    
    # ========== PHASE 1: Fit all models and collect statistics ==========
    model_stats = []
    results_list = []
    diagnostics_list = []
    all_coefficients = []
    
    n_metabolites = len(df_intensities)
    n_excluded_insufficient = 0
    
    for idx, row in df_intensities.iterrows():
        metabolite_id = row[metabolite_id_col] if (metabolite_id_col and metabolite_id_col in df_intensities.columns) else f"Metabolite_{idx}"
        
        # Initialize result row
        result_row = {}
        if metabolite_id_col:
            result_row[metabolite_id_col] = metabolite_id
        else:
            result_row['metabolite_id'] = metabolite_id
        
        # Add non-sample columns
        for col in df_intensities.columns:
            if col not in sample_cols and col != metabolite_id_col:
                result_row[col] = row[col]
        
        # Add sample intensity values
        for sample in sample_cols:
            result_row[sample] = row[sample]
        
        # Check valid data per group
        group_valid_counts = {}
        skip_metabolite = False
        
        if min_samples_type is not None and min_samples_per_group is not None:
            # Standard filtering mode: apply min_samples thresholds
            for group in unique_groups:
                group_samples = [s for s in sample_cols if group_map.get(s) == group]
                valid_values = [row[s] for s in group_samples 
                              if pd.notna(row[s]) and row[s] != 0]
                group_valid_counts[group] = len(valid_values)
            
            # Check minimum samples threshold
            for g1, g2 in pairwise_combos:
                g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
                g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
                
                if min_samples_type == 'percentage':
                    g1_min_required = max(1, int(np.ceil(len(g1_samples) * min_samples_per_group / 100.0)))
                    g2_min_required = max(1, int(np.ceil(len(g2_samples) * min_samples_per_group / 100.0)))
                else:
                    g1_min_required = min_samples_per_group
                    g2_min_required = min_samples_per_group
                
                if group_valid_counts.get(g1, 0) < g1_min_required or group_valid_counts.get(g2, 0) < g2_min_required:
                    skip_metabolite = True
                    break
        else:
            # Imputation mode: filtering already applied, skip here
            for group in unique_groups:
                group_samples = [s for s in sample_cols if group_map.get(s) == group]
                valid_values = [row[s] for s in group_samples 
                              if pd.notna(row[s]) and row[s] != 0]
                group_valid_counts[group] = len(valid_values)
            # Do not skip here: in imputation mode the pre-filter already enforced
            # minimum valid data, so rows should proceed to model fitting.
            skip_metabolite = False
        
        # Store n per group
        for group in unique_groups:
            result_row[f'n_{group}'] = group_valid_counts.get(group, 0)
        
        if skip_metabolite:
            n_excluded_insufficient += 1
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': 'insufficient_valid_data'
            })
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
                result_row[f'{comp_name}_moderated_t'] = np.nan
                result_row[f'{comp_name}_B'] = np.nan
            results_list.append(result_row)
            model_stats.append(None)
            continue
        
        # Extract intensity values
        y = row[sample_cols].values.astype(float)
        
        # Remove missing/zero values for regression (zeros are treated as missing)
        valid_mask = ~(np.isnan(y) | (y == 0) | np.isnan(X).any(axis=1))
        y_clean = y[valid_mask]
        X_clean = X[valid_mask]
        
        min_obs_required = X.shape[1] + 1  # one residual df beyond full rank model
        if len(y_clean) < min_obs_required:
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': 'insufficient_data'
            })
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
                result_row[f'{comp_name}_moderated_t'] = np.nan
                result_row[f'{comp_name}_B'] = np.nan
            results_list.append(result_row)
            model_stats.append(None)
            continue
        
        try:
            # Fit OLS model
            model = sm.OLS(y_clean, X_clean).fit()
            
            # Store statistics for empirical Bayes
            model_stats.append({
                'model': model,
                'metabolite_id': metabolite_id,
                'idx': idx,
                'result_row': result_row,
                'group_valid_counts': group_valid_counts,
                'y_clean': y_clean,
                'X_clean': X_clean,
                'valid_mask': valid_mask
            })
            
            # Store diagnostics
            diagnostics_list.append({
                metabolite_id_col: metabolite_id,
                'r_squared': model.rsquared,
                'adj_r_squared': model.rsquared_adj,
                'f_statistic': model.fvalue,
                'f_pvalue': model.f_pvalue,
                'aic': model.aic,
                'bic': model.bic,
                'log_likelihood': model.llf,
                'n_obs': int(model.nobs),
                'df_model': int(model.df_model),
                'df_resid': int(model.df_resid)
            })
            
        except Exception as e:
            result_row.update({
                'n_samples': sum(group_valid_counts.values()),
                'status': f'error: {str(e)}'
            })
            for g1, g2 in pairwise_combos:
                comp_name = f"{g1}_vs_{g2}"
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_adj_p'] = np.nan
                result_row[f'{g1}_Mean'] = np.nan
                result_row[f'{g2}_Mean'] = np.nan
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
                result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
                result_row[f'{comp_name}_moderated_t'] = np.nan
                result_row[f'{comp_name}_B'] = np.nan
            results_list.append(result_row)
            model_stats.append(None)
    
    # ========== PHASE 2: Empirical Bayes moderation ==========
    valid_models = [m for m in model_stats if m is not None]
    
    if len(valid_models) < 3:
        warnings.warn("Too few valid models for empirical Bayes moderation. Using standard OLS results.")
        d0, s0_squared = 0.0, 1.0
    else:
        # Extract sample variances and degrees of freedom
        s2_array = np.array([m['model'].mse_resid for m in valid_models])
        df_array = np.array([m['model'].df_resid for m in valid_models])
        
        # Fit prior parameters using method of moments
        # Filter valid values
        valid_mask_eb = (s2_array > 0) & np.isfinite(s2_array) & (df_array > 0)
        s2_valid = s2_array[valid_mask_eb]
        df_valid = df_array[valid_mask_eb]
        
        if len(s2_valid) < 3:
            d0, s0_squared = 0.0, np.median(s2_array)
        else:
            # Method of moments estimation for prior parameters
            log_s2 = np.log(s2_valid)
            mean_log_s2 = np.mean(log_s2)
            var_log_s2 = np.var(log_s2, ddof=1)
            
            # Expected variance under the model
            expected_var = np.mean([polygamma(1, d/2) for d in df_valid])
            excess_var = var_log_s2 - expected_var
            
            if excess_var <= 0:
                # No excess variance - strong shrinkage
                d0 = 1e6
            else:
                # Approximate d0 from excess variance
                # polygamma(1, x) ≈ 1/x for small x
                d0 = max(2.0 / excess_var, 0.01)
                d0 = min(d0, 1e6)  # Cap at reasonable value
            
            # Estimate s0² from mean of log(s2)
            correction = np.mean([digamma(d/2) - np.log(d/2) for d in df_valid])
            log_s0_sq = mean_log_s2 - correction
            s0_squared = np.exp(log_s0_sq)
            s0_squared = np.clip(s0_squared, 1e-10, np.max(s2_valid) * 10)
    
    # ========== PHASE 3: Compute moderated statistics ==========
    pairwise_results_dict = {comp: [] for comp in pairwise_combos}
    
    for m_stats in model_stats:
        if m_stats is None:
            continue
        
        model = m_stats['model']
        metabolite_id = m_stats['metabolite_id']
        result_row = m_stats['result_row']
        group_valid_counts = m_stats['group_valid_counts']
        idx = m_stats['idx']
        row = df_intensities.loc[idx]
        
        # Get model statistics
        s2 = model.mse_resid
        df_resid = model.df_resid
        
        # Compute moderated variance (shrinkage estimator)
        s2_mod = (d0 * s0_squared + df_resid * s2) / (d0 + df_resid)
        df_mod = d0 + df_resid
        
        result_row.update({
            'n_samples': sum(group_valid_counts.values()),
            'r_squared': model.rsquared,
            'adj_r_squared': model.rsquared_adj,
            'f_statistic': model.fvalue,
            'f_pvalue': model.f_pvalue,
            'prior_df': d0,
            'prior_var': s0_squared,
            'moderated_var': s2_mod,
            'moderated_df': df_mod,
            'status': 'success'
        })
        
        # Get confidence intervals
        conf_int = model.conf_int(alpha=0.05)
        
        # Calculate pairwise comparisons with moderated statistics
        for g1, g2 in pairwise_combos:
            comp_name = f"{g1}_vs_{g2}"
            
            # Get sample indices for each group
            g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
            g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
            
            # Get valid intensity values
            g1_values = [row[s] for s in g1_samples if pd.notna(row[s]) and row[s] != 0]
            g2_values = [row[s] for s in g2_samples if pd.notna(row[s]) and row[s] != 0]
            
            g1_mean = np.mean(g1_values) if g1_values else np.nan
            g2_mean = np.mean(g2_values) if g2_values else np.nan
            result_row[f'{g1}_Mean'] = g1_mean
            result_row[f'{g2}_Mean'] = g2_mean
            
            # Calculate FC and log2FC
            if not np.isnan(g1_mean) and not np.isnan(g2_mean) and g1_mean > 0:
                fc = g2_mean / g1_mean
                result_row[f'{comp_name}_FC'] = fc
                result_row[f'{comp_name}_log2FC'] = np.log2(fc) if fc > 0 else np.nan
            else:
                result_row[f'{comp_name}_FC'] = np.nan
                result_row[f'{comp_name}_log2FC'] = np.nan
            
            # Compute moderated t-test for pairwise comparison
            g1_col = f'{group_var_name}_{g1}'
            g2_col = f'{group_var_name}_{g2}'
            
            try:
                params = model.params
                if hasattr(params, 'index'):
                    param_names = list(params.index)
                    param_values = np.asarray(params.values, dtype=float)
                else:
                    param_names = list(X_clean.columns)
                    param_values = np.asarray(params, dtype=float)

                name_to_idx = {name: i for i, name in enumerate(param_names)}
                contrast = np.zeros(len(param_values), dtype=float)

                if g1_col in name_to_idx and g2_col in name_to_idx:
                    contrast[name_to_idx[g2_col]] = 1.0
                    contrast[name_to_idx[g1_col]] = -1.0
                elif g2_col in name_to_idx:
                    contrast[name_to_idx[g2_col]] = 1.0
                elif g1_col in name_to_idx:
                    contrast[name_to_idx[g1_col]] = -1.0
                else:
                    contrast = None

                if contrast is not None:
                    coef_diff = float(np.dot(contrast, param_values))
                    normalized_cov = np.asarray(model.normalized_cov_params, dtype=float)
                    # Use normalized covariance directly with moderated variance for stability.
                    var_mod = float(np.dot(contrast, np.dot(normalized_cov, contrast)) * s2_mod)
                    se_mod = np.sqrt(var_mod) if var_mod > 0 else np.nan
                else:
                    coef_diff = np.nan
                    se_mod = np.nan

                ci_lower = np.nan
                ci_upper = np.nan

                if not np.isnan(coef_diff) and not np.isnan(se_mod) and se_mod > 0:
                    t_mod = coef_diff / se_mod
                    
                    # P-value from t-distribution with moderated df
                    pvalue = 2 * scipy_stats.t.sf(np.abs(t_mod), df_mod)

                    # 95% CI from moderated t critical value
                    t_crit_mod = scipy_stats.t.ppf(0.975, df_mod) if df_mod > 0 else np.nan
                    if np.isfinite(t_crit_mod):
                        ci_lower = coef_diff - t_crit_mod * se_mod
                        ci_upper = coef_diff + t_crit_mod * se_mod
                    
                    # B statistic (log-odds of differential expression)
                    # Simplified version: B = log(P(DE)/P(not DE))
                    # Here we use a heuristic based on the moderated t-statistic
                    var_prior = s0_squared if s0_squared > 0 else 1.0
                    B = np.log(1 + t_mod**2 / df_mod) - np.log(var_prior / s2_mod) if s2_mod > 0 else np.nan
                    
                    result_row[f'{comp_name}_moderated_t'] = t_mod
                    result_row[f'{comp_name}_t_statistic'] = t_mod
                    result_row[f'{comp_name}_pvalue'] = pvalue
                    result_row[f'{comp_name}_B'] = B
                    
                    # Store model-based columns (coefficient + SE with moderated variance)
                    result_row[f'{comp_name}_model_effect'] = coef_diff
                    result_row[f'{comp_name}_model_se'] = se_mod
                else:
                    result_row[f'{comp_name}_moderated_t'] = np.nan
                    result_row[f'{comp_name}_t_statistic'] = np.nan
                    result_row[f'{comp_name}_pvalue'] = np.nan
                    result_row[f'{comp_name}_B'] = np.nan
                    
                    # Model-based columns remain NaN if computation failed
                    result_row[f'{comp_name}_model_effect'] = np.nan
                    result_row[f'{comp_name}_model_se'] = np.nan

                # Keep CI schema aligned with OLS output
                result_row[f'{comp_name}_ci_lower_95'] = ci_lower
                result_row[f'{comp_name}_ci_upper_95'] = ci_upper
                    
            except Exception:
                result_row[f'{comp_name}_moderated_t'] = np.nan
                result_row[f'{comp_name}_t_statistic'] = np.nan
                result_row[f'{comp_name}_pvalue'] = np.nan
                result_row[f'{comp_name}_B'] = np.nan
                result_row[f'{comp_name}_ci_lower_95'] = np.nan
                result_row[f'{comp_name}_ci_upper_95'] = np.nan
            
            result_row[f'{comp_name}_adj_p'] = np.nan
            result_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
            
            pairwise_results_dict[(g1, g2)].append(result_row.copy())
        
        results_list.append(result_row)
        
        # Store coefficients
        for coef_name, coef_value in model.params.items():
            ci_lower = conf_int.loc[coef_name, 0] if coef_name in conf_int.index else np.nan
            ci_upper = conf_int.loc[coef_name, 1] if coef_name in conf_int.index else np.nan
            
            all_coefficients.append({
                metabolite_id_col: metabolite_id,
                'coefficient': coef_name,
                'value': coef_value,
                'std_err': model.bse[coef_name] if coef_name in model.bse else np.nan,
                't_statistic': model.tvalues[coef_name] if coef_name in model.tvalues else np.nan,
                'pvalue': model.pvalues[coef_name] if coef_name in model.pvalues else np.nan,
                'ci_lower_95': ci_lower,
                'ci_upper_95': ci_upper
            })
    
    # Convert to DataFrame
    df_results = pd.DataFrame(results_list)
    df_diagnostics = pd.DataFrame(diagnostics_list)
    df_coefficients = pd.DataFrame(all_coefficients)
    
    # Apply FDR correction
    # Always create adj_p columns - if apply_fdr=False, just copy raw p-values
    if len(df_results) > 0:
        for g1, g2 in pairwise_combos:
            comp_name = f"{g1}_vs_{g2}"
            pval_col = f'{comp_name}_pvalue'
            adj_pval_col = f'{comp_name}_adj_p'
            neg_log_col = f'{comp_name}_neg_log10_adj_p'
            
            if pval_col in df_results.columns:
                valid_pvals = df_results[pval_col].dropna()
                if len(valid_pvals) > 0:
                    if not apply_fdr or fdr_method == 'None':
                        # No correction - use raw p-values
                        df_results.loc[valid_pvals.index, adj_pval_col] = valid_pvals.values
                    elif fdr_method == 'BH':
                        n = len(valid_pvals)
                        sorted_idx = valid_pvals.argsort()
                        sorted_pvals = valid_pvals.iloc[sorted_idx]
                        ranks = np.arange(1, n + 1)
                        adjusted = np.minimum.accumulate((sorted_pvals * n / ranks)[::-1])[::-1]
                        adjusted = np.clip(adjusted, 0, 1)
                        adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                        df_results.loc[valid_pvals.index, adj_pval_col] = [adj_dict[i] for i in valid_pvals.index]
                    elif fdr_method == 'Bonferroni':
                        adjusted = np.minimum(valid_pvals.values * len(valid_pvals), 1.0)
                        df_results.loc[valid_pvals.index, adj_pval_col] = adjusted
                    elif fdr_method == 'Holm':
                        n = len(valid_pvals)
                        sorted_idx = valid_pvals.argsort()
                        sorted_pvals = valid_pvals.iloc[sorted_idx]
                        ranks = np.arange(1, n + 1)
                        adjusted = np.maximum.accumulate(np.minimum(sorted_pvals.values * (n - ranks + 1), 1.0))
                        adj_dict = dict(zip(valid_pvals.iloc[sorted_idx].index, adjusted))
                        df_results.loc[valid_pvals.index, adj_pval_col] = [adj_dict[i] for i in valid_pvals.index]
                    
                    # Calculate -log10(adj_p)
                    adj_pvals = df_results[adj_pval_col].replace(0, np.finfo(float).eps)
                    df_results[neg_log_col] = -np.log10(adj_pvals)
    
    # Mark significant metabolites
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        adj_pval_col = f'{comp_name}_adj_p'
        sig_col = f'{comp_name}_significant'
        if adj_pval_col in df_results.columns:
            df_results[sig_col] = df_results[adj_pval_col] < alpha
    
    # Build pairwise DataFrames
    pairwise_dfs = {}
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        if not df_results.empty:
            pairwise_df = df_results.copy()
            pval_col = f'{comp_name}_pvalue'
            if pval_col in pairwise_df.columns:
                pairwise_df = pairwise_df[pairwise_df[pval_col].notna()]
            if 'status' in pairwise_df.columns:
                pairwise_df = pairwise_df[pairwise_df['status'] == 'success']

            if pairwise_df.empty:
                continue

            # Keep only samples and statistics relevant to this pairwise comparison.
            g1_samples = [s for s in sample_cols if group_map.get(s) == g1]
            g2_samples = [s for s in sample_cols if group_map.get(s) == g2]
            comp_samples = g1_samples + g2_samples

            metadata_cols = [c for c in df_results.columns if c not in sample_cols and c != metabolite_id_col]
            keep_cols = [metabolite_id_col] + comp_samples + metadata_cols

            allowed_non_prefixed = {
                f'{g1}_Mean',
                f'{g2}_Mean',
                f'n_{g1}',
                f'n_{g2}',
                'r_squared',
                'adj_r_squared',
                'f_statistic',
                'f_pvalue',
                'status',
                'n_samples',
                'prior_df',
                'prior_var',
                'moderated_var',
                'moderated_df'
            }

            keep_cols = [
                c for c in keep_cols
                if c in pairwise_df.columns and (
                    c == metabolite_id_col or
                    c in comp_samples or
                    comp_name in c or
                    c in allowed_non_prefixed
                )
            ]
            keep_cols = list(dict.fromkeys(keep_cols))

            pairwise_df = pairwise_df[keep_cols].copy()
            rename_map = {
                f'{comp_name}_pvalue': 'pvalue',
                f'{comp_name}_adj_p': 'adj_p',
                f'{comp_name}_FC': 'FC',
                f'{comp_name}_log2FC': 'log2FC',
                f'{comp_name}_neg_log10_adj_p': 'neg_log10_adj_p',
                f'{comp_name}_significant': 'significant',
                f'{comp_name}_moderated_t': 'moderated_t',
                f'{comp_name}_t_statistic': 't_statistic',
                f'{comp_name}_B': 'B'
            }
            pairwise_df = pairwise_df.rename(columns=rename_map)

            if not pairwise_df.empty:
                pairwise_dfs[comp_name] = pairwise_df
    
    # Summary statistics
    sig_counts = {}
    tested_counts = {}
    for g1, g2 in pairwise_combos:
        comp_name = f"{g1}_vs_{g2}"
        sig_col = f'{comp_name}_significant'
        if sig_col in df_results.columns:
            valid_rows = df_results[df_results['status'] == 'success']
            sig_counts[comp_name] = int(valid_rows[sig_col].sum())
            tested_counts[comp_name] = len(valid_rows)
    
    group_sample_counts = {}
    for group in unique_groups:
        group_samples = [s for s in sample_cols if group_map.get(s) == group]
        group_sample_counts[group] = len(group_samples)
    
    n_tested = len([m for m in model_stats if m is not None])
    
    summary = {
        'n_metabolites_total': n_metabolites,
        'n_metabolites_tested': n_tested,
        'n_metabolites_excluded': n_excluded_insufficient,
        'n_covariates_requested': len(covariate_cols) if covariate_cols else 0,
        'n_covariates': len(included_covariates),
        'covariates': ', '.join(included_covariates) if included_covariates else 'None',
        'covariates_requested': ', '.join(covariate_cols) if covariate_cols else 'None',
        'covariates_missing': ', '.join(missing_covariates) if missing_covariates else 'None',
        'n_groups': len(unique_groups),
        'groups': ', '.join(unique_groups),
        'group_sample_counts': group_sample_counts,
        'pairwise_comparisons': pairwise_combos,
        'n_tested_per_comparison': tested_counts,
        'n_significant_per_comparison': sig_counts,
        'alpha': alpha,
        'fdr_method': fdr_method if apply_fdr else 'None',
        'method': 'LIMMA (Empirical Bayes - Pure Python)',
        'prior_df': d0,
        'prior_var': s0_squared,
        'reference_group': effective_reference_group if effective_reference_group is not None else 'None',
        'design_matrix_columns': ', '.join(X.columns.astype(str).tolist())
    }
    
    result = CovariateAnalysisResult(
        metabolite_results=df_results,
        model_diagnostics=df_diagnostics,
        coefficient_table=df_coefficients,
        adjusted_intensities=None,
        summary_stats=summary
    )
    result.pairwise_sheets = pairwise_dfs
    
    return result


def export_covariate_results(
    result: CovariateAnalysisResult,
    output_path: str,
    include_diagnostics: bool = True,
    include_coefficients: bool = True,
    class_results: Optional[CovariateAnalysisResult] = None
) -> str:
    """
    Export covariate analysis results to Excel file(s) with multiple sheets.
    
    If class_results are provided, exports to TWO files:
    - Main file (output_path): Metabolite/lipid results
    - Class file (output_path with _class suffix): Lipid class results
    
    Otherwise exports everything to a single file.
    
    Parameters
    ----------
    result : CovariateAnalysisResult
        Analysis results to export
    output_path : str
        Path for output Excel file
    include_diagnostics : bool
        Whether to include model diagnostics sheet
    include_coefficients : bool
        Whether to include coefficients sheet
    class_results : Optional[CovariateAnalysisResult]
        Optional lipid class analysis results to export (in separate file)
        
    Returns
    -------
    str
        Path to saved file (or tuple if two files created)
    """
    # Construct class file path if needed
    class_file_path = None
    if class_results is not None:
        # Insert "_class" before the file extension
        # E.g., "stat.xlsx" -> "stat_class.xlsx"
        base, ext = os.path.splitext(output_path)
        class_file_path = f"{base}_class{ext}"
    
    try:
        # ========================================
        # MAIN FILE: Metabolite/Lipid Results
        # ========================================
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            sheets_written = []
            
            # Main results - rename to 'Complete Results' for consistency
            if result.metabolite_results is not None and not result.metabolite_results.empty:
                result.metabolite_results.to_excel(
                    writer, sheet_name='Complete Results', index=False
                )
                sheets_written.append('Complete Results')
            
            # Individual pairwise comparison sheets
            if hasattr(result, 'pairwise_sheets') and result.pairwise_sheets:
                for comp_name, comp_df in result.pairwise_sheets.items():
                    if comp_df is not None and not comp_df.empty:
                        # Truncate sheet name if too long (Excel limit is 31 chars)
                        sheet_name = comp_name[:31] if len(comp_name) > 31 else comp_name
                        comp_df.to_excel(writer, sheet_name=sheet_name, index=False)
                        sheets_written.append(sheet_name)
            
            # Diagnostics
            if include_diagnostics and result.model_diagnostics is not None and not result.model_diagnostics.empty:
                result.model_diagnostics.to_excel(
                    writer, sheet_name='Model_Diagnostics', index=False
                )
                sheets_written.append('Model_Diagnostics')
            
            # Coefficients
            if include_coefficients and result.coefficient_table is not None and not result.coefficient_table.empty:
                result.coefficient_table.to_excel(
                    writer, sheet_name='Model_Coefficients', index=False
                )
                sheets_written.append('Model_Coefficients')
            
            # Adjusted intensities
            if result.adjusted_intensities is not None and not result.adjusted_intensities.empty:
                result.adjusted_intensities.to_excel(
                    writer, sheet_name='Adjusted_Intensities', index=True
                )
                sheets_written.append('Adjusted_Intensities')
            
            # Summary - Create a well-formatted summary table
            if result.summary_stats is not None and len(result.summary_stats) > 0:
                try:
                    summary = result.summary_stats
                    summary_rows = []
                    
                    # Section 1: Analysis Parameters
                    summary_rows.append(['ANALYSIS PARAMETERS', ''])
                    summary_rows.append(['Total metabolites in dataset', summary.get('n_metabolites_total', 'N/A')])
                    summary_rows.append(['Metabolites tested', summary.get('n_metabolites_tested', 'N/A')])
                    summary_rows.append(['Metabolites excluded (insufficient data)', summary.get('n_metabolites_excluded', 'N/A')])
                    summary_rows.append(['Covariates adjusted for', summary.get('covariates', 'N/A')])
                    summary_rows.append(['Number of covariates', summary.get('n_covariates', 'N/A')])
                    summary_rows.append(['Significance threshold (alpha)', summary.get('alpha', 0.05)])
                    summary_rows.append(['FDR correction method', summary.get('fdr_method', 'N/A')])
                    summary_rows.append(['', ''])
                    
                    # Section 2: Group Information
                    summary_rows.append(['GROUP INFORMATION', ''])
                    summary_rows.append(['Number of groups', summary.get('n_groups', 'N/A')])
                    summary_rows.append(['Group names', summary.get('groups', 'N/A')])
                    
                    # Add sample counts per group
                    group_counts = summary.get('group_sample_counts', {})
                    for group, count in group_counts.items():
                        summary_rows.append([f'  {group} sample count', count])
                    summary_rows.append(['', ''])
                    
                    # Section 3: Pairwise Comparison Results
                    summary_rows.append(['PAIRWISE COMPARISON RESULTS', '', '', ''])
                    summary_rows.append(['Comparison', 'Metabolites Tested', 'Significant', '% Significant'])
                    
                    pairwise_comps = summary.get('pairwise_comparisons', [])
                    tested_per_comp = summary.get('n_tested_per_comparison', {})
                    sig_per_comp = summary.get('n_significant_per_comparison', {})
                    
                    for g1, g2 in pairwise_comps:
                        comp_name = f"{g1}_vs_{g2}"
                        tested = tested_per_comp.get(comp_name, 0)
                        significant = sig_per_comp.get(comp_name, 0)
                        pct_sig = (significant / tested * 100) if tested > 0 else 0
                        summary_rows.append([comp_name, tested, significant, f'{pct_sig:.1f}%'])
                    
                    # Create DataFrame with 4 columns to match data
                    df_summary = pd.DataFrame(summary_rows, columns=['Parameter', 'Value', 'Column3', 'Column4'])
                    df_summary.to_excel(writer, sheet_name='Summary', index=False)
                    sheets_written.append('Summary')
                    
                except Exception as e:
                    # If summary export fails, continue without it
                    warnings.warn(f"Could not export summary statistics: {e}")
        
        # ========================================
        # CLASS FILE: Lipid Class Results
        # ========================================
        if class_results is not None and class_file_path is not None:
            with pd.ExcelWriter(class_file_path, engine='openpyxl') as writer:
                sheets_written_class = []
                
                # Class Complete Results
                if class_results.metabolite_results is not None and not class_results.metabolite_results.empty:
                    class_results.metabolite_results.to_excel(
                        writer, sheet_name='Complete Results', index=False
                    )
                    sheets_written_class.append('Complete Results')
                
                # Class pairwise comparison sheets
                if hasattr(class_results, 'pairwise_sheets') and class_results.pairwise_sheets:
                    for comp_name, comp_df in class_results.pairwise_sheets.items():
                        if comp_df is not None and not comp_df.empty:
                            # Keep original comp_name without 'Class_' prefix (it's already in class file)
                            sheet_name = comp_name[:31] if len(comp_name) > 31 else comp_name
                            comp_df.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written_class.append(sheet_name)
                
                # Class Summary
                if class_results.summary_stats is not None and len(class_results.summary_stats) > 0:
                    try:
                        summary = class_results.summary_stats
                        summary_rows = []
                        
                        # Section 1: Analysis Parameters
                        summary_rows.append(['CLASS ANALYSIS PARAMETERS', ''])
                        summary_rows.append(['Total lipid classes', summary.get('n_metabolites_total', 'N/A')])
                        summary_rows.append(['Classes tested', summary.get('n_metabolites_tested', 'N/A')])
                        summary_rows.append(['Classes excluded', summary.get('n_metabolites_excluded', 'N/A')])
                        summary_rows.append(['Significance threshold (alpha)', summary.get('alpha', 0.05)])
                        summary_rows.append(['FDR correction method', summary.get('fdr_method', 'N/A')])
                        summary_rows.append(['', ''])
                        
                        # Section 2: Group Information
                        summary_rows.append(['GROUP INFORMATION', ''])
                        summary_rows.append(['Number of groups', summary.get('n_groups', 'N/A')])
                        summary_rows.append(['Group names', summary.get('groups', 'N/A')])
                        summary_rows.append(['', ''])
                        
                        # Section 3: Pairwise Comparison Results
                        summary_rows.append(['CLASS COMPARISON RESULTS', '', '', ''])
                        summary_rows.append(['Comparison', 'Classes Tested', 'Significant', '% Significant'])
                        
                        pairwise_comps = summary.get('pairwise_comparisons', [])
                        tested_per_comp = summary.get('n_tested_per_comparison', {})
                        sig_per_comp = summary.get('n_significant_per_comparison', {})
                        
                        for g1, g2 in pairwise_comps:
                            comp_name = f"{g1}_vs_{g2}"
                            tested = tested_per_comp.get(comp_name, 0)
                            significant = sig_per_comp.get(comp_name, 0)
                            pct_sig = (significant / tested * 100) if tested > 0 else 0
                            summary_rows.append([comp_name, tested, significant, f'{pct_sig:.1f}%'])
                        
                        # Create DataFrame with 4 columns to match data
                        df_summary = pd.DataFrame(summary_rows, columns=['Parameter', 'Value', 'Column3', 'Column4'])
                        df_summary.to_excel(writer, sheet_name='Summary', index=False)
                        sheets_written_class.append('Summary')
                        
                    except Exception as e:
                        # If summary export fails, continue without it
                        warnings.warn(f"Could not export class summary statistics: {e}")

    except ImportError as e:
        raise ImportError(
            "openpyxl is required for Excel export. "
            "Install it with: pip install openpyxl"
        ) from e
    except Exception as e:
        raise RuntimeError(f"Failed to export covariate results: {e}") from e
    
    # Return both files if class file was created, otherwise just main file
    if class_file_path is not None and os.path.exists(class_file_path):
        return {
            'metabolite_file': output_path,
            'class_file': class_file_path
        }
    
    return output_path


# Convenience function for common use case
def covariate_analysis_from_files(
    metabolite_file: str,
    covariate_file: str,
    sample_cols: List[str],
    group_map: Dict[str, str],
    covariate_cols: List[str],
    output_file: str,
    **kwargs
) -> CovariateAnalysisResult:
    """
    Convenience function to run covariate analysis from file paths.
    
    Parameters
    ----------
    metabolite_file : str
        Path to metabolite intensity file (Excel/CSV)
    covariate_file : str
        Path to covariate file
    sample_cols : List[str]
        Sample column names
    group_map : Dict[str, str]
        Sample to group mapping
    covariate_cols : List[str]
        Covariate columns to use
    output_file : str
        Output Excel file path
    **kwargs
        Additional arguments for run_covariate_adjusted_analysis
        
    Returns
    -------
    CovariateAnalysisResult
        Analysis results
    """
    # Load metabolite data
    if metabolite_file.endswith(('.xlsx', '.xls')):
        df_metab = pd.read_excel(metabolite_file)
    else:
        df_metab = pd.read_csv(metabolite_file)
    
    # Load covariate data
    df_cov = load_covariate_file(covariate_file)
    
    # Run analysis
    result = run_covariate_adjusted_analysis(
        df_metab, sample_cols, group_map, df_cov, covariate_cols, **kwargs
    )
    
    # Export results
    export_covariate_results(result, output_file)
    
    return result
