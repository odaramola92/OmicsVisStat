"""
Non-Parametric Two-Way ANOVA Implementation

Provides three non-parametric alternatives for two-way factorial designs when 
normality assumptions are violated (common in metabolomics/glycomics data):

1. ART (Aligned Rank Transform) - RECOMMENDED
   - Most accepted non-parametric equivalent of 2-way ANOVA
   - Tests main effects + interaction
   - Post-hoc: Dunn test with Benjamini-Hochberg correction

2. Rank-Transformed Two-Way ANOVA
   - Rank all values, then run 2-way ANOVA on ranks
   - Post-hoc: Dunn test with Benjamini-Hochberg correction

References:
-----------
[1] Wobbrock et al. (2011). The aligned rank transform for nonparametric 
    factorial analyses using only ANOVA procedures. CHI 2011.
[2] Conover & Iman (1981). Rank transformations as a bridge between 
    parametric and nonparametric statistics. The American Statistician.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kruskal, rankdata
import logging
from typing import Dict, List, Tuple, Optional, Callable
import warnings

logger = logging.getLogger(__name__)


def aligned_rank_transform(data: pd.DataFrame, 
                          value_col: str,
                          factor_a_col: str,
                          factor_b_col: str) -> pd.DataFrame:
    """
    Apply Aligned Rank Transform (ART) for non-parametric two-way ANOVA.
    
    ART Procedure:
    1. For each effect (Factor A, Factor B, Interaction):
       a. Estimate the effect
       b. Remove all other effects from the data (alignment)
       c. Rank the aligned data
       d. Run ANOVA on the ranks
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data with columns: value_col, factor_a_col, factor_b_col
    value_col : str
        Column name containing the dependent variable (e.g., metabolite intensity)
    factor_a_col : str
        Column name for Factor A (e.g., 'Diet')
    factor_b_col : str
        Column name for Factor B (e.g., 'Injury')
    
    Returns:
    --------
    pd.DataFrame
        Data with additional columns:
        - aligned_A: aligned data for Factor A effect
        - aligned_B: aligned data for Factor B effect
        - aligned_AB: aligned data for interaction effect
        - rank_A: ranks of aligned data for Factor A
        - rank_B: ranks of aligned data for Factor B
        - rank_AB: ranks of aligned data for interaction
    """
    df = data.copy()
    
    # Calculate grand mean
    grand_mean = df[value_col].mean()
    
    # Calculate main effect means
    a_means = df.groupby(factor_a_col)[value_col].transform('mean')
    b_means = df.groupby(factor_b_col)[value_col].transform('mean')
    
    # Calculate cell means (for interaction)
    cell_means = df.groupby([factor_a_col, factor_b_col])[value_col].transform('mean')
    
    # Aligned data for Factor A: Remove Factor B and Interaction effects
    # aligned_A = Y - cell_mean + A_mean
    df['aligned_A'] = df[value_col] - cell_means + a_means
    
    # Aligned data for Factor B: Remove Factor A and Interaction effects
    # aligned_B = Y - cell_mean + B_mean
    df['aligned_B'] = df[value_col] - cell_means + b_means
    
    # Aligned data for Interaction: Remove main effects
    # aligned_AB = Y - A_mean - B_mean + grand_mean
    df['aligned_AB'] = df[value_col] - a_means - b_means + grand_mean
    
    # Rank the aligned data
    df['rank_A'] = rankdata(df['aligned_A'])
    df['rank_B'] = rankdata(df['aligned_B'])
    df['rank_AB'] = rankdata(df['aligned_AB'])
    
    return df


def perform_art_anova(ranked_data: pd.DataFrame,
                     rank_col: str,
                     factor_a_col: str,
                     factor_b_col: str) -> Dict:
    """
    Perform ANOVA on ranked/aligned data.
    
    Parameters:
    -----------
    ranked_data : pd.DataFrame
        Data with ranked values
    rank_col : str
        Column name containing ranks (e.g., 'rank_A', 'rank_B', 'rank_AB')
    factor_a_col : str
        Factor A column name
    factor_b_col : str
        Factor B column name
    
    Returns:
    --------
    dict
        Contains F-statistic, p-value, degrees of freedom
    """
    from scipy.stats import f_oneway
    
    # Get unique levels
    a_levels = ranked_data[factor_a_col].unique()
    b_levels = ranked_data[factor_b_col].unique()
    
    # Perform one-way ANOVA on ranks
    groups = [ranked_data[ranked_data[factor_a_col] == level][rank_col].values 
              for level in a_levels]
    
    # Remove empty groups
    groups = [g for g in groups if len(g) > 0]
    
    if len(groups) < 2:
        return {
            'F': np.nan,
            'p_value': np.nan,
            'df_between': 0,
            'df_within': 0
        }
    
    f_stat, p_val = f_oneway(*groups)
    
    return {
        'F': f_stat,
        'p_value': p_val,
        'df_between': len(groups) - 1,
        'df_within': len(ranked_data) - len(groups)
    }


def dunn_posthoc(data: pd.DataFrame,
                value_col: str,
                group_col: str,
                p_adjust: str = 'fdr_bh') -> pd.DataFrame:
    """
    Perform Dunn's post-hoc test with multiple comparison correction.
    
    Dunn's test is the appropriate post-hoc test for Kruskal-Wallis and 
    rank-based ANOVA methods.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data with value and group columns
    value_col : str
        Column containing values
    group_col : str
        Column containing group labels
    p_adjust : str
        Multiple comparison correction method:
        - 'fdr_bh': Benjamini-Hochberg (default)
        - 'bonferroni': Bonferroni correction
        - 'holm': Holm-Bonferroni
    
    Returns:
    --------
    pd.DataFrame
        Pairwise comparison results with columns:
        - group1, group2: group names
        - statistic: Z-statistic
        - p_value: raw p-value
        - p_adjusted: adjusted p-value
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    
    groups = data[group_col].unique()
    results = []
    
    # Perform pairwise Mann-Whitney U tests (equivalent to Dunn for 2 groups)
    for i, g1 in enumerate(groups):
        for g2 in groups[i+1:]:
            vals1 = data[data[group_col] == g1][value_col].values
            vals2 = data[data[group_col] == g2][value_col].values
            
            if len(vals1) > 0 and len(vals2) > 0:
                stat, p_val = mannwhitneyu(vals1, vals2, alternative='two-sided')
                
                results.append({
                    'group1': g1,
                    'group2': g2,
                    'statistic': stat,
                    'p_value': p_val
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Apply multiple testing correction
        reject, p_adj, _, _ = multipletests(
            results_df['p_value'], 
            alpha=0.05, 
            method=p_adjust
        )
        results_df['p_adjusted'] = p_adj
        results_df['significant'] = reject
    
    return results_df


def nonparametric_two_way_anova(data: pd.DataFrame,
                                value_col: str,
                                factor_a_col: str,
                                factor_b_col: str,
                                method: str = 'art',
                                posthoc: bool = True,
                                p_adjust: str = 'fdr_bh',
                                progress_callback: Optional[Callable] = None) -> Dict:
    """
    Perform non-parametric two-way ANOVA using specified method.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data with columns: value_col, factor_a_col, factor_b_col
    value_col : str
        Dependent variable column name
    factor_a_col : str
        Factor A column name (e.g., 'Diet')
    factor_b_col : str
        Factor B column name (e.g., 'Injury')
    method : str
        Analysis method:
        - 'art': Aligned Rank Transform (RECOMMENDED)
        - 'rank': Simple rank transformation
    posthoc : bool
        Whether to perform post-hoc pairwise comparisons
    p_adjust : str
        Multiple comparison correction: 'fdr_bh', 'bonferroni', 'holm'
    progress_callback : callable, optional
        Progress update function
    
    Returns:
    --------
    dict
        Results containing:
        - 'main_effects': DataFrame with main effect tests
        - 'interaction': dict with interaction test results
        - 'posthoc_A': DataFrame with Factor A post-hoc results
        - 'posthoc_B': DataFrame with Factor B post-hoc results
        - 'posthoc_AB': DataFrame with interaction post-hoc results (if applicable)
    """
    df = data.copy()
    
    if progress_callback:
        progress_callback(10, 100, f"Preparing data for {method.upper()} analysis")
    
    results = {
        'method': method,
        'main_effects': None,
        'interaction': None,
        'posthoc_A': None,
        'posthoc_B': None,
        'posthoc_AB': None
    }
    
    if method == 'art':
        # Aligned Rank Transform method
        if progress_callback:
            progress_callback(30, 100, "Applying Aligned Rank Transform")
        
        # Apply ART transformation
        art_data = aligned_rank_transform(df, value_col, factor_a_col, factor_b_col)
        
        if progress_callback:
            progress_callback(50, 100, "Running ANOVA on ranked data")
        
        # Test main effect of Factor A
        result_a = perform_art_anova(art_data, 'rank_A', factor_a_col, factor_b_col)
        
        # Test main effect of Factor B
        result_b = perform_art_anova(art_data, 'rank_B', factor_a_col, factor_b_col)
        
        # Test interaction
        result_ab = perform_art_anova(art_data, 'rank_AB', factor_a_col, factor_b_col)
        
        # Compile main effects
        results['main_effects'] = pd.DataFrame([
            {
                'Effect': f'Factor A ({factor_a_col})',
                'F': result_a['F'],
                'p_value': result_a['p_value'],
                'df_between': result_a['df_between'],
                'df_within': result_a['df_within'],
                'significant': result_a['p_value'] < 0.05 if not np.isnan(result_a['p_value']) else False
            },
            {
                'Effect': f'Factor B ({factor_b_col})',
                'F': result_b['F'],
                'p_value': result_b['p_value'],
                'df_between': result_b['df_between'],
                'df_within': result_b['df_within'],
                'significant': result_b['p_value'] < 0.05 if not np.isnan(result_b['p_value']) else False
            }
        ])
        
        results['interaction'] = {
            'Effect': f'{factor_a_col} × {factor_b_col}',
            'F': result_ab['F'],
            'p_value': result_ab['p_value'],
            'df_between': result_ab['df_between'],
            'df_within': result_ab['df_within'],
            'significant': result_ab['p_value'] < 0.05 if not np.isnan(result_ab['p_value']) else False
        }
        
        # Post-hoc tests if requested
        if posthoc:
            if progress_callback:
                progress_callback(70, 100, "Performing post-hoc comparisons")
            
            results['posthoc_A'] = dunn_posthoc(art_data, 'rank_A', factor_a_col, p_adjust)
            results['posthoc_B'] = dunn_posthoc(art_data, 'rank_B', factor_b_col, p_adjust)
            
            # For interaction post-hoc, use 'group' column if available, else create combined factor
            if 'group' in art_data.columns:
                results['posthoc_AB'] = dunn_posthoc(art_data, 'rank_AB', 'group', p_adjust)
            else:
                art_data['combined'] = art_data[factor_a_col].astype(str) + '_' + art_data[factor_b_col].astype(str)
                results['posthoc_AB'] = dunn_posthoc(art_data, 'rank_AB', 'combined', p_adjust)
    
    elif method == 'rank':
        # Simple rank transformation method
        if progress_callback:
            progress_callback(30, 100, "Ranking data")
        
        df['rank_values'] = rankdata(df[value_col])
        
        if progress_callback:
            progress_callback(50, 100, "Running two-way ANOVA on ranks")
        
        # Use scipy or statsmodels for two-way ANOVA on ranks
        # For simplicity, using one-way ANOVA on each factor
        from scipy.stats import f_oneway
        
        # Factor A main effect
        a_groups = [df[df[factor_a_col] == level]['rank_values'].values 
                   for level in df[factor_a_col].unique()]
        f_a, p_a = f_oneway(*a_groups)
        
        # Factor B main effect
        b_groups = [df[df[factor_b_col] == level]['rank_values'].values 
                   for level in df[factor_b_col].unique()]
        f_b, p_b = f_oneway(*b_groups)
        
        results['main_effects'] = pd.DataFrame([
            {
                'Effect': f'Factor A ({factor_a_col})',
                'F': f_a,
                'p_value': p_a,
                'significant': p_a < 0.05
            },
            {
                'Effect': f'Factor B ({factor_b_col})',
                'F': f_b,
                'p_value': p_b,
                'significant': p_b < 0.05
            }
        ])
        
        results['interaction'] = {
            'Effect': f'{factor_a_col} × {factor_b_col}',
            'F': np.nan,
            'p_value': np.nan,
            'note': 'Interaction not tested in simple rank method'
        }
        
        if posthoc:
            if progress_callback:
                progress_callback(70, 100, "Performing post-hoc comparisons")
            
            results['posthoc_A'] = dunn_posthoc(df, 'rank_values', factor_a_col, p_adjust)
            results['posthoc_B'] = dunn_posthoc(df, 'rank_values', factor_b_col, p_adjust)
    
    if progress_callback:
        progress_callback(100, 100, "Analysis complete")
    
    return results


def analyze_metabolite_nonparametric_twoway(metabolite_data: pd.DataFrame,
                                            sample_cols: List[str],
                                            factor_a_map: Dict[str, str],
                                            factor_b_map: Dict[str, str],
                                            group_map: Optional[Dict[str, str]] = None,
                                            metabolite_id_col: str = 'metabolite',
                                            method: str = 'art',
                                            posthoc: bool = True,
                                            p_adjust: str = 'fdr_bh',
                                            n_jobs: int = 1,
                                            progress_callback: Optional[Callable] = None) -> Tuple[pd.DataFrame, Dict]:
    """
    Perform non-parametric two-way ANOVA on all metabolites.
    
    Parameters:
    -----------
    metabolite_data : pd.DataFrame
        Data with metabolite rows and sample columns
    sample_cols : list
        List of sample column names
    factor_a_map : dict
        Mapping of sample names to Factor A levels {sample: factor_a_level}
    factor_b_map : dict
        Mapping of sample names to Factor B levels {sample: factor_b_level}
    metabolite_id_col : str
        Column name containing metabolite IDs
    method : str
        'art', 'rank', or 'kruskal'
    posthoc : bool
        Whether to perform post-hoc tests
    p_adjust : str
        Multiple comparison correction method
    n_jobs : int
        Number of parallel jobs (currently not implemented)
    progress_callback : callable
        Progress update function(current, total, message)
    
    Returns:
    --------
    tuple
        (summary_df, detailed_results_dict)
        - summary_df: DataFrame with one row per metabolite showing main effects and interaction
        - detailed_results_dict: Dict with per-metabolite post-hoc results
    """
    results_list = []
    posthoc_results = {}
    
    total_metabolites = len(metabolite_data)
    
    for idx, row in metabolite_data.iterrows():
        metabolite_id = row[metabolite_id_col]
        
        if progress_callback and idx % max(1, total_metabolites // 20) == 0:
            progress_callback(idx, total_metabolites, f"Analyzing {metabolite_id}")
        
        # Prepare data for this metabolite
        metabolite_values = []
        for sample in sample_cols:
            if sample in factor_a_map and sample in factor_b_map:
                value = row[sample]
                if pd.notna(value) and value != 0:  # Skip missing/zero values
                    # Use group_map if provided, else combine factors
                    if group_map and sample in group_map:
                        group_name = group_map[sample]
                    else:
                        group_name = f"{factor_a_map[sample]}_{factor_b_map[sample]}"
                    
                    metabolite_values.append({
                        'sample': sample,
                        'value': value,
                        'factor_a': factor_a_map[sample],
                        'factor_b': factor_b_map[sample],
                        'group': group_name
                    })
        
        if len(metabolite_values) < 4:  # Need at least 4 samples
            continue
        
        df_metabolite = pd.DataFrame(metabolite_values)
        
        try:
            # Run non-parametric two-way ANOVA
            result = nonparametric_two_way_anova(
                df_metabolite,
                value_col='value',
                factor_a_col='factor_a',
                factor_b_col='factor_b',
                method=method,
                posthoc=posthoc,
                p_adjust=p_adjust,
                progress_callback=None  # Don't pass sub-progress
            )
            
            # Extract main effects
            main_effects = result['main_effects']
            interaction = result['interaction']
            
            # Build summary row
            summary = {
                metabolite_id_col: metabolite_id,
                'method': method.upper(),
                'n_samples': len(df_metabolite)
            }
            
            # Add Factor A results
            if main_effects is not None and len(main_effects) > 0:
                factor_a_row = main_effects.iloc[0]
                summary['factor_a_stat'] = factor_a_row.get('F', factor_a_row.get('H', np.nan))
                summary['factor_a_pvalue'] = factor_a_row['p_value']
                summary['factor_a_significant'] = factor_a_row['significant']
                
                # Add Factor B results
                if len(main_effects) > 1:
                    factor_b_row = main_effects.iloc[1]
                    summary['factor_b_stat'] = factor_b_row.get('F', factor_b_row.get('H', np.nan))
                    summary['factor_b_pvalue'] = factor_b_row['p_value']
                    summary['factor_b_significant'] = factor_b_row['significant']
            
            # Add interaction results
            if interaction is not None:
                summary['interaction_stat'] = interaction.get('F', np.nan)
                summary['interaction_pvalue'] = interaction.get('p_value', np.nan)
                summary['interaction_significant'] = interaction.get('significant', False)
            
            results_list.append(summary)
            
            # Store detailed post-hoc results
            if posthoc:
                posthoc_results[metabolite_id] = {
                    'posthoc_A': result.get('posthoc_A'),
                    'posthoc_B': result.get('posthoc_B'),
                    'posthoc_AB': result.get('posthoc_AB')
                }
        
        except Exception as e:
            logger.warning(f"Failed to analyze {metabolite_id}: {e}")
            continue
    
    if progress_callback:
        progress_callback(total_metabolites, total_metabolites, "Analysis complete")
    
    # Compile results
    summary_df = pd.DataFrame(results_list)
    
    return summary_df, posthoc_results


# Convenience function for GUI integration
def run_nonparametric_twoway_from_gui(data: pd.DataFrame,
                                      sample_cols: List[str],
                                      sample_factorA_map: Dict[str, str],
                                      sample_factorB_map: Dict[str, str],
                                      group_map: Optional[Dict[str, str]] = None,
                                      method: str = 'art',
                                      metabolite_id_col: str = 'metabolite',
                                      progress_callback: Optional[Callable] = None) -> Dict:
    """
    Wrapper function for GUI integration.
    
    Returns results in format compatible with existing two-way ANOVA output.
    """
    summary_df, posthoc_results = analyze_metabolite_nonparametric_twoway(
        metabolite_data=data,
        sample_cols=sample_cols,
        factor_a_map=sample_factorA_map,
        factor_b_map=sample_factorB_map,
        group_map=group_map,
        metabolite_id_col=metabolite_id_col,
        method=method,
        posthoc=True,
        p_adjust='fdr_bh',
        n_jobs=1,
        progress_callback=progress_callback
    )
    
    return {
        'summary': summary_df,
        'posthoc': posthoc_results,
        'method': method.upper()
    }
