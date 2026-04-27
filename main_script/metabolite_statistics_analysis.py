"""Statistics and normalization utilities for the metabolomics GUI.

This module is intentionally GUI-agnostic. It provides functions to:
 1. Detect candidate numeric sample columns.
 2. Perform several normalization methods per polarity table.
 3. Run row-wise statistical tests across user-defined cohorts.

Design assumptions:
 - Each row corresponds to a metabolite / feature.
 - Each numeric sample column represents one replicate measurement (e.g., subject/sample) for that metabolite.
 - A grouping (cohort) map assigns each sample column to exactly one group.
 - Statistical tests are performed per metabolite across groups using the replicate values.

NOTE: For large datasets the per-row Python loop may be computationally expensive; this
implementation emphasizes clarity and correctness first. Potential optimizations
(vectorization, numba, parallelization) can be explored later if profiling indicates
this function is a bottleneck.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple, Optional, Any, Callable
import numpy as np
import pandas as pd
from scipy import stats
import warnings
from concurrent.futures import ThreadPoolExecutor
import time

# Suppress specific warnings from scipy and statsmodels for constant input data
warnings.filterwarnings('ignore', category=stats.ConstantInputWarning)
warnings.filterwarnings('ignore', message='invalid value encountered in divide')
warnings.filterwarnings('ignore', message='Each of the input arrays is constant')
warnings.filterwarnings('ignore', message='Precision loss occurred in moment calculation due to catastrophic cancellation')

try:
    import statsmodels.stats.multicomp as _sm_mc  # optional for Tukey HSD
except Exception:  # pragma: no cover
    _sm_mc = None

# Optional: statsmodels for two-way ANOVA support
try:  # pragma: no cover - imported lazily when two-way ANOVA is used
    import statsmodels.api as _sm
    import statsmodels.formula.api as _smf
except Exception:  # pragma: no cover
    _sm = None
    _smf = None

# Global floor for any p-values to avoid exact zero (underflow) in outputs/log transforms.
# Using machine epsilon to avoid display issues while preserving actual precision.
PVAL_FLOOR: float = np.finfo(float).eps

def _floor_pval(value: Optional[float]) -> Optional[float]:
    """Return p-value with only a machine-tiny floor or NaN for invalid values.

    Final display flooring is applied per comparison column via
    _apply_dynamic_pvalue_floor_series().
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


def _apply_dynamic_pvalue_floor_series(series: pd.Series, scale: float = 0.001) -> pd.Series:
    """Apply dynamic floor to one p-value series.

    floor = scale * smallest positive p-value in the series.
    Default scale=0.01 means two orders of magnitude smaller.
    """
    s = pd.to_numeric(series, errors='coerce')
    finite_mask = np.isfinite(s)
    positive = s[finite_mask & (s > 0)]
    if positive.empty:
        return s

    dynamic_floor = float(positive.min()) * float(scale)
    if not np.isfinite(dynamic_floor) or dynamic_floor <= 0:
        return s
    s = s.mask(finite_mask & (s <= 0), dynamic_floor)
    s = s.mask(np.isfinite(s) & (s < dynamic_floor), dynamic_floor)
    return s


def _apply_dynamic_pvalue_floor_frame(df: pd.DataFrame, columns: Sequence[str], scale: float = 0.001) -> pd.DataFrame:
    """Apply dynamic p-value floor for each specified column independently."""
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = _apply_dynamic_pvalue_floor_series(out[col], scale=scale)
    return out

def _effect_sizes(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float]:
    """Compute Cohen's d and Cliff's delta between two 1-D numeric arrays.
    Returns (cohen_d, cliffs_delta). NaNs and zeros ignored. If insufficient data returns (nan, nan)."""
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    # Filter out both NaNs and zeros
    v1 = v1[~np.isnan(v1) & (v1 != 0)]
    v2 = v2[~np.isnan(v2) & (v2 != 0)]
    if len(v1) == 0 or len(v2) == 0:
        return (np.nan, np.nan)
    # Cohen's d (pooled SD)
    m1, m2 = v1.mean(), v2.mean()
    s1, s2 = v1.std(ddof=1) if len(v1) > 1 else 0.0, v2.std(ddof=1) if len(v2) > 1 else 0.0
    n1, n2 = len(v1), len(v2)
    if n1 + n2 - 2 > 0:
        pooled = np.sqrt(((n1 - 1)*s1**2 + (n2 - 1)*s2**2) / (n1 + n2 - 2))
    else:
        pooled = np.nan
    cohen_d = (m2 - m1) / pooled if pooled not in (0, np.nan) else np.nan
    # Cliff's delta - vectorized implementation (O(n log n) instead of O(n²))
    # For each value in v1, count how many in v2 are greater/lesser
    v1_reshaped = v1[:, np.newaxis]  # Shape: (n1, 1)
    v2_reshaped = v2[np.newaxis, :]  # Shape: (1, n2)
    comparison = np.sign(v2_reshaped - v1_reshaped)  # Shape: (n1, n2)
    cliffs = np.sum(comparison) / (n1 * n2)
    return cohen_d, float(cliffs)

def _safe_fc(m1: Optional[float], m2: Optional[float]) -> tuple[float, float]:
    """Return (fold_change, log2FC) guarding against zero/NaN divisions."""
    try:
        if m1 is None or m2 is None:
            return (np.nan, np.nan)
        m1 = float(m1)
        m2 = float(m2)
        if np.isnan(m1) or np.isnan(m2) or m1 == 0:
            return (np.nan, np.nan)
        fc = m2 / m1
        if fc <= 0:
            return (fc, np.nan)
        return (fc, np.log2(fc))
    except Exception:
        return (np.nan, np.nan)


def _welch_effect_ci(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float, float, float]:
    """Return mean-difference effect, SE, and 95% CI for (group2 - group1)."""
    try:
        v1 = np.array(v1, dtype=float)
        v2 = np.array(v2, dtype=float)
        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]

        n1 = len(v1)
        n2 = len(v2)
        if n1 < 2 or n2 < 2:
            return (np.nan, np.nan, np.nan, np.nan)

        mean_diff = float(np.mean(v2) - np.mean(v1))
        var1 = float(np.var(v1, ddof=1))
        var2 = float(np.var(v2, ddof=1))
        se2 = (var1 / n1) + (var2 / n2)

        if not np.isfinite(se2) or se2 <= 0:
            return (mean_diff, np.nan, np.nan, np.nan)

        se = float(np.sqrt(se2))
        num = se2 ** 2
        den = ((var1 / n1) ** 2 / (n1 - 1)) + ((var2 / n2) ** 2 / (n2 - 1))
        if not np.isfinite(den) or den <= 0:
            return (mean_diff, se, np.nan, np.nan)

        df = num / den
        if not np.isfinite(df) or df <= 0:
            return (mean_diff, se, np.nan, np.nan)

        t_crit = float(stats.t.ppf(0.975, df))
        if not np.isfinite(t_crit):
            return (mean_diff, se, np.nan, np.nan)

        ci_low = mean_diff - t_crit * se
        ci_high = mean_diff + t_crit * se
        return (mean_diff, se, float(ci_low), float(ci_high))
    except Exception:
        return (np.nan, np.nan, np.nan, np.nan)


def _rank_biserial_effect_ci(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float, float, float]:
    """Return rank-biserial effect, asymptotic SE, and 95% CI for (group2 vs group1)."""
    try:
        v1 = np.array(v1, dtype=float)
        v2 = np.array(v2, dtype=float)
        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]

        n1 = len(v1)
        n2 = len(v2)
        if n1 < 1 or n2 < 1:
            return (np.nan, np.nan, np.nan, np.nan)

        mw = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        u1 = float(mw.statistic)
        total_pairs = float(n1 * n2)
        if total_pairs <= 0:
            return (np.nan, np.nan, np.nan, np.nan)

        u2 = total_pairs - u1
        effect = (2.0 * u2 / total_pairs) - 1.0

        var_u = total_pairs * (n1 + n2 + 1.0) / 12.0
        if not np.isfinite(var_u) or var_u <= 0:
            return (float(effect), np.nan, np.nan, np.nan)

        se = 2.0 * np.sqrt(var_u) / total_pairs
        if not np.isfinite(se) or se <= 0:
            return (float(effect), np.nan, np.nan, np.nan)

        z_crit = 1.959963984540054
        ci_low = max(-1.0, effect - z_crit * se)
        ci_high = min(1.0, effect + z_crit * se)
        return (float(effect), float(se), float(ci_low), float(ci_high))
    except Exception:
        return (np.nan, np.nan, np.nan, np.nan)

def _rots_test(v1: np.ndarray, v2: np.ndarray, *, B: int = 1000, K: int = 100, 
               alpha: float = 0.1, seed: Optional[int] = None,
               rng: Optional[np.random.Generator] = None) -> tuple[float, float, float, float]:
    """Perform ROTS (Reproducibility-Optimized Test Statistic) between two groups.
    
    This is a fast single-metabolite implementation that uses a simplified approach
    with fixed s0 estimation. For batch processing of many metabolites, use 
    _rots_batch() which is much more efficient.
    
    Parameters
    ----------
    v1, v2 : array-like
        Sample values for group 1 and group 2
    B : int
        Number of bootstrap/permutation iterations (default 1000)
    K : int
        Not used in single-metabolite version (kept for API compatibility)
    alpha : float
        Not used in single-metabolite version (kept for API compatibility)
    seed : int or None
        Random seed for reproducibility
        
    Returns
    -------
    statistic : float
        ROTS test statistic (modified t-statistic)
    p_value : float
        Permutation-based p-value
    s0 : float
        Variance stabilization parameter (median SE from bootstrap)
    reproducibility : float
        Always 1.0 for single metabolite (placeholder)
        
    References
    ----------
    Elo et al. (2008). "Reproducibility-optimized test statistic for ranking genes in 
    microarray experiments." IEEE/ACM Trans Comput Biol Bioinform, 5(3):423-431.
    """
    # Initialize independent random generator per call to stay thread-safe
    local_rng = rng
    if local_rng is None:
        local_rng = np.random.default_rng(seed)
    elif seed is not None:
        local_rng = np.random.default_rng(seed)
    
    # Clean data - remove NaN and zeros
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    v1 = v1[~np.isnan(v1) & (v1 != 0)]
    v2 = v2[~np.isnan(v2) & (v2 != 0)]
    
    if len(v1) < 2 or len(v2) < 2:
        return (np.nan, np.nan, np.nan, np.nan)
    
    n1, n2 = len(v1), len(v2)
    
    # Compute observed statistic components
    mean_diff = np.mean(v2) - np.mean(v1)
    var1 = np.var(v1, ddof=1)
    var2 = np.var(v2, ddof=1)
    se = np.sqrt(var1/n1 + var2/n2)
    
    # Estimate s0 via bootstrap - vectorized for speed
    # Generate all bootstrap indices at once
    boot_idx1 = local_rng.integers(0, n1, size=(B, n1))
    boot_idx2 = local_rng.integers(0, n2, size=(B, n2))
    
    # Get bootstrap samples - shape (B, n1) and (B, n2)
    boot_v1 = v1[boot_idx1]
    boot_v2 = v2[boot_idx2]
    
    # Calculate bootstrap SEs - vectorized across all B iterations
    boot_var1 = np.var(boot_v1, axis=1, ddof=1)
    boot_var2 = np.var(boot_v2, axis=1, ddof=1)
    boot_ses = np.sqrt(boot_var1/n1 + boot_var2/n2)
    
    # Use median SE as s0 (variance stabilization parameter)
    valid_ses = boot_ses[boot_ses > 0]
    s0 = np.median(valid_ses) if len(valid_ses) > 0 else se
    if np.isnan(s0) or s0 == 0:
        s0 = se if se > 0 else 1.0
    
    # Calculate ROTS statistic
    rots_stat = mean_diff / (se + s0)
    
    # Permutation test for p-value - fully vectorized
    combined = np.concatenate([v1, v2])
    n_total = n1 + n2
    
    # Generate all permutation indices at once
    perm_indices = np.array([local_rng.permutation(n_total) for _ in range(B)])
    
    # Split permuted data
    perm_v1 = combined[perm_indices[:, :n1]]  # Shape: (B, n1)
    perm_v2 = combined[perm_indices[:, n1:]]  # Shape: (B, n2)
    
    # Calculate null statistics - vectorized
    perm_mean_diff = np.mean(perm_v2, axis=1) - np.mean(perm_v1, axis=1)
    perm_var1 = np.var(perm_v1, axis=1, ddof=1)
    perm_var2 = np.var(perm_v2, axis=1, ddof=1)
    perm_se = np.sqrt(perm_var1/n1 + perm_var2/n2)
    
    null_stats = perm_mean_diff / (perm_se + s0)
    
    # Calculate two-tailed p-value
    p_value = np.mean(np.abs(null_stats) >= np.abs(rots_stat))
    p_value = max(p_value, 1.0 / B)  # Minimum p-value is 1/B
    
    return (float(rots_stat), float(p_value), float(s0), 1.0)


def _rots_batch(data_matrix: np.ndarray, group1_indices: np.ndarray, group2_indices: np.ndarray,
                *, B: int = 500, K: int = None, seed: Optional[int] = None,
                progress_callback: Optional[Callable[[int, int, str], None]] = None
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Perform vectorized ROTS analysis on all metabolites at once.
    
    This is the efficient batch implementation that processes all metabolites
    simultaneously using matrix operations, similar to the R ROTS package.
    
    Parameters
    ----------
    data_matrix : np.ndarray
        2D array of shape (n_metabolites, n_samples) with measurements.
        NaN values are handled per-metabolite.
    group1_indices : np.ndarray
        Column indices for group 1 samples
    group2_indices : np.ndarray
        Column indices for group 2 samples
    B : int
        Number of bootstrap iterations (default 500, sufficient for most cases)
    K : int or None
        Number of top features for reproducibility. If None, uses min(100, n_metabolites/4)
    seed : int or None
        Random seed for reproducibility
    progress_callback : callable or None
        Optional callback(current, total, message) for progress updates
        
    Returns
    -------
    statistics : np.ndarray
        ROTS statistics for each metabolite
    p_values : np.ndarray
        Permutation-based p-values for each metabolite
    s0 : float
        Global optimized variance stabilization parameter
    reproducibility : float
        Reproducibility score at optimal parameters
    """
    local_rng = np.random.default_rng(seed)
    
    n_metabolites = data_matrix.shape[0]
    n1 = len(group1_indices)
    n2 = len(group2_indices)
    
    if K is None:
        K = min(100, max(10, n_metabolites // 4))
    
    # Extract group data - shape (n_metabolites, n_samples_per_group)
    g1_data = data_matrix[:, group1_indices]
    g2_data = data_matrix[:, group2_indices]
    
    # Replace zeros with NaN for proper handling
    g1_data = np.where(g1_data == 0, np.nan, g1_data)
    g2_data = np.where(g2_data == 0, np.nan, g2_data)
    
    # Calculate observed statistics for all metabolites at once
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        # Nanmean handles missing values
        mean1 = np.nanmean(g1_data, axis=1)
        mean2 = np.nanmean(g2_data, axis=1)
        mean_diff = mean2 - mean1
        
        # Count valid (non-NaN) samples per metabolite per group
        n1_valid = np.sum(~np.isnan(g1_data), axis=1)
        n2_valid = np.sum(~np.isnan(g2_data), axis=1)
        
        # Variance with ddof=1, handling NaN
        var1 = np.nanvar(g1_data, axis=1, ddof=1)
        var2 = np.nanvar(g2_data, axis=1, ddof=1)
        
        # Standard error - use valid counts
        se = np.sqrt(var1/np.maximum(n1_valid, 1) + var2/np.maximum(n2_valid, 1))
    
    # Mark invalid metabolites (not enough samples)
    valid_mask = (n1_valid >= 2) & (n2_valid >= 2)
    
    if progress_callback:
        progress_callback(0, 100, "Estimating s0 parameter...")
    
    # Estimate s0 via bootstrap across all metabolites
    # Use a subset of metabolites for s0 estimation if dataset is large
    n_for_s0 = min(n_metabolites, 1000)
    s0_indices = local_rng.choice(n_metabolites, size=n_for_s0, replace=False) if n_metabolites > 1000 else np.arange(n_metabolites)
    
    boot_ses_all = []
    n_boot_for_s0 = min(B, 200)  # Fewer bootstrap samples needed for s0 estimation
    
    for b in range(n_boot_for_s0):
        # Bootstrap indices for each group
        boot_idx1 = local_rng.integers(0, n1, size=n1)
        boot_idx2 = local_rng.integers(0, n2, size=n2)
        
        # Get bootstrap data for s0 estimation subset
        boot_g1 = g1_data[s0_indices][:, boot_idx1]
        boot_g2 = g2_data[s0_indices][:, boot_idx2]
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            boot_var1 = np.nanvar(boot_g1, axis=1, ddof=1)
            boot_var2 = np.nanvar(boot_g2, axis=1, ddof=1)
            boot_n1 = np.sum(~np.isnan(boot_g1), axis=1)
            boot_n2 = np.sum(~np.isnan(boot_g2), axis=1)
            boot_se = np.sqrt(boot_var1/np.maximum(boot_n1, 1) + boot_var2/np.maximum(boot_n2, 1))
        
        boot_ses_all.append(boot_se)
    
    # Stack and compute median SE across bootstraps for each metabolite
    boot_ses_matrix = np.vstack(boot_ses_all)  # Shape: (n_boot, n_for_s0)
    
    # Use percentile-based s0 estimation (similar to R ROTS)
    # s0 is chosen to stabilize variance across different expression levels
    median_ses = np.nanmedian(boot_ses_matrix, axis=0)
    s0_candidates = np.nanpercentile(median_ses[~np.isnan(median_ses)], [25, 50, 75])
    
    if progress_callback:
        progress_callback(20, 100, "Optimizing reproducibility...")
    
    # Simple s0 optimization: use median of standard errors
    s0 = s0_candidates[1] if not np.isnan(s0_candidates[1]) else np.nanmedian(se[valid_mask])
    if np.isnan(s0) or s0 == 0:
        s0 = 0.01  # Fallback
    
    # Calculate ROTS statistics with optimized s0
    rots_stats = mean_diff / (se + s0)
    rots_stats[~valid_mask] = np.nan
    
    if progress_callback:
        progress_callback(40, 100, "Computing permutation p-values...")
    
    # Permutation test for p-values - vectorized across all metabolites
    combined_data = np.hstack([g1_data, g2_data])  # Shape: (n_metabolites, n1+n2)
    n_total = n1 + n2
    
    # Count how many permutation statistics exceed observed
    exceed_count = np.zeros(n_metabolites)
    
    # Process permutations in batches to manage memory
    batch_size = min(100, B)
    n_batches = (B + batch_size - 1) // batch_size
    
    for batch_idx in range(n_batches):
        batch_start = batch_idx * batch_size
        batch_end = min((batch_idx + 1) * batch_size, B)
        current_batch_size = batch_end - batch_start
        
        # Generate permutation indices for this batch
        perm_null_stats = np.zeros((current_batch_size, n_metabolites))
        
        for p in range(current_batch_size):
            perm_idx = local_rng.permutation(n_total)
            perm_g1 = combined_data[:, perm_idx[:n1]]
            perm_g2 = combined_data[:, perm_idx[n1:]]
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                perm_mean_diff = np.nanmean(perm_g2, axis=1) - np.nanmean(perm_g1, axis=1)
                perm_n1 = np.sum(~np.isnan(perm_g1), axis=1)
                perm_n2 = np.sum(~np.isnan(perm_g2), axis=1)
                perm_var1 = np.nanvar(perm_g1, axis=1, ddof=1)
                perm_var2 = np.nanvar(perm_g2, axis=1, ddof=1)
                perm_se = np.sqrt(perm_var1/np.maximum(perm_n1, 1) + perm_var2/np.maximum(perm_n2, 1))
            
            perm_null_stats[p] = perm_mean_diff / (perm_se + s0)
        
        # Count exceedances
        exceed_count += np.sum(np.abs(perm_null_stats) >= np.abs(rots_stats), axis=0)
        
        if progress_callback:
            progress = 40 + int(60 * (batch_idx + 1) / n_batches)
            progress_callback(progress, 100, f"Permutation batch {batch_idx + 1}/{n_batches}...")
    
    # Calculate p-values
    p_values = exceed_count / B
    p_values = np.maximum(p_values, 1.0 / B)  # Minimum p-value
    p_values[~valid_mask] = np.nan
    
    # Estimate reproducibility (simplified)
    reproducibility = 1.0  # Placeholder - could implement proper reproducibility scoring
    
    if progress_callback:
        progress_callback(100, 100, "ROTS analysis complete")
    
    return rots_stats, p_values, float(s0), reproducibility


def _limma_batch(data_matrix: np.ndarray, group1_indices: np.ndarray, group2_indices: np.ndarray,
                 *, seed: Optional[int] = None,
                 progress_callback: Optional[Callable[[int, int, str], None]] = None
                 ) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Perform vectorized limma-style analysis on all metabolites at once.
    
    This is a pure Python implementation of limma's empirical Bayes moderated t-test
    for simple two-group comparisons without covariates. It borrows variance information
    across all metabolites to produce more stable variance estimates, which is particularly
    useful for small sample sizes.
    
    Parameters
    ----------
    data_matrix : np.ndarray
        2D array of shape (n_metabolites, n_samples) with measurements.
        NaN values are handled per-metabolite.
    group1_indices : np.ndarray
        Column indices for group 1 samples
    group2_indices : np.ndarray
        Column indices for group 2 samples
    seed : int or None
        Random seed (not used, kept for API compatibility)
    progress_callback : callable or None
        Optional callback(current, total, message) for progress updates
        
    Returns
    -------
    statistics : np.ndarray
        Moderated t-statistics for each metabolite
    p_values : np.ndarray
        P-values from t-distribution with moderated degrees of freedom
    d0 : float
        Prior degrees of freedom (estimated from data)
    s0_squared : float
        Prior variance (estimated from data)
        
    References
    ----------
    Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing 
    differential expression in microarray experiments. Statistical Applications 
    in Genetics and Molecular Biology, 3(1), Article 3.
    """
    from scipy import stats as scipy_stats
    from scipy.special import digamma, polygamma
    
    n_metabolites = data_matrix.shape[0]
    n1 = len(group1_indices)
    n2 = len(group2_indices)
    
    if progress_callback:
        progress_callback(0, 100, "Computing group statistics...")
    
    # Extract group data - shape (n_metabolites, n_samples_per_group)
    g1_data = data_matrix[:, group1_indices].astype(float)
    g2_data = data_matrix[:, group2_indices].astype(float)
    
    # Replace zeros with NaN for proper handling
    g1_data = np.where(g1_data == 0, np.nan, g1_data)
    g2_data = np.where(g2_data == 0, np.nan, g2_data)
    
    # Calculate observed statistics for all metabolites at once
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        # Group means
        mean1 = np.nanmean(g1_data, axis=1)
        mean2 = np.nanmean(g2_data, axis=1)
        mean_diff = mean2 - mean1
        
        # Count valid (non-NaN) samples per metabolite per group
        n1_valid = np.sum(~np.isnan(g1_data), axis=1)
        n2_valid = np.sum(~np.isnan(g2_data), axis=1)
        
        # Variance with ddof=1
        var1 = np.nanvar(g1_data, axis=1, ddof=1)
        var2 = np.nanvar(g2_data, axis=1, ddof=1)
        
        # Pooled variance estimate (assuming equal variances for limma)
        # s^2_pooled = ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2)
        df_resid = n1_valid + n2_valid - 2
        s2_pooled = ((n1_valid - 1) * var1 + (n2_valid - 1) * var2) / np.maximum(df_resid, 1)
        
        # Standard error of the difference
        se = np.sqrt(s2_pooled * (1.0/np.maximum(n1_valid, 1) + 1.0/np.maximum(n2_valid, 1)))
    
    # Mark invalid metabolites (not enough samples)
    valid_mask = (n1_valid >= 2) & (n2_valid >= 2) & (df_resid > 0) & (s2_pooled > 0)
    
    if progress_callback:
        progress_callback(30, 100, "Estimating prior parameters...")
    
    # ========== Empirical Bayes: Estimate prior parameters ==========
    # Filter valid values for prior estimation
    s2_valid = s2_pooled[valid_mask]
    df_valid = df_resid[valid_mask]
    
    if len(s2_valid) < 3:
        # Too few valid metabolites - use standard t-test (no shrinkage)
        d0, s0_squared = 0.0, 1.0
    else:
        # Method of moments estimation for prior parameters
        log_s2 = np.log(s2_valid)
        mean_log_s2 = np.mean(log_s2)
        var_log_s2 = np.var(log_s2, ddof=1)
        
        # Expected variance of log(s2) under the model
        expected_var = np.mean([polygamma(1, d/2) for d in df_valid])
        excess_var = var_log_s2 - expected_var
        
        if excess_var <= 0:
            # No excess variance - strong shrinkage toward prior
            d0 = 1e6
        else:
            # Approximate d0 from excess variance
            d0 = max(2.0 / excess_var, 0.01)
            d0 = min(d0, 1e6)  # Cap at reasonable value
        
        # Estimate s0² from mean of log(s2)
        correction = np.mean([digamma(d/2) - np.log(d/2) for d in df_valid])
        log_s0_sq = mean_log_s2 - correction
        s0_squared = np.exp(log_s0_sq)
        s0_squared = np.clip(s0_squared, 1e-10, np.max(s2_valid) * 10)
    
    if progress_callback:
        progress_callback(60, 100, "Computing moderated statistics...")
    
    # ========== Compute moderated statistics ==========
    # Moderated variance: shrink toward prior
    s2_mod = (d0 * s0_squared + df_resid * s2_pooled) / (d0 + df_resid)
    df_mod = d0 + df_resid
    
    # Moderated standard error
    se_mod = np.sqrt(s2_mod * (1.0/np.maximum(n1_valid, 1) + 1.0/np.maximum(n2_valid, 1)))
    
    # Moderated t-statistic
    t_mod = np.zeros(n_metabolites)
    t_mod[:] = np.nan
    nonzero_se = (se_mod > 0) & valid_mask
    t_mod[nonzero_se] = mean_diff[nonzero_se] / se_mod[nonzero_se]
    
    # P-values from t-distribution with moderated df
    p_values = np.zeros(n_metabolites)
    p_values[:] = np.nan
    
    if progress_callback:
        progress_callback(80, 100, "Computing p-values...")
    
    # Vectorized p-value computation
    valid_for_pval = valid_mask & np.isfinite(t_mod) & (df_mod > 0)
    if np.any(valid_for_pval):
        p_values[valid_for_pval] = 2 * scipy_stats.t.sf(np.abs(t_mod[valid_for_pval]), df_mod[valid_for_pval])
    
    # Floor p-values
    p_values = np.maximum(p_values, np.finfo(float).eps)
    
    if progress_callback:
        progress_callback(100, 100, "Limma analysis complete")
    
    return t_mod, p_values, float(d0), float(s0_squared)


def _limma_test(v1: np.ndarray, v2: np.ndarray, *, seed: Optional[int] = None,
                s0_squared: float = None, d0: float = None) -> tuple[float, float, float, float]:
    """Perform limma-style moderated t-test between two groups (single metabolite).
    
    This is a simplified single-metabolite version. For batch processing, use _limma_batch().
    If s0_squared and d0 are not provided, falls back to standard Welch t-test.
    
    Parameters
    ----------
    v1, v2 : array-like
        Sample values for group 1 and group 2
    seed : int or None
        Not used (kept for API compatibility)
    s0_squared : float or None
        Prior variance from batch estimation. If None, uses standard t-test.
    d0 : float or None
        Prior degrees of freedom. If None, uses standard t-test.
        
    Returns
    -------
    statistic : float
        Moderated t-statistic
    p_value : float
        P-value from moderated t-distribution
    d0 : float
        Prior degrees of freedom used
    s0_squared : float
        Prior variance used
    """
    from scipy import stats as scipy_stats
    
    # Clean data - remove NaN and zeros
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    v1 = v1[~np.isnan(v1) & (v1 != 0)]
    v2 = v2[~np.isnan(v2) & (v2 != 0)]
    
    if len(v1) < 2 or len(v2) < 2:
        return (np.nan, np.nan, np.nan, np.nan)
    
    n1, n2 = len(v1), len(v2)
    
    # Group statistics
    mean_diff = np.mean(v2) - np.mean(v1)
    var1 = np.var(v1, ddof=1)
    var2 = np.var(v2, ddof=1)
    
    # Pooled variance
    df_resid = n1 + n2 - 2
    s2_pooled = ((n1 - 1) * var1 + (n2 - 1) * var2) / df_resid
    
    if s0_squared is None or d0 is None:
        # No prior information - use standard pooled t-test
        se = np.sqrt(s2_pooled * (1.0/n1 + 1.0/n2))
        t_stat = mean_diff / se if se > 0 else np.nan
        p_value = 2 * scipy_stats.t.sf(np.abs(t_stat), df_resid) if not np.isnan(t_stat) else np.nan
        return (float(t_stat), float(p_value), 0.0, float(s2_pooled))
    
    # Moderated variance
    s2_mod = (d0 * s0_squared + df_resid * s2_pooled) / (d0 + df_resid)
    df_mod = d0 + df_resid
    
    # Moderated standard error and t-statistic
    se_mod = np.sqrt(s2_mod * (1.0/n1 + 1.0/n2))
    t_mod = mean_diff / se_mod if se_mod > 0 else np.nan
    
    # P-value from moderated t-distribution
    p_value = 2 * scipy_stats.t.sf(np.abs(t_mod), df_mod) if not np.isnan(t_mod) else np.nan
    
    return (float(t_mod), float(p_value), float(d0), float(s0_squared))

# ---------------------- Detection ----------------------

COMMON_ID_COLUMNS = {
    'Name','HMDB_ID','KEGG_ID','PubChem_CID','ChEBI_ID','LipidMaps_ID',
    'InChIKey','InChI','SMILES','CAS','Formula','Endogenous_Source'
}

FEATURE_COLUMNS_CANONICAL = [
    'Name','Name_Key','Formula','Molecular_Formula','Molecular Formula','MW','Molecular_Weight','ppm',
    'Reference Ion','Reference_Ion','MS2','m/z','RT','RT [min]','Area (Max.)','Polarity','MS2_Purity','MS2 Purity [%]',
    'LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID','ChEBI_ID','CAS','SMILES','InChI','InChIKey','IUPAC_Name',
    'Super_Class','Class','Sub_Class','Endogenous_Source','Metabolika Pathways','BioCyc Pathways',
    # Additional from user examples
    'LipidID', 'Class_name', 'CalcMz', 'BaseRt', 'AdductIon', 'LipidMaps_ID_Match_Type', 'Systematic_Name', 'Preferred_Name',
    'Abbreviation', 'KEGG_Match_Type', 'match_source', 'annotation_sources', 'Endogenous', 'metabolite_id',
    # Mass spectrometry feature columns (commonly seen variations)
    'Calc_mass', 'Calc. MW', 'Calc_mass_(M+H)', 'Calc Mass', 'CalcMass', 'Calculated Mass', 'Theoretical Mass',
    'Delta', 'Delta_Mod', 'DeltaMass', 'Delta_Mass', 'Delta Mass', 'Annot. DeltaMass [ppm]', 'Delta_ppm', 'Delta_PPM',
    'Obs_mass', 'Obs. Mass', 'Observed Mass', 'ObsMass', 'Exact Mass', 'ExactMass', 'Monoisotopic Mass',
    'Charge', 'Adduct', 'Adduct_Ion', 'Ion_Mode', 'IonMode', 'Ion Mode', 'Ion_Formula', 'IonFormula',
    'Score', 'Match_Score', 'Confidence', 'Annotation_Score', 'Ann. Score', 'Annot_Score',
    'Compound', 'Compound_Name', 'Metabolite', 'Metabolite_Name', 'Feature', 'Feature_ID', 'ID',
    'DB_ID', 'Database_ID', 'Ref_ID', 'Reference_ID', 'Index', 'Row_ID', 'RowID',
    'Modification', 'Mod', 'Mod_Type', 'PTM', 'MassError', 'Mass_Error', 'PPM_Error', 'ppm_error',
    'Isotope', 'Isotope_Pattern', 'MS1_Score', 'MS2_Score', 'Fragmentation_Score',
    'Area', 'Peak_Area', 'Intensity', 'Height', 'Peak_Height', 'Signal', 'Response',
    'Ion_Intensity', 'Precursor_Intensity', 'Fragment_Intensity',
]

try:
    from gui.shared.utils import LIPID_FEATURE_CANONICAL, is_statistics_metadata_col as _is_stat_metadata_col
except Exception:
    LIPID_FEATURE_CANONICAL = [
        'lipidid', 'class', 'lipidgroup', 'charge', 'calcmz', 'basert', 'subclass',
        'adduction', 'ionformula', 'molstructure', 'obsmz', 'obsrt', 'ppmdiff', 'polarity'
    ]

    def _is_stat_metadata_col(col_name: str) -> bool:
        try:
            lowered = str(col_name).lower()
        except Exception:
            lowered = str(col_name)
        if lowered.startswith('stat'):
            return True
        if lowered.startswith('adj') and ('p' in lowered or lowered.startswith('adjp')):
            return True
        if lowered.startswith('neg_log10'):
            return True
        for token in ('stat_', '_stat', 'statistic', 'adj_p', 'adjp', 'padj', '_pvalue', '_p_value', '_p', 'neg_log10'):
            if token in lowered:
                return True
        return False

# Ion order constants removed - data cleaner has its own implementation

def detect_sample_columns(df: pd.DataFrame, min_numeric_ratio: float = 0.95) -> List[str]:
    """Heuristically detect sample (intensity/abundance) columns.

    Rules:
      * Column must be numeric dtype OR > min_numeric_ratio of non-null values convertible to float.
      * Column not in the known ID/metadata set.
    """
    sample_cols: List[str] = []
    for col in df.columns:
        if any(col.lower() == c.lower() for c in COMMON_ID_COLUMNS):
            continue
        if _is_stat_metadata_col(col):
            continue
        s = df[col]
        if pd.api.types.is_numeric_dtype(s):
            sample_cols.append(col)
        else:
            # attempt soft numeric detection
            non_null = s.dropna()
            if len(non_null) == 0:
                continue
            convertible = 0
            for v in non_null.head(200):  # sample up to 200
                try:
                    float(str(v).strip())
                    convertible += 1
                except Exception:
                    pass
            if convertible / max(1, len(non_null.head(200))) >= min_numeric_ratio:
                sample_cols.append(col)
    return sample_cols

def detect_feature_and_sample_columns(df: pd.DataFrame, *, min_numeric_ratio: float = 0.95) -> Tuple[List[str], List[str]]:
    """Return (feature_columns_present, sample_columns) using a canonical feature list and numeric detection.

    Any column NOT in the detected feature set and passing numeric heuristics is treated as a sample column.
    """
    feature_present = [col for col in df.columns if any(col.lower() == c.lower() for c in FEATURE_COLUMNS_CANONICAL)]
    sample_cols: List[str] = []
    for col in df.columns:
        if col in feature_present:
            continue
        if _is_stat_metadata_col(col):
            continue
        s = df[col]
        if pd.api.types.is_numeric_dtype(s):
            sample_cols.append(col)
        else:
            non_null = s.dropna()
            if len(non_null) == 0:
                continue
            convertible = 0
            for v in non_null.head(200):
                try:
                    float(str(v).strip())
                    convertible += 1
                except Exception:
                    pass
            if convertible / max(1, len(non_null.head(200))) >= min_numeric_ratio:
                sample_cols.append(col)
    return feature_present, sample_cols

# calculate_ion_order function removed - data cleaner has its own implementation

def clean_sample_column_names(sample_cols: List[str], polarity: str) -> Dict[str, str]:
    """Clean sample column names by removing prefixes and polarity suffixes.
    
    Enhanced to handle multiple polarity suffix variations for robust merging.
    
    Parameters:
    -----------
    sample_cols : List[str]
        List of sample column names
    polarity : str  
        'positive' or 'negative'
        
    Returns:
    --------
    Dict[str, str]
        Mapping of original_name -> cleaned_name
    """
    cleaned_mapping = {}
    
    # Define all possible polarity suffixes to remove (case-insensitive)
    positive_suffixes = ['_pos', '_positive', ' pos', ' positive', '-pos', '-positive']
    negative_suffixes = ['_neg', '_negative', ' neg', ' negative', '-neg', '-negative']
    
    # Choose which suffixes to remove based on polarity
    if polarity.lower() == 'positive':
        suffixes_to_check = positive_suffixes
    else:
        suffixes_to_check = negative_suffixes
    
    for col in sample_cols:
        cleaned = col
        
        # Remove "Group Area: " prefix if present
        if cleaned.startswith('Group Area: '):
            cleaned = cleaned[len('Group Area: '):]
        
        # Remove polarity suffixes (case-insensitive, try all variations)
        col_lower = cleaned.lower()
        for suffix in suffixes_to_check:
            if col_lower.endswith(suffix.lower()):
                # Remove the suffix (preserve original case for the base name)
                cleaned = cleaned[:len(cleaned) - len(suffix)]
                break  # Only remove first matching suffix
        
        # Also strip any trailing/leading whitespace after suffix removal
        cleaned = cleaned.strip()
        
        cleaned_mapping[col] = cleaned
    
    return cleaned_mapping

# ---------------------- Group-Based Data Filtering ----------------------

def apply_min_group_size_filter(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Dict[str, str],
    min_group_size: int = 2,
    min_group_size_type: str = 'absolute',
    min_group_size_percent: float = 50.0
) -> tuple[pd.DataFrame, Dict[str, Dict[str, int]]]:
    """Apply minimum group size filter by setting values to 0 for insufficient groups.
    
    This is applied BEFORE normalization. For each metabolite, independently checks each group.
    If a group has fewer valid values than the threshold, ALL values for that group are zeroed.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with metabolite rows and sample columns
    sample_cols : Sequence[str]
        List of sample column names
    group_map : Dict[str, str]
        Mapping of sample column name to group label
    min_group_size : int
        Minimum number of valid (non-zero, non-NaN) values required per group (for 'absolute' mode)
    min_group_size_type : str
        'absolute' for fixed number, 'percentage' for percent of group size
    min_group_size_percent : float
        Minimum percentage (0-100) required when type is 'percentage'
        
    Returns
    -------
    tuple[pd.DataFrame, Dict[str, Dict[str, int]]]
        Filtered dataframe and dictionary of metabolite counts per group:
        {group_name: {'total_samples': int, 'before': count, 'after': count, 'removed': count}}
    """
    filtered_df = df.copy()
    
    # Build group to columns mapping
    group_cols: Dict[str, List[str]] = {}
    for col in sample_cols:
        grp = group_map.get(col)
        if grp:
            if grp not in group_cols:
                group_cols[grp] = []
            group_cols[grp].append(col)
    
    groups = list(group_cols.keys())
    
    # Calculate threshold per group based on type
    group_thresholds: Dict[str, int] = {}
    for grp, cols in group_cols.items():
        group_size = len(cols)
        if min_group_size_type == 'percentage':
            # Calculate required samples based on percentage
            threshold = int(np.ceil(group_size * min_group_size_percent / 100.0))
        else:
            # Use absolute number
            threshold = min_group_size
        group_thresholds[grp] = threshold
    
    # Track statistics per group
    group_stats: Dict[str, Dict[str, int]] = {}
    for grp in groups:
        group_stats[grp] = {
            'total_samples': len(group_cols[grp]),
            'threshold': group_thresholds[grp],
            'before': 0,
            'after': 0,
            'removed': 0
        }
    
    # Process each metabolite row
    total_rows = len(filtered_df)
    rows_with_changes = 0
    total_cells_zeroed = 0
    
    for idx, row in filtered_df.iterrows():
        row_had_changes = False
        
        for grp, cols in group_cols.items():
            # Get values for this group
            vals = row[cols].apply(pd.to_numeric, errors='coerce').values
            # Count non-zero, non-NaN values
            valid_vals = vals[~np.isnan(vals) & (vals != 0)]
            valid_count = len(valid_vals)
            
            # Track "before" counts (metabolites with ANY valid value in this group)
            if valid_count > 0:
                group_stats[grp]['before'] += 1
            
            # Check against threshold for this specific group
            threshold = group_thresholds[grp]
            
            if valid_count < threshold:
                # Insufficient valid values - zero out ALL values for this group in this metabolite
                if valid_count > 0:  # Only track if we're actually changing something
                    row_had_changes = True
                    total_cells_zeroed += len(cols)
                    group_stats[grp]['removed'] += 1
                # Set all columns for this group to 0
                filtered_df.loc[idx, cols] = 0
            else:
                # Group passes threshold for this metabolite
                group_stats[grp]['after'] += 1
        
        if row_had_changes:
            rows_with_changes += 1
    
    # Remove rows where ALL sample columns are now zero or NaN
    # (i.e., metabolites that failed the filter in ALL groups)
    def is_all_zero_or_nan(row):
        vals = row[sample_cols].apply(pd.to_numeric, errors='coerce').values
        return np.all(np.isnan(vals) | (vals == 0))
    
    rows_to_drop = []
    for idx, row in filtered_df.iterrows():
        if is_all_zero_or_nan(row):
            rows_to_drop.append(idx)
    
    if rows_to_drop:
        filtered_df = filtered_df.drop(rows_to_drop)
    
    return filtered_df, group_stats

def get_group_metabolite_counts(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Dict[str, str],
    min_group_size: int = 2
) -> Dict[str, int]:
    """Count how many metabolites have sufficient data per group.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with metabolite rows and sample columns
    sample_cols : Sequence[str]
        List of sample column names
    group_map : Dict[str, str]
        Mapping of sample column name to group label
    min_group_size : int
        Minimum number of valid (non-zero, non-NaN) values required
        
    Returns
    -------
    Dict[str, int]
        Dictionary mapping group name to count of metabolites with sufficient data
    """
    # Build group to columns mapping
    group_cols: Dict[str, List[str]] = {}
    for col in sample_cols:
        grp = group_map.get(col)
        if grp:
            if grp not in group_cols:
                group_cols[grp] = []
            group_cols[grp].append(col)
    
    # Count metabolites per group
    group_counts: Dict[str, int] = {grp: 0 for grp in group_cols.keys()}
    
    for idx, row in df.iterrows():
        for grp, cols in group_cols.items():
            # Get values for this group
            vals = row[cols].apply(pd.to_numeric, errors='coerce').values
            # Count non-zero, non-NaN values
            valid_vals = vals[~np.isnan(vals) & (vals != 0)]
            
            if len(valid_vals) >= min_group_size:
                group_counts[grp] += 1
    
    return group_counts


def apply_variability_filter(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Optional[Dict[str, str]] = None,
    variance_percentile: float = 10.0,
    compute_anova: bool = True,
    require_testable_rows: bool = False,
) -> tuple[pd.DataFrame, Dict[str, Any]]:
    """Filter low-variability features using a variance percentile threshold.

    ANOVA can be used for diagnostics; when require_testable_rows=True,
    rows that cannot be evaluated across >=2 groups are removed.
    """
    out = df.copy()
    valid_cols = [c for c in sample_cols if c in out.columns]
    if not valid_cols:
        return out, {
            'applied': False,
            'reason': 'no_sample_columns',
            'removed_features': 0,
            'kept_features': len(out)
        }

    var_pct = float(np.clip(variance_percentile, 0.0, 100.0))
    matrix = out[valid_cols].apply(pd.to_numeric, errors='coerce')
    matrix = matrix.replace(0, np.nan)

    row_var = matrix.var(axis=1, ddof=1).fillna(0.0)
    threshold = float(np.nanpercentile(row_var.values, var_pct)) if len(row_var) else 0.0
    keep_mask_low_variance = row_var > threshold

    if not keep_mask_low_variance.any() and len(row_var):
        keep_mask_low_variance = row_var == row_var.max()

    keep_mask = keep_mask_low_variance.copy()
    removed_low_variance = int((~keep_mask_low_variance).sum())

    testable_mask = pd.Series(False, index=out.index)
    removed_untestable = 0

    anova_summary = None
    if compute_anova and group_map:
        groups = {}
        for c in valid_cols:
            grp = group_map.get(c)
            if grp:
                groups.setdefault(grp, []).append(c)

        if len(groups) >= 2:
            f_vals = []
            p_vals = []
            for idx, row in matrix.iterrows():
                arrays = []
                for cols in groups.values():
                    vals = pd.to_numeric(row[cols], errors='coerce').dropna().values
                    vals = vals[vals != 0]
                    if len(vals) > 0:
                        arrays.append(vals)
                if len(arrays) >= 2:
                    testable_mask.loc[idx] = True
                    try:
                        f_stat, p_val = stats.f_oneway(*arrays)
                        if np.isfinite(f_stat):
                            f_vals.append(float(f_stat))
                        if np.isfinite(p_val):
                            p_vals.append(float(p_val))
                    except Exception:
                        continue

            if require_testable_rows:
                keep_mask = keep_mask & testable_mask
                if not keep_mask.any() and testable_mask.any():
                    max_var = row_var[testable_mask].max()
                    keep_mask = testable_mask & (row_var == max_var)
                removed_untestable = int((keep_mask_low_variance & ~testable_mask).sum())

            anova_summary = {
                'rows_total': int(len(out)),
                'rows_testable': int(testable_mask.sum()),
                'rows_tested': len(f_vals),
                'rows_untested': int(len(out) - testable_mask.sum()),
                'mean_f': float(np.nanmean(f_vals)) if f_vals else np.nan,
                'median_f': float(np.nanmedian(f_vals)) if f_vals else np.nan,
                'median_p': float(np.nanmedian(p_vals)) if p_vals else np.nan
            }

    removed = int((~keep_mask).sum())
    filtered = out.loc[keep_mask].copy()

    return filtered, {
        'applied': True,
        'variance_percentile': var_pct,
        'variance_threshold': threshold,
        'testable_rows_required': bool(require_testable_rows),
        'removed_low_variance_features': removed_low_variance,
        'removed_untestable_features': int(removed_untestable),
        'removed_features': removed,
        'kept_features': int(len(filtered)),
        'anova_summary': anova_summary
    }


def apply_imputation(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    method: str = 'half_min',
    knn_neighbors: int = 5,
) -> tuple[pd.DataFrame, Dict[str, Any]]:
    """Impute missing values (NaN/0) in sample columns.

    Supported methods: 'half_min', 'median_per_group', 'median_global', and 'knn'.
    """
    out = df.copy()
    valid_cols = [c for c in sample_cols if c in out.columns]
    if not valid_cols:
        return out, {
            'applied': False,
            'reason': 'no_sample_columns',
            'method': method,
            'imputed_cells': 0
        }

    method_norm = str(method).strip().lower()
    sub = out[valid_cols].apply(pd.to_numeric, errors='coerce')
    sub = sub.replace(0, np.nan)

    missing_before = int(sub.isna().sum().sum())

    unresolved_before_fallback = 0

    if method_norm == 'half_min':
        imputed = sub.copy()
        all_pos_vals = sub.values[np.isfinite(sub.values) & (sub.values > 0)]
        global_half_min = float(np.min(all_pos_vals) / 2.0) if len(all_pos_vals) else 1e-12
        for idx, row in imputed.iterrows():
            vals = row.dropna().values
            vals = vals[vals > 0]
            half_min = float(np.min(vals) / 2.0) if len(vals) else global_half_min
            row = row.fillna(half_min)
            imputed.loc[idx] = row
    elif method_norm == 'median_per_group':
        imputed = sub.copy()
        all_pos_vals = sub.values[np.isfinite(sub.values) & (sub.values > 0)]
        global_median = float(np.median(all_pos_vals)) if len(all_pos_vals) else 1e-12
        for idx, row in imputed.iterrows():
            vals = row.dropna().values
            vals = vals[vals > 0]
            row_median = float(np.median(vals)) if len(vals) else global_median
            row = row.fillna(row_median)
            imputed.loc[idx] = row
    elif method_norm == 'median_global':
        imputed = sub.copy()
        all_pos_vals = sub.values[np.isfinite(sub.values) & (sub.values > 0)]
        global_median = float(np.median(all_pos_vals)) if len(all_pos_vals) else 1e-12
        imputed = imputed.fillna(global_median)
    elif method_norm == 'knn':
        try:
            from sklearn.impute import KNNImputer
        except Exception as e:
            raise ImportError('KNN imputation requires scikit-learn.') from e

        n_neighbors = max(1, int(knn_neighbors))
        # Some columns can become entirely NaN after 0->NaN conversion.
        # KNNImputer may drop these depending on sklearn version/config,
        # so impute only columns with at least one observed value and then
        # restore full shape in original column order.
        knn_cols = [c for c in sub.columns if sub[c].notna().any()]
        dropped_all_nan_cols = [c for c in sub.columns if c not in knn_cols]

        if knn_cols:
            imputer = KNNImputer(n_neighbors=n_neighbors)
            imputed_arr = imputer.fit_transform(sub[knn_cols])
            imputed_knn = pd.DataFrame(imputed_arr, index=sub.index, columns=knn_cols)
        else:
            imputed_knn = pd.DataFrame(index=sub.index)

        imputed = pd.DataFrame(index=sub.index, columns=sub.columns, dtype=float)
        if knn_cols:
            imputed.loc[:, knn_cols] = imputed_knn
        if dropped_all_nan_cols:
            # No information exists for these columns after preprocessing.
            # Keep as NaN so zeros are never reintroduced as pseudo-observations.
            imputed.loc[:, dropped_all_nan_cols] = np.nan

        # Keep unresolved values as NaN for now; we'll apply a deterministic fallback below.
        imputed = imputed.fillna(np.nan)

        # Fallback for unresolved KNN cells (can happen with sparse patterns):
        # use row-wise half-min of observed positive values; if unavailable, use global half-min.
        unresolved_before_fallback = int(imputed.isna().sum().sum())
        if unresolved_before_fallback > 0:
            all_pos_vals = sub.values[np.isfinite(sub.values) & (sub.values > 0)]
            global_half_min = float(np.min(all_pos_vals) / 2.0) if len(all_pos_vals) else 1e-12

            for idx in imputed.index:
                row = imputed.loc[idx]
                if row.isna().any():
                    src_row = sub.loc[idx]
                    row_pos = src_row[(~src_row.isna()) & (src_row > 0)].values
                    row_half_min = float(np.min(row_pos) / 2.0) if len(row_pos) else global_half_min
                    imputed.loc[idx] = row.fillna(row_half_min)
    else:
        raise ValueError(f"Unknown imputation method: {method}")

    out.loc[:, valid_cols] = imputed
    missing_after = int(pd.DataFrame(imputed).isna().sum().sum())
    zeros_after = int((pd.DataFrame(imputed) == 0).sum().sum())
    unresolved_after_fallback = int(pd.DataFrame(imputed).isna().sum().sum()) if method_norm == 'knn' else 0

    return out, {
        'applied': True,
        'method': method_norm,
        'knn_neighbors': int(knn_neighbors) if method_norm == 'knn' else None,
        'imputed_cells': max(0, missing_before - missing_after),
        'missing_before': missing_before,
        'missing_after': missing_after,
        'zeros_after': zeros_after,
        'unresolved_before_fallback': unresolved_before_fallback if method_norm == 'knn' else 0,
        'unresolved_after_fallback': unresolved_after_fallback
    }


def apply_pca_outlier_filter(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Optional[Dict[str, str]] = None,
    threshold_sd: float = 3.0,
) -> tuple[pd.DataFrame, Dict[str, Any]]:
    """Remove sample outliers using PCA scores (PC1-PC3 distance from center)."""
    out = df.copy()
    valid_cols = [c for c in sample_cols if c in out.columns]
    if len(valid_cols) < 3:
        return out, {
            'applied': False,
            'reason': 'insufficient_samples',
            'removed_samples': [],
            'candidate_outliers': []
        }

    matrix = out[valid_cols].apply(pd.to_numeric, errors='coerce').replace(0, np.nan)
    if matrix.isna().any().any():
        return out, {
            'applied': False,
            'reason': 'missing_values_present',
            'removed_samples': [],
            'candidate_outliers': []
        }

    # Scale each feature (row-wise z-score), then transpose to sample x feature for PCA.
    arr = matrix.values.astype(float)
    row_mean = np.nanmean(arr, axis=1, keepdims=True)
    row_std = np.nanstd(arr, axis=1, ddof=1, keepdims=True)
    row_std[(~np.isfinite(row_std)) | (row_std == 0)] = 1.0
    scaled = (arr - row_mean) / row_std
    scaled = np.nan_to_num(scaled, nan=0.0)
    sample_matrix = scaled.T

    n_components = int(min(3, sample_matrix.shape[0], sample_matrix.shape[1]))
    if n_components < 2:
        return out, {
            'applied': False,
            'reason': 'insufficient_dimensions',
            'removed_samples': [],
            'candidate_outliers': []
        }

    try:
        from sklearn.decomposition import PCA
    except Exception as e:
        raise ImportError('PCA outlier detection requires scikit-learn.') from e

    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(sample_matrix)

    center = np.mean(scores[:, :n_components], axis=0)
    dists = np.sqrt(np.sum((scores[:, :n_components] - center) ** 2, axis=1))
    dist_mean = float(np.mean(dists))
    dist_sd = float(np.std(dists, ddof=1)) if len(dists) > 1 else 0.0
    thresh = dist_mean + float(threshold_sd) * dist_sd

    candidate_idx = np.where(dists > thresh)[0].tolist()
    candidate_samples = [valid_cols[i] for i in candidate_idx]

    removed_samples = candidate_samples.copy()
    protected_samples: List[str] = []

    # Conservative guard: keep candidates that are not outliers within their own group.
    if group_map and candidate_samples:
        kept = []
        for sample in candidate_samples:
            grp = group_map.get(sample)
            if not grp:
                kept.append(sample)
                continue
            grp_samples = [s for s in valid_cols if group_map.get(s) == grp]
            if len(grp_samples) < 3:
                kept.append(sample)
                continue
            grp_idx = [valid_cols.index(s) for s in grp_samples]
            grp_scores = scores[grp_idx, :n_components]
            grp_center = np.mean(grp_scores, axis=0)
            grp_d = np.sqrt(np.sum((grp_scores - grp_center) ** 2, axis=1))
            grp_thresh = float(np.mean(grp_d) + float(threshold_sd) * np.std(grp_d, ddof=1)) if len(grp_d) > 1 else np.inf
            sample_d = float(np.sqrt(np.sum((scores[valid_cols.index(sample), :n_components] - grp_center) ** 2)))
            if sample_d <= grp_thresh:
                kept.append(sample)

        protected_samples = kept
        removed_samples = [s for s in candidate_samples if s not in protected_samples]

    if removed_samples:
        out = out.drop(columns=removed_samples, errors='ignore')

    scores_df = pd.DataFrame(
        {
            'sample': valid_cols,
            'distance': dists,
            'is_candidate_outlier': [s in candidate_samples for s in valid_cols],
            'is_removed_outlier': [s in removed_samples for s in valid_cols],
            'group': [group_map.get(s) if group_map else None for s in valid_cols],
        }
    )
    for i in range(n_components):
        scores_df[f'PC{i + 1}'] = scores[:, i]

    return out, {
        'applied': True,
        'threshold_sd': float(threshold_sd),
        'threshold_distance': float(thresh),
        'distance_mean': dist_mean,
        'distance_sd': dist_sd,
        'candidate_outliers': candidate_samples,
        'protected_samples': protected_samples,
        'removed_samples': removed_samples,
        'scores': scores_df
    }

# ---------------------- Normalization ----------------------

def _quantile_normalize(values: pd.DataFrame) -> pd.DataFrame:
    """Perform quantile normalization across columns.

    Standard algorithm:
      1. Sort each column.
      2. Compute mean for each rank across columns.
      3. Map mean values back to original ranks per column.
    """
    if values.shape[1] <= 1:
        return values.copy()
    # Sort values per column
    # Force numpy float arrays to avoid ExtensionArray typing issues
    sorted_dict = {}
    for c in values.columns:
        arr = pd.to_numeric(values[c], errors='coerce').to_numpy(dtype=float, copy=True)
        sorted_dict[c] = np.sort(arr)
    sorted_df = pd.DataFrame(sorted_dict)
    rank_means = sorted_df.mean(axis=1).values
    # Assign back
    result = values.copy()
    for c in values.columns:
        col_arr = pd.to_numeric(values[c], errors='coerce').to_numpy(dtype=float, copy=False)
        order = np.argsort(col_arr)
        reverse_order = np.empty_like(order)
        reverse_order[order] = np.arange(len(order))
        result[c] = rank_means[reverse_order]
    return result


def _pqn_normalize(values: pd.DataFrame, reference_sample_idx: Optional[int] = None) -> pd.DataFrame:
    """Perform Probabilistic Quotient Normalization (PQN).
    
    PQN is widely used in metabolomics to correct for dilution effects.
    Algorithm (Dieterle et al. 2006):
      1. Compute reference spectrum (median across all samples by default, or use a specific QC sample)
      2. For each sample, calculate quotients of all variables relative to reference
      3. Calculate median of quotients (excluding zeros and NaNs)
      4. Divide sample by its median quotient
    
    Args:
        values: DataFrame with numeric intensity values
        reference_sample_idx: Optional index of reference sample column (default: use median spectrum)
    
    Returns:
        PQN-normalized DataFrame
    """
    if values.shape[1] <= 1:
        return values.copy()
    
    result = values.copy()
    
    # Step 1: Create reference spectrum (median of all samples or specific reference)
    if reference_sample_idx is not None and 0 <= reference_sample_idx < len(values.columns):
        reference = values.iloc[:, reference_sample_idx].values
    else:
        # Use median across all samples for each metabolite
        reference = values.median(axis=1).values
    
    # Step 2-4: For each sample, calculate quotients and normalize
    for col in values.columns:
        sample_values = values[col].values
        
        # Calculate quotients (sample / reference)
        with np.errstate(divide='ignore', invalid='ignore'):
            quotients = sample_values / reference
        
        # Filter out zeros, infinities, and NaNs for median calculation
        valid_quotients = quotients[(np.isfinite(quotients)) & (quotients > 0) & (reference > 0)]
        
        if len(valid_quotients) > 0:
            # Median quotient for this sample
            median_quotient = np.median(valid_quotients)
            
            # Normalize sample by dividing by median quotient
            if median_quotient > 0 and np.isfinite(median_quotient):
                result[col] = sample_values / median_quotient
            else:
                result[col] = sample_values
        else:
            result[col] = sample_values
    
    return result


def _vsn_normalize(values: pd.DataFrame) -> pd.DataFrame:
    """Perform Variance Stabilizing Normalization (VSN).
    
    VSN applies an arcsinh transformation that stabilizes variance across the intensity range.
    This is particularly useful for mass spectrometry data where variance increases with intensity.
    
    The arcsinh transformation: arcsinh(x) = log(x + sqrt(x^2 + 1))
    
    VSN is similar to log transformation but:
    - Better handles low intensities (no need for pseudo-count)
    - Stabilizes variance more effectively across full intensity range
    - Commonly used in proteomics and metabolomics
    
    Reference: Huber et al. (2002) "Variance stabilization applied to microarray data calibration"
    
    Args:
        values: DataFrame with numeric intensity values
    
    Returns:
        VSN-transformed DataFrame
    """
    if values.shape[1] == 0:
        return values.copy()
    
    result = values.copy()
    
    # Apply arcsinh transformation: arcsinh(x) = log(x + sqrt(x^2 + 1))
    # NumPy has arcsinh built-in
    result = pd.DataFrame(
        np.arcsinh(values.values),
        index=values.index,
        columns=values.columns
    )
    
    return result


def _median_normalize(values: pd.DataFrame) -> pd.DataFrame:
    """Perform median normalization.
    
    Each sample column is divided by its median value, then optionally scaled.
    This is robust to outliers and commonly used in metabolomics.
    
    Algorithm:
      1. Calculate median of each sample (column)
      2. Divide each sample by its median
      3. Multiply by global median of all medians (to preserve original scale)
    
    Args:
        values: DataFrame with numeric intensity values
    
    Returns:
        Median-normalized DataFrame
    """
    if values.shape[1] <= 1:
        return values.copy()
    
    result = values.copy()
    
    # Calculate median for each sample
    sample_medians = values.median(axis=0)
    
    # Global median (median of all sample medians) for scaling
    global_median = sample_medians.median()
    
    # Normalize each sample
    for col in values.columns:
        col_median = sample_medians[col]
        if col_median > 0 and np.isfinite(col_median):
            result[col] = (values[col] / col_median) * global_median
        else:
            result[col] = values[col]
    
    return result


def _tic_normalize(values: pd.DataFrame) -> pd.DataFrame:
    """Perform Total Ion Current (TIC) normalization.
    
    TIC normalization divides each sample by its total ion count (sum of all features),
    then multiplies by the mean TIC to preserve scale. This is standard in mass spectrometry.
    
    Algorithm:
      1. Calculate sum (total ion current) for each sample
      2. Divide each sample by its TIC
      3. Multiply by mean TIC across all samples
    
    Args:
        values: DataFrame with numeric intensity values
    
    Returns:
        TIC-normalized DataFrame
    """
    if values.shape[1] <= 1:
        return values.copy()
    
    result = values.copy()
    
    # Calculate TIC (sum) for each sample
    sample_sums = values.sum(axis=0)
    
    # Mean TIC for scaling
    mean_tic = sample_sums.mean()
    
    # Normalize each sample
    for col in values.columns:
        col_sum = sample_sums[col]
        if col_sum > 0 and np.isfinite(col_sum):
            result[col] = (values[col] / col_sum) * mean_tic
        else:
            result[col] = values[col]
    
    return result


def _clr_transform(values: pd.DataFrame, pseudo_count: float = 1.0) -> pd.DataFrame:
    """Perform Centered Log-Ratio (CLR) transformation.
    
    CLR is designed for compositional data (e.g., relative abundances that sum to a constant).
    It addresses the issue of compositional data being constrained to a simplex.
    
    Algorithm (Aitchison 1986):
      1. Add pseudo-count to avoid log(0)
      2. Calculate geometric mean for each sample
      3. Take log-ratio of each feature to geometric mean
    
    Formula: CLR(x_i) = log(x_i / geometric_mean(x))
    
    Args:
        values: DataFrame with numeric intensity values
        pseudo_count: Small value added to avoid log(0) (default: 1.0)
    
    Returns:
        CLR-transformed DataFrame
    """
    if values.shape[1] <= 1:
        return values.copy()
    
    result = values.copy()
    
    # Add pseudo-count to avoid log(0)
    data_with_pseudo = values + pseudo_count
    
    # For each sample (column), calculate CLR
    for col in values.columns:
        sample = data_with_pseudo[col].values
        
        # Calculate geometric mean (using log to avoid overflow)
        with np.errstate(divide='ignore', invalid='ignore'):
            log_sample = np.log(sample)
            valid_logs = log_sample[np.isfinite(log_sample)]
            
            if len(valid_logs) > 0:
                geometric_mean = np.exp(np.mean(valid_logs))
                # CLR transformation
                result[col] = np.log(sample / geometric_mean)
            else:
                result[col] = np.nan
    
    return result


def _loess_qc_correction(values: pd.DataFrame, qc_indices: Optional[List[int]] = None, 
                         frac: float = 0.3) -> pd.DataFrame:
    """Perform LOESS-based QC drift correction.
    
    This method corrects for systematic drift in MS signal over the analytical run order.
    It uses QC (Quality Control) samples run throughout the batch to model and correct drift.
    
    Algorithm:
      1. Identify QC samples (or assume evenly spaced if not provided)
      2. For each metabolite, fit LOESS curve to QC samples vs. run order
      3. Use LOESS model to predict drift at each sample position
      4. Correct all samples by dividing by predicted drift
    
    Args:
        values: DataFrame with numeric intensity values (rows=features, cols=samples in run order)
        qc_indices: List of column indices that are QC samples (if None, auto-detect or use all)
        frac: Fraction of data used for LOESS smoothing (default: 0.3)
    
    Returns:
        Drift-corrected DataFrame
    
    Note: Requires scikit-learn for lowess. If not available, returns original data.
    """
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
    except ImportError:
        warnings.warn("statsmodels not available for LOESS correction. Returning uncorrected data.")
        return values.copy()
    
    if values.shape[1] <= 3:  # Need at least 3 samples for LOESS
        return values.copy()
    
    result = values.copy()
    n_samples = values.shape[1]
    
    # If no QC indices provided, use all samples (or detect by name pattern)
    if qc_indices is None:
        # Try to detect QC samples by column name (look for 'qc', 'quality', 'pool', etc.)
        qc_indices = []
        qc_keywords = ['qc', 'quality', 'pool', 'blank', 'pooled']
        for idx, col_name in enumerate(values.columns):
            if any(keyword in str(col_name).lower() for keyword in qc_keywords):
                qc_indices.append(idx)
        
        # If insufficient QC samples, supplement with evenly spaced pseudo-QC points
        if len(qc_indices) < 3:
            existing_qc = set(qc_indices)
            step = max(1, n_samples // 10)  # Use ~10% of samples as pseudo-QC
            pseudo_qc = [i for i in range(0, n_samples, step) if i not in existing_qc]
            qc_indices = sorted(list(existing_qc) + pseudo_qc)
    
    if len(qc_indices) < 3:
        warnings.warn("Insufficient QC samples for LOESS correction. Need at least 3 samples.")
        return values.copy()
    
    # Run order (sample index)
    run_order = np.arange(n_samples)
    qc_run_order = np.array(qc_indices)
    
    # For each metabolite (row), fit LOESS and correct
    for row_idx in range(len(values)):
        metabolite_values = values.iloc[row_idx, :].values
        qc_values = metabolite_values[qc_indices]
        
        # Filter out NaN and zero values from QC samples
        valid_mask = np.isfinite(qc_values) & (qc_values > 0)
        if valid_mask.sum() < 3:  # Need at least 3 valid QC points
            continue
        
        valid_qc_order = qc_run_order[valid_mask]
        valid_qc_values = qc_values[valid_mask]
        
        try:
            # Fit LOESS to QC samples
            # lowess returns (x, y_smoothed) sorted by x
            smoothed = lowess(valid_qc_values, valid_qc_order, frac=frac, return_sorted=True)
            
            # Interpolate LOESS curve to all sample positions
            # Use linear interpolation for positions between QC samples
            drift_curve = np.interp(run_order, smoothed[:, 0], smoothed[:, 1])
            
            # Prevent division by zero or negative drift
            drift_curve = np.maximum(drift_curve, 1e-10)
            
            # Normalize to mean drift = 1 (preserve overall intensity)
            drift_curve = drift_curve / np.mean(drift_curve)
            
            # Correct all samples by dividing by drift
            result.iloc[row_idx, :] = metabolite_values / drift_curve
            
        except Exception as e:
            # If LOESS fails for this metabolite, leave it uncorrected
            warnings.warn(f"LOESS correction failed for row {row_idx}: {e}")
            continue
    
    return result


def _internal_standard_normalize(values: pd.DataFrame, is_column: Optional[str] = None) -> pd.DataFrame:
    """Perform Internal Standard (IS) normalization.
    
    Internal standard normalization divides each sample by the intensity of a known
    internal standard added to all samples at a fixed concentration. This is the gold
    standard for targeted metabolomics and lipidomics.
    
    Algorithm:
      1. Identify the internal standard feature (row)
      2. For each sample, divide all features by the IS intensity
      3. Multiply by mean IS intensity to preserve scale
    
    Args:
        values: DataFrame with numeric intensity values (rows=features, cols=samples)
        is_column: Name/identifier of the internal standard feature (row). If None, uses first row.
    
    Returns:
        IS-normalized DataFrame
    
    Note: Internal standard should be present in 'Name' or first column of parent dataframe
    """
    if values.shape[0] <= 1:  # Need at least 2 features (IS + analytes)
        warnings.warn("Internal Standard normalization requires at least 2 features.")
        return values.copy()
    
    result = values.copy()
    
    # If no IS specified, issue warning and return unchanged
    if is_column is None:
        warnings.warn("No internal standard specified for IS normalization. Data returned unchanged. "
                     "Please specify internal standard in column verification step.")
        return values.copy()
    
    # Internal standard should be identified by row name/index
    # This assumes the is_column parameter points to a row identifier
    # The actual IS row lookup will be handled in normalize_dataframe
    # For now, this is a placeholder that will be called with pre-extracted IS values
    
    return result


def normalize_dataframe(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    method: str,
    *,
    log_pseudo: float = 1.0,
    preserve_zeros_for_log: bool = True,
    is_feature_name: Optional[str] = None,
    qc_sample_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Return a new DataFrame with specified normalization applied to sample columns.

    Supported methods:
      - none: return df copy unchanged
      - Rel_Abundance(%) (Percent Relative Abundance): col_value / column_sum * 100
      - log2: log2(x + log_pseudo)
      - quantile: quantile normalization across sample columns
      - log2_quantile: log2 after quantile normalization
      - zscore: feature-wise z-score across samples (per row)
      - PQN: Probabilistic Quotient Normalization (Dieterle et al. 2006) - metabolomics standard
      - median: Median normalization - robust to outliers
      - TIC: Total Ion Current normalization - standard for mass spectrometry
      - VSN: Variance Stabilizing Normalization - arcsinh transformation for MS data
      - CLR: Centered Log-Ratio transformation - for compositional data
      - LOESS_QC: LOESS-based QC drift correction - corrects run-order effects
      - IS: Internal Standard normalization - gold standard for targeted metabolomics/lipidomics

    Chaining methods:
      - Use '+' or '->' to chain multiple normalizations sequentially
      - Examples: 'median+log2', 'TIC->log2', 'IS+quantile+log2'
      - Each method is applied in order from left to right
      - 'log2_quantile' is equivalent to 'quantile+log2'

    Behavior regarding zeros:
        - If preserve_zeros_for_log is True (default), positions that were exactly 0 in the raw
          data (within the provided sample_cols) are set back to 0 after applying 'quantile',
          'log2', 'log2_quantile', or 'zscore'. This ensures raw zeros remain zeros so that
          downstream filtering/statistics that treat 0 as "absent" continue to exclude them
          consistently across methods.

    Args:
        df: DataFrame with metabolite/lipid data
        sample_cols: List of sample column names
        method: Normalization method name or chain (e.g., 'median+log2')
        log_pseudo: Pseudo-count for log transformations (default: 1.0)
        preserve_zeros_for_log: Whether to preserve zeros after normalization (default: True)
        is_feature_name: Name of internal standard feature (for IS normalization)
        qc_sample_cols: List of QC sample column names (for LOESS_QC)
    """
    original_method = method or 'none'
    
    # Check if chaining is requested (using '+' or '->' as separators)
    if '+' in original_method or '->' in original_method:
        # Split by both separators
        chain = original_method.replace('->', '+').split('+')
        chain = [m.strip() for m in chain if m.strip()]
        
        # Apply each normalization in sequence
        current_df = df.copy()
        for step_method in chain:
            current_df = normalize_dataframe(
                current_df, 
                sample_cols, 
                step_method,
                log_pseudo=log_pseudo,
                preserve_zeros_for_log=preserve_zeros_for_log,
                is_feature_name=is_feature_name,
                qc_sample_cols=qc_sample_cols
            )
        return current_df
    
    method = original_method.lower()
    new_df = df.copy()
    if not sample_cols:
        return new_df

    sub = new_df[list(sample_cols)].apply(pd.to_numeric, errors='coerce')
    # Treat exact zeros as missing values for all normalization workflows.
    zero_mask = (sub == 0)
    sub = sub.mask(zero_mask, np.nan)

    if method == 'none':
        new_df[sample_cols] = sub
        return new_df

    # Accept both exact case display name and lowercase variant
    if original_method == 'Rel_Abundance(%)' or method in ('rel_abundance(%)','rel_abundance','relative','relative_abundance'):
        col_sums = sub.sum(axis=0)
        col_sums_replaced = col_sums.replace({0: np.nan})
        norm = (sub / col_sums_replaced) * 100.0
        new_df[sample_cols] = norm
        return new_df

    if method == 'quantile':
        norm = _quantile_normalize(sub)
        if preserve_zeros_for_log:
            # Keep original zeros as missing after normalization.
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'log2':
        norm = pd.DataFrame(np.log2(sub + float(log_pseudo)), index=sub.index, columns=sub.columns)
        if preserve_zeros_for_log:
            # Keep original zeros as missing after normalization.
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method in ('log2_quantile','log2quantile'):
        qn = _quantile_normalize(sub)
        norm = pd.DataFrame(np.log2(qn + float(log_pseudo)), index=qn.index, columns=qn.columns)
        if preserve_zeros_for_log:
            # Keep original zeros as missing after normalization.
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'zscore':
        # z = (x - mean_row)/std_row across sample columns (per row)
        row_means = sub.mean(axis=1)
        row_stds = sub.std(axis=1).replace({0: np.nan})
        norm = (sub.sub(row_means, axis=0)).div(row_stds, axis=0)
        if preserve_zeros_for_log:
            # Keep original zeros as missing after normalization.
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'pqn':
        # Probabilistic Quotient Normalization
        norm = _pqn_normalize(sub)
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'median':
        # Median normalization
        norm = _median_normalize(sub)
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'tic':
        # Total Ion Current normalization
        norm = _tic_normalize(sub)
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'vsn':
        # Variance Stabilizing Normalization
        norm = _vsn_normalize(sub)
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method == 'clr':
        # Centered Log-Ratio transformation
        norm = _clr_transform(sub, pseudo_count=log_pseudo)
        # CLR already involves log, so preserve zeros differently
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method in ('loess_qc', 'loess', 'qc_correction'):
        # LOESS QC drift correction
        norm = _loess_qc_correction(sub)
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=np.nan)
        new_df[sample_cols] = norm
        return new_df

    if method in ('is', 'internal_standard', 'istd'):
        # Internal Standard normalization
        if is_feature_name is None:
            warnings.warn("Internal Standard normalization selected but no IS feature specified. "
                         "Please identify Internal Standard in Column Verification step. "
                         "Returning data unchanged.")
            return new_df
        
        # Find the IS feature row
        id_col = None
        for possible_id in ['Name', 'Metabolite', 'Feature ID', 'metabolite_id', 'LipidID', 'Lipid_ID']:
            if possible_id in df.columns:
                id_col = possible_id
                break
        
        if id_col is None:
            warnings.warn("Cannot identify feature name column for IS normalization. "
                         "Returning data unchanged.")
            return new_df
        
        # Find IS row
        is_row_mask = df[id_col].astype(str).str.strip().str.lower() == str(is_feature_name).strip().lower()
        if not is_row_mask.any():
            warnings.warn(f"Internal Standard '{is_feature_name}' not found in data. "
                         "Returning data unchanged.")
            return new_df
        
        # Extract IS intensities for normalization
        is_row_data = sub[is_row_mask]
        if len(is_row_data) == 0:
            warnings.warn(f"No valid IS data found for '{is_feature_name}'. Returning data unchanged.")
            return new_df
        
        is_intensities = is_row_data.iloc[0]  # Get first matching row
        mean_is = is_intensities.mean()
        
        # Normalize: divide each sample by its IS intensity, multiply by mean IS
        norm = sub.copy()
        for col in sample_cols:
            is_val = is_intensities[col]
            if is_val > 0 and np.isfinite(is_val):
                norm[col] = (sub[col] / is_val) * mean_is
            else:
                # If IS is zero or invalid, leave sample unchanged
                norm[col] = sub[col]
        
        if preserve_zeros_for_log:
            norm = norm.where(~zero_mask, other=0.0)
        new_df[sample_cols] = norm
        return new_df

    raise ValueError(f"Unknown normalization method: {method}")

# ---------------------- Normality Testing ----------------------

def perform_normality_test(df: pd.DataFrame, sample_cols: Sequence[str], group_definitions: Optional[Dict[str, List[str]]] = None, progress_callback: Optional[Callable[[int, int, str], None]] = None) -> pd.DataFrame:
    """Test normality of data distribution for each metabolite across samples using Shapiro-Wilk test.
    
    IMPORTANT: Tests normality PER METABOLITE across samples, not per sample across metabolites.
    Each metabolite has multiple sample values, and we test if those values are normally distributed.
    
    Args:
        df: DataFrame containing metabolite data with an identifier column
        sample_cols: List of sample column names to test
        group_definitions: Optional dict mapping group names to lists of sample columns.
                          If provided, tests normality per metabolite within each group.
                          Otherwise tests each metabolite across all samples.
        progress_callback: Optional callback(current, total, metabolite_name) for progress updates.
    
    Returns:
        DataFrame with columns: ['Metabolite', 'Group', 'N', 'Shapiro_p', 'Is_Normal']
        Is_Normal is 'Yes' if p-value > 0.05, 'No' if p-value <= 0.05
    """
    from scipy import stats
    
    results = []
    
    # Find identifier column
    id_col = None
    for possible_id in ['Name', 'Metabolite', 'Feature ID', 'metabolite_id', 'LipidID']:
        if possible_id in df.columns:
            id_col = possible_id
            break
    
    if id_col is None:
        # Use index as identifier
        df_with_id = df.copy()
        df_with_id['_temp_id'] = df_with_id.index
        id_col = '_temp_id'
    else:
        df_with_id = df
    
    # Get numeric data for sample columns
    sample_data = df_with_id[list(sample_cols)].apply(pd.to_numeric, errors='coerce')
    metabolite_ids = df_with_id[id_col]
    
    if group_definitions:
        # Test each metabolite within each group
        total_tests = len(metabolite_ids) * len(group_definitions)
        callback_freq = max(1, total_tests // 100)
        processed = 0
        
        for group_name, group_samples in group_definitions.items():
            # Get valid samples that exist in the data
            valid_samples = [s for s in group_samples if s in sample_data.columns]
            if not valid_samples:
                continue
            
            # Test each metabolite
            for idx, metabolite_id in enumerate(metabolite_ids):
                processed += 1
                if progress_callback and (processed % callback_freq == 0 or processed == total_tests):
                    try:
                        progress_callback(processed, total_tests, f"{metabolite_id} ({group_name})")
                    except Exception:
                        pass
                
                # Get values for this metabolite across samples in this group
                metabolite_values = sample_data.iloc[idx][valid_samples].values
                # Remove NaN values
                metabolite_values = metabolite_values[~np.isnan(metabolite_values)]
                
                n = len(metabolite_values)
                
                if n < 3:
                    # Not enough data for normality test
                    results.append({
                        'Metabolite': metabolite_id,
                        'Sample/Group': group_name,
                        'Kind': 'GROUP',
                        'N': n,
                        'Shapiro_p': np.nan,
                        'Is_Normal': 'Insufficient data'
                    })
                    continue
                
                # Shapiro-Wilk test (best for small to medium samples)
                shapiro_p = np.nan
                if n < 5000:  # Shapiro-Wilk has limitations for very large samples
                    try:
                        _, shapiro_p = stats.shapiro(metabolite_values)
                    except Exception:
                        shapiro_p = np.nan
                
                # Determine normality: use Shapiro-Wilk p-value
                is_normal = 'Yes' if (not np.isnan(shapiro_p) and shapiro_p > 0.05) else 'No' if not np.isnan(shapiro_p) else 'Unknown'
                
                results.append({
                    'Metabolite': metabolite_id,
                    'Sample/Group': group_name,
                    'Kind': 'GROUP',
                    'N': n,
                    'Shapiro_p': shapiro_p,
                    'Is_Normal': is_normal
                })
    else:
        # Test each metabolite across all samples
        total_metabolites = len(metabolite_ids)
        callback_freq = max(1, total_metabolites // 100)
        
        for idx, metabolite_id in enumerate(metabolite_ids):
            if progress_callback and ((idx + 1) % callback_freq == 0 or (idx + 1) == total_metabolites):
                try:
                    progress_callback(idx + 1, total_metabolites, str(metabolite_id))
                except Exception:
                    pass
            
            # Get values for this metabolite across all samples
            metabolite_values = sample_data.iloc[idx].values
            # Remove NaN values
            metabolite_values = metabolite_values[~np.isnan(metabolite_values)]
            
            n = len(metabolite_values)
            
            if n < 3:
                results.append({
                    'Metabolite': metabolite_id,
                    'Sample/Group': 'All Samples',
                    'Kind': 'SAMPLE',
                    'N': n,
                    'Shapiro_p': np.nan,
                    'Is_Normal': 'Insufficient data'
                })
                continue
            
            # Shapiro-Wilk test (best for small to medium samples)
            shapiro_p = np.nan
            if n < 5000:  # Shapiro-Wilk has limitations for very large samples
                try:
                    _, shapiro_p = stats.shapiro(metabolite_values)
                except Exception:
                    shapiro_p = np.nan
            
            # Determine normality: use Shapiro-Wilk p-value
            is_normal = 'Yes' if (not np.isnan(shapiro_p) and shapiro_p > 0.05) else 'No' if not np.isnan(shapiro_p) else 'Unknown'
            
            results.append({
                'Metabolite': metabolite_id,
                'Sample/Group': 'All Samples',
                'Kind': 'SAMPLE',
                'N': n,
                'Shapiro_p': shapiro_p,
                'Is_Normal': is_normal
            })
    
    return pd.DataFrame(results)

# ---------------------- Statistics ----------------------

@dataclass
class OverallTestResult:
    metabolite: Any
    statistic: float
    p_value: float

@dataclass
class PairwiseTestResult:
    metabolite: Any
    group1: str
    group2: str
    statistic: float
    p_value: float

def _bh_adjust(pvals: List[float]) -> List[float]:
    """Apply Benjamini-Hochberg FDR correction, preserving NaN values.
    
    NaN p-values are excluded from adjustment and returned as NaN in output.
    """
    # Identify valid (non-NaN) p-values
    valid_mask = np.array([p is not None and not np.isnan(p) for p in pvals])
    valid_pvals = np.array([_floor_pval(p) if p is not None else 1.0 for p, v in zip(pvals, valid_mask) if v], dtype=float)
    
    if len(valid_pvals) == 0:
        # All NaN - return all NaN
        return [np.nan] * len(pvals)
    
    n = len(valid_pvals)
    order = np.argsort(valid_pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n+1)
    q = valid_pvals * n / ranks
    # cumulative minimum from largest to smallest
    q_corrected = np.minimum.accumulate(q[order][::-1])[::-1]
    adjusted_valid = np.empty_like(q_corrected)
    adjusted_valid[order] = np.clip(q_corrected, 0, 1)
    
    # Map back to original positions, inserting NaN where original was NaN
    result = []
    valid_idx = 0
    for i, is_valid in enumerate(valid_mask):
        if is_valid:
            result.append(adjusted_valid[valid_idx])
            valid_idx += 1
        else:
            result.append(np.nan)
    
    return result

def perform_statistical_analysis(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    group_map: Dict[str, str],
    *,
    group_order: Optional[Sequence[str]] = None,
    overall_test: Optional[str] = None,
    pairwise_test: Optional[str] = None,
    fdr: bool = True,
    base_group: Optional[str] = None,
    custom_comparisons: Optional[List[Tuple[str, str]]] = None,
    fdr_scope: str = 'per-metabolite',  # 'per-comparison', 'per-metabolite', or 'global'
    alpha: float = 0.05,
    use_adjusted_pvalues: bool = True,
    min_group_size: int = 2,
    min_group_size_percent: Optional[float] = None,
    gate_posthoc_by_overall: bool = True,
    n_jobs: int = 3,
    id_column_name: str = 'metabolite',
    rots_B: int = 1000,
    rots_K: int = 100,
    rots_alpha: float = 0.1,
    rots_seed: Optional[int] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> Dict[str, pd.DataFrame]:
    """Run per-metabolite statistical analysis with comprehensive output.

    Parameters
    ----------
    df : DataFrame
        Input (already normalized) data.
    sample_cols : list
        Sample columns considered.
    group_map : dict
        Mapping sample_col -> group name (must cover all sample_cols used).
    overall_test : str | None
        One of {'anova','kruskal'} for >2 groups. Ignored for 2 groups.
    pairwise_test : str | None
        One of {'welch','mannwhitney','rots','limma'} or None. For >2 groups does all pairwise; for 2 groups overrides overall_test.
    fdr : bool
        Apply Benjamini-Hochberg to p-values.
    fdr_scope : str
        FDR correction scope: 'per-comparison' (adjust within each pairwise comparison across metabolites), 
        'per-metabolite' (adjust within each metabolite across comparisons, R-like), or 'global' (adjust all p-values together).
    base_group : str | None
        Reference group for comparisons. If provided, only comparisons against this base group will be performed.
        Takes precedence over custom_comparisons.
    custom_comparisons : List[Tuple[str, str]] | None
        List of specific group pairs to compare, e.g., [('Group1', 'Group2'), ('Group3', 'Group4')].
        Only used if base_group is None. If both are None, all pairwise comparisons are performed.
    alpha : float
        Significance threshold for classification.
    use_adjusted_pvalues : bool
        If True (default), adjusted p-values are saved in the combined result columns.
        If False, raw p-values are saved instead (affects column naming: *_pvalue vs *_adj_p).
    min_group_size : int
        Minimum number of replicates required per group. Groups with fewer samples are excluded
        from all statistical tests and downstream visualizations.
        Used when min_group_size_percent is None (absolute count mode).
    min_group_size_percent : float | None
        If provided, minimum percentage (0-100) of group size required per comparison.
        For example, 70.0 means at least 70% of samples in each group must have non-zero values.
        When set, this takes precedence over min_group_size for per-metabolite filtering.
    rots_B : int
        Number of bootstrap iterations for ROTS (default 1000). Only used when pairwise_test='rots'.
    rots_K : int
        Number of top features for ROTS reproducibility optimization (default 100). Only used when pairwise_test='rots'.
    rots_alpha : float
        Top proportion parameter for ROTS optimization (0-1, default 0.1). Only used when pairwise_test='rots'.
    rots_seed : int or None
        Random seed for ROTS reproducibility. If None, results may vary between runs. Only used when pairwise_test='rots'.
    
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary containing:
        - 'overall': Overall test results with p-values, adj p-values, and visualization metrics
        - 'pairwise': Pairwise test results with fold changes and comprehensive metrics
        - 'group_means': Group means for each metabolite
        - 'group_counts': Sample counts per group
        - 'enhanced_metabolites': Original metabolite data enhanced with statistical results
    """
    results: Dict[str, pd.DataFrame] = {}
    stats_timer_start = time.time()
    
    if not sample_cols:
        return results

    # Validate group map
    missing = [c for c in sample_cols if c not in group_map]
    if missing:
        raise ValueError(f"Group mapping missing for columns: {missing}")

    # Normalize minimum group size requirement
    try:
        min_group_size = max(1, int(min_group_size))
    except Exception:
        min_group_size = 1

    # Determine group ordering:
    # 1. If group_order provided, keep only those present (preserving given order).
    # 2. Otherwise preserve first-seen order from sample_cols.
    if group_order:
        present = []
        seen_set = set()
        for g in group_order:
            if g in group_map.values() and g not in seen_set:
                present.append(g)
                seen_set.add(g)
        groups = present
    else:
        raw_groups = [group_map[c] for c in sample_cols]
        groups = []
        seen_local = set()
        for g in raw_groups:
            if g not in seen_local:
                seen_local.add(g)
                groups.append(g)
    # Build column lists per group (pre-filter)
    group_cols_all: Dict[str, List[str]] = {g: [c for c in sample_cols if group_map[c] == g] for g in groups}

    # Create group counts summary and filter by minimum size
    group_counts_data = []
    included_groups: List[str] = []
    group_cols: Dict[str, List[str]] = {}
    for g in groups:
        cols = group_cols_all[g]
        count = len(cols)
        include = count >= min_group_size
        group_counts_data.append({
            'group': g,
            'sample_count': count,
            'sample_columns': ', '.join(cols),
            'included_in_analysis': include,
            'min_required': min_group_size
        })
        if include:
            included_groups.append(g)
            group_cols[g] = cols
    results['group_counts'] = pd.DataFrame(group_counts_data)
    groups = included_groups
    sample_cols = [c for c in sample_cols if group_map[c] in groups]
    k = len(groups)
    group_sizes = {g: len(group_cols[g]) for g in groups}

    # Determine metabolite identifier
    metabolite_id_col = None
    # First try standard candidates including Protein
    for candidate in ['Name', 'Protein', 'Metabolite', 'Molecule', 'Compound', 'LipidID', 'Lipid_ID', 'Lipid_Class']:
        if candidate in df.columns:
            metabolite_id_col = candidate
            break
            
    # If no standard candidate found, use first feature column
    if metabolite_id_col is None and len(df.columns) > 0:
        metabolite_id_col = df.columns[0]  # Fall back to first column

    overall_rows: List[Dict] = []
    pairwise_rows: List[Dict] = []
    dunn_rows: List[Dict] = []
    tukey_rows: List[Dict] = []
    enhanced_metabolite_rows: List[Dict] = []

    # Decide test usage
    use_pairwise_only_two = (k == 2 and pairwise_test in ('welch','mannwhitney','rots','limma'))

    use_dunn = (overall_test == 'kruskal')  # If user chose Kruskal, provide Dunn pairwise comparisons
    use_tukey = (overall_test == 'anova')   # Provide Tukey HSD when ANOVA selected

    # Normalize base group usage
    if base_group and base_group not in groups:
        base_group = None  # ignore invalid

    # Pre-compute desired pair ordering (triangular in group order) for consistent output
    desired_pairs: List[tuple] = []
    if base_group:
        # If restricted to base group comparisons, put those only in defined group order
        desired_pairs = [(base_group, g) for g in groups if g != base_group]
    else:
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                desired_pairs.append((groups[i], groups[j]))
    order_index_map = {p: idx for idx, p in enumerate(desired_pairs)}

    # Enable parallel processing whenever user requests multiple jobs; ROTS is now thread-safe
    can_parallel = bool(n_jobs and n_jobs > 1)

    # ===== BATCH ROTS PRE-COMPUTATION =====
    # If ROTS is selected, pre-compute all statistics using efficient vectorized batch processing
    # This replaces the slow per-metabolite ROTS computation with a fast matrix-based approach
    rots_precomputed: Dict[Tuple[str, str, Any], Tuple[float, float, float]] = {}  # (g1, g2, metabolite_id) -> (stat, pval, s0)
    
    if pairwise_test == 'rots' and ((k > 2 and not use_dunn and not use_tukey) or use_pairwise_only_two):
        print(f"[STATS_TIMER] Using batch ROTS computation for {len(df)} metabolites...")
        batch_rots_start = time.time()
        
        # Build data matrix from dataframe - shape (n_metabolites, n_samples)
        sample_cols_array = np.array(sample_cols)
        data_matrix = df[sample_cols].values.astype(float)
        
        # Get metabolite IDs for lookup
        if metabolite_id_col and metabolite_id_col in df.columns:
            metabolite_ids = df[metabolite_id_col].values
        else:
            metabolite_ids = np.arange(len(df))
        
        # Determine comparisons to run
        if base_group:
            batch_comparisons = [(base_group, g) for g in groups if g != base_group]
        elif custom_comparisons and len(custom_comparisons) > 0:
            batch_comparisons = [(g1, g2) for g1, g2 in custom_comparisons if g1 in groups and g2 in groups]
        else:
            batch_comparisons = [(groups[i], groups[j]) for i in range(k) for j in range(i+1, k)]
        
        # Run batch ROTS for each comparison pair
        for comp_idx, (g1, g2) in enumerate(batch_comparisons):
            # Get column indices for each group
            g1_cols = group_cols[g1]
            g2_cols = group_cols[g2]
            g1_indices = np.array([sample_cols.index(c) for c in g1_cols])
            g2_indices = np.array([sample_cols.index(c) for c in g2_cols])
            
            if progress_callback:
                progress_callback(0, 100, f"ROTS batch: {g1} vs {g2} ({comp_idx+1}/{len(batch_comparisons)})")
            
            # Run vectorized ROTS on all metabolites at once
            try:
                stats_arr, pvals_arr, s0, _ = _rots_batch(
                    data_matrix, g1_indices, g2_indices,
                    B=rots_B, K=rots_K, seed=rots_seed,
                    progress_callback=progress_callback
                )
                
                # Store results in lookup dictionary
                for i, met_id in enumerate(metabolite_ids):
                    if not np.isnan(stats_arr[i]):
                        rots_precomputed[(g1, g2, met_id)] = (stats_arr[i], pvals_arr[i], s0)
                
            except Exception as e:
                print(f"[STATS_TIMER] Batch ROTS failed for {g1} vs {g2}: {e}, falling back to per-metabolite")
        
        batch_rots_elapsed = time.time() - batch_rots_start
        print(f"[STATS_TIMER] Batch ROTS completed in {batch_rots_elapsed:.1f}s for {len(rots_precomputed)} valid comparisons")

    # ===== BATCH LIMMA PRE-COMPUTATION =====
    # If limma is selected, pre-compute all statistics using efficient vectorized batch processing
    limma_precomputed: Dict[Tuple[str, str, Any], Tuple[float, float, float, float]] = {}  # (g1, g2, metabolite_id) -> (stat, pval, d0, s0_sq)
    
    if pairwise_test == 'limma' and ((k > 2 and not use_dunn and not use_tukey) or use_pairwise_only_two):
        print(f"[STATS_TIMER] Using batch limma computation for {len(df)} metabolites...")
        batch_limma_start = time.time()
        
        # Build data matrix from dataframe - shape (n_metabolites, n_samples)
        sample_cols_array = np.array(sample_cols)
        data_matrix = df[sample_cols].values.astype(float)
        
        # Get metabolite IDs for lookup
        if metabolite_id_col and metabolite_id_col in df.columns:
            metabolite_ids = df[metabolite_id_col].values
        else:
            metabolite_ids = np.arange(len(df))
        
        # Determine comparisons to run
        if base_group:
            batch_comparisons = [(base_group, g) for g in groups if g != base_group]
        elif custom_comparisons and len(custom_comparisons) > 0:
            batch_comparisons = [(g1, g2) for g1, g2 in custom_comparisons if g1 in groups and g2 in groups]
        else:
            batch_comparisons = [(groups[i], groups[j]) for i in range(k) for j in range(i+1, k)]
        
        # Run batch limma for each comparison pair
        for comp_idx, (g1, g2) in enumerate(batch_comparisons):
            # Get column indices for each group
            g1_cols = group_cols[g1]
            g2_cols = group_cols[g2]
            g1_indices = np.array([sample_cols.index(c) for c in g1_cols])
            g2_indices = np.array([sample_cols.index(c) for c in g2_cols])
            
            if progress_callback:
                progress_callback(0, 100, f"Limma batch: {g1} vs {g2} ({comp_idx+1}/{len(batch_comparisons)})")
            
            # Run vectorized limma on all metabolites at once
            try:
                stats_arr, pvals_arr, d0, s0_sq = _limma_batch(
                    data_matrix, g1_indices, g2_indices,
                    progress_callback=progress_callback
                )
                
                # Store results in lookup dictionary
                for i, met_id in enumerate(metabolite_ids):
                    if not np.isnan(stats_arr[i]):
                        limma_precomputed[(g1, g2, met_id)] = (stats_arr[i], pvals_arr[i], d0, s0_sq)
                
            except Exception as e:
                print(f"[STATS_TIMER] Batch limma failed for {g1} vs {g2}: {e}, falling back to per-metabolite")
        
        batch_limma_elapsed = time.time() - batch_limma_start
        print(f"[STATS_TIMER] Batch limma completed in {batch_limma_elapsed:.1f}s for {len(limma_precomputed)} valid comparisons")

    def _process_one_row(item):
        idx, row = item
        local_overall_rows = []
        local_pairwise_rows = []
        local_dunn_rows = []
        local_tukey_rows = []
        local_enhanced_entry = None
        row_rng = None
        if pairwise_test == 'rots' and rots_seed is None:
            row_rng = np.random.default_rng()
        
        # Timing tracking for this metabolite
        test_times = {'welch': 0, 'mannwhitney': 0, 'rots': 0, 'limma': 0, 'overall': 0}

        group_values: Dict[str, np.ndarray] = {}
        group_means: Dict[str, float] = {}
        group_total_samples: Dict[str, int] = {}
        valid_for_overall = True
        for g, cols in group_cols.items():
            vals = row[cols].astype(float).values
            group_total_samples[g] = len(cols)
            vals = vals[~np.isnan(vals) & (vals != 0)]
            if len(vals) == 0:
                valid_for_overall = False
                group_means[g] = np.nan
            else:
                group_means[g] = float(vals.mean())
            group_values[g] = vals

        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        enhanced_entry = row.to_dict()
        enhanced_entry['metabolite_id'] = metabolite_id
        for g in groups:
            enhanced_entry[f'mean_{g}'] = group_means[g]
            enhanced_entry[f'n_{g}'] = len(group_values[g])
        local_enhanced_entry = enhanced_entry

        overall_stat, overall_p = np.nan, np.nan
        has_insufficient_samples = False
        sufficient_groups = 0
        if min_group_size_percent is not None:
            for g in groups:
                total_samples = group_total_samples[g]
                non_zero_count = len(group_values[g])
                required_count = int(np.ceil(total_samples * min_group_size_percent / 100.0))
                if non_zero_count >= required_count:
                    sufficient_groups += 1
        else:
            for g in groups:
                non_zero_count = len(group_values[g])
                if non_zero_count >= min_group_size:
                    sufficient_groups += 1
        # Require at least 2 groups with sufficient data for overall test
        if sufficient_groups < 2:
            has_insufficient_samples = True
            has_insufficient_samples = any(len(group_values[g]) < min_group_size for g in groups)

        if k > 2 and overall_test in ('anova','kruskal') and valid_for_overall:
            if has_insufficient_samples:
                overall_stat, overall_p = np.nan, np.nan
            else:
                try:
                    t_start = time.time()
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=stats.ConstantInputWarning)
                        warnings.filterwarnings('ignore', category=RuntimeWarning)
                        if overall_test == 'anova':
                            overall_stat, overall_p = stats.f_oneway(*[group_values[g] for g in groups])
                        else:
                            overall_stat, overall_p = stats.kruskal(*[group_values[g] for g in groups])
                    test_times['overall'] += time.time() - t_start
                    overall_stat, overall_p = float(overall_stat), float(overall_p)
                except Exception:
                    overall_stat, overall_p = np.nan, np.nan
        elif k == 2 and overall_test in ('anova','kruskal'):
            g1, g2 = groups
            if has_insufficient_samples:
                overall_stat, overall_p = np.nan, np.nan
            else:
                try:
                    t_start = time.time()
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=stats.ConstantInputWarning)
                        warnings.filterwarnings('ignore', category=RuntimeWarning)
                        if overall_test == 'anova':
                            overall_stat, overall_p = stats.f_oneway(group_values[g1], group_values[g2])
                        else:
                            overall_stat, overall_p = stats.kruskal(group_values[g1], group_values[g2])
                    test_times['overall'] += time.time() - t_start
                    overall_stat, overall_p = float(overall_stat), float(overall_p)
                except Exception:
                    overall_stat, overall_p = np.nan, np.nan

        if not np.isnan(overall_p):
            overall_p = _floor_pval(overall_p)
            local_overall_rows.append({id_column_name: metabolite_id, 'statistic': overall_stat, 'p_value': overall_p})

        overall_sig = (not np.isnan(overall_p)) and (overall_p <= alpha)
        allow_posthoc = not has_insufficient_samples

        if use_dunn and k >= 2 and allow_posthoc:
            pooled_vals = []
            group_labels = []
            for g in groups:
                vals = group_values[g]
                pooled_vals.extend(vals.tolist())
                group_labels.extend([g] * len(vals))
            pooled_vals = np.array(pooled_vals, dtype=float)
            N = len(pooled_vals)
            if N > 0:
                ranks = stats.rankdata(pooled_vals, method='average')
                avg_ranks = {}
                for g in groups:
                    mask = np.array(group_labels) == g
                    if mask.sum():
                        avg_ranks[g] = ranks[mask].mean()
                    else:
                        avg_ranks[g] = np.nan
                unique_vals, counts = np.unique(pooled_vals, return_counts=True)
                tie_term = 0.0
                if N > 1:
                    tie_term = (np.sum(counts**3 - counts)) / (12.0 * N * (N - 1))
                base = N * (N + 1) / 12.0
                if base_group:
                    comp_iter = [(base_group, g) for g in groups if g != base_group]
                elif custom_comparisons and len(custom_comparisons) > 0:
                    comp_iter = [(g1, g2) for g1, g2 in custom_comparisons if g1 in groups and g2 in groups]
                else:
                    # Generate all pairwise combinations
                    comp_iter = []
                    for i in range(k):
                        for j in range(i + 1, k):
                            comp_iter.append((groups[i], groups[j]))
                for g1, g2 in comp_iter:
                    v1, v2 = group_values[g1], group_values[g2]
                    mean1 = group_means[g1]
                    mean2 = group_means[g2]
                    if not np.isnan(mean1) and not np.isnan(mean2) and mean1 != 0:
                        fc = mean2 / mean1
                        log2_fc = np.log2(fc) if fc > 0 else np.nan
                    else:
                        fc = np.nan
                        log2_fc = np.nan
                    # For Dunn post-hoc: Metabolites already passed overall Kruskal filter.
                    # Skip comparisons where groups have insufficient samples.
                    threshold1 = int(np.ceil(group_sizes[g1] * min_group_size_percent / 100.0)) if min_group_size_percent else min_group_size
                    threshold2 = int(np.ceil(group_sizes[g2] * min_group_size_percent / 100.0)) if min_group_size_percent else min_group_size
                    
                    # Skip if insufficient data - don't add to results at all
                    if len(v1) == 0 and len(v2) == 0:
                        continue
                    if len(v1) < threshold1 or len(v2) < threshold2:
                        continue
                    
                    var_ij = (base - tie_term) * (1.0 / len(v1) + 1.0 / len(v2))
                    if var_ij <= 0:
                        continue
                    
                    z_stat = (avg_ranks[g1] - avg_ranks[g2]) / np.sqrt(var_ij)
                    p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
                    cohen_d, cliffs = _effect_sizes(v1, v2)
                    eff, se_eff, ci_low, ci_high = _rank_biserial_effect_ci(v1, v2)
                    local_dunn_rows.append({
                        id_column_name: metabolite_id,
                        'group1': g1,
                        'group2': g2,
                        'mean_group1': mean1,
                        'mean_group2': mean2,
                        'fold_change': fc,
                        'log2_fold_change': log2_fc,
                        'dunn_p_value': _floor_pval(p_val),
                        'model_effect': eff,
                        'model_se': se_eff,
                        'ci_lower_95': ci_low,
                        'ci_upper_95': ci_high,
                        'cohen_d': cohen_d,
                        'cliffs_delta': cliffs,
                        'n_group1': len(v1),
                        'n_group2': len(v2),
                        'order_idx': order_index_map.get((g1, g2), 999999)
                    })
        elif ((k > 2 and pairwise_test in ('welch','mannwhitney','rots','limma') and allow_posthoc)
              or use_pairwise_only_two):
            if base_group:
                comp_iter = [(base_group, g) for g in groups if g != base_group]
            elif custom_comparisons and len(custom_comparisons) > 0:
                comp_iter = [(g1, g2) for g1, g2 in custom_comparisons if g1 in groups and g2 in groups]
            else:
                # Generate all pairwise combinations
                comp_iter = [(groups[i], groups[j]) for i in range(k) for j in range(i+1, k)]
            for g1, g2 in comp_iter:
                v1, v2 = group_values[g1], group_values[g2]
                mean1 = group_means[g1]
                mean2 = group_means[g2]
                if not np.isnan(mean1) and not np.isnan(mean2) and mean1 != 0:
                    fc = mean2 / mean1
                    log2_fc = np.log2(fc) if fc > 0 else np.nan
                else:
                    fc = np.nan
                    log2_fc = np.nan
                insufficient = False
                rots_s = np.nan
                if min_group_size_percent is not None:
                    for g in [g1, g2]:
                        total_samples = group_total_samples[g]
                        non_zero_count = len(group_values[g])
                        required_count = int(np.ceil(total_samples * min_group_size_percent / 100.0))
                        if non_zero_count < required_count:
                            insufficient = True
                            break
                else:
                    insufficient = (len(v1) < min_group_size or len(v2) < min_group_size)
                if insufficient or len(v1) == 0 or len(v2) == 0:
                    # Skip this comparison entirely - don't add to results
                    continue
                else:
                    try:
                        eff = np.nan
                        se_eff = np.nan
                        ci_low = np.nan
                        ci_high = np.nan
                        if pairwise_test == 'welch':
                            t_start = time.time()
                            tt = stats.ttest_ind(v1, v2, equal_var=False)
                            stat, p = float(tt.statistic), float(tt.pvalue)
                            eff, se_eff, ci_low, ci_high = _welch_effect_ci(v1, v2)
                            test_times['welch'] += time.time() - t_start
                        elif pairwise_test == 'mannwhitney':
                            t_start = time.time()
                            mw_result = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                            stat, p = float(mw_result.statistic), float(mw_result.pvalue)
                            eff, se_eff, ci_low, ci_high = _rank_biserial_effect_ci(v1, v2)
                            test_times['mannwhitney'] += time.time() - t_start
                        elif pairwise_test == 'rots':
                            t_start = time.time()
                            # Check if we have precomputed batch results
                            precomputed_key = (g1, g2, metabolite_id)
                            if precomputed_key in rots_precomputed:
                                stat, p, rots_s = rots_precomputed[precomputed_key]
                            else:
                                # Fallback to single-metabolite computation (should rarely happen)
                                rots_stat, rots_pval, rots_s0, _ = _rots_test(
                                    v1, v2, B=rots_B, K=rots_K, alpha=rots_alpha, seed=rots_seed
                                )
                                stat, p = rots_stat, rots_pval
                                rots_s = rots_s0
                            # Report interpretable pairwise effect/CI alongside ROTS p-values.
                            eff, se_eff, ci_low, ci_high = _welch_effect_ci(v1, v2)
                            test_times['rots'] += time.time() - t_start
                        elif pairwise_test == 'limma':
                            t_start = time.time()
                            # Check if we have precomputed batch results
                            precomputed_key = (g1, g2, metabolite_id)
                            if precomputed_key in limma_precomputed:
                                stat, p, limma_d0, limma_s0 = limma_precomputed[precomputed_key]
                            else:
                                # Fallback to single-metabolite computation (no prior info)
                                stat, p, limma_d0, limma_s0 = _limma_test(v1, v2)
                            test_times['limma'] += time.time() - t_start
                        else:
                            # Unknown test - skip
                            continue
                    except Exception:
                        # Test failed - skip this comparison
                        continue
                
                cohen_d, cliffs = _effect_sizes(v1, v2)
                result_dict = {
                    id_column_name: metabolite_id,
                    'group1': g1,
                    'group2': g2,
                    'mean_group1': mean1,
                    'mean_group2': mean2,
                    'fold_change': fc,
                    'log2_fold_change': log2_fc,
                    'statistic': stat,
                    'p_value': _floor_pval(p),
                    'model_effect': eff,
                    'model_se': se_eff,
                    'ci_lower_95': ci_low,
                    'ci_upper_95': ci_high,
                    'cohen_d': cohen_d,
                    'cliffs_delta': cliffs,
                    'n_group1': len(v1),
                    'n_group2': len(v2),
                    'order_idx': order_index_map.get((g1, g2), 999999)
                }
                if pairwise_test == 'rots':
                    result_dict['rots_s0'] = rots_s
                local_pairwise_rows.append(result_dict)
        
        # Only run Tukey HSD when ANOVA overall test is selected
        if use_tukey and k > 2 and allow_posthoc:
            try:
                vals_all = []
                labels_all = []
                for g in groups:
                    vals = group_values[g]
                    if len(vals):
                        vals_all.extend(vals.tolist())
                        labels_all.extend([g] * len(vals))
                
                if len(vals_all) and len(set(labels_all)) >= 2 and _sm_mc is not None:
                    # Run Tukey HSD (only when ANOVA selected and allow_posthoc is True)
                    mc = _sm_mc.MultiComparison(vals_all, labels_all)
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
                        tuk = mc.tukeyhsd(alpha=alpha)
                    summary_df = pd.DataFrame(data=tuk._results_table.data[1:], columns=tuk._results_table.data[0])
                    for _, r in summary_df.iterrows():
                        g1 = str(r['group1']); g2 = str(r['group2'])
                        if base_group and base_group not in (g1, g2):
                            continue
                        v1 = group_values.get(g1, np.array([])); v2 = group_values.get(g2, np.array([]))
                        
                        # Skip if insufficient data - don't add to results at all
                        if len(v1) < min_group_size or len(v2) < min_group_size:
                            continue
                        
                        mean1 = group_means.get(g1, np.nan); mean2 = group_means.get(g2, np.nan)
                        if not np.isnan(mean1) and not np.isnan(mean2) and mean1 != 0:
                            fc = mean2 / mean1; log2_fc = np.log2(fc) if fc > 0 else np.nan
                        else:
                            fc = np.nan; log2_fc = np.nan
                        cohen_d, cliffs = _effect_sizes(v1, v2)
                        p_value_adj = _floor_pval(r.get('p-adj', np.nan))
                        reject = r.get('reject', False)
                        # Tukey HSD already provides 95% CI bounds in the result table.
                        ci_low = r.get('lower', np.nan)
                        ci_high = r.get('upper', np.nan)
                        local_tukey_rows.append({
                            id_column_name: metabolite_id,
                            'group1': g1,
                            'group2': g2,
                            'mean_group1': mean1,
                            'mean_group2': mean2,
                            'fold_change': fc,
                            'log2_fold_change': log2_fc,
                            'meandiff': r.get('meandiff', np.nan),
                            'p_value_adj': p_value_adj,
                            'model_effect': r.get('meandiff', np.nan),
                            'model_se': np.nan,
                            'ci_lower_95': ci_low,
                            'ci_upper_95': ci_high,
                            'reject': reject,
                            'cohen_d': cohen_d,
                            'cliffs_delta': cliffs,
                            'n_group1': len(v1),
                            'n_group2': len(v2),
                            'order_idx': order_index_map.get((g1, g2), 999999)
                        })
            except Exception:
                pass

        return (local_overall_rows, local_pairwise_rows, local_dunn_rows, local_tukey_rows, local_enhanced_entry)

    total_metabolites = len(df)
    callback_frequency_stats = max(1, total_metabolites // 100)  # Call every 1%
    
    processing_start = time.time()
    print(f"[STATS_TIMER] Starting statistical analysis for {total_metabolites} metabolites...")
    print(f"[STATS_TIMER] Configuration: overall_test={overall_test}, pairwise_test={pairwise_test}, n_jobs={n_jobs}")
    
    if can_parallel:
        completed = 0
        for local in ThreadPoolExecutor(max_workers=max(1, int(n_jobs))).map(_process_one_row, df.iterrows()):
            o, p, d, t, e = local
            if o: overall_rows.extend(o)
            if p: pairwise_rows.extend(p)
            if d: dunn_rows.extend(d)
            if t: tukey_rows.extend(t)
            if e is not None: enhanced_metabolite_rows.append(e)
            completed += 1
            if progress_callback and (completed % callback_frequency_stats == 0 or completed == total_metabolites):
                try:
                    metab_name = e.get('metabolite_id', 'Unknown') if e else 'Unknown'
                    elapsed = time.time() - processing_start
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta_remaining = (total_metabolites - completed) / rate if rate > 0 else 0
                    print(f"[STATS_TIMER] Progress: {completed}/{total_metabolites} ({100*completed/total_metabolites:.1f}%) | Elapsed: {elapsed:.1f}s | Rate: {rate:.1f} met/s | ETA: {eta_remaining:.1f}s")
                    progress_callback(completed, total_metabolites, str(metab_name))
                except Exception:
                    pass
    else:
        completed = 0
        for item in df.iterrows():
            o, p, d, t, e = _process_one_row(item)
            if o:
                overall_rows.extend(o)
            if p:
                pairwise_rows.extend(p)
            if d:
                dunn_rows.extend(d)
            if t:
                tukey_rows.extend(t)
            if e is not None:
                enhanced_metabolite_rows.append(e)
            completed += 1
            if progress_callback and (completed % callback_frequency_stats == 0 or completed == total_metabolites):
                try:
                    metab_name = e.get('metabolite_id', 'Unknown') if e else 'Unknown'
                    elapsed = time.time() - processing_start
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta_remaining = (total_metabolites - completed) / rate if rate > 0 else 0
                    print(f"[STATS_TIMER] Progress: {completed}/{total_metabolites} ({100*completed/total_metabolites:.1f}%) | Elapsed: {elapsed:.1f}s | Rate: {rate:.1f} met/s | ETA: {eta_remaining:.1f}s")
                    progress_callback(completed, total_metabolites, str(metab_name))
                except Exception:
                    pass

    # Post-processing phase with timing
    postproc_start = time.time()
    print(f"[STATS_TIMER] Post-processing results (FDR correction, visualization metrics)...")
    
    # Process overall results (ROTS finalization code removed - now computed directly)
    if overall_rows:
        overall_df = pd.DataFrame(overall_rows)
        if fdr and 'p_value' in overall_df.columns:
            p_vals = [float(_floor_pval(p) if p is not None else 1.0) for p in overall_df['p_value'].tolist()]
            overall_df['p_value_adj'] = _bh_adjust(p_vals)
        overall_df = _apply_dynamic_pvalue_floor_frame(overall_df, ['p_value', 'p_value_adj'])
        # Add visualization metrics using vectorized operations
        if 'p_value' in overall_df.columns:
            overall_df['neg_log10_p'] = -np.log10(overall_df['p_value'].replace(0, np.finfo(float).eps))
        if 'p_value_adj' in overall_df.columns:
            overall_df['neg_log10_p_adj'] = -np.log10(overall_df['p_value_adj'].replace(0, np.finfo(float).eps))
            
        results['overall'] = overall_df

    # Process pairwise results
    # Apply FDR across ALL comparisons globally if requested (collect raw p-values first)
    global_p_map = None
    if fdr and fdr_scope == 'global':
        raw_records = []
        if dunn_rows:
            raw_records.extend([(r['group1'], r['group2'], r[id_column_name], r['dunn_p_value']) for r in dunn_rows])
        if pairwise_rows:
            raw_records.extend([(r['group1'], r['group2'], r[id_column_name], r['p_value']) for r in pairwise_rows])
        if raw_records:
            vals = [p if p is not None and not np.isnan(p) else 1.0 for (_, _, _, p) in raw_records]
            adj_all = _bh_adjust(vals)
            global_p_map = {}
            for (rec, adj) in zip(raw_records, adj_all):
                g1, g2, metab, _ = rec
                global_p_map[(g1, g2, metab)] = adj

    # Prefer Tukey rows if available for ANOVA; they already carry adjusted p-values
    if 'use_tukey' in locals() and tukey_rows:
        tuk_df = pd.DataFrame(tukey_rows)
        if 'order_idx' in tuk_df.columns:
            tuk_df.sort_values('order_idx', inplace=True)
    # Expression classification based on adjusted p-value (already adjusted)
        tuk_df['Expression'] = 'Not Significant'
        mask_sig = (tuk_df['p_value_adj'] <= alpha) & tuk_df['log2_fold_change'].notna()
        tuk_df.loc[mask_sig & (tuk_df['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
        tuk_df.loc[mask_sig & (tuk_df['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        tuk_df.rename(columns={'p_value_adj':'tukey_p_adj'}, inplace=True)
        tuk_df = _apply_dynamic_pvalue_floor_frame(tuk_df, ['tukey_p_adj'])
        tuk_df['neg_log10_adj_p'] = -np.log10(tuk_df['tukey_p_adj'].replace(0, np.finfo(float).eps))
        tuk_df['p_value_adj'] = tuk_df['tukey_p_adj']
        # Provide a consistent metabolite_id column for downstream merging/exports
        if id_column_name in tuk_df.columns and 'metabolite_id' not in tuk_df.columns:
            tuk_df['metabolite_id'] = tuk_df[id_column_name]
        results['pairwise'] = tuk_df
    elif dunn_rows:
        dunn_df = pd.DataFrame(dunn_rows)
        if 'order_idx' in dunn_df.columns:
            dunn_df.sort_values('order_idx', inplace=True)
        # BH adjust based on selected scope - ALWAYS compute adjusted p-values
        if 'dunn_p_value' in dunn_df.columns:
            # Default to per-comparison if fdr not explicitly set
            effective_fdr = fdr if fdr is not None else True
            effective_scope = fdr_scope if effective_fdr else 'per-comparison'
            
            if effective_fdr:
                if effective_scope == 'per-comparison':
                    # Adjust within each pair (group1, group2) across metabolites
                    dunn_df['dunn_p_adj'] = np.nan
                    for key, sub_idx in dunn_df.groupby(['group1', 'group2']).groups.items():
                        g1, g2 = key
                        sub = dunn_df.loc[sub_idx]
                        adj = _bh_adjust([float(_floor_pval(p) if p is not None else 1.0) for p in sub['dunn_p_value'].tolist()])
                        dunn_df.loc[sub_idx, 'dunn_p_adj'] = adj
                elif effective_scope == 'per-metabolite':
                    # Adjust within each metabolite across all pairwise comparisons (R-like approach)
                    dunn_df['dunn_p_adj'] = np.nan
                    for metabolite, sub_idx in dunn_df.groupby(id_column_name).groups.items():
                        sub = dunn_df.loc[sub_idx]
                        adj = _bh_adjust([float(_floor_pval(p) if p is not None else 1.0) for p in sub['dunn_p_value'].tolist()])
                        dunn_df.loc[sub_idx, 'dunn_p_adj'] = adj
                else:  # global
                    dunn_df['dunn_p_adj'] = [(global_p_map or {}).get((g1, g2, m), np.nan) for g1, g2, m in zip(dunn_df['group1'], dunn_df['group2'], dunn_df[id_column_name])]
            else:
                # Even if fdr=False, still provide adjusted p as copy for consistency
                dunn_df['dunn_p_adj'] = dunn_df['dunn_p_value']
        else:
            dunn_df['dunn_p_adj'] = dunn_df.get('dunn_p_value')
        # Dynamic per-column floor
        dunn_df = _apply_dynamic_pvalue_floor_frame(dunn_df, ['dunn_p_value', 'dunn_p_adj'])
        # Compute -log10 of adjusted only
        dunn_df['neg_log10_dunn_p_adj'] = -np.log10(dunn_df['dunn_p_adj'].replace(0, np.finfo(float).eps))
        # Expression classification
        dunn_df['Expression'] = 'Not Significant'
        mask_sig = (dunn_df['dunn_p_adj'] <= alpha) & dunn_df['log2_fold_change'].notna()
        dunn_df.loc[mask_sig & (dunn_df['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
        dunn_df.loc[mask_sig & (dunn_df['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        # Drop raw internal columns user does not want
        drop_cols = [c for c in ['dunn_z', 'dunn_p_value', 'dunn_p_value_adj', 'neg_log10_dunn_p', 'neg_log10_dunn_p_adj'] if c in dunn_df.columns]
        if drop_cols:
            dunn_df.drop(columns=drop_cols, inplace=True)
        # For compatibility still provide generic names expected by exporter if needed
        dunn_df['p_value_adj'] = dunn_df['dunn_p_adj']
        # Provide a consistent metabolite_id column for downstream merging/exports
        if id_column_name in dunn_df.columns and 'metabolite_id' not in dunn_df.columns:
            dunn_df['metabolite_id'] = dunn_df[id_column_name]
        results['pairwise'] = dunn_df
    elif pairwise_rows:
        pairwise_df = pd.DataFrame(pairwise_rows)
        if 'order_idx' in pairwise_df.columns:
            pairwise_df.sort_values('order_idx', inplace=True)
        # ALWAYS compute adjusted p-values for pairwise comparisons, using per-comparison scope by default
        if 'p_value' in pairwise_df.columns:
            # Default to per-comparison if fdr not explicitly set
            effective_fdr = fdr if fdr is not None else True
            effective_scope = fdr_scope if effective_fdr else 'per-comparison'
            
            if effective_fdr:
                if effective_scope == 'per-comparison':
                    # Optimize: use apply with vectorized operations
                    pairwise_df['p_value'] = pairwise_df['p_value'].apply(lambda x: float(_floor_pval(x)) if pd.notna(x) else 1.0)
                    pairwise_df['p_value_adj'] = np.nan
                    for key, sub_idx in pairwise_df.groupby(['group1', 'group2'], sort=False).groups.items():
                        adj = _bh_adjust(pairwise_df.loc[sub_idx, 'p_value'].tolist())
                        pairwise_df.loc[sub_idx, 'p_value_adj'] = adj
                elif effective_scope == 'per-metabolite':
                    # Adjust within each metabolite across its pairwise comparisons (R-like)
                    pairwise_df['p_value'] = pairwise_df['p_value'].apply(lambda x: float(_floor_pval(x)) if pd.notna(x) else 1.0)
                    pairwise_df['p_value_adj'] = np.nan
                    for metabolite, sub_idx in pairwise_df.groupby(id_column_name, sort=False).groups.items():
                        adj = _bh_adjust(pairwise_df.loc[sub_idx, 'p_value'].tolist())
                        pairwise_df.loc[sub_idx, 'p_value_adj'] = adj
                elif effective_scope == 'global':
                    # Apply global BH across all pairwise p-values
                    pairwise_df['p_value'] = pairwise_df['p_value'].apply(lambda x: float(_floor_pval(x)) if pd.notna(x) else 1.0)
                    adj_all = _bh_adjust(pairwise_df['p_value'].tolist())
                    pairwise_df['p_value_adj'] = adj_all
            else:
                # Even if fdr=False, still provide p_value_adj as copy of p_value for consistency
                pairwise_df['p_value'] = pairwise_df['p_value'].apply(lambda x: float(_floor_pval(x)) if pd.notna(x) else 1.0)
                pairwise_df['p_value_adj'] = pairwise_df['p_value']

        pairwise_df = _apply_dynamic_pvalue_floor_frame(pairwise_df, ['p_value', 'p_value_adj'])
        
        # Vectorized operations instead of .apply() for speed
        if 'p_value' in pairwise_df.columns:
            pairwise_df['neg_log10_p'] = -np.log10(pairwise_df['p_value'].replace(0, np.finfo(float).eps))
        if 'p_value_adj' in pairwise_df.columns:
            pairwise_df['neg_log10_p_adj'] = -np.log10(pairwise_df['p_value_adj'].replace(0, np.finfo(float).eps))
        
        # Expression classification: prefer adjusted p-values if present, otherwise raw p-values
        if 'log2_fold_change' in pairwise_df.columns:
            pairwise_df['Expression'] = 'Not Significant'
            if 'p_value_adj' in pairwise_df.columns:
                mask_sig = (pairwise_df['p_value_adj'] <= alpha) & pairwise_df['log2_fold_change'].notna()
            else:
                mask_sig = (pairwise_df['p_value'] <= alpha) & pairwise_df['log2_fold_change'].notna()
            pairwise_df.loc[mask_sig & (pairwise_df['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
            pairwise_df.loc[mask_sig & (pairwise_df['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        # Ensure a metabolite_id column exists for exporters to merge feature metadata reliably
        if id_column_name in pairwise_df.columns and 'metabolite_id' not in pairwise_df.columns:
            pairwise_df['metabolite_id'] = pairwise_df[id_column_name]
        results['pairwise'] = pairwise_df

    # Enhanced metabolite data (original data + group means + statistics)
    if enhanced_metabolite_rows:
        enhanced_df = pd.DataFrame(enhanced_metabolite_rows)
        
        # Ensure metabolite_id column exists for downstream merging/exports
        if 'metabolite_id' not in enhanced_df.columns:
            # Fallback: try to get from a known identifier column, then first available feature column
            if 'metabolite' in enhanced_df.columns:
                enhanced_df['metabolite_id'] = enhanced_df['metabolite']
            else:
                # Try known identifier columns first
                for cand in ['Name', 'Metabolite', 'Molecule', 'Compound', 'LipidID', 'Lipid_ID', 'Lipid_Class']:
                    if cand in enhanced_df.columns:
                        enhanced_df['metabolite_id'] = enhanced_df[cand]
                        break
                # If still missing, use first non-numeric column as fallback
                if 'metabolite_id' not in enhanced_df.columns:
                    for col in enhanced_df.columns:
                        if not pd.api.types.is_numeric_dtype(enhanced_df[col]):
                            enhanced_df['metabolite_id'] = enhanced_df[col]
                            break
        
        # Merge overall statistics if available
        if overall_rows:
            overall_stats = pd.DataFrame(overall_rows)[[id_column_name, 'statistic', 'p_value']].rename(columns={
                'statistic': 'overall_statistic',
                'p_value': 'overall_p_value'
            })
            # Always include -log10 raw overall p for completeness
            overall_stats['overall_neg_log10_p'] = -np.log10(overall_stats['overall_p_value'].replace(0, np.finfo(float).eps))
            if fdr:
                overall_stats['overall_p_value_adj'] = _bh_adjust([float(_floor_pval(p) if p is not None else 1.0) for p in overall_stats['overall_p_value'].tolist()])
                overall_stats['overall_neg_log10_p_adj'] = -np.log10(overall_stats['overall_p_value_adj'].replace(0, np.finfo(float).eps))

            overall_stats = _apply_dynamic_pvalue_floor_frame(overall_stats, ['overall_p_value', 'overall_p_value_adj'])
            overall_stats['overall_neg_log10_p'] = -np.log10(overall_stats['overall_p_value'].replace(0, np.finfo(float).eps))
            if 'overall_p_value_adj' in overall_stats.columns:
                overall_stats['overall_neg_log10_p_adj'] = -np.log10(overall_stats['overall_p_value_adj'].replace(0, np.finfo(float).eps))
            
            # Merge with enhanced data
            enhanced_df = enhanced_df.merge(
                overall_stats,
                left_on='metabolite_id',
                right_on=id_column_name,
                how='left',
                suffixes=("", "_dup")
            )
            enhanced_df.drop(id_column_name, axis=1, inplace=True)

            # If earlier runs (cached) produced legacy columns, collapse any duplicates safely
            # Priority order: adjusted values if present; keep single base names.
            def _collapse(col_base):
                c_main = col_base
                c_dup1 = f"{col_base}_x"
                c_dup2 = f"{col_base}_y"
                c_dup_generic = f"{col_base}_dup"
                existing = [c for c in [c_main, c_dup1, c_dup2, c_dup_generic] if c in enhanced_df.columns]
                if len(existing) <= 1:
                    return
                # Choose first non-null series preference order existing list
                for c in existing:
                    series = enhanced_df[c]
                    if series.notna().any():
                        enhanced_df[c_main] = series
                        break
                # Drop all others except canonical name
                for c in existing:
                    if c != c_main:
                        enhanced_df.drop(columns=c, inplace=True, errors='ignore')

            for base in ['overall_statistic', 'overall_p_value', 'overall_p_value_adj', 'overall_neg_log10_p_adj']:
                _collapse(base)
        
        # Wide aggregation of pairwise metrics for Complete Results
        if 'pairwise' in results:
            # Skip if already has pairwise columns (e.g., Two-Way ANOVA)
            if any('_vs_' in col for col in enhanced_df.columns):
                pass  # Already merged
            else:
                pw = results['pairwise']
                # Handle different formats: dict of lists (ANOVA) or DataFrame (Kruskal)
                if isinstance(pw, dict):
                    # Convert dict of lists to DataFrame
                    rows = []
                    for comp, items in pw.items():
                        if isinstance(items, list):
                            for item in items:
                                rows.append({**item, 'comparison': comp})
                        elif isinstance(items, dict):
                            rows.append({**items, 'comparison': comp})
                        elif isinstance(items, pd.DataFrame):
                            for _, row in items.iterrows():
                                rows.append({**row.to_dict(), 'comparison': comp})
                        else:
                            continue
                    pw = pd.DataFrame(rows)
                # Now pw is DataFrame
                wide_frames = []
                # Choose unified adjusted p column
                if 'dunn_p_adj' in pw.columns:
                    adj_col = 'dunn_p_adj'
                    raw_col = 'dunn_p_value'
                elif 'p_value_adj' in pw.columns:
                    adj_col = 'p_value_adj'
                    raw_col = 'p_value'
                else:
                    adj_col = None
                    raw_col = 'p_value'  # fallback
            
            # Determine which p-value column to use based on use_adjusted_pvalues parameter
            # For Dunn results, always use adjusted p-values since they are always computed
            if 'dunn_p_adj' in pw.columns:
                p_col_to_use = 'dunn_p_adj'
                p_col_suffix = '_adj_p'
            elif use_adjusted_pvalues and adj_col:
                p_col_to_use = adj_col
                p_col_suffix = '_adj_p'
            else:
                p_col_to_use = raw_col if raw_col in pw.columns else adj_col
                p_col_suffix = '_pvalue'
            
            for (g1, g2), sub in pw.groupby(['group1','group2']):
                comp_label = f'{g1}_vs_{g2}'
                sub_part = sub[[id_column_name,'fold_change','log2_fold_change']].copy()
                sub_part.rename(columns={
                    'fold_change': f'{comp_label}_FC',
                    'log2_fold_change': f'{comp_label}_log2FC'
                }, inplace=True)
                # Carry model-based effect/uncertainty terms into Complete Results when available.
                for src_col, dst_suffix in [
                    ('model_effect', 'model_effect'),
                    ('model_se', 'model_se'),
                    ('ci_lower_95', 'ci_lower_95'),
                    ('ci_upper_95', 'ci_upper_95')
                ]:
                    if src_col in sub.columns:
                        sub_part[f'{comp_label}_{dst_suffix}'] = sub[src_col].values
                if p_col_to_use and p_col_to_use in sub.columns:
                    # Add the selected p-value (adjusted or raw)
                    sub_part[f'{comp_label}{p_col_suffix}'] = sub[p_col_to_use].values
                    # Also add the corresponding -log10 column for volcano/visualization needs
                    if p_col_suffix == '_adj_p':
                        neg_name = f'{comp_label}_neg_log10_adj_p'
                    else:
                        neg_name = f'{comp_label}_neg_log10_p'
                    series_vals = pd.to_numeric(sub[p_col_to_use], errors='coerce').replace(0, np.finfo(float).eps)
                    sub_part[neg_name] = -np.log10(series_vals.to_numpy())
                wide_frames.append(sub_part)
            if wide_frames:
                # Use join for better performance
                enhanced_df = enhanced_df.set_index('metabolite_id')
                for wf in wide_frames:
                    wf = wf.set_index(id_column_name)
                    enhanced_df = enhanced_df.join(wf, how='left')
                enhanced_df = enhanced_df.reset_index()
        
        # Final safeguard: ensure metabolite_id column exists before export
        if 'metabolite_id' not in enhanced_df.columns:
            # Try to recover from available identifier columns
            for cand in ['Name', 'Metabolite', 'Molecule', 'Compound', 'LipidID', 'Lipid_Class', 'HMDB_ID']:
                if cand in enhanced_df.columns:
                    enhanced_df['metabolite_id'] = enhanced_df[cand]
                    break
            # If still missing, use first non-numeric column as fallback (not index!)
            if 'metabolite_id' not in enhanced_df.columns:
                for col in enhanced_df.columns:
                    if not pd.api.types.is_numeric_dtype(enhanced_df[col]):
                        enhanced_df['metabolite_id'] = enhanced_df[col]
                        break
        
        results['enhanced_metabolites'] = enhanced_df

    # Group means for convenience (separate table)
    means_rows = []
    for idx, row in df.iterrows():
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        entry = {id_column_name: metabolite_id}
        for g, cols in group_cols.items():
            vals = row[cols].astype(float).values
            # Exclude both NaNs and zeros
            vals = vals[~np.isnan(vals) & (vals != 0)]
            entry[f'mean_{g}'] = float(vals.mean()) if len(vals) else np.nan
            entry[f'n_{g}'] = len(vals)
        means_rows.append(entry)
    results['group_means'] = pd.DataFrame(means_rows)
    
    # Add pairwise comparison summary showing actual metabolite counts tested per comparison
    if 'pairwise' in results and not results['pairwise'].empty:
        pairwise_summary = []
        pairwise_df = results['pairwise']
        summary_p_cols = [c for c in ['p_value', 'p_value_adj', 'dunn_p_adj', 'tukey_p_adj'] if c in pairwise_df.columns]
        
        # Group by comparison pairs
        for (g1, g2), sub_df in pairwise_df.groupby(['group1', 'group2'], sort=False):
            # Count valid comparisons using any available p-value column.
            # Some class export paths retain only adjusted p-value columns.
            total = len(sub_df)
            if summary_p_cols:
                tested = sub_df[summary_p_cols].notna().any(axis=1).sum()
            else:
                tested = 0
            skipped = total - tested
            
            pairwise_summary.append({
                'comparison': f'{g1}_vs_{g2}',
                'group1': g1,
                'group2': g2,
                'total_metabolites': total,
                'tested_metabolites': tested,
                'skipped_insufficient_n': skipped
            })
        
        results['pairwise_summary'] = pd.DataFrame(pairwise_summary)
    
    # Final safeguard: enforce ordering if order_idx present
    if 'pairwise' in results and 'order_idx' in results['pairwise'].columns:
        # Use the appropriate ID column name for sorting (metabolite, LipidID, or Class)
        sort_cols = ['order_idx']
        if id_column_name in results['pairwise'].columns:
            sort_cols.append(id_column_name)
        results['pairwise'] = results['pairwise'].sort_values(sort_cols).reset_index(drop=True)
    
    # Final timing report
    total_time = time.time() - stats_timer_start
    postproc_time = time.time() - postproc_start
    processing_time = time.time() - processing_start
    
    print(f"\n[STATS_TIMER] ========== STATISTICAL ANALYSIS COMPLETE ==========")
    print(f"[STATS_TIMER] Total metabolites processed: {total_metabolites}")
    print(f"[STATS_TIMER] Processing time: {processing_time:.2f}s ({100*processing_time/total_time:.1f}%)")
    print(f"[STATS_TIMER] Post-processing time: {postproc_time:.2f}s ({100*postproc_time/total_time:.1f}%)")
    print(f"[STATS_TIMER] Total time: {total_time:.2f}s")
    print(f"[STATS_TIMER] Average per metabolite: {total_time/total_metabolites*1000:.2f}ms" if total_metabolites > 0 else "[STATS_TIMER] Average per metabolite: N/A")
    print(f"[STATS_TIMER] Throughput: {total_metabolites/total_time:.2f} metabolites/s" if total_time > 0 else "[STATS_TIMER] Throughput: N/A")
    print(f"[STATS_TIMER] Overall results: {len(results.get('overall', []))} rows")
    print(f"[STATS_TIMER] Pairwise results: {len(results.get('pairwise', []))} rows")
    print(f"[STATS_TIMER] =====================================================\n")
    
    return results

def get_group_summary(group_map: Dict[str, str]) -> str:
    """Generate a formatted summary of groups with sample counts.
    
    Parameters
    ----------
    group_map : Dict[str, str]
        Mapping of sample_col -> group name
        
    Returns
    -------
    str
        Formatted string like "HFD_TBI:6, HFD:7, TBI:5, NC:4"
    """
    from collections import Counter
    group_counts = Counter(group_map.values())
    sorted_groups = sorted(group_counts.items())
    return ', '.join([f'{group}:{count}' for group, count in sorted_groups])

def perform_two_way_anova(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    factor_a_map: Dict[str, str],
    factor_b_map: Dict[str, str],
    *,
    drop_zeros: bool = True,
    min_per_cell: int = 1,
    fdr: bool = True,
    typ: int = 2,
    n_jobs: int = 3,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> pd.DataFrame:
    """Run two-way (factorial) ANOVA per metabolite with main effects and interaction.

    Parameters
    ----------
    df : DataFrame
        Input data with metabolite rows and sample columns.
    sample_cols : Sequence[str]
        Sample columns to include.
    factor_a_map : Dict[str, str]
        Mapping sample_col -> level for factor A.
    factor_b_map : Dict[str, str]
        Mapping sample_col -> level for factor B.
    drop_zeros : bool
        If True, treat zeros as missing and drop them.
    min_per_cell : int
        Minimum observations required per cell (A,B) to fit model.
    fdr : bool
        Apply Benjamini-Hochberg per-effect across metabolites.
    typ : int
        Type for ANOVA table (2 or 3). Default 2.
    n_jobs : int
        Number of parallel workers (default 3).
    progress_callback : Callable[[int, int, str], None] | None
        Optional callback(current, total, metabolite_name) for progress updates.

    Returns
    -------
    DataFrame
        Columns: metabolite, F_A, p_A, F_B, p_B, F_interaction, p_interaction,
        mean_A_<level> for each A level, mean_B_<level> for each B level,
        all sample columns with their values, and adjusted p-values if fdr=True.
    """
    # Ensure statsmodels formula API is available
    if _sm is None or _smf is None:
        raise ImportError("statsmodels is required for two-way ANOVA (install statsmodels)")

    # Validate mapping coverage
    missing_a = [c for c in sample_cols if c not in factor_a_map]
    missing_b = [c for c in sample_cols if c not in factor_b_map]
    if missing_a:
        raise ValueError(f"Factor A mapping missing for columns: {missing_a}")
    if missing_b:
        raise ValueError(f"Factor B mapping missing for columns: {missing_b}")

    # Determine metabolite identifier column (reuse logic from overall analysis)
    metabolite_id_col: Optional[str] = None
    for candidate in ['Name', 'Protein', 'Metabolite', 'Molecule', 'Compound', 'LipidID', 'Lipid_ID', 'Lipid_Class']:
        if candidate in df.columns:
            metabolite_id_col = candidate
            break
    if metabolite_id_col is None and len(df.columns) > 0:
        metabolite_id_col = df.columns[0]

    # Get unique levels for means reporting
    unique_a_levels = sorted(set(factor_a_map.values()))
    unique_b_levels = sorted(set(factor_b_map.values()))

    def _compute_two_way_entry(item):
        idx, row = item
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        values: List[float] = []
        A_levels: List[str] = []
        B_levels: List[str] = []
        sample_values: Dict[str, float] = {}
        
        for s in sample_cols:
            try:
                v = float(pd.to_numeric(row[s], errors='coerce'))
            except Exception:
                v = np.nan
            sample_values[s] = v  # Store all sample values
            
            if np.isnan(v):
                continue
            if drop_zeros and v == 0:
                continue
            values.append(v)
            A_levels.append(factor_a_map[s])
            B_levels.append(factor_b_map[s])

        # Check minimal structure
        if len(values) < 3 or len(set(A_levels)) < 2 or len(set(B_levels)) < 2:
            result = {
                'metabolite': metabolite_id,
                'F_A': np.nan, 'p_A': np.nan,
                'F_B': np.nan, 'p_B': np.nan,
                'F_interaction': np.nan, 'p_interaction': np.nan
            }
            # Add mean columns for all levels
            for al in unique_a_levels:
                result[f'mean_A_{al}'] = np.nan
            for bl in unique_b_levels:
                result[f'mean_B_{bl}'] = np.nan
            # Add sample columns
            result.update(sample_values)
            return result

        if min_per_cell > 1:
            from collections import Counter
            cell_counts = Counter(zip(A_levels, B_levels))
            if any(ct < min_per_cell for ct in cell_counts.values()):
                result = {
                    'metabolite': metabolite_id,
                    'F_A': np.nan, 'p_A': np.nan,
                    'F_B': np.nan, 'p_B': np.nan,
                    'F_interaction': np.nan, 'p_interaction': np.nan
                }
                # Add mean columns
                for al in unique_a_levels:
                    result[f'mean_A_{al}'] = np.nan
                for bl in unique_b_levels:
                    result[f'mean_B_{bl}'] = np.nan
                # Add sample columns
                result.update(sample_values)
                return result

        long_df_local = pd.DataFrame({'value': values, 'A': A_levels, 'B': B_levels})
        try:
            model = _smf.ols('value ~ C(A) + C(B) + C(A):C(B)', data=long_df_local).fit()
            anova_tbl = _sm.stats.anova_lm(model, typ=typ)
            def _get(row_key: str, col: str) -> float:
                try:
                    return float(anova_tbl.loc[row_key, col])
                except Exception:
                    return np.nan
            F_A = _get('C(A)', 'F')
            p_A = _get('C(A)', 'PR(>F)')
            F_B = _get('C(B)', 'F')
            p_B = _get('C(B)', 'PR(>F)')
            F_I = _get('C(A):C(B)', 'F')
            p_I = _get('C(A):C(B)', 'PR(>F)')
        except Exception:
            F_A = F_B = F_I = np.nan
            p_A = p_B = p_I = np.nan

        p_A = _floor_pval(p_A) if not np.isnan(p_A) else np.nan
        p_B = _floor_pval(p_B) if not np.isnan(p_B) else np.nan
        p_I = _floor_pval(p_I) if not np.isnan(p_I) else np.nan
        
        # Compute means for each level
        means_a = {}
        for al in unique_a_levels:
            mask_a = [i for i, a in enumerate(A_levels) if a == al]
            if mask_a:
                means_a[al] = float(np.mean([values[i] for i in mask_a]))
            else:
                means_a[al] = np.nan
        
        means_b = {}
        for bl in unique_b_levels:
            mask_b = [i for i, b in enumerate(B_levels) if b == bl]
            if mask_b:
                means_b[bl] = float(np.mean([values[i] for i in mask_b]))
            else:
                means_b[bl] = np.nan
        
        result = {
            'metabolite': metabolite_id,
            'F_A': F_A, 'p_A': p_A,
            'F_B': F_B, 'p_B': p_B,
            'F_interaction': F_I, 'p_interaction': p_I
        }
        # Add mean columns
        for al in unique_a_levels:
            result[f'mean_A_{al}'] = means_a.get(al, np.nan)
        for bl in unique_b_levels:
            result[f'mean_B_{bl}'] = means_b.get(bl, np.nan)
        # Add all sample columns
        result.update(sample_values)
        return result

    items = list(df.iterrows())
    total = len(items)
    callback_frequency = max(1, total // 100)  # Call every 1% of items, or at least every item if < 100 items
    
    if n_jobs and n_jobs > 1:
        with ThreadPoolExecutor(max_workers=max(1, int(n_jobs))) as ex:
            rows = []
            for i, result in enumerate(ex.map(_compute_two_way_entry, items), 1):
                rows.append(result)
                # Call callback more frequently - every 1% or every item if small dataset
                if progress_callback and (i % callback_frequency == 0 or i == total):
                    try:
                        metab_name = result.get('metabolite', 'Unknown')
                        progress_callback(i, total, str(metab_name))
                    except Exception:
                        pass
    else:
        rows = []
        for i, it in enumerate(items, 1):
            result = _compute_two_way_entry(it)
            rows.append(result)
            # Call callback more frequently - every 1% or every item if small dataset
            if progress_callback and (i % callback_frequency == 0 or i == total):
                try:
                    metab_name = result.get('metabolite', 'Unknown')
                    progress_callback(i, total, str(metab_name))
                except Exception:
                    pass

    out = pd.DataFrame(rows)
    
    # Rename 'metabolite' to the actual identifier column name if it's not 'metabolite'
    if metabolite_id_col and metabolite_id_col != 'metabolite' and 'metabolite' in out.columns:
        out = out.rename(columns={'metabolite': metabolite_id_col})

    if fdr and not out.empty:
        out['p_A_adj'] = _bh_adjust(out['p_A'].tolist())
        out['p_B_adj'] = _bh_adjust(out['p_B'].tolist())
        out['p_interaction_adj'] = _bh_adjust(out['p_interaction'].tolist())

    out = _apply_dynamic_pvalue_floor_frame(
        out,
        ['p_A', 'p_B', 'p_interaction', 'p_A_adj', 'p_B_adj', 'p_interaction_adj']
    )

    return out

def perform_two_way_anova_posthoc(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    factor_a_map: Dict[str, str],
    factor_b_map: Dict[str, str],
    *,
    anova_results: Optional[pd.DataFrame] = None,
    alpha: float = 0.05,
    drop_zeros: bool = True,
    min_per_level: int = 1,
    min_group_size: Optional[int] = None,
    min_group_size_percent: Optional[float] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> Dict[str, pd.DataFrame]:
    """Compute Tukey HSD posthoc pairwise comparisons for Two-Way ANOVA factors.
    
    Only computes Tukey HSD when effects are significant:
    - Factor A posthoc only if p_A_adj < alpha
    - Factor B posthoc only if p_B_adj < alpha  
    - All pairwise only if p_interaction_adj < alpha
    
    For non-significant effects, returns entries with ANOVA p-values but computed stats.

    Parameters
    ----------
    min_group_size : int, optional
        Absolute minimum number of non-zero samples required per group for pairwise comparison.
    min_group_size_percent : float, optional
        Percentage-based minimum (0-100) of non-zero samples required per group.
        If provided, takes precedence over min_group_size.

    Returns separate long-format DataFrames for factor A and factor B comparisons.

        Output columns (per row):
            <id_col>, factor, group1, group2, mean_group1, mean_group2,
            fold_change, log2_fold_change, meandiff, p_value_adj, reject,
            cohen_d, cliffs_delta, n_group1, n_group2, factor_significant
        Pairwise interaction comparisons also include `interaction_significant`.
    """
    results: Dict[str, pd.DataFrame] = {}
    if _sm_mc is None:
        return results

    # Validate mapping coverage
    missing_a = [c for c in sample_cols if c not in factor_a_map]
    missing_b = [c for c in sample_cols if c not in factor_b_map]
    if missing_a or missing_b:
        return results

    # Determine metabolite identifier column (preserve original name, don't rename to 'metabolite')
    metabolite_id_col: Optional[str] = None
    for candidate in ['LipidID', 'Name', 'Class', 'Protein', 'Metabolite', 'Molecule', 'Compound', 'Lipid_ID', 'Lipid_Class']:
        if candidate in df.columns:
            metabolite_id_col = candidate
            break
    if metabolite_id_col is None and len(df.columns) > 0:
        metabolite_id_col = df.columns[0]
    
    # Compute total samples per factor level for percentage-based thresholds
    total_samples_per_A = {}
    total_samples_per_B = {}
    for s in sample_cols:
        a = factor_a_map.get(s)
        b = factor_b_map.get(s)
        if a:
            total_samples_per_A[a] = total_samples_per_A.get(a, 0) + 1
        if b:
            total_samples_per_B[b] = total_samples_per_B.get(b, 0) + 1
    
    # Build lookup for ANOVA p-values if provided
    anova_pvals: Dict[Any, Dict[str, float]] = {}
    if anova_results is not None and metabolite_id_col in anova_results.columns:
        for _, row in anova_results.iterrows():
            met_id = row[metabolite_id_col]
            anova_pvals[met_id] = {
                'p_A_adj': row.get('p_A_adj', row.get('p_A', 1.0)),
                'p_B_adj': row.get('p_B_adj', row.get('p_B', 1.0)),
                'p_interaction_adj': row.get('p_interaction_adj', row.get('p_interaction', 1.0))
            }

    total_metabolites = len(df)
    processed = 0
    callback_frequency_posthoc = max(1, total_metabolites // 100)  # Call every 1%
    
    rows_A = []
    rows_B = []
    
    for idx, row in df.iterrows():
        processed += 1
        if progress_callback and (processed % callback_frequency_posthoc == 0 or processed == total_metabolites):
            try:
                metabolite_id_temp = row[metabolite_id_col] if metabolite_id_col else idx
                progress_callback(processed, total_metabolites, str(metabolite_id_temp))
            except Exception:
                pass
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        
        # Get ANOVA p-values for this metabolite
        met_pvals = anova_pvals.get(metabolite_id, {})
        p_A_adj = met_pvals.get('p_A_adj', 1.0)
        p_B_adj = met_pvals.get('p_B_adj', 1.0)
        p_int_adj = met_pvals.get('p_interaction_adj', 1.0)
        
        # Record significance against alpha for reporting, but always run Tukey posthoc
        sig_A = p_A_adj < alpha if anova_results is not None else False
        sig_B = p_B_adj < alpha if anova_results is not None else False
        sig_interaction = p_int_adj < alpha if anova_results is not None else False

        # Build long-form arrays
        vals: List[float] = []
        A_levels: List[str] = []
        B_levels: List[str] = []
        for s in sample_cols:
            try:
                v = float(pd.to_numeric(row[s], errors='coerce'))
            except Exception:
                v = np.nan
            if np.isnan(v):
                continue
            if drop_zeros and v == 0:
                continue
            vals.append(v)
            A_levels.append(factor_a_map[s])
            B_levels.append(factor_b_map[s])

        if len(vals) < 2:
            continue

        # Factor A posthoc across all B (main effect)
        try:
            unique_A = set(A_levels)
            if len(unique_A) >= 2:
                counts_A = pd.Series(A_levels).value_counts().to_dict()
                # Calculate thresholds for each A level
                thresholds_A = {}
                for a in unique_A:
                    if min_group_size_percent is not None:
                        total_a = total_samples_per_A.get(a, 0)
                        thresholds_A[a] = max(1, int(np.ceil(total_a * min_group_size_percent / 100.0)))
                    elif min_group_size is not None:
                        thresholds_A[a] = max(1, min_group_size)
                    else:
                        thresholds_A[a] = min_per_level
                valid_idx = [i for i, a in enumerate(A_levels) if counts_A.get(a, 0) >= thresholds_A.get(a, min_per_level)]
                if len(valid_idx) >= 2:
                    vals_A = [vals[i] for i in valid_idx]
                    labs_A = [A_levels[i] for i in valid_idx]
                    df_A = pd.DataFrame({'v': vals_A, 'A': labs_A})
                    means_A = df_A.groupby('A')['v'].mean().to_dict()
                    n_A = df_A.groupby('A')['v'].size().to_dict()

                    mcA = _sm_mc.MultiComparison(vals_A, labs_A)
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
                        tukA = mcA.tukeyhsd(alpha=alpha)
                    sumA = pd.DataFrame(data=tukA._results_table.data[1:], columns=tukA._results_table.data[0])
                    for _, rA in sumA.iterrows():
                        g1 = str(rA['group1']); g2 = str(rA['group2'])
                        mean1 = float(means_A.get(g1, np.nan)); mean2 = float(means_A.get(g2, np.nan))
                        fc, l2fc = _safe_fc(mean1, mean2)
                        v1 = df_A[df_A['A'] == g1]['v'].values
                        v2 = df_A[df_A['A'] == g2]['v'].values
                        cohen_d, cliffs = _effect_sizes(v1, v2)
                        rows_A.append({
                            metabolite_id_col: metabolite_id,
                            'factor': 'A',
                            'group1': g1,
                            'group2': g2,
                            'mean_group1': mean1,
                            'mean_group2': mean2,
                            'fold_change': fc,
                            'log2_fold_change': l2fc,
                            'meandiff': rA.get('meandiff', np.nan),
                            'p_value_adj': _floor_pval(rA.get('p-adj', np.nan)),
                            'reject': bool(rA.get('reject', False)),
                            'cohen_d': cohen_d,
                            'cliffs_delta': cliffs,
                            'n_group1': int(n_A.get(g1, 0)),
                            'n_group2': int(n_A.get(g2, 0)),
                            'factor_significant': sig_A
                        })
        except Exception:
            pass

        # Factor B posthoc across all A (main effect)
        try:
            unique_B = set(B_levels)
            if len(unique_B) >= 2:
                counts_B = pd.Series(B_levels).value_counts().to_dict()
                # Calculate thresholds for each B level
                thresholds_B = {}
                for b in unique_B:
                    if min_group_size_percent is not None:
                        total_b = total_samples_per_B.get(b, 0)
                        thresholds_B[b] = max(1, int(np.ceil(total_b * min_group_size_percent / 100.0)))
                    elif min_group_size is not None:
                        thresholds_B[b] = max(1, min_group_size)
                    else:
                        thresholds_B[b] = min_per_level
                valid_idx_B = [i for i, b in enumerate(B_levels) if counts_B.get(b, 0) >= thresholds_B.get(b, min_per_level)]
                if len(valid_idx_B) >= 2:
                    vals_B = [vals[i] for i in valid_idx_B]
                    labs_B = [B_levels[i] for i in valid_idx_B]
                    df_B = pd.DataFrame({'v': vals_B, 'B': labs_B})
                    means_B = df_B.groupby('B')['v'].mean().to_dict()
                    n_B = df_B.groupby('B')['v'].size().to_dict()

                    mcB = _sm_mc.MultiComparison(vals_B, labs_B)
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
                        tukB = mcB.tukeyhsd(alpha=alpha)
                    sumB = pd.DataFrame(data=tukB._results_table.data[1:], columns=tukB._results_table.data[0])
                    for _, rB in sumB.iterrows():
                        g1 = str(rB['group1']); g2 = str(rB['group2'])
                        mean1 = float(means_B.get(g1, np.nan)); mean2 = float(means_B.get(g2, np.nan))
                        fc, l2fc = _safe_fc(mean1, mean2)
                        v1 = df_B[df_B['B'] == g1]['v'].values
                        v2 = df_B[df_B['B'] == g2]['v'].values
                        cohen_d, cliffs = _effect_sizes(v1, v2)
                        rows_B.append({
                            metabolite_id_col: metabolite_id,
                            'factor': 'B',
                            'group1': g1,
                            'group2': g2,
                            'mean_group1': mean1,
                            'mean_group2': mean2,
                            'fold_change': fc,
                            'log2_fold_change': l2fc,
                            'meandiff': rB.get('meandiff', np.nan),
                            'p_value_adj': _floor_pval(rB.get('p-adj', np.nan)),
                            'reject': bool(rB.get('reject', False)),
                            'cohen_d': cohen_d,
                            'cliffs_delta': cliffs,
                            'n_group1': int(n_B.get(g1, 0)),
                            'n_group2': int(n_B.get(g2, 0)),
                            'factor_significant': sig_B
                        })
        except Exception:
            pass

    # Build DataFrames and add derived columns
    if rows_A:
        dfA = pd.DataFrame(rows_A)
        dfA = _apply_dynamic_pvalue_floor_frame(dfA, ['p_value_adj'])
        dfA['Expression'] = 'Not Significant'
        maskA = dfA['p_value_adj'].le(alpha) & dfA['log2_fold_change'].notna()
        dfA.loc[maskA & (dfA['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
        dfA.loc[maskA & (dfA['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        dfA['neg_log10_adj_p'] = -np.log10(dfA['p_value_adj'].replace(0, np.finfo(float).eps))
        results['A'] = dfA

    if rows_B:
        dfB = pd.DataFrame(rows_B)
        dfB = _apply_dynamic_pvalue_floor_frame(dfB, ['p_value_adj'])
        dfB['Expression'] = 'Not Significant'
        maskB = dfB['p_value_adj'].le(alpha) & dfB['log2_fold_change'].notna()
        dfB.loc[maskB & (dfB['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
        dfB.loc[maskB & (dfB['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        dfB['neg_log10_adj_p'] = -np.log10(dfB['p_value_adj'].replace(0, np.finfo(float).eps))
        results['B'] = dfB
    
    # Add ALL pairwise comparisons (A×B combinations)
    # Use individual level names (not combined labels) for comparison names
    rows_pairwise: List[Dict[str, Any]] = []
    processed = 0
    callback_frequency_pairwise = max(1, total_metabolites // 100)  # Call every 1%
    
    for idx, row in df.iterrows():
        processed += 1
        if progress_callback and (processed % callback_frequency_pairwise == 0 or processed == total_metabolites):
            try:
                metabolite_id_temp = row[metabolite_id_col] if metabolite_id_col else idx
                progress_callback(processed, total_metabolites, str(metabolite_id_temp))
            except Exception:
                pass
        metabolite_id = row[metabolite_id_col] if metabolite_id_col else idx
        
        # Group samples by A×B combination, track which samples are in each group
        combined_groups: Dict[str, List[float]] = {}
        combined_group_samples: Dict[str, List[str]] = {}  # Track sample columns per group
        for s in sample_cols:
            try:
                v = float(pd.to_numeric(row[s], errors='coerce'))
            except Exception:
                v = np.nan
            if np.isnan(v):
                continue
            if drop_zeros and v == 0:
                continue
            
            # Create combined group label (e.g., "Control_HFD")
            factor_a = factor_a_map[s]
            factor_b = factor_b_map[s]
            combined_label = f"{factor_a}_{factor_b}" if factor_a and factor_b else factor_a or factor_b
            
            if combined_label not in combined_groups:
                combined_groups[combined_label] = []
                combined_group_samples[combined_label] = []
            combined_groups[combined_label].append(v)
            combined_group_samples[combined_label].append(s)
        
        # Filter groups with insufficient samples (apply threshold settings)
        valid_groups = {}
        valid_group_samples = {}
        for g, vals in combined_groups.items():
            # Calculate required minimum based on settings
            if min_group_size_percent is not None:
                # Percentage-based threshold: need X% of group size to be non-zero
                # Count total samples in this group across all metabolites
                group_samples = [s for s in sample_cols if f"{factor_a_map[s]}_{factor_b_map[s]}" == g]
                total_group_size = len(group_samples)
                required_min = max(1, int(np.ceil(total_group_size * min_group_size_percent / 100.0)))
            elif min_group_size is not None:
                required_min = max(1, min_group_size)
            else:
                required_min = min_per_level  # Fallback to min_per_level
            
            # Only include group if it has enough non-zero values
            if len(vals) >= required_min:
                valid_groups[g] = vals
                valid_group_samples[g] = combined_group_samples[g]
        
        if len(valid_groups) < 2:
            continue
        
        # Get interaction p-value for this metabolite (for reporting only)
        p_interaction_adj_val = anova_pvals.get(metabolite_id, {}).get('p_interaction_adj', 1.0)
        interaction_sig = p_interaction_adj_val < alpha if anova_results is not None else False

        # Perform Tukey HSD on all combined groups
        try:
            all_vals = []
            all_labels = []
            for g, vals in valid_groups.items():
                all_vals.extend(vals)
                all_labels.extend([g] * len(vals))

            if len(all_vals) < 2:
                continue

            mc_result = _sm_mc.MultiComparison(all_vals, all_labels)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
                tukey_result = mc_result.tukeyhsd(alpha=alpha)
            sumPW = pd.DataFrame(data=tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])

            for _, row_pw in sumPW.iterrows():
                g1 = str(row_pw['group1']); g2 = str(row_pw['group2'])
                if g1 not in valid_groups or g2 not in valid_groups:
                    continue
                mean1 = float(np.mean(valid_groups[g1])); mean2 = float(np.mean(valid_groups[g2]))
                fc, l2fc = _safe_fc(mean1, mean2)
                cohen_d, cliffs = _effect_sizes(valid_groups[g1], valid_groups[g2])
                
                # Create comparison label using level names (not factor names)
                # e.g., "Control_vs_HFD" (NOT "Injury_vs_Diet")
                comparison_name = f"{g1}_vs_{g2}"
                
                # Build a base row with all sample columns set to empty/blank
                base_row = {
                    metabolite_id_col: metabolite_id,
                    'comparison': comparison_name,
                    'group1': g1,
                    'group2': g2,
                    'mean_group1': mean1,
                    'mean_group2': mean2,
                    'fold_change': fc,
                    'log2_fold_change': l2fc,
                    'meandiff': row_pw.get('meandiff', np.nan),
                    'p_value_adj': _floor_pval(row_pw.get('p-adj', np.nan)),
                    'reject': bool(row_pw.get('reject', False)),
                    'cohen_d': cohen_d,
                    'cliffs_delta': cliffs,
                    'n_group1': len(valid_groups[g1]),
                    'n_group2': len(valid_groups[g2]),
                    'interaction_significant': interaction_sig
                }
                
                # Add sample columns - mark which samples belong to each group
                for s in sample_cols:
                    base_row[s] = ""  # Initialize all sample columns
                
                # Mark group membership for group1 and group2
                for s in valid_group_samples.get(g1, []):
                    base_row[s] = g1
                for s in valid_group_samples.get(g2, []):
                    base_row[s] = g2
                
                rows_pairwise.append(base_row)
        except Exception:
            pass
    
    if rows_pairwise:
        dfPW = pd.DataFrame(rows_pairwise)
        dfPW = _apply_dynamic_pvalue_floor_frame(dfPW, ['p_value_adj'])
        dfPW['Expression'] = 'Not Significant'
        maskPW = dfPW['p_value_adj'].le(alpha) & dfPW['log2_fold_change'].notna()
        dfPW.loc[maskPW & (dfPW['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
        dfPW.loc[maskPW & (dfPW['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
        dfPW['neg_log10_adj_p'] = -np.log10(dfPW['p_value_adj'].replace(0, np.finfo(float).eps))
        results['pairwise'] = dfPW

    return results

__all__ = [
    'detect_sample_columns',
    'detect_feature_and_sample_columns',
    'normalize_dataframe',
    'perform_statistical_analysis',
    'get_group_summary',
    'apply_variability_filter',
    'apply_imputation',
    'apply_pca_outlier_filter',
    'perform_two_way_anova',
    'perform_two_way_anova_posthoc'
]
