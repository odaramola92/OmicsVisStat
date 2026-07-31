"""Refactored visualization services for GUI integration.

This module provides clean APIs for generating plots from statistical results,
extracted from the standalone metabolite_visualization.py CLI tool.
All functions are designed to be called from a GUI thread with proper context.
"""
from __future__ import annotations

import os
import itertools
from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Optional, Set, Union

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

# Set matplotlib to non-interactive backend before importing pyplot
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for thread safety
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
from scipy.stats import chi2
from matplotlib.patches import Ellipse, Circle
from matplotlib.ticker import MaxNLocator
import importlib
import importlib.util
from itertools import combinations
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


venn2 = None
venn3 = None
venn4 = None

try:
    if importlib.util.find_spec("matplotlib_venn") is not None:
        mv = importlib.import_module("matplotlib_venn")
        venn2 = getattr(mv, "venn2", None)
        venn3 = getattr(mv, "venn3", None)
        # venn4 may be provided by optional extensions or newer versions
        venn4 = getattr(mv, "venn4", None)
    else:
        venn2 = None
        venn3 = None
        venn4 = None
except Exception:
    # Any import/runtime error - treat as not available
    venn2 = None
    venn3 = None
    venn4 = None

# Optional 3D interactive support
try:
    import plotly.express as px
    import plotly.graph_objects as go
except ImportError:
    px = None
    go = None

# Constants from original module
FEATURE_COLUMNS_CANONICAL = [
    'Name','Name_Key','Formula','Molecular_Formula','Molecular Formula','MW','Molecular_Weight','ppm',
    'Reference Ion','Reference_Ion','MS2','m/z','RT','RT [min]','RT_mean','Area (Max.)','Polarity','MS2_Purity','MS2 Purity [%]',
    'LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID','ChEBI_ID','CAS','SMILES','InChI','InChIKey','IUPAC_Name',
    'Super_Class','Class','Sub_Class','Endogenous_Source','Metabolika Pathways','BioCyc Pathways',
    # Additional from user examples
    'LipidID', 'Class_name', 'CalcMz', 'BaseRt', 'AdductIon', 'LipidMaps_ID_Match_Type', 'Systematic_Name', 'Preferred_Name',
    'Abbreviation', 'KEGG_Match_Type', 'match_source', 'annotation_sources', 'Endogenous', 'metabolite_id'
]

try:
    from gui.shared.utils import LIPID_FEATURE_CANONICAL as LIPID_FEATURE_CANONICAL
except Exception:
    # Fallback (should rarely be used)
    LIPID_FEATURE_CANONICAL = [
        'lipidid', 'class', 'lipidgroup', 'charge', 'calcmz', 'basert', 'subclass',
        'adduction', 'ionformula', 'molstructure', 'obsmz', 'obsrt', 'ppmdiff', 'polarity'
    ]

STAT_EXCLUDE_PREFIXES = [
    'metabolite_id', 'mean_', 'n_', 'overall_',
    # Treat statistical result columns as non-sample so they are not mis-assigned
    'F_', 'p_',
]

# Exact match column names to exclude from samples (case-insensitive)
STAT_EXCLUDE_EXACT = [
    'z',  # Z-score or other non-sample column
]

# Terms that indicate a statistical column (used with substring matching)
# NOTE: Patterns without leading underscore will match anywhere in column name
# Patterns WITH leading underscore require underscore before the term
STAT_EXCLUDE_CONTAINS = [
    '_vs_',         # Pairwise comparison columns (e.g., Control_vs_PD_log2FC)
    '_statistic',   # Statistical test result columns
    '_p_value',     # P-value columns (with underscore)
    '_pvalue',      # P-value columns (without space)
    '_neg_log10',   # Negative log10 p-values
    '_log2FC',      # Log2 fold change (with underscore)
    '_log2fc',      # Log2 fold change (lowercase)
    '_adj_p',       # Adjusted p-values
    '_p_adj',       # Adjusted p-values (alternate format)
    '_adj',         # Generic adjusted suffix (captures p_Injury_adj, etc.)
    'foldchange',   # Fold change (no underscore, matches FoldChange columns)
    'fold_change',  # Fold change with underscore
]

# Parameter dataclasses
@dataclass
class CommonVizContext:
    """Shared context for all visualization functions.

    Added width/height/dpi defaults so GUI can override figure sizing.
    Includes lipid mode support and class-level data.
    """
    complete_df: pd.DataFrame
    groups: List[str]
    sample_cols: List[str]
    sample_to_group: Dict[str, str]
    output_dir: str
    color_map: Dict[str, Any]
    fig_width: float = 8.0
    fig_height: float = 6.0
    fig_dpi: int = 220
    preferred_group_order: Optional[List[str]] = None
    use_adj_p: bool = True
    is_lipid_mode: bool = False
    verified_assignments: Optional[Dict[str, str]] = None
    stat_column_assignments: Optional[Dict] = None  # From Configure Stat Columns dialog
    id_column: Optional[str] = None  # Explicit ID column propagated from GUI config
    lipid_class_df: Optional[pd.DataFrame] = None
    _lipid_feature_columns_removed: Optional[List[str]] = None
    
@dataclass
class PCAParams:
    """Parameters for PCA visualization."""
    components: int = 10
    plot_3d: bool = False
    interactive_3d: bool = False
    scree: bool = False
    loadings: bool = True
    loadings_top_k: int = 30
    pairwise: bool = True
    specific_pairs: Optional[List[Tuple[str, str]]] = None
    fig_width: float = 8.0
    fig_height: float = 6.0
    fig_dpi: int = 220
    # Custom group selection for "combined" PCA - if specified, only these groups plotted
    custom_groups: Optional[List[str]] = None
    point_size_2d: float = 35
    point_size_3d: float = 30
    # Font size controls for labels and titles
    xlabel_fontsize: int = 11
    ylabel_fontsize: int = 11
    title_fontsize: int = 12
    tick_fontsize: int = 11
    legend_fontsize: int = 10
    show_legend: bool = True  # Show legend on all PCA plots (2D and 3D)
    # Save options - allow users to choose which plots to generate
    save_2d: bool = True
    save_3d: bool = False
    save_excel: bool = True
    # 3D viewing angles (azimuth and elevation)
    view_azim: float = -60  # Azimuth angle for 3D plot rotation (degrees)
    view_elev: float = 30   # Elevation angle for 3D plot rotation (degrees)
    # Comparison selection for pairwise - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None
    # Lipid class PCA options
    include_lipid_class: bool = False
    class_subdir_name: Optional[str] = 'pca_class'

@dataclass
class VolcanoParams:
    """Parameters for volcano plots."""
    p_threshold: float = 0.05
    fc_threshold: float = 2.0
    annotate_top_n: int = 0
    annotate: bool = False
    interactive: bool = False
    fig_width: float = 5.0
    fig_height: float = 4.0
    fig_dpi: int = 230
    point_size_sig: float = 28
    point_size_nonsig: float = 18
    # Font size controls
    xlabel_fontsize: int = 11
    ylabel_fontsize: int = 11
    title_fontsize: int = 18
    tick_fontsize: int = 11
    legend_fontsize: int = 10
    # Count annotation controls
    count_fontsize: int = 9
    total_fontsize: int = 8
    count_background: str = 'colored'  # 'colored' or 'transparent'
    save_excel: bool = True
    # Comparison selection - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None

@dataclass
class BoxplotParams:
    """Parameters for boxplot generation."""
    top_n: int = 10
    no_limit: bool = False
    annotate: bool = True
    alpha: float = 0.05
    include_metabolites: Optional[List[str]] = None
    verbose: bool = False
    fig_width: float = 3.0
    fig_height: float = 3.0
    fig_dpi: int = 240
    p_threshold: float = 0.05
    fc_threshold: float = 2.0
    filter_mode: str = 'any'  # any | all | specific
    filter_pairs: Optional[List[Tuple[str, str]]] = None  # only used if filter_mode == 'specific'
    use_custom_only: bool = False  # if True, ignore statistical filtering and show ONLY custom list
    # Font size controls
    xlabel_fontsize: int = 11
    ylabel_fontsize: int = 11
    title_fontsize: int = 12
    tick_fontsize: int = 11
    legend_fontsize: int = 10
    save_excel: bool = True
    # Annotation filtering - None means annotate all, List means only annotate specified pairs
    annotate_comparisons: Optional[List[Tuple[str, str]]] = None
    # Comparison selection - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None
    # Group selection - None means all groups, List means only display specified groups
    selected_groups: Optional[List[str]] = None
    # Filter comparison for All/Specific modes
    filter_comparison: Optional[Tuple[str, str]] = None
    # Tick label rotation controls
    rotate_xticks: bool = True
    xtick_rotation: int = 45
    # Y-axis label customization
    ylabel_text: str = 'Relative Abundance (%)'
    # Title wrapping control
    title_wrap_width: int = 25  # Character limit before wrapping title to multiple lines

@dataclass
class BargraphParams:
    """Parameters for bar graph generation.

    Behavior intentionally mirrors BoxplotParams, with a display-mode toggle.
    """
    top_n: int = 10
    no_limit: bool = False
    annotate: bool = True
    alpha: float = 0.05
    include_metabolites: Optional[List[str]] = None
    verbose: bool = False
    fig_width: float = 3.0
    fig_height: float = 3.0
    fig_dpi: int = 240
    p_threshold: float = 0.05
    fc_threshold: float = 2.0
    filter_mode: str = 'any'  # any | all | specific
    filter_pairs: Optional[List[Tuple[str, str]]] = None  # only used if filter_mode == 'specific'
    use_custom_only: bool = False  # if True, ignore statistical filtering and show ONLY custom list
    # Font size controls
    xlabel_fontsize: int = 11
    ylabel_fontsize: int = 11
    title_fontsize: int = 12
    tick_fontsize: int = 11
    legend_fontsize: int = 10
    save_excel: bool = True
    # Annotation filtering - None means annotate all, List means only annotate specified pairs
    annotate_comparisons: Optional[List[Tuple[str, str]]] = None
    # Comparison selection - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None
    # Group selection - None means all groups, List means only display specified groups
    selected_groups: Optional[List[str]] = None
    # Filter comparison for All/Specific modes
    filter_comparison: Optional[Tuple[str, str]] = None
    # Tick label rotation controls
    rotate_xticks: bool = True
    xtick_rotation: int = 45
    # Y-axis label customization
    ylabel_text: str = 'Relative Abundance (%)'
    # Title wrapping control
    title_wrap_width: int = 25
    # Bar graph specific mode
    display_mode: str = 'separate'  # separate | grouped
    grouped_title: str = ''
    # Grouped-mode low-value visual boost (plot-only scaling)
    low_value_boost_enabled: bool = False
    low_value_boost_threshold: float = 0.25
    low_value_boost_factor: float = 2.0

@dataclass
class HeatmapParams:
    """Parameters for heatmap generation."""
    p_threshold: float = 0.05
    fc_threshold: float = 2.0
    max_metabolites: int = 0
    cluster: bool = True
    scale: str = 'row'
    column_split_counts: Optional[str] = None
    row_split: int = 0
    row_split_mode: str = 'cluster'
    show_fc_divider: bool = True  # Show dividing line between up/down regulated
    # Layout controls
    show_colorbar: bool = True  # Whether to display a colorbar (top row)
    dendrogram_width_ratio: float = 0.18  # Relative width for dendrogram column vs heatmap column (GridSpec width_ratios)
    colorbar_height_ratio: float = 0.12   # Relative height for colorbar row vs heatmap row (GridSpec height_ratios) - DEPRECATED, use colorbar_height_inches
    colorbar_height_inches: float = 0.6   # Fixed height for colorbar in inches (replaces ratio-based sizing)
    font_size: int = 14
    legend_font_size: int = 12
    title_align: str = 'center'
    no_col_split: bool = False
    combined: bool = False
    combined_mode: str = 'union'  # union | intersection | concatenate
    order_by_mean: bool = True
    export_clusters: bool = False
    show_dendrogram: bool = True
    include_metabolites: Optional[List[str]] = None
    metabolite_list: Optional[List[str]] = None
    fig_width: float = 8.0
    fig_height: float = 6.0
    fig_dpi: int = 190
    auto_size: bool = True  # Auto-calculate dimensions based on data
    # Added dynamic filters & scaling
    use_fixed_scale: bool = True  # Use default -3 to 3 scale (checked by default)
    auto_scale: bool = False      # Auto scale to 5th-95th percentile
    vmin: float = -3.0           # Manual min (used when both above are False)
    vmax: float = 3.0            # Manual max (used when both above are False)
    # Extended filtering controls (not yet fully wired for heatmap but reserved)
    filter_mode: str = 'any'
    filter_pairs: Optional[List[Tuple[str, str]]] = None
    use_custom_only: bool = False
    # Font size controls - simplified for heatmaps (rows = features/metabolites, columns = samples)
    feature_fontsize: int = 10  # For y-axis labels (metabolite names on rows)
    sample_fontsize: int = 10   # For x-axis labels (sample names on columns)
    title_fontsize: int = 14    # For plot title
    save_excel: bool = True
    # Per-comparison metabolite lists: key is tuple(g1, g2) or 'all' for global
    metabolite_lists: Dict[Union[str, Tuple[str, str]], List[str]] = field(default_factory=dict)
    skip_unlisted_comparisons: bool = False
    # Comparison selection - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None
    # Metabolite scope for combined heatmaps: 'selected' or 'all'
    metabolite_scope: str = 'selected'  # 'selected' for selected comparisons, 'all' for all metabolites

@dataclass
class ROCParams:
    """Parameters for ROC curve generation."""
    target_positive: str = ''
    target_negative: str = ''
    metabolites: Optional[List[str]] = None
    max_metabolites: int = 20
    min_auc: float = 0.5
    include_combined: bool = False  # Combined ROC from all plotted metabolites
    force_combined: bool = False
    all_pairs: bool = False
    top_n_per_pair: int = 5
    fig_width: float = 8.0
    fig_height: float = 6.8
    fig_dpi: int = 260
    p_threshold: float = 0.05
    fc_threshold: float = 2.0
    filter_mode: str = 'any'
    filter_pairs: Optional[List[Tuple[str, str]]] = None
    use_custom_only: bool = False
    # Font size controls
    xlabel_fontsize: int = 11
    ylabel_fontsize: int = 11
    title_fontsize: int = 12
    tick_fontsize: int = 11
    legend_fontsize: int = 10
    save_excel: bool = True
    excel_only: bool = False
    # Per-comparison metabolite lists: key is tuple(g1, g2) or 'all' for global
    metabolite_lists: Dict[Union[str, Tuple[str, str]], List[str]] = field(default_factory=dict)
    skip_unlisted_comparisons: bool = False
    # Comparison selection - None means all comparisons
    selected_comparisons: Optional[List[Tuple[str, str]]] = None

@dataclass
class PCAResult:
    """Results from PCA analysis."""
    scores: pd.DataFrame
    explained_variance: Tuple[float, float]
    pca_model: Optional[PCA] = None
    feature_ids: Optional[List[str]] = None

@dataclass
@dataclass
class VizResults:
    """Results from visualization generation."""
    files_created: List[str]
    errors: List[str]
    summary: str
    venn_summaries: Optional[List[str]] = None  # Per-venn text summaries

# ------------------ Venn Diagram Support ------------------
from dataclasses import field

@dataclass
class VennSpec:
    """Definition of a single Venn plot: a name and a list of pairwise comparisons.
    comparisons: list of (group1, group2) tuples; order matters for label naming only.
    """
    name: str
    comparisons: List[Tuple[str, str]]

# ------------------ Multi-Column Metabolite Matching ------------------

def match_metabolites_multi_column(df: pd.DataFrame, metabolite_list, id_col: str = 'metabolite_id') -> pd.Series:
    """Match metabolites using multiple columns for better accuracy.
    
    Supports both legacy format (simple list of names) and new dict format with IDs.
    
    Args:
        df: DataFrame with metabolite data
        metabolite_list: Either a list of metabolite names (legacy) or a dict with:
            {
                'names': [list of metabolite names],
                'pubchem_ids': [list of PubChem IDs] (optional),
                'hmdb_ids': [list of HMDB IDs] (optional),
                'cas_ids': [list of CAS IDs] (optional),
                'gene_symbols': [list of gene symbols] (optional),
                'accessions': [list of protein accessions e.g., UniProt] (optional)
            }
        id_col: Primary ID column in df (usually 'metabolite_id' or 'Name')
    
    Returns:
        Boolean Series indicating which rows match the metabolite list
    """
    if metabolite_list is None:
        return pd.Series([False] * len(df), index=df.index)
    
    # Handle legacy format (simple list of names)
    if isinstance(metabolite_list, list):
        name_set = {str(n).lower().strip() for n in metabolite_list}
        return df[id_col].astype(str).str.lower().str.strip().isin(name_set)
    
    # Handle new dict format with multi-column support
    if not isinstance(metabolite_list, dict):
        return pd.Series([False] * len(df), index=df.index)
    
    # Start with all False
    match_mask = pd.Series([False] * len(df), index=df.index)
    matched_by_type = {}
    
    # Match by Name (case-insensitive)
    if 'names' in metabolite_list:
        name_set = {str(n).lower().strip() for n in metabolite_list['names']}
        name_matches = df[id_col].astype(str).str.lower().str.strip().isin(name_set)
        match_mask |= name_matches
        matched_by_type['name'] = name_matches.sum()
    
    # Match by PubChem ID if available in both list and dataframe
    if 'pubchem_ids' in metabolite_list:
        pubchem_col = None
        for col in ['PubChem_CID', 'PubChem', 'pubchem_cid', 'pubchem']:
            if col in df.columns:
                pubchem_col = col
                break
        
        if pubchem_col:
            pubchem_set = {str(pid).strip() for pid in metabolite_list['pubchem_ids']}
            pubchem_matches = df[pubchem_col].astype(str).str.strip().isin(pubchem_set)
            before_count = match_mask.sum()
            match_mask |= pubchem_matches
            matched_by_type['pubchem'] = match_mask.sum() - before_count
    
    # Match by HMDB ID if available
    if 'hmdb_ids' in metabolite_list:
        hmdb_col = None
        for col in ['HMDB_ID', 'HMDB', 'hmdb_id', 'hmdb']:
            if col in df.columns:
                hmdb_col = col
                break
        
        if hmdb_col:
            hmdb_set = {str(hid).strip() for hid in metabolite_list['hmdb_ids']}
            hmdb_matches = df[hmdb_col].astype(str).str.strip().isin(hmdb_set)
            before_count = match_mask.sum()
            match_mask |= hmdb_matches
            matched_by_type['hmdb'] = match_mask.sum() - before_count
    
    # Match by CAS ID if available
    if 'cas_ids' in metabolite_list:
        cas_col = None
        for col in ['CAS', 'cas', 'CAS_ID', 'cas_id']:
            if col in df.columns:
                cas_col = col
                break
        
        if cas_col:
            cas_set = {str(cid).strip() for cid in metabolite_list['cas_ids']}
            cas_matches = df[cas_col].astype(str).str.strip().isin(cas_set)
            before_count = match_mask.sum()
            match_mask |= cas_matches
            matched_by_type['cas'] = match_mask.sum() - before_count

    # Match by Gene symbol if available
    if 'gene_symbols' in metabolite_list:
        gene_col = None
        for col in ['Gene', 'gene', 'Symbol', 'symbol', 'Gene_Symbol', 'GeneSymbol', 'gene_symbol', 'genesymbol']:
            if col in df.columns:
                gene_col = col
                break
        if gene_col:
            gene_set = {str(gs).strip().upper() for gs in metabolite_list['gene_symbols']}
            # Normalize df gene column to uppercase no surrounding spaces
            gene_matches = df[gene_col].astype(str).str.strip().str.upper().isin(gene_set)
            before_count = match_mask.sum()
            match_mask |= gene_matches
            matched_by_type['gene'] = match_mask.sum() - before_count

    # Match by protein accession if available
    if 'accessions' in metabolite_list:
        acc_col = None
        for col in [
            'Accession', 'accession', 'Protein_Accession', 'Protein Accession',
            'UniProt', 'Uniprot', 'UniProt_ID', 'Uniprot_ID', 'UniProt Accession', 'Uniprot Accession',
            'Swiss-Prot', 'SwissProt', 'GenPept', 'Genpept', 'ProteinID', 'Protein_ID'
        ]:
            if col in df.columns:
                acc_col = col
                break
        if acc_col:
            acc_set = {str(a).strip().upper() for a in metabolite_list['accessions']}
            acc_matches = df[acc_col].astype(str).str.strip().str.upper().isin(acc_set)
            before_count = match_mask.sum()
            match_mask |= acc_matches
            matched_by_type['accession'] = match_mask.sum() - before_count
    
    # Debug output with logger
    total_matched = match_mask.sum()
    logger.info(f"🔍 Multi-column matching: {total_matched}/{len(df)} metabolites matched")
    for match_type, count in matched_by_type.items():
        if count > 0:
            logger.info(f"   ✓ {match_type}: {count} additional matches")
    
    # Show which columns were checked
    checked_cols = [id_col]
    if 'pubchem_ids' in metabolite_list:
        checked_cols.append(f"PubChem ({matched_by_type.get('pubchem', 0)} matches)")
    if 'hmdb_ids' in metabolite_list:
        checked_cols.append(f"HMDB ({matched_by_type.get('hmdb', 0)} matches)")
    if 'cas_ids' in metabolite_list:
        checked_cols.append(f"CAS ({matched_by_type.get('cas', 0)} matches)")
    if 'gene_symbols' in metabolite_list:
        checked_cols.append(f"Gene ({matched_by_type.get('gene', 0)} matches)")
    if 'accessions' in metabolite_list:
        checked_cols.append(f"Accession ({matched_by_type.get('accession', 0)} matches)")
    logger.info(f"   Columns checked: {', '.join(checked_cols)}")
    
    return match_mask

@dataclass
class VennParams:
    """Parameters for generating Venn diagrams from pairwise stats."""
    p_threshold: float = 0.05
    fc_threshold: float = 2.0  # set to 0 to ignore FC cutoff
    venn_specs: List[VennSpec] = field(default_factory=list)
    fig_width: float = 6.0
    fig_height: float = 6.0
    fig_dpi: int = 220
    # If True (default), prefer adjusted p-value columns when present
    prefer_adj_p: bool = True
    # Venn diagram number font size controls
    venn_number_fontsize: int = 16  # Default to bigger, bold font
    venn_label_fontsize: int = 11
    # Skip all cutoffs and use presence/absence instead
    skip_all_cutoffs: bool = False
    min_presence_type: str = 'count'  # 'count' or 'percentage'
    min_presence_count: int = 3
    min_presence_percent: float = 50.0
    # Output format ('png' or 'svg')
    output_format: str = 'png'

def _get_metabolites_for_comparison(params, g1: str, g2: str) -> Optional[List[str]]:
    """Get metabolite list for a specific comparison from params.
    
    Returns the metabolite list for the given comparison pair, trying both
    orderings. Falls back to 'all' if no specific list found.
    
    Parameters
    ----------
    params : HeatmapParams or ROCParams
        Parameters object with metabolite_lists attribute
    g1, g2 : str
        Group names for the comparison
        
    Returns
    -------
    Optional[List[str]]
        List of metabolite IDs/names, or None if no list configured
    """
    if not hasattr(params, 'metabolite_lists'):
        return None
    
    # Try exact match
    if (g1, g2) in params.metabolite_lists:
        return params.metabolite_lists[(g1, g2)]
    
    # Try reversed
    if (g2, g1) in params.metabolite_lists:
        return params.metabolite_lists[(g2, g1)]
    
    # Fall back to global list
    if 'all' in params.metabolite_lists:
        return params.metabolite_lists['all']
    
    # No list configured
    return None

def _should_skip_comparison(params, g1: str, g2: str) -> bool:
    """Check if a comparison should be skipped based on params.selected_comparisons.
    
    Parameters
    ----------
    params : Any params object with selected_comparisons attribute
        Parameters object
    g1, g2 : str
        Group names for the comparison
        
    Returns
    -------
    bool
        True if comparison should be skipped, False otherwise
    """
    if not hasattr(params, 'selected_comparisons') or params.selected_comparisons is None:
        return False  # No filtering, include all
    
    # Check if this comparison is in the selected list (either order)
    return not ((g1, g2) in params.selected_comparisons or (g2, g1) in params.selected_comparisons)


def _normalize_feature_id_for_venn(value: Any) -> Optional[str]:
    """Normalize feature IDs for Venn export while preserving textual identifiers.

    This keeps plain text untouched, but converts integer-like numeric values such as
    44101.0 to "44101" so Excel sheets do not show trailing .0 for IDs.
    """
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    if isinstance(value, (int, np.integer)):
        return str(int(value))

    if isinstance(value, (float, np.floating)):
        fvalue = float(value)
        if fvalue.is_integer():
            return str(int(fvalue))
        return str(value)

    text = str(value).strip()
    if not text:
        return None

    # Only collapse decimal text with an all-zero fractional part (e.g., "44101.0").
    if "." in text:
        try:
            whole, frac = text.split(".", 1)
            if frac and set(frac) == {"0"}:
                nvalue = float(text)
                if nvalue.is_integer():
                    return str(int(nvalue))
        except Exception:
            pass

    return text

def _significant_ids_for_pair(df: pd.DataFrame, g1: str, g2: str, *, p_thresh: float, fc_thresh: float | int,
                              prefer_adj: bool = True, stat_column_assignments: Optional[Dict] = None, id_column: Optional[str] = None) -> Tuple[Optional[str], Set[str]]:
    found = _locate_pair_columns(df, g1, g2, prefer_adj=prefer_adj, stat_column_assignments=stat_column_assignments)
    if not found:
        logger.warning(f"   Venn: Could not find stats columns for {g1} vs {g2}. Ensure statistics analysis was run first.")
        return None, set()
    log2fc_col, p_col, used_base = found
    
    # Use provided id_column (must be configured via Configure Stat Columns dialog)
    if not id_column or id_column not in df.columns:
        logger.warning(f"   Venn: No ID column configured for {g1} vs {g2}, returning empty set. Configure via 'Configure Stat Columns' dialog.")
        return None, set()
    
    id_col = id_column
    
    logger.debug(f"   Venn: Using ID column '{id_col}' for {g1} vs {g2}, p_col={p_col}, fc_col={log2fc_col}")
    
    try:
        series_p = df[p_col].astype(float)
    except Exception:
        series_p = pd.to_numeric(df[p_col], errors='coerce')
    try:
        series_fc = df[log2fc_col].astype(float)
    except Exception:
        series_fc = pd.to_numeric(df[log2fc_col], errors='coerce')

    sig_mask = series_p < p_thresh
    try:
        use_fc = (float(fc_thresh) > 0)
    except Exception:
        use_fc = False
    if use_fc:
        try:
            log2th = float(np.log2(float(fc_thresh))) if float(fc_thresh) > 1.0 else 0.0
        except Exception:
            log2th = 0.0
        sig_mask = sig_mask & (series_fc.abs() >= log2th)

    selected_ids = df.loc[sig_mask & series_p.notna() & series_fc.notna(), id_col]
    ids = {
        norm
        for norm in selected_ids.map(_normalize_feature_id_for_venn)
        if norm is not None
    }
    
    # Log the filtering results
    logger.info(f"   Venn: {g1} vs {g2}: p<{p_thresh}, fc_thresh={fc_thresh} (use_fc={use_fc}) -> {len(ids)} significant metabolites")
    if len(ids) == 0:
        logger.warning(f"   Venn: No significant metabolites found for {g1} vs {g2}. Check if statistical analysis was run and columns exist.")
        logger.debug(f"   Venn: Available columns: {list(df.columns)[:20]}...")
    
    # Use g1_vs_g2 format (baseline-first) for label instead of column name
    label = f"{g1}_vs_{g2}"
    return label, ids

def _present_ids_for_group(df: pd.DataFrame, group: str, sample_cols: List[str], sample_to_group: Dict[str, str], 
                           min_presence_type: str, min_presence_count: int, min_presence_percent: float, 
                           id_column: Optional[str] = None) -> Tuple[str, Set[str]]:
    """Get metabolite IDs present in a group based on non-zero/non-NaN sample count.
    
    Parameters
    ----------
    df : pd.DataFrame
        Complete dataframe with sample columns
    group : str
        Group name to check
    sample_cols : List[str]
        All sample column names
    sample_to_group : Dict[str, str]
        Mapping of sample columns to group names
    min_presence_type : str
        'count' or 'percentage'
    min_presence_count : int
        Minimum number of non-zero samples required
    min_presence_percent : float
        Minimum percentage of non-zero samples required (0-100)
    id_column : Optional[str]
        Explicit ID column name to use (e.g., 'Class', 'Class_name' for lipid class data)
        
    Returns
    -------
    Tuple[str, Set[str]]
        (group_name, set of metabolite IDs present in group)
    """
    # Use explicit id_column if provided, otherwise fallback to standard column names
    if id_column and id_column in df.columns:
        id_col = id_column
    else:
        id_col = 'metabolite_id' if 'metabolite_id' in df.columns else 'Name' if 'Name' in df.columns else None
    
    if not id_col:
        return group, set()
    
    # Get sample columns belonging to this group
    group_samples = [col for col in sample_cols if sample_to_group.get(col) == group]
    if not group_samples:
        return group, set()
    
    n_samples = len(group_samples)
    
    # Determine threshold
    if min_presence_type == 'percentage':
        threshold = int(np.ceil(n_samples * min_presence_percent / 100.0))
    else:
        threshold = min(min_presence_count, n_samples)  # Can't require more than available samples
    
    # Count non-zero/non-NaN values per metabolite
    present_ids = set()
    for idx, row in df.iterrows():
        non_zero_count = 0
        for col in group_samples:
            try:
                val = pd.to_numeric(row[col], errors='coerce')
                if not pd.isna(val) and val != 0:
                    non_zero_count += 1
            except Exception:
                pass
        
        if non_zero_count >= threshold:
            try:
                normalized_id = _normalize_feature_id_for_venn(row[id_col])
                if normalized_id is not None:
                    present_ids.add(normalized_id)
            except Exception:
                pass
    
    return group, present_ids

# ----------------- Equal-size Venn rendering -----------------
def _draw_equal_venn2(ax, name1: str, s1: Set[str], name2: str, s2: Set[str], *, labels_on_top: bool = False, number_fontsize: int = 16, label_fontsize: int = 11) -> None:
    """Equal-radius two-set Venn with truthful counts. Overlap is a fixed pleasant amount.
    
    Args:
        labels_on_top: If True, place labels above circles; if False, below circles (for alternating)
        number_fontsize: Font size for the count numbers (default 16, bold)
        label_fontsize: Font size for the set labels (default 11, bold)
    """
    # geometry (why: constant sizes → no visual bias)
    r = 1.0
    d = 1.2  # center distance (0<d<2r); tweak for how much overlap you want
    c1 = (-d/2, 0.0)
    c2 = ( d/2, 0.0)

    # circles
    circ1 = Circle(c1, r, linewidth=1.5, edgecolor="black", facecolor="#FF6B6B", alpha=0.5)
    circ2 = Circle(c2, r, linewidth=1.5, edgecolor="black", facecolor="#51CF66", alpha=0.5)
    ax.add_patch(circ1)
    ax.add_patch(circ2)

    # counts
    only1 = len(s1 - s2)
    only2 = len(s2 - s1)
    both  = len(s1 & s2)

    # labels (make counts bold and use custom font size)
    ax.text(c1[0] - 0.45, 0.0, str(only1), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')
    ax.text(c2[0] + 0.45, 0.0, str(only2), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')
    ax.text(0.0, 0.0, str(both), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')

    # set titles (bold) - alternate position for multiple 2-way Venns
    if labels_on_top:
        # Labels on top
        ax.text(c1[0], 1.35, name1, ha="center", va="bottom", fontsize=label_fontsize, fontweight='bold')
        ax.text(c2[0], 1.35, name2, ha="center", va="bottom", fontsize=label_fontsize, fontweight='bold')
    else:
        # Labels on bottom (default)
        ax.text(c1[0], -1.35, name1, ha="center", va="top", fontsize=label_fontsize, fontweight='bold')
        ax.text(c2[0], -1.35, name2, ha="center", va="top", fontsize=label_fontsize, fontweight='bold')

    # frame/layout
    ax.set_aspect("equal")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-1.8, 1.8)
    ax.axis("off")

def _draw_equal_venn3(ax, n1: str, s1: Set[str], n2: str, s2: Set[str], n3: str, s3: Set[str], *, number_fontsize: int = 16, label_fontsize: int = 11) -> None:
    """Equal-radius three-set Venn with proper overlap; triangle layout matching R's VennDiagram.
    
    Args:
        number_fontsize: Font size for the count numbers (default 16, bold)
        label_fontsize: Font size for the set labels (default 11, bold)
    """
    r = 1.0
    # Adjusted positions for better 3-way overlap (similar to R's VennDiagram package)
    cA = (-0.7,  0.4)   # Top-left
    cB = ( 0.7,  0.4)   # Top-right
    cC = ( 0.0, -0.6)   # Bottom center

    circA = Circle(cA, r, linewidth=1.5, edgecolor="black", facecolor="#74C0FC", alpha=0.5)
    circB = Circle(cB, r, linewidth=1.5, edgecolor="black", facecolor="#FFD43B", alpha=0.5)
    circC = Circle(cC, r, linewidth=1.5, edgecolor="black", facecolor="#FF8787", alpha=0.5)
    for c in (circA, circB, circC):
        ax.add_patch(c)

    # region counts
    Aonly = len(s1 - s2 - s3)
    Bonly = len(s2 - s1 - s3)
    Conly = len(s3 - s1 - s2)
    AB    = len((s1 & s2) - s3)
    AC    = len((s1 & s3) - s2)
    BC    = len((s2 & s3) - s1)
    ABC   = len(s1 & s2 & s3)

    # text positions (adjusted for better overlap geometry)
    # region counts (bold, custom font size)
    ax.text(cA[0] - 0.5, cA[1] + 0.15, str(Aonly), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')
    ax.text(cB[0] + 0.5, cB[1] + 0.15, str(Bonly), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')
    ax.text(cC[0],        cC[1] - 0.5, str(Conly), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')

    ax.text(0.0,   0.75, str(AB),  ha="center", va="center", fontsize=number_fontsize, fontweight='bold')  # AB overlap
    ax.text(-0.55, -0.15, str(AC), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')  # AC overlap
    ax.text( 0.55, -0.15, str(BC), ha="center", va="center", fontsize=number_fontsize, fontweight='bold')  # BC overlap
    ax.text(0.0,   0.0, str(ABC),  ha="center", va="center", fontsize=number_fontsize+1, fontweight='bold')  # Center intersection

    # set titles (bold) - positioned outside circles like R
    ax.text(cA[0], cA[1] + 1.2, n1, ha="center", va="bottom", fontsize=label_fontsize, fontweight='bold')
    ax.text(cB[0], cB[1] + 1.2, n2, ha="center", va="bottom", fontsize=label_fontsize, fontweight='bold')
    ax.text(cC[0], cC[1] - 1.2, n3, ha="center", va="top", fontsize=label_fontsize, fontweight='bold')

    ax.set_aspect("equal")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-2.0, 2.0)
    ax.axis("off")

def draw_venn4(ax, names, sets, number_fontsize=16, label_fontsize=12):
    """
    Clean 4-way Venn diagram using 4 overlapping ellipses.
    Creates a cleaner layout with better visibility of all regions.

    Args:
        ax: Matplotlib axis
        names: list of 4 set names
        sets: list of 4 Python sets
        number_fontsize: font size for counts
        label_fontsize: font size for set labels
    """
    from matplotlib.patches import Ellipse

    # =========================
    # ELLIPSE COLORS (softer, more visible)
    # =========================
    colors = [
        "#A6CEE3",  # light blue (top-left)
        "#B2DF8A",  # light green (top-right)
        "#FB9A99",  # light pink (bottom-right)
        "#FDBF6F",  # light orange (bottom-left)
    ]

    # =========================
    # ELLIPSE POSITIONS
    # Larger ellipses with better overlap for cleaner appearance
    # =========================
    ellipse_params = [
        (-0.40,  0.40, 2.4, 1.6, 0),    # Top-left (blue)
        ( 0.40,  0.40, 2.4, 1.6, 0),    # Top-right (green)
        ( 0.40, -0.40, 2.4, 1.6, 0),    # Bottom-right (pink)
        (-0.40, -0.40, 2.4, 1.6, 0),    # Bottom-left (orange)
    ]

    # Draw ellipses with cleaner styling
    for (x, y, w, h, angle), color in zip(ellipse_params, colors):
        ax.add_patch(
            Ellipse((x, y), w, h, angle=angle,
                    edgecolor="black", facecolor=color, alpha=0.4, lw=1.5)
        )

    # Unpack sets
    A, B, C, D = sets

    # =========================
    # REGION COUNTS
    # =========================
    onlyA = len(A - B - C - D)
    onlyB = len(B - A - C - D)
    onlyC = len(C - A - B - D)
    onlyD = len(D - A - B - C)

    AB  = len((A & B) - C - D)
    AC  = len((A & C) - B - D)
    AD  = len((A & D) - B - C)
    BC  = len((B & C) - A - D)
    BD  = len((B & D) - A - C)
    CD  = len((C & D) - A - B)

    ABC = len((A & B & C) - D)
    ABD = len((A & B & D) - C)
    ACD = len((A & C & D) - B)
    BCD = len((B & C & D) - A)

    ABCD = len(A & B & C & D)

    fs = number_fontsize

    # Helper function to only show non-zero counts
    def show_count(x, y, count, fontsize=fs, bold=False):
        if count > 0:
            weight = "bold" if bold else "normal"
            ax.text(x, y, str(count), ha="center", va="center", 
                   fontsize=fontsize, fontweight=weight)

    # =========================
    # TEXT LABEL POSITIONS (optimized for cleaner layout)
    # =========================

    # Single-set outer regions (corners)
    show_count(-0.95,  0.95, onlyA, bold=True)
    show_count( 0.95,  0.95, onlyB, bold=True)
    show_count( 0.95, -0.95, onlyC, bold=True)
    show_count(-0.95, -0.95, onlyD, bold=True)

    # Pairwise overlaps
    show_count( 0.00,  0.85, AB)  # Top center (A&B)
    show_count(-0.60,  0.00, AD)  # Left center (A&D)
    show_count( 0.60,  0.00, BC)  # Right center (B&C)
    show_count( 0.00, -0.85, CD)  # Bottom center (C&D)
    show_count(-0.35,  0.35, AC)  # Top-left inner (A&C)
    show_count( 0.35, -0.35, BD)  # Bottom-right inner (B&D)

    # Triple overlaps (smaller font)
    show_count( 0.00,  0.35, ABC, fontsize=fs-2)  # Top inner
    show_count( 0.35,  0.00, ABD, fontsize=fs-2)  # Right inner
    show_count(-0.35,  0.00, ACD, fontsize=fs-2)  # Left inner
    show_count( 0.00, -0.35, BCD, fontsize=fs-2)  # Bottom inner

    # All four center (largest, bold)
    show_count(0.00, 0.00, ABCD, fontsize=fs+4, bold=True)

    # =========================
    # Set Labels (positioned outside ellipses)
    # =========================
    ax.text(-1.00,  1.30, names[0], ha="center", fontsize=label_fontsize+2, fontweight="bold")
    ax.text( 1.00,  1.30, names[1], ha="center", fontsize=label_fontsize+2, fontweight="bold")
    ax.text( 1.00, -1.30, names[2], ha="center", fontsize=label_fontsize+2, fontweight="bold")
    ax.text(-1.00, -1.30, names[3], ha="center", fontsize=label_fontsize+2, fontweight="bold")

    # Final axis settings
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.8, 1.8)
    ax.set_aspect("equal")
    ax.axis("off")

def plot_venn(
    sets_named: List[Tuple[str, Set[str]]],
    out_png: str,
    *,
    title: str,
    width: float,
    height: float,
    dpi: int,
    equal_size: bool = True,  # <- key switch (default: equal circles)
    venn_index: int = 0,  # <- for alternating 2-way Venn label positions
    number_fontsize: int = 16,  # <- font size for numbers (default bigger and bold)
    label_fontsize: int = 11,  # <- font size for labels
    output_format: str = 'png',  # <- 'png' or 'svg'
) -> Optional[str]:
    """
    Draw Venn diagrams for 2-4 sets.
    - 2 sets: Side-by-side circles with alternating legend positions (odd on top, even on bottom)
    - 3 sets: Triangle layout with R-style legend positioning
    - 4 sets: Custom ellipse layout (R VennDiagram style)
    
    If equal_size=True: draw equal-radius Venns (geometry constant; counts truthful).
    If equal_size=False: fall back to matplotlib-venn with area-proportional geometry.
    
    Args:
        venn_index: Index of this Venn (0-based). Used to alternate label positions for 2-way Venns.
        number_fontsize: Font size for the count numbers (default 16, bold)
        label_fontsize: Font size for the set labels (default 11, bold)
    """
    try:
        n = len(sets_named)
        fig = plt.figure(figsize=(width, height))
        ax = fig.add_subplot(111)

        if n == 2:
            (name1, s1), (name2, s2) = sets_named
            if equal_size:
                # Alternate: odd indices (1st, 3rd, 5th...) on top, even (2nd, 4th...) on bottom
                labels_on_top = (venn_index % 2 == 0)
                _draw_equal_venn2(ax, name1, s1, name2, s2, labels_on_top=labels_on_top, number_fontsize=number_fontsize, label_fontsize=label_fontsize)
            else:
                if venn2 is None:
                    raise RuntimeError("matplotlib-venn not installed")
                v = venn2(subsets=(s1, s2), set_labels=(name1, name2), ax=ax)
                if v:
                    for pid in ("10", "01", "11"):
                        p = v.get_patch_by_id(pid)
                        if p: p.set_alpha(0.5); p.set_linewidth(1.2)

        elif n == 3:
            (n1, s1), (n2, s2), (n3, s3) = sets_named
            if equal_size:
                _draw_equal_venn3(ax, n1, s1, n2, s2, n3, s3, number_fontsize=number_fontsize, label_fontsize=label_fontsize)
            else:
                if venn3 is None:
                    raise RuntimeError("matplotlib-venn not installed")
                v = venn3(subsets=(s1, s2, s3), set_labels=(n1, n2, n3), ax=ax)
                if v:
                    for pid in ("100","010","001","110","101","011","111"):
                        p = v.get_patch_by_id(pid)
                        if p: p.set_alpha(0.5); p.set_linewidth(1.2)
        
        elif n == 4:
            # 4-way Venn with custom ellipse layout
            names = [name for name, _ in sets_named]
            sets = [s for _, s in sets_named]
            draw_venn4(ax, names, sets, number_fontsize=number_fontsize, label_fontsize=label_fontsize)
        
        elif n > 4:
            ax.axis("off")
            ax.text(0.5, 0.5,
                    f"Venn diagrams support 2–4 sets.\nReceived {n} sets for '{title}'.\n\n"
                    f"For {n} sets, consider using UpSet plots instead.",
                    ha="center", va="center", fontsize=12, wrap=True)
            plt.suptitle(title, y=0.98, fontsize=12, fontweight="bold")
            plt.tight_layout()
            save_kwargs: Dict[str, Any] = {'bbox_inches': 'tight', 'format': output_format}
            if output_format == 'png':
                save_kwargs['dpi'] = dpi
            plt.savefig(out_png, **save_kwargs)
            plt.close()
            return out_png

        plt.suptitle(title, y=0.98, fontsize=13, fontweight="bold")
        plt.tight_layout()
        save_kwargs: Dict[str, Any] = {'bbox_inches': 'tight', 'format': output_format}
        if output_format == 'png':
            save_kwargs['dpi'] = dpi
        plt.savefig(out_png, **save_kwargs)
        plt.close()
        return out_png

    except Exception as e:
        print(f"Error plotting Venn: {e}")
        import traceback
        traceback.print_exc()
        try:
            plt.close()
        except Exception:
            pass
        return None

def write_venn_excel(sets_named: List[Tuple[str, Set[str]]], excel_path: str) -> Optional[str]:
    """Write Excel with sheets: summary, overlap_all, union_all, unique_<set>, pairwise overlaps."""
    try:
        with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
            summary_rows = [{"Comparison": name, "Count": len(s)} for name, s in sets_named]
            if sets_named:
                all_sets = [s for _, s in sets_named]
                union_all = set().union(*all_sets) if all_sets else set()
                inter_all = set(all_sets[0]) if len(all_sets) == 1 else set.intersection(*all_sets)
            else:
                union_all, inter_all = set(), set()

            pd.DataFrame(summary_rows + [
                {"Comparison": "UNION (all)", "Count": len(union_all)},
                {"Comparison": "INTERSECTION (all)", "Count": len(inter_all)},
            ]).to_excel(writer, sheet_name="summary", index=False)

            pd.DataFrame({"Molecule": sorted(inter_all)}).to_excel(writer, sheet_name="overlap_all", index=False)
            pd.DataFrame({"Molecule": sorted(union_all)}).to_excel(writer, sheet_name="union_all", index=False)

            for i, (name_i, s_i) in enumerate(sets_named):
                others = set().union(*[s for j, (_, s) in enumerate(sets_named) if j != i]) if len(sets_named) > 1 else set()
                uniq = s_i - others
                safe_name = name_i.replace("_vs_", "_").replace(" ", "_")[:25]
                sheet_name = f"unique_{safe_name}" if len(sets_named) > 1 else f"unique_{i+1}"
                pd.DataFrame({"Molecule": sorted(uniq)}).to_excel(writer, sheet_name=sheet_name, index=False)

            if len(sets_named) > 2:
                for (i, (n1, s1)), (j, (n2, s2)) in combinations(list(enumerate(sets_named)), 2):
                    inter = s1 & s2
                    # Use full comparison names for clarity, sanitize for sheet name
                    sheet_name = f"overlap_{n1}_and_{n2}"
                    # Clean sheet name: replace special characters and limit to 31 chars (Excel limit)
                    sheet_name = sheet_name.replace("_vs_", "-").replace(" ", "_").replace("|", "-")[:31]
                    pd.DataFrame({"Molecule": sorted(inter)}).to_excel(writer, sheet_name=sheet_name, index=False)
        return excel_path
    except Exception as e:
        print(f"[ERROR] write_venn_excel: {e}")
        return None


def _excel_symbol_for_index(idx: int) -> str:
    """Return Excel-style symbols: A..Z, AA..AZ, BA... for set indices."""
    n = idx + 1
    chars = []
    while n > 0:
        n, rem = divmod(n - 1, 26)
        chars.append(chr(ord("A") + rem))
    return "".join(reversed(chars))


def _excel_safe_sheet_name(base_name: str, used_names: Set[str]) -> str:
    """Sanitize and deduplicate sheet name within Excel's 31-char limit."""
    cleaned = (
        base_name.replace("[", "_")
        .replace("]", "_")
        .replace(":", "_")
        .replace("*", "_")
        .replace("?", "_")
        .replace("/", "_")
        .replace("\\", "_")
    )
    cleaned = cleaned.strip("'") or "Sheet"
    candidate = cleaned[:31]
    if candidate not in used_names:
        used_names.add(candidate)
        return candidate

    i = 2
    while True:
        suffix = f"_{i}"
        truncated = cleaned[: max(1, 31 - len(suffix))]
        candidate = f"{truncated}{suffix}"
        if candidate not in used_names:
            used_names.add(candidate)
            return candidate
        i += 1


def _exclusive_pair_overlap(all_sets: List[Set[str]], i: int, j: int) -> Set[str]:
    """Return region in i∩j excluding any other set (Venn pair-only overlap)."""
    pair_inter = all_sets[i] & all_sets[j]
    if len(all_sets) <= 2:
        return pair_inter
    others_union = set().union(*[s for k, s in enumerate(all_sets) if k not in (i, j)])
    return pair_inter - others_union

def write_venn_excel_consolidated(all_venn_data: List[Tuple[str, List[Tuple[str, Set[str]]]]], excel_path: str) -> Optional[str]:
    """
    Write consolidated Excel with all Venn diagrams as separate sheets.
    
    Args:
        all_venn_data: List of (venn_name, sets_named) tuples
        excel_path: Output Excel file path
        
    Returns:
        Path to Excel file if successful, None otherwise
        
    Sheet naming pattern (compact):
        - V{n}_desc: Symbol map (A/B/C...) to original circle labels
        - V{n}_summary: Counts summary
        - V{n}_A_B_union: Pairwise union
        - V{n}_A_B_ovlp: Pairwise exclusive overlap region
        - V{n}_A_unique: Unique to one circle
        - V{n}_A_B_C_union / V{n}_A_B_C_common: All-set union/common
    """
    try:
        with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
            used_sheet_names: Set[str] = set()

            for venn_idx, (venn_name, sets_named) in enumerate(all_venn_data, start=1):
                block_prefix = f"V{venn_idx}"

                symbols = [_excel_symbol_for_index(i) for i in range(len(sets_named))]
                all_sets = [s for _, s in sets_named]

                if sets_named:
                    union_all = set().union(*all_sets)
                    inter_all = set(all_sets[0]) if len(all_sets) == 1 else set.intersection(*all_sets)
                else:
                    union_all, inter_all = set(), set()

                # Description sheet: maps short symbols to full circle labels.
                desc_rows = [
                    {
                        "Venn": venn_name,
                        "Symbol": sym,
                        "Circle": circle_name,
                        "Count": len(circle_set),
                    }
                    for sym, (circle_name, circle_set) in zip(symbols, sets_named)
                ]
                desc_df = pd.DataFrame(desc_rows)
                desc_sheet = _excel_safe_sheet_name(f"{block_prefix}_desc", used_sheet_names)
                desc_df.to_excel(writer, sheet_name=desc_sheet, index=False)

                # Summary sheet with key region counts.
                summary_rows = [
                    {"Metric": "Venn Name", "Value": venn_name},
                    {"Metric": "Set Count", "Value": len(sets_named)},
                    {"Metric": f"{'_'.join(symbols)}_union", "Value": len(union_all)} if symbols else {"Metric": "union", "Value": 0},
                    {"Metric": f"{'_'.join(symbols)}_common", "Value": len(inter_all)} if symbols else {"Metric": "common", "Value": 0},
                ]
                for sym, (_, s) in zip(symbols, sets_named):
                    summary_rows.append({"Metric": f"{sym}_count", "Value": len(s)})
                summary_df = pd.DataFrame(summary_rows)
                summary_sheet = _excel_safe_sheet_name(f"{block_prefix}_summary", used_sheet_names)
                summary_df.to_excel(writer, sheet_name=summary_sheet, index=False)

                # All-set union/common sheets (short names, symbol-based).
                if symbols:
                    symbol_join = "_".join(symbols)
                    union_sheet = _excel_safe_sheet_name(f"{block_prefix}_{symbol_join}_union", used_sheet_names)
                    common_sheet = _excel_safe_sheet_name(f"{block_prefix}_{symbol_join}_common", used_sheet_names)
                    pd.DataFrame({"Molecule": sorted(union_all)}).to_excel(writer, sheet_name=union_sheet, index=False)
                    pd.DataFrame({"Molecule": sorted(inter_all)}).to_excel(writer, sheet_name=common_sheet, index=False)

                # Unique per set: A_unique, B_unique, ...
                for i, (sym, (_, s_i)) in enumerate(zip(symbols, sets_named)):
                    others = set().union(*[s for j, s in enumerate(all_sets) if j != i]) if len(all_sets) > 1 else set()
                    uniq = s_i - others
                    uniq_sheet = _excel_safe_sheet_name(f"{block_prefix}_{sym}_unique", used_sheet_names)
                    pd.DataFrame({"Molecule": sorted(uniq)}).to_excel(writer, sheet_name=uniq_sheet, index=False)

                # Pairwise sheets:
                # - A_B_union = union of A and B
                # - A_B_ovlp = exclusive pairwise overlap region (excludes all-common with other sets)
                if len(all_sets) >= 2:
                    for i, j in combinations(range(len(all_sets)), 2):
                        sym_i, sym_j = symbols[i], symbols[j]
                        pair_union = all_sets[i] | all_sets[j]
                        pair_ovlp_exclusive = _exclusive_pair_overlap(all_sets, i, j)

                        union_sheet = _excel_safe_sheet_name(f"{block_prefix}_{sym_i}_{sym_j}_union", used_sheet_names)
                        ovlp_sheet = _excel_safe_sheet_name(f"{block_prefix}_{sym_i}_{sym_j}_ovlp", used_sheet_names)

                        pd.DataFrame({"Molecule": sorted(pair_union)}).to_excel(writer, sheet_name=union_sheet, index=False)
                        pd.DataFrame({"Molecule": sorted(pair_ovlp_exclusive)}).to_excel(writer, sheet_name=ovlp_sheet, index=False)

        return excel_path
    except Exception as e:
        print(f"[ERROR] write_venn_excel_consolidated: {e}")
        return None

def run_venn_analysis(ctx: CommonVizContext, params: VennParams) -> VizResults:
    """
    Generate Venn diagrams per user-defined specs based on pairwise statistics or presence/absence.
    
    Creates PNG files for each Venn and a single consolidated Excel file with multiple sheets.
    """
    files_created: List[str] = []
    errors: List[str] = []
    venn_summaries: List[str] = []

    try:
        venn_dir = ctx.output_dir if os.path.basename(ctx.output_dir.rstrip(os.sep)).lower() == "venn" else os.path.join(ctx.output_dir, "venn")
        os.makedirs(venn_dir, exist_ok=True)

        prefer_adj = params.prefer_adj_p
        if hasattr(ctx, "use_adj_p") and ctx.use_adj_p is not None:
            prefer_adj = bool(ctx.use_adj_p)

        # Collect all sets for consolidated Excel export
        all_venn_data = []  # List of (spec_name, sets_named) tuples
        
        for idx, spec in enumerate(params.venn_specs):
            sets_named: List[Tuple[str, Set[str]]] = []
            set_sizes: List[int] = []

            if params.skip_all_cutoffs:
                # Use presence/absence logic - extract groups from comparisons
                # For each comparison (g1, g2), create sets for g1 and g2 individually
                groups_seen = set()
                for g1, g2 in spec.comparisons:
                    groups_seen.add(g1)
                    groups_seen.add(g2)
                
                # Determine ID column for reporting
                id_col = ctx.id_column if ctx.id_column else ('metabolite_id' if 'metabolite_id' in ctx.complete_df.columns else 'Name')
                total_rows = len(ctx.complete_df)
                total_unique = ctx.complete_df[id_col].nunique() if id_col and id_col in ctx.complete_df.columns else total_rows
                
                # Log clear summary of what Venn received
                logger.info(f"  ════════════════════════════════════════════════════════════")
                if total_rows != total_unique:
                    logger.info(f"  VENN INPUT: {total_rows} rows → {total_unique} unique features (ID column: '{id_col}')")
                else:
                    logger.info(f"  VENN INPUT: Received {total_unique} features (ID column: '{id_col}')")
                
                # Describe threshold being used
                if params.min_presence_type == 'count':
                    logger.info(f"  THRESHOLD: Absolute count ≥ {params.min_presence_count} samples (user-selected)")
                else:
                    logger.info(f"  THRESHOLD: Percentage ≥ {params.min_presence_percent}% of samples (user-selected)")
                logger.info(f"  ────────────────────────────────────────────────────────────")
                
                # Generate presence sets for each unique group
                for group in sorted(groups_seen):
                    label, ids = _present_ids_for_group(
                        ctx.complete_df,
                        group,
                        ctx.sample_cols,
                        ctx.sample_to_group,
                        params.min_presence_type,
                        params.min_presence_count,
                        params.min_presence_percent,
                        id_column=ctx.id_column
                    )
                    # Count samples in this group for debug clarity
                    group_sample_count = len([col for col in ctx.sample_cols if ctx.sample_to_group.get(col) == group])
                    
                    # Calculate threshold based on type
                    if params.min_presence_type == 'percentage':
                        threshold = int(np.ceil(group_sample_count * params.min_presence_percent / 100.0))
                        logger.info(f"  {group}: n={group_sample_count} samples, {params.min_presence_percent}% = ≥{threshold} → keeps {len(ids)} features")
                    else:
                        threshold = min(params.min_presence_count, group_sample_count)
                        logger.info(f"  {group}: n={group_sample_count} samples, threshold ≥{threshold} → keeps {len(ids)} features")
                    
                    sets_named.append((label, ids))
                    set_sizes.append(len(ids))
            else:
                # Use statistical significance logic (standard)
                for g1, g2 in spec.comparisons:
                    label, ids = _significant_ids_for_pair(
                        ctx.complete_df,
                        g1,
                        g2,
                        p_thresh=params.p_threshold,
                        fc_thresh=params.fc_threshold,
                        prefer_adj=prefer_adj,
                        stat_column_assignments=ctx.stat_column_assignments,
                        id_column=ctx.id_column
                    )
                    sig_label = label if label is not None else f"{g1}_vs_{g2}"

                    sets_named.append((sig_label, ids))
                    set_sizes.append(len(ids))

            set_info = ", ".join([f"{name}={len(s)}" for name, s in sets_named])

            if sets_named:
                all_sets = [s for _, s in sets_named]
                union_size = len(set().union(*all_sets)) if all_sets else 0
                intersection_size = len(all_sets[0]) if len(all_sets) == 1 else len(set.intersection(*all_sets))
                
                # Calculate total UNIQUE features (not rows) for accurate excluded count
                id_col = ctx.id_column if ctx.id_column else ('metabolite_id' if 'metabolite_id' in ctx.complete_df.columns else 'Name')
                total_unique_features = ctx.complete_df[id_col].nunique() if id_col and id_col in ctx.complete_df.columns else len(ctx.complete_df)
                total_rows = len(ctx.complete_df)
                
                venn_summary = f"{spec.name}: {len(sets_named)} sets ({set_info}), Union={union_size}, Intersection={intersection_size}"
                
                # Add clarifying debug info about presence threshold vs total analyzed
                if params.skip_all_cutoffs:
                    logger.info(f"  ────────────────────────────────────────────────────────────")
                    logger.info(f"  AFTER FILTERING: Union={union_size}, Intersection={intersection_size}")
                    
                    # Check for duplicate IDs in the data
                    if total_rows != total_unique_features:
                        logger.info(f"  ℹ️  Note: {total_rows} rows contain {total_unique_features} unique '{id_col}' values ({total_rows - total_unique_features} duplicates)")
                    
                    excluded = total_unique_features - union_size
                    if excluded > 0:
                        logger.info(f"  ⚠️  {excluded} features EXCLUDED (below threshold in ALL groups)")
                        
                        # Show detailed counts for excluded features
                        logger.info(f"  ")
                        logger.info(f"  📋 EXCLUDED FEATURES DETAIL:")
                        logger.info(f"  ────────────────────────────────────────────────────────────")
                        
                        # Find excluded IDs
                        all_present_ids = set().union(*all_sets) if all_sets else set()
                        
                        # Get all UNIQUE IDs from dataframe as strings
                        all_df_ids = {
                            norm
                            for norm in pd.Series(ctx.complete_df[id_col].unique()).map(_normalize_feature_id_for_venn)
                            if norm is not None
                        }
                        excluded_ids = all_df_ids - all_present_ids
                        
                        if not excluded_ids:
                            logger.info(f"    ⚠️ No excluded IDs found (calculation mismatch)")
                        else:
                            # Get counts for each excluded feature in each group  
                            excluded_details = []
                            normalized_col = ctx.complete_df[id_col].map(_normalize_feature_id_for_venn)
                            for feat_id in excluded_ids:
                                # More robust row matching
                                mask = normalized_col == feat_id
                                if not mask.any():
                                    continue
                                feat_row = ctx.complete_df[mask].iloc[0]
                                
                                group_counts = {}
                                for group in sorted(groups_seen):
                                    group_samples = [col for col in ctx.sample_cols if ctx.sample_to_group.get(col) == group]
                                    non_zero_count = 0
                                    for col in group_samples:
                                        try:
                                            val = pd.to_numeric(feat_row[col], errors='coerce')
                                            if not pd.isna(val) and val != 0:
                                                non_zero_count += 1
                                        except Exception:
                                            pass
                                    group_counts[group] = non_zero_count
                                
                                excluded_details.append((feat_id, group_counts))
                            
                            # Show first 20 excluded features
                            show_limit = 20
                            for i, (feat_id, group_counts) in enumerate(sorted(excluded_details)[:show_limit]):
                                counts_str = ", ".join([f"{g}={n}" for g, n in group_counts.items()])
                                logger.info(f"    {feat_id}: {counts_str}")
                            
                            if len(excluded_details) > show_limit:
                                logger.info(f"    ... and {len(excluded_details) - show_limit} more")
                        
                        logger.info(f"  ────────────────────────────────────────────────────────────")
                    logger.info(f"  ════════════════════════════════════════════════════════════")
            else:
                venn_summary = f"{spec.name}: No sets generated"
            venn_summaries.append(venn_summary)

            # Store for consolidated Excel
            all_venn_data.append((spec.name, sets_named))
            
            # Generate Venn plot with appropriate format and extension
            safe_name = spec.name.replace(" ", "_")
            # Add suffix for All Molecules to prevent overwriting filtered Venns
            name_suffix = "_AllMolecules" if params.skip_all_cutoffs else ""
            file_ext = '.svg' if params.output_format == 'svg' else '.png'
            out_path = os.path.join(venn_dir, f"venn_{safe_name}{name_suffix}{file_ext}")
            plotted = plot_venn(
                sets_named,
                out_path,
                title=spec.name,
                width=params.fig_width,
                height=params.fig_height,
                dpi=params.fig_dpi,
                venn_index=idx,  # Pass index for alternating 2-way Venn labels
                number_fontsize=params.venn_number_fontsize,
                label_fontsize=params.venn_label_fontsize,
                output_format=params.output_format,
            )
            if plotted:
                files_created.append(plotted)
            else:
                errors.append(f"Failed to plot venn for '{spec.name}'")

        # Write single consolidated Excel file with all Venns as separate sheets
        if all_venn_data:
            suffix = "_AllMolecules" if params.skip_all_cutoffs else ""
            consolidated_xlsx = os.path.join(venn_dir, f"venn_all_comparisons{suffix}.xlsx")
            written = write_venn_excel_consolidated(all_venn_data, consolidated_xlsx)
            if written:
                files_created.append(written)
            else:
                errors.append("Failed to write consolidated Excel file")

        summary = f"Venn analysis complete: {len(params.venn_specs)} diagram(s), {len(files_created)} files generated"
        return VizResults(files_created=files_created, errors=errors, summary=summary, venn_summaries=venn_summaries)
    except Exception as e:
        errors.append(f"Venn analysis failed: {e}")
        return VizResults(files_created=files_created, errors=errors, summary="Venn analysis failed", venn_summaries=venn_summaries)

# Utility Functions
def detect_lipid_mode(df: pd.DataFrame) -> bool:
    """Detect if the data is in lipid mode based on presence of lipid-specific columns."""
    lipid_indicators = ['lipidid', 'class', 'lipidgroup', 'subclass']
    present = sum(1 for col in lipid_indicators if col in df.columns)
    return present >= 2  # At least 2 lipid columns present

def aggregate_by_lipid_class(df: pd.DataFrame, sample_cols: List[str], *, 
                              class_col: str = 'class',
                              original_data: pd.DataFrame = None) -> pd.DataFrame:
    """Aggregate lipid data by class (sum of sample values per class).
    
    Parameters
    ----------
    df : pd.DataFrame
        Complete lipid dataframe with individual lipid rows (may be imputed)
    sample_cols : List[str]
        Sample column names to aggregate
    class_col : str
        Column name containing lipid class (default: 'class')
    original_data : pd.DataFrame, optional
        Original data before imputation. Must have the same index as df.
        If provided, only values that were non-zero in the original data are 
        included in the aggregation. This prevents imputed values from creating 
        classes that don't exist in negative polarity or other data subsets.
        Completely excludes classes if all their lipids were imputed-only.
        
    Returns
    -------
    pd.DataFrame
        Aggregated dataframe with one row per lipid class.
        Classes with all-imputed lipids are excluded if original_data is provided.
    """
    if class_col not in df.columns:
        raise ValueError(f"Class column '{class_col}' not found in dataframe")
    
    # Ensure original_data has same index as df if provided
    if original_data is not None and not original_data.empty:
        if len(original_data) != len(df):
            raise ValueError(f"original_data must have same length as df: {len(original_data)} vs {len(df)}")
        original_data = original_data.reset_index(drop=True)
    
    # Group by class and compute sum for sample columns
    # Exclude zeros and NaNs before aggregating (matching stats behavior)
    def safe_sum(series):
        vals = pd.to_numeric(series, errors='coerce')
        vals = vals[~vals.isna() & (vals != 0)]
        return vals.sum() if len(vals) > 0 else np.nan
    
    def safe_sum_with_original(series):
        """Aggregation function that filters based on original data."""
        col_name = series.name
        idx = series.index
        
        # Get values from current series (may be imputed)
        vals = pd.to_numeric(series, errors='coerce')
        
        # Get original values for the same rows
        orig_vals = pd.to_numeric(original_data.loc[idx, col_name], errors='coerce')
        
        # Create mask: only include where original was non-zero/non-NaN
        # (This means we only include values that actually existed in the original data)
        mask = (orig_vals > 0) | (orig_vals.notna() & (orig_vals != 0))
        vals_filtered = vals[mask]
        
        # Also exclude zeros and NaNs from the current values
        vals_filtered = vals_filtered[~vals_filtered.isna() & (vals_filtered != 0)]
        
        return vals_filtered.sum() if len(vals_filtered) > 0 else np.nan
    
    # Choose aggregation function
    if original_data is not None and not original_data.empty:
        agg_dict = {col: safe_sum_with_original for col in sample_cols}
    else:
        agg_dict = {col: safe_sum for col in sample_cols}
    
    df_reset = df.reset_index(drop=True)
    class_df = df_reset.groupby(class_col, as_index=False).agg(agg_dict)
    
    # Add lipid_class_id column for identification
    class_df.insert(0, 'lipid_class_id', class_df[class_col])
    
    # Add count column
    class_counts = df_reset.groupby(class_col).size()
    class_df['n_lipids'] = class_df[class_col].map(class_counts)
    
    # If filtering by original data, remove classes that ended up with all NaN values
    # This prevents imputed-only classes from appearing in the results
    if original_data is not None and not original_data.empty:
        # Check which classes have at least one non-NaN value in the sample columns
        has_values = class_df[sample_cols].notna().any(axis=1)
        class_df = class_df[has_values].reset_index(drop=True)
    
    return class_df

def detect_groups_from_complete(df: pd.DataFrame) -> List[str]:
    """Detect group names from n_<Group> columns."""
    groups = []
    for col in df.columns:
        if col.startswith('n_') and len(col) > 2:
            g = col[2:]
            if g not in groups:
                groups.append(g)
    return groups

def identify_sample_columns(df: pd.DataFrame, groups: List[str], pattern_map: Dict[str, List[str]]) -> Tuple[List[str], Dict[str, str]]:
    """Identify sample columns and map them to groups."""
    sample_cols: List[str] = []
    sample_to_group: Dict[str, str] = {}

    def is_feature(col: str) -> bool:
        return col in FEATURE_COLUMNS_CANONICAL
    
    def is_stat(col: str) -> bool:
        lc = col.lower()
        # Check exact matches (case-insensitive)
        if lc in [s.lower() for s in STAT_EXCLUDE_EXACT]:
            return True
        # Check prefixes
        if any(col.startswith(p) for p in STAT_EXCLUDE_PREFIXES):
            return True
        # Check substring patterns
        if any(tok in lc for tok in [s.lower() for s in STAT_EXCLUDE_CONTAINS]):
            return True
        # Special handling: _FC pattern (with underscore before FC)
        # This matches "Control_vs_PD_FC" but NOT "FC1", "FC2"
        if '_fc' in lc or lc.endswith('_fc'):
            return True
        # Mean columns (e.g., Control_Mean, PD_Mean)
        if 'mean' in lc and ('_mean' in lc or lc.endswith('mean')):
            return True
        return False

    for col in df.columns:
        if is_feature(col) or is_stat(col):
            continue
        
        # Try to assign to a group based on patterns
        assigned_group = None
        for g, pats in pattern_map.items():
            for pat in pats:
                if pat and pat.lower() in col.lower():
                    assigned_group = g
                    break
            if assigned_group:
                break
        
        # Consider only numeric columns
        if assigned_group is not None:
            series = pd.to_numeric(df[col], errors='coerce')
            if series.notna().sum() > 0:  # at least one numeric value
                sample_cols.append(col)
                sample_to_group[col] = assigned_group
    
    return sample_cols, sample_to_group

def _get_feature_ids(df: pd.DataFrame) -> List[str]:
    """Extract feature identifiers from DataFrame."""
    # Priority list with Protein first, then other feature identifiers
    for candidate in ['Protein', 'Name', 'Gene', 'metabolite_id', 'lipid_class_id', 'lipidid', 'class', 'Compound', 'Feature_ID', 'Metabolite']:
        if candidate in df.columns:
            return df[candidate].astype(str).tolist()
    # If no standard identifiers found, use the first column instead of numeric indices
    return df.iloc[:, 0].astype(str).tolist() if not df.empty else []

def prepare_lipid_class_stats(lipid_class_df: pd.DataFrame, groups: List[str], 
                               sample_cols: List[str], sample_to_group: Dict[str, str]) -> pd.DataFrame:
    """Compute group means and counts for lipid class data.
    
    This mimics the statistical summary structure expected by visualization functions.
    
    Parameters
    ----------
    lipid_class_df : pd.DataFrame
        Aggregated lipid class dataframe
    groups : List[str]
        List of group names
    sample_cols : List[str]
        Sample column names
    sample_to_group : Dict[str, str]
        Mapping of sample columns to groups
        
    Returns
    -------
    pd.DataFrame
        Enhanced dataframe with mean_<group> and n_<group> columns
    """
    result = lipid_class_df.copy()
    
    # Build group-wise column lists
    group_cols = {g: [c for c in sample_cols if sample_to_group.get(c) == g] for g in groups}
    
    # Compute means and counts per group
    for g, cols in group_cols.items():
        if not cols:
            result[f'mean_{g}'] = np.nan
            result[f'n_{g}'] = 0
            continue
        
        # Extract values and compute mean (excluding zeros and NaNs)
        means = []
        counts = []
        for idx, row in result.iterrows():
            vals = row[cols].apply(pd.to_numeric, errors='coerce').values
            vals = vals[~np.isnan(vals) & (vals != 0)]
            means.append(vals.mean() if len(vals) > 0 else np.nan)
            counts.append(len(vals))
        
        result[f'mean_{g}'] = means
        result[f'n_{g}'] = counts
    
    return result
    """Extract feature identifiers from DataFrame."""
    for candidate in ['metabolite_id','Name','Compound','Feature_ID','Metabolite']:
        if candidate in df.columns:
            return df[candidate].astype(str).tolist()
    return [str(i) for i in range(len(df))]

def _style_axes(ax):
    """Apply consistent styling to matplotlib axes."""
    for spine in ['top','right']:
        if spine in ax.spines:
            ax.spines[spine].set_visible(False)
    for spine in ['left','bottom']:
        if spine in ax.spines:
            ax.spines[spine].set_linewidth(2.2)
    ax.tick_params(axis='both', labelsize=11, width=1.5)
    for lbl in ax.get_xticklabels()+ax.get_yticklabels():
        lbl.set_fontweight('bold')

def _confidence_ellipse(ax, x, y, color, alpha=0.15):
    """Add 95% confidence ellipse to 2D plot."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if np.linalg.det(cov) <= 0:
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    chi2_val = chi2.ppf(0.95, df=2)
    width, height = 2 * np.sqrt(vals * chi2_val)
    angle = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    mean_x, mean_y = np.mean(x), np.mean(y)
    ell = Ellipse((mean_x, mean_y), width=width, height=height, angle=angle,
                  facecolor=color, edgecolor='none', alpha=alpha, zorder=0)
    ax.add_patch(ell)

# PCA Functions
def run_pca(df: pd.DataFrame, sample_cols: List[str], sample_to_group: Dict[str, str], *, n_components: int = 2) -> PCAResult:
    """Run PCA analysis on sample data."""
    if not sample_cols:
        raise ValueError("No sample columns provided for PCA")
    
    sub = df[sample_cols].apply(pd.to_numeric, errors='coerce')
    sub_imputed = sub.copy()
    col_medians = sub_imputed.median(axis=0, skipna=True)
    sub_imputed = sub_imputed.fillna(col_medians)
    
    scaler = StandardScaler()
    X = scaler.fit_transform(sub_imputed.T)
    
    max_comp = min(n_components, min(X.shape))
    if max_comp < 2:
        raise ValueError("Insufficient samples/features for PCA")
    
    pca = PCA(n_components=max_comp, random_state=42)
    comps = pca.fit_transform(X)
    
    samples = sub_imputed.columns.tolist()
    groups = [sample_to_group.get(s, 'NA') for s in samples]
    
    scores_dict = {
        'sample': samples,
        'group': groups,
        'PC1': comps[:,0],
        'PC2': comps[:,1]
    }
    if max_comp >= 3:
        scores_dict['PC3'] = comps[:,2]
    
    scores = pd.DataFrame(scores_dict)
    feat_ids = _get_feature_ids(df)
    
    return PCAResult(
        scores=scores, 
        explained_variance=(pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]),
        pca_model=pca, 
        feature_ids=feat_ids
    )

def plot_pca(result: PCAResult, title: str, output_path: str, *, color_map: Optional[Dict[str, Any]] = None, group_order: Optional[List[str]] = None, fig_size: Tuple[float,float]=(8,6), dpi: int = 220, point_size: float = 35, xlabel_fontsize: int = 11, ylabel_fontsize: int = 11, title_fontsize: int = 12, tick_fontsize: int = 11, legend_fontsize: int = 10, show_legend: bool = True):
    """Generate 2D PCA scatter plot preserving supplied group order.

    If group_order is provided we keep that order (filtering missing groups). Otherwise we
    preserve first appearance order (no alphabetical re-sorting) to match GUI baseline intent.
    """
    scores = result.scores
    var1, var2 = result.explained_variance

    plt.figure(figsize=fig_size)
    ax = plt.gca()
    try:
        ax.set_facecolor('white')
        ax.figure.set_facecolor('white')
    except Exception:
        pass

    if group_order:
        groups_seq = [g for g in group_order if g in scores['group'].unique()]
    else:
        groups_seq = list(dict.fromkeys(scores['group'].tolist()))  # first appearance order

    resolved_color_map: Dict[str, Any]
    if color_map is None:
        palette = sns.color_palette('tab10')
        resolved_color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups_seq)}
    else:
        resolved_color_map = color_map

    for g in groups_seq:
        sub = scores[scores['group']==g]
        if sub.empty:
            continue
        ax.scatter(sub['PC1'], sub['PC2'], s=point_size, edgecolor='k', linewidth=0.4,
               color=resolved_color_map.get(g, '#999999'), label=f"{g} (n={len(sub)})")
        _confidence_ellipse(ax, sub['PC1'].values, sub['PC2'].values, resolved_color_map.get(g, '#999999'))

    ax.set_xlabel(f'PC1 ({var1*100:.2f}%)', fontweight='bold', fontsize=xlabel_fontsize)
    ax.set_ylabel(f'PC2 ({var2*100:.2f}%)', fontweight='bold', fontsize=ylabel_fontsize)
    ax.set_title(title, fontweight='bold', fontsize=title_fontsize)
    # Explicitly ensure axis label font sizes are applied
    ax.xaxis.label.set_fontsize(xlabel_fontsize)
    ax.yaxis.label.set_fontsize(ylabel_fontsize)
    ax.tick_params(axis='both', labelsize=tick_fontsize, width=1.5)
    try:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    except Exception:
        pass
    for spine in ['top','right']:
        if spine in ax.spines:
            ax.spines[spine].set_visible(False)
    for spine in ['left','bottom']:
        if spine in ax.spines:
            ax.spines[spine].set_linewidth(2.2)
    for lbl in ax.get_xticklabels()+ax.get_yticklabels():
        try:
            lbl.set_fontweight('bold')
        except Exception:
            pass
    
    # Handle legend display based on show_legend parameter
    if show_legend:
        # Place legend outside the plot to the right to avoid overlap with data points
        try:
            leg = ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=legend_fontsize)
            for t in leg.get_texts():
                t.set_fontweight('bold')
                t.set_fontsize(legend_fontsize)  # Ensure fontsize is applied to each text element
            # Make room on the right for the legend
            plt.tight_layout(rect=(0, 0, 0.85, 1.0))
        except Exception:
            try:
                leg = ax.legend(frameon=False, loc='upper right', fontsize=legend_fontsize)
                for t in leg.get_texts():
                    t.set_fontweight('bold')
                    t.set_fontsize(legend_fontsize)  # Ensure fontsize is applied to each text element
                plt.tight_layout()
            except Exception:
                plt.tight_layout()
    else:
        # No legend - use full plot area
        plt.tight_layout()
        
    try:
        plt.savefig(output_path, dpi=dpi)
        return True
    except Exception as e:
        print(f"⚠️ Failed saving plot {output_path}: {e}")
        return False
    finally:
        plt.close()

def _confidence_ellipsoid_3d(ax, data: np.ndarray, color, alpha=0.12):
    if data.shape[0] < 4:
        return
    cov = np.cov(data.T)
    if np.linalg.det(cov) <= 0:
        return
    center = data.mean(axis=0)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    chi2_val = chi2.ppf(0.95, df=3)
    radii = np.sqrt(vals * chi2_val)
    u = np.linspace(0, 2*np.pi, 30)
    v = np.linspace(0, np.pi, 22)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    sphere = np.stack([xs, ys, zs], axis=-1)
    ellipsoid = sphere @ np.diag(radii) @ vecs.T + center
    ax.plot_surface(ellipsoid[...,0], ellipsoid[...,1], ellipsoid[...,2], color=color, alpha=alpha, linewidth=0, rstride=1, cstride=1)

def plot_pca_3d(result: PCAResult, title: str, output_path: str, *, color_map: Optional[Dict[str, Any]] = None, group_order: Optional[List[str]] = None, fig_size: Tuple[float,float]=(9,8), dpi: int = 230, point_size: float = 30, xlabel_fontsize: int = 11, ylabel_fontsize: int = 11, zlabel_fontsize: Optional[int] = None, title_fontsize: int = 12, tick_fontsize: int = 11, legend_fontsize: int = 10, show_legend: bool = True, view_azim: float = -60, view_elev: float = 30):
    scores = result.scores
    if 'PC3' not in scores.columns:
        return False
    
    # Use ylabel_fontsize for zlabel if not explicitly provided (backward compatibility)
    if zlabel_fontsize is None:
        zlabel_fontsize = ylabel_fontsize
    
    pca_model = result.pca_model
    if pca_model is not None:
        var = list(pca_model.explained_variance_ratio_[:3]) + [np.nan]*(3-len(pca_model.explained_variance_ratio_[:3]))
    else:
        var = [np.nan, np.nan, np.nan]
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111, projection='3d')

    # Standardize to white background for export-ready figures
    try:
        fig.patch.set_facecolor('white')
    except Exception:
        pass
    try:
        ax.set_facecolor('white')
    except Exception:
        pass
    try:
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):  # type: ignore[attr-defined]
            pane.set_facecolor('white')
            pane.set_edgecolor('#d0d0d0')
    except Exception:
        # Older Matplotlib uses w_*axis helpers
        w_xaxis = getattr(ax, 'w_xaxis', None)
        w_yaxis = getattr(ax, 'w_yaxis', None)
        w_zaxis = getattr(ax, 'w_zaxis', None)
        if w_xaxis is not None:
            try:
                w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
                if w_yaxis is not None:
                    w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
                if w_zaxis is not None:
                    w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
            except Exception:
                pass

    # Light-touch grid styling keeps focus on points/ellipsoids
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis._axinfo['grid']['linewidth'] = 0.35  # type: ignore[attr-defined]
            axis._axinfo['grid']['color'] = (0.88, 0.88, 0.88, 0.7)  # type: ignore[attr-defined]
        except Exception:
            pass

    try:
        ax.set_box_aspect((1.0, 1.0, 1.0))
    except Exception:
        pass
    groups_seq = ([g for g in group_order if g in scores['group'].unique()] if group_order
                  else sorted(scores['group'].unique()))
    resolved_color_map: Dict[str, Any]
    if color_map is None:
        palette = sns.color_palette('tab10')
        resolved_color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups_seq)}
    else:
        resolved_color_map = color_map
    for g in groups_seq:
        sub = scores[scores['group']==g]
        if sub.empty:
            continue
        zs_values: Any = sub['PC3'].to_numpy(dtype=float)
        ax_any: Any = ax
        ax_any.scatter(
            sub['PC1'].to_numpy(dtype=float),
            sub['PC2'].to_numpy(dtype=float),
            zs=zs_values,
            s=int(point_size),
            depthshade=False,
            color=resolved_color_map.get(g,'#999999'),
            edgecolor='k',
            linewidth=0.3,
            alpha=0.95,
            label=f"{g} (n={len(sub)})",
        )  # pyright: ignore[reportArgumentType]
        _confidence_ellipsoid_3d(ax, sub[['PC1','PC2','PC3']].values, resolved_color_map.get(g,'#999999'))
    
    # Set axis labels OUTSIDE the loop (so they're only set once)
    # Further tighten labelpad so text hugs the plot box per user feedback
    ax.set_xlabel(f"PC1 ({var[0]*100:.2f}%)" if not np.isnan(var[0]) else 'PC1', fontweight='bold', fontsize=xlabel_fontsize, labelpad=2)
    ax.set_ylabel(f"PC2 ({var[1]*100:.2f}%)" if not np.isnan(var[1]) else 'PC2', fontweight='bold', fontsize=ylabel_fontsize, labelpad=4)
    
    # For PC3 (z-axis), use set_zlabel with proper positioning
    pc3_label = f"PC3 ({var[2]*100:.2f}%)" if not np.isnan(var[2]) else 'PC3'
    ax.set_zlabel(pc3_label, fontweight='bold', fontsize=zlabel_fontsize, labelpad=4)
    
    # Explicitly ensure axis label font sizes are applied and centered along the axes
    ax.xaxis.label.set_fontsize(xlabel_fontsize)
    ax.yaxis.label.set_fontsize(ylabel_fontsize)
    ax.zaxis.label.set_fontsize(zlabel_fontsize)
    for label in (ax.xaxis.label, ax.yaxis.label, ax.zaxis.label):
        try:
            label.set_horizontalalignment('center')
            label.set_verticalalignment('center')
        except Exception:
            pass
    try:
        ax.xaxis.set_label_coords(0.5, -0.03)
        ax.yaxis.set_label_coords(-0.06, 0.5)
        ax.zaxis.set_label_coords(0.08, 0.5)
    except Exception:
        pass
    ax.set_title(title + "\n(95% Confidence Ellipsoids)", fontweight='bold', pad=10, fontsize=title_fontsize)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    ax.tick_params(axis='z', labelsize=tick_fontsize)  # pyright: ignore[reportArgumentType]
    try:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.zaxis.set_major_locator(MaxNLocator(nbins=5))
    except Exception:
        pass
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        for lbl in axis.get_ticklabels():
            try:
                lbl.set_fontweight('bold')
            except Exception:
                pass

    # Provide consistent padding so labels stay visible when exporting
    # Increase margins to ensure all axis labels are visible
    try:
        fig.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.92)
    except Exception:
        pass
    
    # Set viewing angle for 3D plot (user-configurable for optimal group separation)
    ax.view_init(elev=view_elev, azim=view_azim)
        
    if show_legend:
        try:
            box = ax.get_position()
            ax.set_position((box.x0, box.y0, box.width * 0.72, box.height))
        except Exception:
            pass
        try:
            leg = ax.legend(
                loc='center left',
                bbox_to_anchor=(1.08, 0.99),
                frameon=False,
                borderaxespad=0,
                fontsize=legend_fontsize
            )
            try:
                for text in leg.get_texts():
                    text.set_fontsize(legend_fontsize)
                    text.set_fontweight('bold')
                if leg.get_title() is not None:
                    leg.get_title().set_fontweight('bold')
            except Exception:
                pass
        except Exception:
            try:
                ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.05, 0.3), fontsize=legend_fontsize)
            except Exception:
                pass
    else:
        # Expand plot back to full width when legend is hidden
        try:
            box = ax.get_position()
            ax.set_position((box.x0, box.y0, box.width, box.height))
        except Exception:
            pass

    try:
        # Save with bbox_inches='tight' and appropriate padding to prevent label cutoff
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', pad_inches=0.3)
        return True
    except Exception as e:
        print(f"⚠️ Failed saving 3D PCA plot {output_path}: {e}")
        return False
    finally:
        plt.close()

def save_interactive_pca_3d(result: PCAResult, title: str, out_path: str, *, color_map: Optional[Dict[str, Any]] = None, point_size: float = 30):
    if px is None or go is None or 'PC3' not in result.scores.columns:
        return False
    try:
        scores = result.scores.copy()
        counts = scores['group'].value_counts()
        # Use group order from color_map if provided, else unique order
        group_order = list(color_map.keys()) if color_map else list(scores['group'].unique())
        scores['group_label'] = scores['group'].map(lambda g: f"{g} (n={counts[g]})")
        cmap = None
        if color_map:
            cmap = {f"{g} (n={counts[g]})": color_map[g] for g in group_order if g in counts and g in color_map}
        
        # Create the scatter plot
        fig = px.scatter_3d(
            scores,
            x='PC1', y='PC2', z='PC3',
            color='group_label',
            hover_data=['sample'],
            title=title + "<br>(95% Confidence Ellipsoids)",
            color_discrete_map=cmap,
            category_orders={'group_label': [f"{g} (n={counts[g]})" for g in group_order if g in counts]}
        )
        
        # Update scatter markers - match static 3D point size
        for tr in fig.data:
            tr_any: Any = tr
            if getattr(tr_any, 'type', None) == 'scatter3d':
                trace_update = getattr(tr_any, 'update', None)
                if callable(trace_update):
                    trace_update(marker=dict(size=point_size/5, line=dict(color='black', width=1)))
        
        # Add 95% confidence ellipsoids for each group
        for g in group_order:
            if g not in scores['group'].values:
                continue
            sub = scores[scores['group'] == g]
            if len(sub) < 4:
                continue
            
            # Get data for this group
            data = sub[['PC1', 'PC2', 'PC3']].values
            
            # Compute covariance and ellipsoid
            cov = np.cov(data.T)
            if np.linalg.det(cov) <= 0:
                continue
            
            center = data.mean(axis=0)
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            chi2_val = chi2.ppf(0.95, df=3)
            radii = np.sqrt(vals * chi2_val)
            
            # Generate ellipsoid surface
            u = np.linspace(0, 2*np.pi, 30)
            v = np.linspace(0, np.pi, 22)
            xs = np.outer(np.cos(u), np.sin(v))
            ys = np.outer(np.sin(u), np.sin(v))
            zs = np.outer(np.ones_like(u), np.cos(v))
            sphere = np.stack([xs, ys, zs], axis=-1)
            ellipsoid = sphere @ np.diag(radii) @ vecs.T + center
            
            # Get color for this group
            group_color = color_map.get(g, '#999999') if color_map else '#999999'
            
            # Convert color to rgba with transparency
            if isinstance(group_color, str):
                if group_color.startswith('#'):
                    # Convert hex to rgb
                    h = group_color.lstrip('#')
                    rgb = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
                    rgba = f'rgba({rgb[0]},{rgb[1]},{rgb[2]},0.15)'
                else:
                    rgba = group_color
            else:
                # Assume it's already rgb tuple
                rgba = f'rgba({int(group_color[0]*255)},{int(group_color[1]*255)},{int(group_color[2]*255)},0.15)'
            
            # Add ellipsoid as Mesh3d
            fig.add_trace(go.Mesh3d(
                x=ellipsoid[..., 0].flatten(),
                y=ellipsoid[..., 1].flatten(),
                z=ellipsoid[..., 2].flatten(),
                alphahull=0,
                opacity=0.15,
                color=group_color,
                name=f'{g} 95% CI',
                showlegend=False,
                hoverinfo='skip'
            ))
        
        # Move plotly legend to the right outside the scene - match static 3D positioning
        try:
            fig.update_layout(
                legend=dict(
                    x=1.05,
                    y=0.2,
                    xanchor='left',
                    yanchor='middle',
                    traceorder='normal',
                    # Choose a bold font family so legend appears bold
                    font=dict(family="Arial Black, Arial, sans-serif", size=11)
                ),
                # Match static plot appearance: white background, gridlines
                scene=dict(
                    bgcolor='white',
                    xaxis=dict(
                        backgroundcolor='white',
                        gridcolor='lightgray',
                        showbackground=True,
                        gridwidth=1
                    ),
                    yaxis=dict(
                        backgroundcolor='white',
                        gridcolor='lightgray',
                        showbackground=True,
                        gridwidth=1
                    ),
                    zaxis=dict(
                        backgroundcolor='white',
                        gridcolor='lightgray',
                        showbackground=True,
                        gridwidth=1
                    )
                ),
                # Set default camera angle to match static 3D plot (elev=20, azim=45)
                scene_camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.2)
                ),
                paper_bgcolor='white',
                plot_bgcolor='white'
            )
        except Exception:
            pass

        fig.write_html(out_path, include_plotlyjs='cdn')
        return True
    except Exception as e:
        print(f"⚠️ Failed interactive 3D PCA save {out_path}: {e}")
        return False

def plot_scree(result: PCAResult, title: str, output_path: str, *, fig_size: Tuple[float,float]=(8,5), dpi: int = 220,
              xlabel_fontsize: int = 11, ylabel_fontsize: int = 11, title_fontsize: int = 13, tick_fontsize: int = 10):
    """Generate scree plot showing explained variance per principal component."""
    if result.pca_model is None:
        return False
    
    try:
        n_components = len(result.pca_model.explained_variance_ratio_)
        variance_explained = result.pca_model.explained_variance_ratio_ * 100
        cumulative_variance = np.cumsum(variance_explained)
        
        fig, ax1 = plt.subplots(figsize=fig_size, dpi=dpi)
        
        # Bar plot for individual variance
        x = np.arange(1, n_components + 1)
        ax1.bar(x, variance_explained, alpha=0.7, color='steelblue', label='Individual')
        ax1.set_xlabel('Principal Component', fontsize=xlabel_fontsize, fontweight='bold')
        ax1.set_ylabel('Variance Explained (%)', fontsize=ylabel_fontsize, color='steelblue', fontweight='bold')
        # Explicitly ensure axis label font sizes are applied
        ax1.xaxis.label.set_fontsize(xlabel_fontsize)
        ax1.yaxis.label.set_fontsize(ylabel_fontsize)
        ax1.tick_params(axis='y', labelcolor='steelblue', labelsize=tick_fontsize)
        ax1.tick_params(axis='x', labelsize=tick_fontsize)
        ax1.set_xticks(x)
        
        # Line plot for cumulative variance on secondary y-axis
        ax2 = ax1.twinx()
        ax2.plot(x, cumulative_variance, color='darkorange', marker='o', linewidth=2, markersize=6, label='Cumulative')
        ax2.set_ylabel('Cumulative Variance (%)', fontsize=ylabel_fontsize, color='darkorange', fontweight='bold')
        # Explicitly ensure axis label font sizes are applied
        ax2.yaxis.label.set_fontsize(ylabel_fontsize)
        ax2.tick_params(axis='y', labelcolor='darkorange', labelsize=tick_fontsize)
        ax2.set_ylim(0, 105)
        
        # Add grid for readability
        ax1.grid(True, alpha=0.3, linestyle='--')
        
        # Title
        ax1.set_title(title, fontsize=title_fontsize, fontweight='bold', pad=10)
        
        # Add legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True, fancybox=True)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        plt.close()
        return True
    except Exception as e:
        print(f"⚠️ Failed saving scree plot {output_path}: {e}")
        plt.close()
        return False

def run_pca_analysis(ctx: CommonVizContext, params: PCAParams) -> VizResults:
    """Main PCA analysis entry point for GUI."""
    files_created = []
    errors = []
    
    try:
        # Defensive filtering: if in lipid mode, ensure sample_cols exclude canonical lipid feature metadata columns
        if ctx.is_lipid_mode and ctx.sample_cols:
            # Build normalized set for quick membership
            lipid_feat_norm = {str(c).strip().lower() for c in LIPID_FEATURE_CANONICAL}
            filtered = []
            removed = []
            for c in ctx.sample_cols:
                norm = str(c).strip().lower()
                if norm in lipid_feat_norm:
                    removed.append(c)
                    continue
                filtered.append(c)
            if len(filtered) != len(ctx.sample_cols) and len(filtered) >= 2:
                # Only replace if we still have enough columns for PCA
                ctx.sample_cols = filtered
                try:
                    # Attach note so GUI can optionally log later
                    ctx._lipid_feature_columns_removed = removed
                except Exception:
                    pass
            # Also drop from mapping any removed columns to prevent key errors later
            ctx.sample_to_group = {k: v for k, v in ctx.sample_to_group.items() if k in ctx.sample_cols}

        # Enforce preferred group order
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining
        # Create PCA output directory (avoid double 'pca/pca')
        if os.path.basename(ctx.output_dir.lower()) == 'pca':
            pca_dir = ctx.output_dir
        else:
            pca_dir = os.path.join(ctx.output_dir, 'pca')
        os.makedirs(pca_dir, exist_ok=True)

        def _pca3d_fig_size() -> Tuple[float, float]:
            """Expand GUI-configured dimensions for 3D plots by +2 inches each (per requirement)."""
            width = max(params.fig_width + 2.0, 1.0)
            height = max(params.fig_height + 2.0, 1.0)
            return (width, height)

        def _pca3d_font(value: int) -> int:
            """Use fonts 2pt smaller for 3D output, keeping a readable minimum."""
            return max(value - 2, 6)
        
        # Determine if we should skip "All Groups" PCA when we have only 2 groups
        # (since it would be identical to a single comparison)
        skip_all_groups_pca = len(ctx.groups) == 2
        
        # Generate default PCA with all groups (skip if only 2 groups to avoid duplication)
        if not skip_all_groups_pca:
            pca_all = run_pca(ctx.complete_df, ctx.sample_cols, ctx.sample_to_group, n_components=params.components)
            
            # Save scores for all groups (controlled by save_excel parameter)
            if params.save_excel:
                scores_path = os.path.join(pca_dir, 'pca_all_scores.csv')
                pca_all.scores.to_csv(scores_path, index=False)
                files_created.append(scores_path)
            
            # Generate 2D plot for all groups (controlled by save_2d parameter)
            if params.save_2d:
                plot_path = os.path.join(pca_dir, 'pca_all.png')
                if plot_pca(pca_all, 'PCA - All Groups', plot_path, color_map=ctx.color_map, group_order=ctx.groups, 
                           fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                           xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize, 
                           title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                           legend_fontsize=params.legend_fontsize, show_legend=params.show_legend):
                    files_created.append(plot_path)
            
            # Generate scree plot if requested (for all groups)
            if params.scree and pca_all.pca_model:
                scree_path = os.path.join(pca_dir, 'pca_all_scree.png')
                if plot_scree(pca_all, 'Scree Plot - All Groups', scree_path, 
                             fig_size=(params.fig_width, params.fig_height*0.7), dpi=params.fig_dpi,
                             xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize,
                             title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize):
                    files_created.append(scree_path)
        
        # Export loadings if requested (for all groups) - only if we generated pca_all
        if not skip_all_groups_pca and params.loadings and pca_all.pca_model is not None:
            try:
                # Get feature names (metabolites/lipids) - these are the ROWS of the DataFrame
                # The PCA components correspond to features (rows), not sample columns
                feature_ids = _get_feature_ids(ctx.complete_df)
                
                if len(feature_ids) > 0 and hasattr(pca_all.pca_model, 'components_'):
                    # Verify shapes match
                    n_features = pca_all.pca_model.components_.shape[1]
                    if len(feature_ids) != n_features:
                        # Mismatch - use generic feature IDs
                        feature_ids = [f'Feature_{i}' for i in range(n_features)]
                    
                    loadings_df = pd.DataFrame(
                        pca_all.pca_model.components_.T,
                        columns=[f'PC{i+1}' for i in range(pca_all.pca_model.n_components_)],
                        index=feature_ids[:n_features]  # Use only what we need
                    )
                    
                    # Calculate absolute contribution for each metabolite (sum of squared loadings)
                    loadings_df['Total_Contribution'] = np.sqrt((loadings_df ** 2).sum(axis=1))
                    
                    # Sort by total contribution and take top k
                    loadings_df = loadings_df.sort_values('Total_Contribution', ascending=False)
                    top_loadings = loadings_df.head(params.loadings_top_k)
                    
                    # Save to Excel - ONLY if loadings export is enabled
                    loadings_path = os.path.join(pca_dir, f'pca_all_loadings_top{params.loadings_top_k}.xlsx')
                    top_loadings.to_excel(loadings_path)
                    files_created.append(loadings_path)
            except Exception as e:
                errors.append(f"Failed to export PCA loadings: {str(e)}")
        
        # 3D plot if requested (for all groups) - controlled by save_3d parameter
        # Note: plot_3d controls interactive plots, save_3d controls static 3D PNG generation
        # Only generate if we created pca_all (i.e., more than 2 groups)
        if not skip_all_groups_pca and params.save_3d and pca_all.pca_model and pca_all.pca_model.n_components_ >= 3:
            pca_all_3d = pca_all if 'PC3' in pca_all.scores.columns else run_pca(ctx.complete_df, ctx.sample_cols, ctx.sample_to_group, n_components=max(3, params.components))
            plot_3d_path = os.path.join(pca_dir, 'pca_all_3d.png')
            if plot_pca_3d(pca_all_3d, 'PCA (3D) - All Groups', plot_3d_path, color_map=ctx.color_map, group_order=ctx.groups, 
                          fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                          xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                          title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                          legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                          view_azim=params.view_azim, view_elev=params.view_elev):
                files_created.append(plot_3d_path)
            if params.interactive_3d:
                html_path = os.path.join(pca_dir, 'pca_all_3d.html')
                if save_interactive_pca_3d(pca_all_3d, 'PCA (3D) - All Groups', html_path, color_map=ctx.color_map, point_size=params.point_size_3d):
                    files_created.append(html_path)
        
        # ADDITIONAL CUSTOM GROUP PCA if user selected specific groups
        # For 2 groups, we generate comparison PCA with clean naming (not "custom")
        if params.custom_groups and not skip_all_groups_pca:
            # Only process custom groups if we have MORE than 2 groups
            # Filter to only requested groups that exist
            custom_groups_to_plot = [g for g in params.custom_groups if g in ctx.groups]
            
            # Only generate if selection is different from all groups
            if custom_groups_to_plot and len(custom_groups_to_plot) != len(ctx.groups):
                custom_sample_cols = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) in custom_groups_to_plot]
                
                if len(custom_sample_cols) >= 3:
                    # Run custom PCA
                    pca_custom = run_pca(ctx.complete_df, custom_sample_cols, ctx.sample_to_group, n_components=params.components)
                    
                    # Build custom title
                    custom_title_suffix = ' - ' + ' vs '.join(custom_groups_to_plot)
                    custom_filename_suffix = '_' + '_vs_'.join(custom_groups_to_plot)
                    
                    # Save custom scores (controlled by save_excel parameter)
                    if params.save_excel:
                        custom_scores_path = os.path.join(pca_dir, f'pca_custom{custom_filename_suffix}_scores.csv')
                        pca_custom.scores.to_csv(custom_scores_path, index=False)
                        files_created.append(custom_scores_path)
                    
                    # Generate custom 2D plot (controlled by save_2d parameter)
                    if params.save_2d:
                        custom_plot_path = os.path.join(pca_dir, f'pca_custom{custom_filename_suffix}.png')
                        if plot_pca(pca_custom, f'PCA{custom_title_suffix}', custom_plot_path, color_map=ctx.color_map, group_order=custom_groups_to_plot, 
                                   fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                                   xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize, 
                                   title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                                   legend_fontsize=params.legend_fontsize, show_legend=params.show_legend):
                            files_created.append(custom_plot_path)
                    
                    # Custom 3D plot if requested (controlled by save_3d parameter)
                    if params.save_3d and pca_custom.pca_model and pca_custom.pca_model.n_components_ >= 3:
                        pca_custom_3d = pca_custom if 'PC3' in pca_custom.scores.columns else run_pca(ctx.complete_df, custom_sample_cols, ctx.sample_to_group, n_components=max(3, params.components))
                        custom_3d_path = os.path.join(pca_dir, f'pca_custom{custom_filename_suffix}_3d.png')
                        if plot_pca_3d(pca_custom_3d, f'PCA (3D){custom_title_suffix}', custom_3d_path, color_map=ctx.color_map, group_order=custom_groups_to_plot, 
                                      fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                                      xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                                      title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                                      legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                                      view_azim=params.view_azim, view_elev=params.view_elev):
                            files_created.append(custom_3d_path)
                        if params.interactive_3d:
                            custom_html_path = os.path.join(pca_dir, f'pca_custom{custom_filename_suffix}_3d.html')
                            if save_interactive_pca_3d(pca_custom_3d, f'PCA (3D){custom_title_suffix}', custom_html_path, color_map=ctx.color_map, point_size=params.point_size_3d):
                                files_created.append(custom_html_path)
                else:
                    errors.append(f'Insufficient samples for custom PCA with selected groups: {len(custom_sample_cols)}')
        
        # For 2-group case, generate a single comparison PCA with clean naming
        if skip_all_groups_pca and len(ctx.groups) == 2:
            pca_comparison = run_pca(ctx.complete_df, ctx.sample_cols, ctx.sample_to_group, n_components=params.components)
            
            # Build title with group names
            comparison_title = f'PCA - {" vs ".join(ctx.groups)}'
            comparison_filename = f'pca_{"_vs_".join(ctx.groups)}'
            
            # Save scores (controlled by save_excel parameter)
            if params.save_excel:
                scores_path = os.path.join(pca_dir, f'{comparison_filename}_scores.csv')
                pca_comparison.scores.to_csv(scores_path, index=False)
                files_created.append(scores_path)
            
            # Generate 2D plot (controlled by save_2d parameter)
            if params.save_2d:
                plot_path = os.path.join(pca_dir, f'{comparison_filename}.png')
                if plot_pca(pca_comparison, comparison_title, plot_path, color_map=ctx.color_map, group_order=ctx.groups, 
                           fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                           xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize, 
                           title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                           legend_fontsize=params.legend_fontsize, show_legend=params.show_legend):
                    files_created.append(plot_path)
            
            # 3D plot if requested (controlled by save_3d parameter)
            if params.save_3d and pca_comparison.pca_model and pca_comparison.pca_model.n_components_ >= 3:
                pca_comparison_3d = pca_comparison if 'PC3' in pca_comparison.scores.columns else run_pca(ctx.complete_df, ctx.sample_cols, ctx.sample_to_group, n_components=max(3, params.components))
                plot_3d_path = os.path.join(pca_dir, f'{comparison_filename}_3d.png')
                if plot_pca_3d(pca_comparison_3d, f'{comparison_title} (3D)', plot_3d_path, color_map=ctx.color_map, group_order=ctx.groups, 
                              fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                              xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                              title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                              legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                              view_azim=params.view_azim, view_elev=params.view_elev):
                    files_created.append(plot_3d_path)
                if params.interactive_3d:
                    html_path = os.path.join(pca_dir, f'{comparison_filename}_3d.html')
                    if save_interactive_pca_3d(pca_comparison_3d, f'{comparison_title} (3D)', html_path, color_map=ctx.color_map, point_size=params.point_size_3d):
                        files_created.append(html_path)
        
        # Lipid class PCA if requested and available
        if params.include_lipid_class and ctx.is_lipid_mode and ctx.lipid_class_df is not None:
            try:
                # Create dedicated subdirectory for class-level PCA outputs (improves organization)
                class_dir = os.path.join(pca_dir, params.class_subdir_name) if params.class_subdir_name else pca_dir
                os.makedirs(class_dir, exist_ok=True)

                # Determine if custom-only mode should be used (match behavior of normal PCA)
                custom_groups_to_plot = [g for g in (params.custom_groups or []) if g in ctx.groups]
                use_custom_only = bool(custom_groups_to_plot and len(custom_groups_to_plot) != len(ctx.groups))
                
                if not use_custom_only:
                    pca_class = run_pca(ctx.lipid_class_df, ctx.sample_cols, ctx.sample_to_group, n_components=params.components)

                # Save lipid class scores (controlled by save_excel parameter)
                if not use_custom_only and params.save_excel:
                    class_scores_path = os.path.join(class_dir, 'pca_lipid_class_scores.csv')
                    pca_class.scores.to_csv(class_scores_path, index=False)
                    files_created.append(class_scores_path)
                
                # Generate lipid class 2D plot (controlled by save_2d parameter)
                if not use_custom_only and params.save_2d:
                    class_plot_path = os.path.join(class_dir, 'pca_lipid_class.png')
                    if plot_pca(pca_class, 'PCA - Lipid Classes', class_plot_path, color_map=ctx.color_map, group_order=ctx.groups, 
                               fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                               xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize, 
                               title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                               legend_fontsize=params.legend_fontsize, show_legend=params.show_legend):
                        files_created.append(class_plot_path)
                
                # 3D lipid class plot if requested (controlled by save_3d parameter)
                if not use_custom_only and params.save_3d and pca_class.pca_model and pca_class.pca_model.n_components_ >= 3:
                    pca_class_3d = pca_class if 'PC3' in pca_class.scores.columns else run_pca(ctx.lipid_class_df, ctx.sample_cols, ctx.sample_to_group, n_components=max(3, params.components))
                    class_3d_path = os.path.join(class_dir, 'pca_lipid_class_3d.png')
                    if plot_pca_3d(pca_class_3d, 'PCA (3D) - Lipid Classes', class_3d_path, color_map=ctx.color_map, group_order=ctx.groups, 
                                  fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                                  xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                                  title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                                  legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                                  view_azim=params.view_azim, view_elev=params.view_elev):
                        files_created.append(class_3d_path)
                    if params.interactive_3d:
                        class_html_path = os.path.join(class_dir, 'pca_lipid_class_3d.html')
                        if save_interactive_pca_3d(pca_class_3d, 'PCA (3D) - Lipid Classes', class_html_path, color_map=ctx.color_map, point_size=params.point_size_3d):
                            files_created.append(class_html_path)
            except Exception as e:
                errors.append(f"Lipid class PCA failed: {str(e)}")
            
            # Respect custom group selection for lipid class as well
            try:
                custom_groups_to_plot = [g for g in (params.custom_groups or []) if g in ctx.groups]
                if custom_groups_to_plot and len(custom_groups_to_plot) != len(ctx.groups):
                    custom_sample_cols = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) in custom_groups_to_plot]
                    if len(custom_sample_cols) >= 3:
                        pca_class_custom = run_pca(ctx.lipid_class_df, custom_sample_cols, ctx.sample_to_group, n_components=params.components)
                        custom_title_suffix = ' - ' + ' vs '.join(custom_groups_to_plot)
                        custom_filename_suffix = '_' + '_vs_'.join(custom_groups_to_plot)
                        if params.save_excel:
                            class_scores_path = os.path.join(class_dir, f'pca_lipid_class_custom{custom_filename_suffix}_scores.csv')
                            pca_class_custom.scores.to_csv(class_scores_path, index=False)
                            files_created.append(class_scores_path)
                        if params.save_2d:
                            class_plot_path = os.path.join(class_dir, f'pca_lipid_class_custom{custom_filename_suffix}.png')
                            if plot_pca(pca_class_custom, f'PCA - Lipid Classes{custom_title_suffix}', class_plot_path, color_map=ctx.color_map, group_order=custom_groups_to_plot, 
                                       fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                                       xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize, 
                                       title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                                       legend_fontsize=params.legend_fontsize, show_legend=params.show_legend):
                                files_created.append(class_plot_path)
                        if params.save_3d and pca_class_custom.pca_model and pca_class_custom.pca_model.n_components_ >= 3:
                            class_3d_path = os.path.join(class_dir, f'pca_lipid_class_custom{custom_filename_suffix}_3d.png')
                            if plot_pca_3d(pca_class_custom, f'PCA (3D) - Lipid Classes{custom_title_suffix}', class_3d_path, color_map=ctx.color_map, group_order=custom_groups_to_plot, 
                                          fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                                          xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                                          title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                                          legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                                          view_azim=params.view_azim, view_elev=params.view_elev):
                                files_created.append(class_3d_path)
            except Exception:
                pass
            
        # Pairwise PCAs if requested
        if params.pairwise:
            pairs = params.specific_pairs if params.specific_pairs else list(itertools.combinations(ctx.groups, 2))
            
            # Apply selected_comparisons filter
            if params.selected_comparisons is not None:
                pairs = [
                    (g1, g2) for (g1, g2) in pairs
                    if (g1, g2) in params.selected_comparisons or (g2, g1) in params.selected_comparisons
                ]
            
            for g1, g2 in pairs:
                cols_pair = [c for c in ctx.sample_cols if ctx.sample_to_group[c] in (g1, g2)]
                if len(cols_pair) < 3:
                    errors.append(f"Insufficient samples for PCA {g1} vs {g2}: {len(cols_pair)}")
                    continue
                
                pca_pair = run_pca(ctx.complete_df, cols_pair, ctx.sample_to_group, n_components=params.components)
                
                # Save pair scores (controlled by save_excel parameter)
                tag = f"{g1}_vs_{g2}".replace(' ', '_')
                if params.save_excel:
                    pair_scores_path = os.path.join(pca_dir, f'pca_{tag}_scores.csv')
                    pca_pair.scores.to_csv(pair_scores_path, index=False)
                    files_created.append(pair_scores_path)
                
                # Generate pair plot (controlled by save_2d parameter)
                if params.save_2d:
                    pair_plot_path = os.path.join(pca_dir, f'pca_{tag}.png')
                    if plot_pca(
                        pca_pair, f'PCA - {g1} vs {g2}', pair_plot_path, color_map=ctx.color_map, group_order=[g1, g2],
                        fig_size=(params.fig_width, params.fig_height), dpi=params.fig_dpi, point_size=params.point_size_2d,
                        xlabel_fontsize=params.xlabel_fontsize, ylabel_fontsize=params.ylabel_fontsize,
                        title_fontsize=params.title_fontsize, tick_fontsize=params.tick_fontsize,
                        legend_fontsize=params.legend_fontsize, show_legend=params.show_legend
                    ):
                        files_created.append(pair_plot_path)
                
                # 3D pair plot if requested (controlled by save_3d parameter)
                if params.save_3d and pca_pair.pca_model and pca_pair.pca_model.n_components_ >= 3:
                    pca_pair3d = pca_pair if 'PC3' in pca_pair.scores.columns else run_pca(ctx.complete_df, cols_pair, ctx.sample_to_group, n_components=max(3, params.components))
                    pair_3d_path = os.path.join(pca_dir, f'pca_{tag}_3d.png')
                    if plot_pca_3d(pca_pair3d, f'PCA (3D) - {g1} vs {g2}', pair_3d_path, color_map=ctx.color_map, group_order=[g1,g2], 
                                  fig_size=_pca3d_fig_size(), dpi=params.fig_dpi, point_size=params.point_size_3d,
                                  xlabel_fontsize=_pca3d_font(params.xlabel_fontsize), ylabel_fontsize=_pca3d_font(params.ylabel_fontsize), 
                                  title_fontsize=_pca3d_font(params.title_fontsize), tick_fontsize=_pca3d_font(params.tick_fontsize),
                                  legend_fontsize=_pca3d_font(params.legend_fontsize), show_legend=params.show_legend,
                                  view_azim=params.view_azim, view_elev=params.view_elev):
                        files_created.append(pair_3d_path)
                    if params.interactive_3d:
                        pair_html = os.path.join(pca_dir, f'pca_{tag}_3d.html')
                        if save_interactive_pca_3d(pca_pair3d, f'PCA (3D) - {g1} vs {g2}', pair_html, color_map=ctx.color_map, point_size=params.point_size_3d):
                            files_created.append(pair_html)
        
        # Conditional summary reflecting success/partial/failure
        if files_created and not errors:
            summary = f"PCA analysis complete: {len(files_created)} files generated"
        elif files_created and errors:
            summary = f"PCA analysis partial: {len(files_created)} files, {len(errors)} errors"
        else:
            summary = f"PCA analysis failed: {len(errors)} errors"
        
    except Exception as e:
        errors.append(f"PCA analysis failed: {str(e)}")
        summary = "PCA analysis failed"
    
    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def _find_p_value_column(df: pd.DataFrame, base_pref: str, prefer_adj: bool = True) -> Optional[str]:
    """Given a base prefix like 'B_vs_A', return the best p-value column name.

    If prefer_adj is True, try *_adj_p first, then fall back to *_p_value or *_pvalue.
    If prefer_adj is False, try *_p_value / *_pvalue first, then fall back to *_adj_p.
    Returns the column name or None.
    
    IMPORTANT: Excludes *_neg_log10* columns (those are -log10 transforms, not raw p-values).
    """
    candidates_adj = [f"{base_pref}_adj_p", f"{base_pref}_adj_p_value", f"{base_pref}_p_adj"]
    candidates_raw = [f"{base_pref}_p_value", f"{base_pref}_pvalue", f"{base_pref}_p"]
    if prefer_adj:
        for c in candidates_adj + candidates_raw:
            # Skip neg_log10 columns - they are NOT p-values
            if c in df.columns and 'neg_log10' not in c.lower():
                return c
    else:
        for c in candidates_raw + candidates_adj:
            # Skip neg_log10 columns - they are NOT p-values
            if c in df.columns and 'neg_log10' not in c.lower():
                return c
    return None

def _locate_pair_columns(df: pd.DataFrame, g1: str, g2: str, prefer_adj: bool = True, verified_assignments: Optional[Dict[str, str]] = None, stat_column_assignments: Optional[Dict] = None) -> Optional[Tuple[str,str,str]]:
    """Locate pairwise stats columns for two groups, tolerant to orientation.

    Returns (log2fc_col, p_col, base_prefix_used). Searches all available *_vs_* columns 
    and returns the first match where BOTH group names appear (in any order).
    The prefer_adj flag controls whether adjusted p columns are preferred.
    
    If stat_column_assignments is provided (from Configure Stat Columns dialog), use those first.
    If verified_assignments is provided, uses those column mappings second.
    """
    # PRIORITY 1: Use stat_column_assignments from Configure Stat Columns dialog
    if stat_column_assignments and 'comparisons' in stat_column_assignments:
        comparisons = stat_column_assignments['comparisons']
        if not isinstance(comparisons, dict):
            comparisons = {}

        def _get_comp_mapping(a: str, b: str) -> Optional[Dict[str, Any]]:
            """Resolve comparison mapping by tuple or serialized string key."""
            # Native tuple keys
            if (a, b) in comparisons:
                return comparisons.get((a, b))
            if (b, a) in comparisons:
                return comparisons.get((b, a))

            # Serialized keys from JSON persistence
            key_ab = f"{a}|{b}"
            key_ba = f"{b}|{a}"
            if key_ab in comparisons:
                return comparisons.get(key_ab)
            if key_ba in comparisons:
                return comparisons.get(key_ba)
            return None

        comp = _get_comp_mapping(g1, g2)
        if comp:
            log2fc_col = comp.get('log2fc')
            p_col = comp.get('pvalue')
            fc_col = comp.get('fc')
                
            # If only FC provided (not log2FC), calculate log2FC
            if not log2fc_col and fc_col and fc_col in df.columns:
                import numpy as np
                new_col = f"{g1}_vs_{g2}_log2FC_calc"
                if new_col not in df.columns:
                    df[new_col] = np.log2(df[fc_col].replace(0, np.nan).astype(float))
                    logger.info(f"   Calculated log2FC from FC column: {fc_col} -> {new_col}")
                log2fc_col = new_col
                
            if log2fc_col and p_col and log2fc_col in df.columns and p_col in df.columns:
                logger.info(f"   Using configured stat columns for {g1} vs {g2}: {log2fc_col}, {p_col}")
                return log2fc_col, p_col, f"{g1}_vs_{g2}"
            elif log2fc_col and log2fc_col in df.columns:
                logger.debug(f"   Configured log2FC found ({log2fc_col}) but p-value missing for {g1} vs {g2}")
            elif p_col and p_col in df.columns:
                logger.debug(f"   Configured p-value found ({p_col}) but log2FC missing for {g1} vs {g2}")
    
    # PRIORITY 2: If verified assignments provided, use them (they map to column names)
    if verified_assignments:
        log2fc_col = verified_assignments.get('log2FC')
        p_col = verified_assignments.get('pvalue')
        if log2fc_col and p_col and log2fc_col in df.columns and p_col in df.columns:
            logger.info(f"   Using verified column assignments: {log2fc_col}, {p_col}")
            return log2fc_col, p_col, f"{g1}_vs_{g2}"  # prefix is informational
        # If verified assignments incomplete, fall through to standard search
    
    # PRIORITY 3: Standard column name search
    # Build potential prefixes (exact order + reversed)
    base_pref = f"{g1}_vs_{g2}"
    alt_pref = f"{g2}_vs_{g1}"
    
    # Try exact order first
    for pref in [base_pref, alt_pref]:
        fc = f"{pref}_log2FC"
        if fc in df.columns:
            # Found FC column; now find matching p-value
            p = _find_p_value_column(df, pref, prefer_adj=prefer_adj)
            if p and p in df.columns:
                logger.info(f"   Found stats for {g1} vs {g2}: {fc}, {p}")
                return fc, p, pref
    
    # Fallback: Extract all unique *_vs_* prefixes from columns and match groups
    all_prefixes = set()
    for col in df.columns:
        if '_vs_' in col:
            # Extract prefix (everything before _log2FC, _adj_p, _p_value, etc.)
            for suffix in ['_log2FC', '_FC', '_log2fc', '_fc', '_adj_p', '_p_value', '_pvalue', '_p', '_neg_log10_adj_p']:
                if col.endswith(suffix):
                    prefix = col[:-len(suffix)]
                    all_prefixes.add(prefix)
                    break
    
    # Now search all prefixes for one that contains both g1 and g2
    for prefix in sorted(all_prefixes):  # Sorted for deterministic order
        parts = prefix.split('_vs_')
        if len(parts) == 2:
            pg1, pg2 = parts[0].strip(), parts[1].strip()
            # Match if BOTH groups appear (case-insensitive, in either order)
            g1_lower, g2_lower = g1.lower(), g2.lower()
            pg1_lower, pg2_lower = pg1.lower(), pg2.lower()
            
            if (pg1_lower == g1_lower and pg2_lower == g2_lower) or (pg1_lower == g2_lower and pg2_lower == g1_lower):
                # Found a prefix that matches both groups!
                fc = f"{prefix}_log2FC"
                if fc in df.columns:
                    p = _find_p_value_column(df, prefix, prefer_adj=prefer_adj)
                    if p and p in df.columns:
                        logger.info(f"   Found stats for {g1} vs {g2} via fallback: {fc}, {p}")
                        return fc, p, prefix
    
    # Still not found; log all available stats columns to help debug
    logger.warning(f"Could not locate pairwise stats columns for {g1} vs {g2}")
    available_pairs = sorted(set(
        c.replace('_log2FC', '').replace('_FC', '').replace('_adj_p', '').replace('_p_value', '').replace('_pvalue', '').replace('_p', '').replace('_neg_log10_adj_p', '')
        for c in df.columns if '_vs_' in c and any(x in c for x in ['_log2FC', '_adj_p', '_p_value', '_pvalue', '_p'])
    ))
    if available_pairs:
        logger.info(f"  Available pairwise comparisons: {', '.join(available_pairs[:10])}")
        logger.info(f"  Groups requested: {g1}, {g2}")
        logger.info(f"  If your data uses different group names, check the visualization settings.")
    else:
        logger.warning(f"  No pairwise comparison columns found. Run Statistics tab to generate them.")
    
    return None

def generate_volcano_plots(complete_df: pd.DataFrame, groups: List[str], outdir: str, *, p_thresh: float, fc_thresh: float, top_n: int = 0, annotate: bool = False,
                           fig_width: float = 8.5, fig_height: float = 6.5, dpi: int = 230, point_size_sig: float = 28, point_size_nonsig: float = 18,
                           xlabel_fontsize: int = 11, ylabel_fontsize: int = 11, title_fontsize: int = 18, tick_fontsize: int = 11,
                           count_fontsize: int = 9, total_fontsize: int = 8, count_background: str = 'colored', legend_fontsize: int = 10, save_excel: bool = True,
                           selected_comparisons: Optional[List[Tuple[str, str]]] = None, verified_assignments: Optional[Dict[str, str]] = None,
                           stat_column_assignments: Optional[Dict] = None, id_column: Optional[str] = None, annot_fontsize: int = 8, output_format: str = 'png'):
    """Generate volcano plots ensuring baseline-first naming (g1 vs g2) with size and font controls.
    
    Parameters
    ----------
    selected_comparisons : Optional[List[Tuple[str, str]]]
        If provided, only generate plots for these specific comparisons. None means generate all.
    stat_column_assignments : Optional[Dict]
        User-configured explicit column assignments for each comparison.
    """
    from itertools import combinations
    
    # DEBUG: Log entry point with annotation params
    logger.info(f"🌋 generate_volcano_plots ENTRY:")
    logger.info(f"   └─ top_n = {top_n}")
    logger.info(f"   └─ annotate = {annotate}")
    logger.info(f"   └─ p_thresh = {p_thresh}")
    logger.info(f"   └─ fc_thresh = {fc_thresh}")
    logger.info(f"   └─ groups = {groups}")

    volcano_dir = outdir if os.path.basename(outdir.lower()) == 'volcano' else os.path.join(outdir, 'volcano')
    os.makedirs(volcano_dir, exist_ok=True)

    files_created: List[str] = []
    errors: List[str] = []

    palette = {
        'Upregulated': '#d62728',
        'Downregulated': '#2ca02c',
        # Muted sand to fade mid-band points so they read closer to the gray tier
        'SigOnly': '#d8c68a',
        'Not Significant': 'lightgray'
    }

    # Interpret fc_thresh == 0 as "no fold-change cutoff" (use p-value only).
    # Otherwise accept any positive fc_thresh (including 1.0). If fc_thresh <= 0 or invalid,
    # disable fold-change filtering.
    use_fc = True
    log2fc_thresh = None
    try:
        # Allow numeric-like inputs
        if fc_thresh == 0:
            use_fc = False
            log2fc_thresh = None
        else:
            fc_val = float(fc_thresh)
            if fc_val <= 0:
                use_fc = False
                log2fc_thresh = None
            else:
                # Accept fc_val == 1.0 (log2 threshold = 0); caller may want this behavior
                log2fc_thresh = float(np.log2(fc_val))
    except Exception:
        use_fc = False
        log2fc_thresh = None

    # Determine preference for adjusted p-values; default to True
    prefer_adj = True
    if hasattr(complete_df, '_viz_prefer_adj_p'):
        prefer_adj = bool(getattr(complete_df, '_viz_prefer_adj_p'))
    
    # Helper: establish canonical ordering for group pairs based on ctx.groups
    group_order = {g: idx for idx, g in enumerate(groups)} if groups else {}

    def canonical_pair(a: str, b: str) -> Optional[Tuple[str, str]]:
        if not a or not b or a == b:
            return None
        idx_a = group_order.get(a)
        idx_b = group_order.get(b)
        if idx_a is not None and idx_b is not None:
            return (a, b) if idx_a <= idx_b else (b, a)
        if idx_a is not None:
            return (a, b)
        if idx_b is not None:
            return (b, a)
        # Fallback: alphabetical to keep deterministic ordering
        pair = sorted([a, b], key=lambda x: str(x).lower())
        return (pair[0], pair[1]) if len(pair) == 2 else None

    # Detect available comparisons and normalize their ordering
    available_pairs: Set[Tuple[str, str]] = set()
    for col in complete_df.columns:
        if '_vs_' in col and any(sfx in col for sfx in ('_log2FC', '_log2fc', '_FC', '_fc')):
            if col.endswith('_log2FC'):
                prefix = col[:-len('_log2FC')]
            elif col.endswith('_log2fc'):
                prefix = col[:-len('_log2fc')]
            elif col.endswith('_FC'):
                prefix = col[:-len('_FC')]
            elif col.endswith('_fc'):
                prefix = col[:-len('_fc')]
            else:
                continue
            parts = prefix.split('_vs_')
            if len(parts) != 2:
                continue
            canon = canonical_pair(parts[0].strip(), parts[1].strip())
            if canon:
                available_pairs.add(canon)

    # Also include comparisons explicitly configured in Configure Stat Columns.
    if stat_column_assignments and isinstance(stat_column_assignments, dict):
        comparisons = stat_column_assignments.get('comparisons', {})
        if isinstance(comparisons, dict):
            for key in comparisons.keys():
                if isinstance(key, tuple) and len(key) == 2:
                    canon = canonical_pair(str(key[0]).strip(), str(key[1]).strip())
                    if canon:
                        available_pairs.add(canon)
                elif isinstance(key, str) and '|' in key:
                    parts = key.split('|', 1)
                    canon = canonical_pair(parts[0].strip(), parts[1].strip())
                    if canon:
                        available_pairs.add(canon)

    # Primary candidate list is all pairwise combinations from configured groups.
    # This guarantees every expected comparison is attempted and any skip reason is reported.
    expected_pairs: Set[Tuple[str, str]] = set()
    if groups and len(groups) >= 2:
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                canon = canonical_pair(groups[i], groups[j])
                if canon:
                    expected_pairs.add(canon)

    if expected_pairs:
        candidate_pairs = expected_pairs
    elif available_pairs:
        candidate_pairs = available_pairs
    else:
        candidate_pairs = set()

    logger.info(
        f"Volcano comparison discovery: expected_from_groups={len(expected_pairs)}, "
        f"detected_from_data={len(available_pairs)}, candidates={len(candidate_pairs)}"
    )

    def pair_sort_key(pair: Tuple[str, str]):
        g1, g2 = pair
        idx1 = group_order.get(g1, 10**6)
        idx2 = group_order.get(g2, 10**6)
        return (idx1, idx2, str(g1).lower(), str(g2).lower())

    selected_filter_count = 0 if selected_comparisons is None else len(selected_comparisons)
    if selected_comparisons is not None:
        logger.info(f"Volcano comparison filter active: {selected_filter_count} selected pair(s)")

    skipped_by_selection = 0
    skipped_missing_columns = 0

    # Now loop over comparisons in preferred order
    for g1, g2 in sorted(candidate_pairs, key=pair_sort_key):
        # Skip if not in selected comparisons
        if selected_comparisons is not None:
            if not ((g1, g2) in selected_comparisons or (g2, g1) in selected_comparisons):
                skipped_by_selection += 1
                errors.append(
                    f"Skipping volcano {g1} vs {g2} (not selected in Comparison Selection filter)"
                )
                continue

        found = _locate_pair_columns(complete_df, g1, g2, prefer_adj=prefer_adj, stat_column_assignments=stat_column_assignments)
        if not found:
            skipped_missing_columns += 1
            # Determine what is missing
            log2fc_col = f"{g1}_vs_{g2}_log2FC"
            p_col = _find_p_value_column(complete_df, f"{g1}_vs_{g2}", prefer_adj=prefer_adj)
            missing = []
            if log2fc_col not in complete_df.columns:
                missing.append(f"log2FC column '{log2fc_col}'")
            if not p_col or p_col not in complete_df.columns:
                missing.append(f"p-value column (expected formats like '{g1}_vs_{g2}_adj_p')")
            if missing:
                errors.append(f"Skipping volcano {g1} vs {g2} (missing: {', '.join(missing)})")
            else:
                errors.append(f"Skipping volcano {g1} vs {g2} (pairwise columns exist but incompatible)")
            continue
        log2fc_col, p_col, used_base = found

        # Use provided ID column or fallback to detection
        id_col = id_column  # Use the passed parameter first
        if not id_col:
            # Fallback: try stat_column_assignments
            if stat_column_assignments and 'id_column' in stat_column_assignments:
                id_col = stat_column_assignments['id_column']
            # Fallback: try common column names
            elif 'metabolite_id' in complete_df.columns:
                id_col = 'metabolite_id'
            elif 'Name' in complete_df.columns:
                id_col = 'Name'
            elif 'GP' in complete_df.columns:
                id_col = 'GP'
        
        cols = [id_col] if id_col and id_col in complete_df.columns else []
        logger.info(f"   🔍 ID column for annotation: {id_col if cols else 'NOT FOUND (id_column param=' + str(id_column) + ')'}")
        
        df = complete_df[cols + [log2fc_col, p_col]].copy()
        df.rename(columns={log2fc_col: 'log2FC', p_col: 'p_value'}, inplace=True)
        
        # Track total in dataset BEFORE filtering
        total_in_dataset = len(df)
        
        # Filter out metabolites with NaN/empty p-values (these are the ones skipped - not common in both groups)
        df = df[df['p_value'].notna()].copy()
        
        # Calculate tested count AFTER filtering out metabolites with missing p-values
        total_tested = len(df)
        skipped_count = total_in_dataset - total_tested
        
        df['neg_log10_p'] = -np.log10(df['p_value'].replace(0, np.finfo(float).eps))

        df['Expression'] = 'Not Significant'
        sig_mask = df['p_value'] <= p_thresh
        
        # Debug: Log significance calculation details
        total_sig_by_pvalue = sig_mask.sum()
        logger.info(f"  🔍 {g1} vs {g2}: FC filter = {'ENABLED' if use_fc else 'SKIPPED'}")
        logger.info(f"     └─ P-value threshold: {p_thresh}, {total_sig_by_pvalue}/{len(df)} features pass")
        if use_fc and log2fc_thresh is not None:
            logger.info(f"     └─ FC threshold: {fc_thresh} (log2={log2fc_thresh:.3f})")
        
        if use_fc and log2fc_thresh is not None:
            up_mask = sig_mask & (df['log2FC'] >= log2fc_thresh)
            down_mask = sig_mask & (df['log2FC'] <= -log2fc_thresh)
            df.loc[up_mask, 'Expression'] = 'Upregulated'
            df.loc[down_mask, 'Expression'] = 'Downregulated'
            # Highlight significant-but-low-FC points in a dedicated band
            mid_mask = sig_mask & ~(up_mask | down_mask)
            df.loc[mid_mask, 'Expression'] = 'SigOnly'
            logger.info(f"     └─ Significant+FC: Up={up_mask.sum()}, Down={down_mask.sum()}, SigOnly={mid_mask.sum()}")
        else:
            # No fold-change cutoff: classify by p-value and sign of log2FC
            df.loc[sig_mask & (df['log2FC'] > 0), 'Expression'] = 'Upregulated'
            df.loc[sig_mask & (df['log2FC'] < 0), 'Expression'] = 'Downregulated'
            up_by_sign = (sig_mask & (df['log2FC'] > 0)).sum()
            down_by_sign = (sig_mask & (df['log2FC'] < 0)).sum()
            logger.info(f"     └─ Significant (no FC filter): Up={up_by_sign}, Down={down_by_sign}")

        up_count = (df['Expression']=='Upregulated').sum()
        down_count = (df['Expression']=='Downregulated').sum()
        total_sig = up_count + down_count
        logger.info(f"     └─ FINAL: {total_sig} significant ({up_count} up, {down_count} down)")

        tag = f"{g1}_vs_{g2}".replace(' ', '_')
        
        # Save CSV only if save_excel is enabled
        if save_excel:
            csv_path = os.path.join(volcano_dir, f'volcano_{tag}.csv')
            df.to_csv(csv_path, index=False)
            files_created.append(csv_path)

        plt.figure(figsize=(fig_width, fig_height))
        ax = plt.gca()
        ax.scatter(df.loc[df.Expression=='Not Significant','log2FC'], 
                   df.loc[df.Expression=='Not Significant','neg_log10_p'], 
                   c=palette['Not Significant'], s=point_size_nonsig, alpha=0.6, edgecolor='none')
        if use_fc and log2fc_thresh is not None:
            sig_only_mask = df['Expression']=='SigOnly'
            if sig_only_mask.any():
                ax.scatter(df.loc[sig_only_mask, 'log2FC'],
                           df.loc[sig_only_mask, 'neg_log10_p'],
                           c=palette['SigOnly'], s=point_size_sig * 0.9, alpha=0.45,
                           edgecolor='none', linewidth=0.0)
        ax.scatter(df.loc[df.Expression=='Upregulated','log2FC'], 
                   df.loc[df.Expression=='Upregulated','neg_log10_p'], 
                   c=palette['Upregulated'], s=point_size_sig, alpha=0.85, edgecolor='k', linewidth=0.3)
        ax.scatter(df.loc[df.Expression=='Downregulated','log2FC'], 
                   df.loc[df.Expression=='Downregulated','neg_log10_p'], 
                   c=palette['Downregulated'], s=point_size_sig, alpha=0.85, edgecolor='k', linewidth=0.3)
        # Add explicit legend entries matching the expression classes
        try:
            from matplotlib.lines import Line2D
            present_expr = set(df['Expression'].unique())
            sig_only_label = "Below FC threshold"
            not_sig_label = "Not Significant"

            # Build handles
            def make_handle(edgecolor, facecolor, marker_size):
                return Line2D([0], [0], marker='o', color='w', markerfacecolor=facecolor,
                              markeredgecolor=edgecolor if edgecolor != 'none' else 'none', markersize=marker_size)

            down_handle = make_handle('k', palette['Downregulated'], 8) if 'Downregulated' in present_expr else None
            up_handle = make_handle('k', palette['Upregulated'], 8) if 'Upregulated' in present_expr else None
            not_sig_handle = make_handle('none', palette['Not Significant'], 6) if 'Not Significant' in present_expr else None
            sig_only_handle = make_handle('none', palette['SigOnly'], 7) if (use_fc and log2fc_thresh is not None and 'SigOnly' in present_expr) else None

            legend_handles = []
            legend_labels = []
            
            # Determine if we have 4 legends (FC threshold used) or 3 (FC threshold skipped)
            has_fc_legend = sig_only_handle is not None

            if has_fc_legend:
                # 4 legends: vertical arrangement on the right side
                # Order: Upregulated, Downregulated, Below FC threshold, Not Significant
                if up_handle:
                    legend_handles.append(up_handle)
                    legend_labels.append('Upregulated')
                if down_handle:
                    legend_handles.append(down_handle)
                    legend_labels.append('Downregulated')
                if sig_only_handle:
                    legend_handles.append(sig_only_handle)
                    legend_labels.append(sig_only_label)
                if not_sig_handle:
                    legend_handles.append(not_sig_handle)
                    legend_labels.append(not_sig_label)
                
                leg = None
                if legend_handles:
                    leg = ax.legend(
                        legend_handles,
                        legend_labels,
                        frameon=False,
                        loc='center left',
                        ncol=1,
                        bbox_to_anchor=(1.02, 0.5),
                        prop={'size': legend_fontsize},
                        markerscale=0.9,
                        labelspacing=0.8
                    )
            else:
                # 3 legends: horizontal arrangement at bottom center
                # Order: Downregulated, Not Significant, Upregulated
                if down_handle:
                    legend_handles.append(down_handle)
                    legend_labels.append('Downregulated')
                if not_sig_handle:
                    legend_handles.append(not_sig_handle)
                    legend_labels.append(not_sig_label)
                if up_handle:
                    legend_handles.append(up_handle)
                    legend_labels.append('Upregulated')
                
                leg = None
                if legend_handles:
                    leg = ax.legend(
                        legend_handles,
                        legend_labels,
                        frameon=False,
                        loc='lower center',
                        ncol=len(legend_handles),
                        bbox_to_anchor=(0.5, -0.40),
                        prop={'size': legend_fontsize},
                        markerscale=0.9,
                        columnspacing=1.0,
                        handletextpad=0.4
                    )

            if leg:
                for t in leg.get_texts():
                    t.set_fontweight('bold')
                    try:
                        t.set_fontsize(legend_fontsize)
                    except Exception:
                        pass
        except Exception:
            pass
        # Horizontal p-value threshold line
        try:
            ax.axhline(-np.log10(p_thresh), color='k', linestyle=(0,(5,5)), linewidth=1.4)
        except Exception:
            pass
        # Vertical fold-change cutoff lines only when fc cutoff is enabled
        if use_fc and log2fc_thresh is not None:
            try:
                ax.axvline(log2fc_thresh, color='k', linestyle=(0,(5,5)), linewidth=1.4)
                ax.axvline(-log2fc_thresh, color='k', linestyle=(0,(5,5)), linewidth=1.4)
            except Exception:
                pass
        ax.set_xlabel('log₂ fold change', fontweight='bold', fontsize=xlabel_fontsize)
        ax.set_ylabel('- log₁₀(p)', fontweight='bold', fontsize=ylabel_fontsize)
        ax.set_title(f'{g1} vs {g2}', fontweight='bold', fontsize=title_fontsize, pad=12)
        # Explicitly ensure axis label font sizes are applied
        ax.xaxis.label.set_fontsize(xlabel_fontsize)
        ax.yaxis.label.set_fontsize(ylabel_fontsize)
        ax.tick_params(axis='both', labelsize=tick_fontsize, width=1.5)
        for spine in ['top','right']:
            if spine in ax.spines:
                ax.spines[spine].set_visible(False)
        for spine in ['left','bottom']:
            if spine in ax.spines:
                ax.spines[spine].set_linewidth(2.2)
        for lbl in ax.get_xticklabels()+ax.get_yticklabels():
            lbl.set_fontweight('bold')
        # Expand upper y-limit so we can place the summary/count boxes inside the axes
        try:
            y_min, y_max = ax.get_ylim()
            y_range = (y_max - y_min) if y_max > y_min else 1.0
            # Reserve a larger top margin so stacked labels don't collide with points
            extra_top = 0.30 * y_range
            ax.set_ylim(y_min, y_max + extra_top)
        except Exception:
            pass

        # Place the summary text and counts INSIDE the axes and stacked to avoid overlap.
        # Use axes-relative coordinates (transAxes) so placement is consistent regardless of data scale.
        # Use count_background parameter to control background color
        try:
            # Determine background colors based on count_background parameter
            if count_background == 'transparent':
                down_bg = 'white'
                up_bg = 'white'
                down_alpha = 0.5
                up_alpha = 0.5
            else:  # 'colored'
                down_bg = palette['Downregulated']
                up_bg = palette['Upregulated']
                down_alpha = 1.0
                up_alpha = 1.0
            
            # Total on the top-left inside the axes
            # Show both total in dataset and tested count with skipped info
            total_text = f'Total: {total_tested}'
            
            ax.text(0.02, 0.96, total_text, transform=ax.transAxes, va='top', ha='left',
                    fontweight='bold', fontsize=total_fontsize, bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=0.8))

            # Counts: stacked below the Total. Place Down on left, Up on right.
            counts_y = 0.84
            ax.text(0.02, counts_y, str(down_count), transform=ax.transAxes, fontweight='bold', fontsize=count_fontsize,
                    bbox=dict(boxstyle='round,pad=0.18', fc=down_bg, ec='black', lw=0.8, alpha=down_alpha), ha='left', va='top')
            ax.text(0.98, counts_y, str(up_count), transform=ax.transAxes, fontweight='bold', fontsize=count_fontsize,
                    bbox=dict(boxstyle='round,pad=0.18', fc=up_bg, ec='black', lw=0.8, alpha=up_alpha), ha='right', va='top')
        except Exception:
            # Fallback to previous single-line placement if anything goes wrong
            try:
                if skipped_count > 0:
                    total_text = f'Tested: {total_tested}/{total_in_dataset} ({skipped_count} skipped)'
                else:
                    total_text = f'Total: {total_tested}'
                
                ax.text(0.02, 0.97, total_text, transform=ax.transAxes, va='top', ha='left',
                        fontweight='bold', fontsize=total_fontsize, bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=0.8))
                ax.text(0.1, 0.95, f'{down_count}', transform=ax.transAxes, fontweight='bold', fontsize=count_fontsize,
                        bbox=dict(boxstyle='round,pad=0.18', fc=palette['Downregulated'], ec='black', lw=0.8), ha='left', va='top')
                ax.text(0.9, 0.95, f'{up_count}', transform=ax.transAxes, fontweight='bold', fontsize=count_fontsize,
                        bbox=dict(boxstyle='round,pad=0.18', fc=palette['Upregulated'], ec='black', lw=0.8), ha='right', va='top')
            except Exception:
                pass
        
        # Annotate top N significant metabolites (only if annotate checkbox is checked)
        if annotate and top_n > 0 and cols:
            sig_df = df[df['Expression'] != 'Not Significant'].copy()
            
            # Debug annotation process
            logger.info(f"  📍 ANNOTATION DEBUG for {g1} vs {g2}:")
            logger.info(f"     └─ top_n parameter = {top_n}")
            logger.info(f"     └─ Total features in volcano = {len(df)}")
            logger.info(f"     └─ Significant features = {len(sig_df)} (Expression != 'Not Significant')")
            logger.info(f"        ├─ Upregulated = {(sig_df['Expression'] == 'Upregulated').sum()}")
            logger.info(f"        ├─ Downregulated = {(sig_df['Expression'] == 'Downregulated').sum()}")
            logger.info(f"        └─ SigOnly = {(sig_df['Expression'] == 'SigOnly').sum()}")
            
            if len(sig_df) == 0:
                logger.info(f"     └─ ⚠️  No significant metabolites to annotate (all p > {p_thresh})")
            elif len(sig_df) < top_n:
                logger.info(f"     └─ ⚠️  Only {len(sig_df)} significant metabolites available (requested top {top_n})")
            
            if len(sig_df) > 0:
                # Sort by combined score: -log10(p-value) * |log2FC|
                # This prioritizes metabolites with both high significance AND large fold changes
                sig_df['abs_fc'] = sig_df['log2FC'].abs()
                sig_df['combined_score'] = sig_df['neg_log10_p'] * sig_df['abs_fc']
                top_rows = sig_df.sort_values('combined_score', ascending=False).head(top_n)
                
                # Log the top features being annotated
                logger.info(f"     └─ ✅ Annotating {len(top_rows)} features (sorted by combined score: -log10(p) × |log2FC|):")
                for idx, (_, r) in enumerate(top_rows.iterrows(), 1):
                    feature_name = str(r[cols[0]])[:50] if cols else "Unknown"
                    logger.info(f"        {idx}. {feature_name}: log2FC={r['log2FC']:.3f}, p={r['p_value']:.2e}, score={r['combined_score']:.2f}")
                
                # Draw annotations on plot with connecting dotted lines
                for _, r in top_rows.iterrows():
                    ax.annotate(str(r[cols[0]]), 
                               xy=(r['log2FC'], r['neg_log10_p']),
                               xytext=(0, 10), textcoords='offset points',
                               fontsize=annot_fontsize, fontweight='bold', ha='center', va='bottom',
                               arrowprops=dict(arrowstyle='-', linestyle=':', linewidth=0.8, color='black', alpha=0.6))
        else:
            if top_n <= 0:
                logger.info(f"  📍 ANNOTATION SKIPPED: top_n = {top_n} (must be > 0)")
            elif not cols:
                logger.info(f"  📍 ANNOTATION SKIPPED: No ID column found (cols = {cols})")
        # Tight layout then expand bottom margin so legend placed below x-axis label
        plt.tight_layout()
        try:
            # Increase bottom margin to accommodate legend below plot
            plt.subplots_adjust(bottom=0.32)
        except Exception:
            pass
        # Determine file extension and path based on output format
        file_ext = '.svg' if output_format == 'svg' else '.png'
        plot_path = os.path.join(volcano_dir, f'volcano_{tag}{file_ext}')
        try:
            # Save with bbox_extra_artists to include the legend in bbox_inches='tight'
            # Get the legend object to pass to savefig
            legend_obj = ax.get_legend()
            save_kwargs: Dict[str, Any] = {'bbox_inches': 'tight', 'format': output_format}
            if output_format == 'png':
                save_kwargs['dpi'] = dpi
            if legend_obj:
                save_kwargs['bbox_extra_artists'] = [legend_obj]
            plt.savefig(plot_path, **save_kwargs)
            files_created.append(plot_path)
        except Exception as e:
            errors.append(f"Failed saving volcano plot {plot_path}: {e}")
        finally:
            plt.close()

    logger.info(
        f"Volcano generation summary: generated={len(files_created)}, "
        f"skipped_by_selection={skipped_by_selection}, skipped_missing_columns={skipped_missing_columns}"
    )
    
    return files_created, errors

def run_volcano_analysis(ctx: CommonVizContext, params: VolcanoParams) -> VizResults:
    """Volcano plot analysis."""
    files_created = []
    errors = []
    
    try:
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining
        # Propagate user preference for adjusted vs raw p-values. The GUI
        # sets ctx.use_adj_p when available; default to True.
        prefer_adj = True
        if hasattr(ctx, 'use_adj_p'):
            prefer_adj = bool(ctx.use_adj_p)
        # Attach temporary attribute to DataFrame for downstream helper
        try:
            setattr(ctx.complete_df, '_viz_prefer_adj_p', prefer_adj)
        except Exception:
            pass

        # Write a small summary file documenting the p-value / fold-change settings
        try:
            volcano_dir = ctx.output_dir if os.path.basename(ctx.output_dir.rstrip(os.sep)).lower() == 'volcano' else os.path.join(ctx.output_dir, 'volcano')
            os.makedirs(volcano_dir, exist_ok=True)
            use_fc = bool(params.fc_threshold and params.fc_threshold > 0)
            if use_fc:
                log2fc_thresh = float(np.log2(params.fc_threshold)) if params.fc_threshold > 1 else 0.0
            else:
                log2fc_thresh = None
            summary_lines = [
                f"Total metabolites: {len(ctx.complete_df)}",
                f"Pairs considered: {len(ctx.groups) if isinstance(ctx.groups, list) else 0}",
                f"p_threshold: {params.p_threshold}",
                f"fc_threshold: {params.fc_threshold}",
                f"log2FC threshold used: {log2fc_thresh if log2fc_thresh is not None else 'SKIPPED'}",
            ]
            with open(os.path.join(volcano_dir, 'volcano_filter_summary.txt'), 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(summary_lines))
        except Exception:
            pass

        files, errs = generate_volcano_plots(
            ctx.complete_df, ctx.groups, ctx.output_dir,
            p_thresh=params.p_threshold,
            fc_thresh=params.fc_threshold,
            top_n=params.annotate_top_n,
            annotate=params.annotate,
            fig_width=params.fig_width,
            fig_height=params.fig_height,
            dpi=params.fig_dpi,
            point_size_sig=params.point_size_sig,
            point_size_nonsig=params.point_size_nonsig,
            xlabel_fontsize=params.xlabel_fontsize,
            ylabel_fontsize=params.ylabel_fontsize,
            title_fontsize=params.title_fontsize,
            tick_fontsize=params.tick_fontsize,
            count_fontsize=params.count_fontsize,
            total_fontsize=params.total_fontsize,
            count_background=params.count_background,
            annot_fontsize=getattr(params, 'annot_fontsize', 8),
            output_format=getattr(params, 'output_format', 'png'),
            legend_fontsize=params.legend_fontsize,
            save_excel=params.save_excel,
            selected_comparisons=params.selected_comparisons,
            verified_assignments=ctx.verified_assignments,
            stat_column_assignments=ctx.stat_column_assignments,
            id_column=ctx.id_column  # Pass the configured ID column
        )
        files_created.extend(files)
        errors.extend(errs)

        if files_created and errors:
            summary = f"Volcano analysis partial: {len(files_created)} files generated, {len(errors)} skipped/warning item(s)"
        elif files_created:
            summary = f"Volcano analysis complete: {len(files_created)} files generated"
        elif errors:
            summary = f"Volcano analysis completed with no files: {len(errors)} skipped/warning item(s)"
        else:
            summary = "Volcano analysis completed with no output"
        
    except Exception as e:
        errors.append(f"Volcano analysis failed: {str(e)}")
        summary = "Volcano analysis failed"
    
    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def generate_boxplots(complete_df: pd.DataFrame, *, groups: List[str], sample_cols: List[str],
                     sample_to_group: Dict[str, str], outdir: str, top_n: int, no_limit: bool,
                     annotate: bool, alpha: float, group_color_map: Dict[str, str],
                     include: Optional[List[str]] = None,
                     fig_width: float = 3.0, fig_height: float = 3.0, dpi: int = 240,
                     filename_prefix: str = 'boxplot_',
                     xlabel_fontsize: int = 11, ylabel_fontsize: int = 8, title_fontsize: int = 12, tick_fontsize: int = 11,
                     save_excel: bool = True, annotate_comparisons: Optional[List[Tuple[str, str]]] = None,
                     using_custom_list: bool = False,
                     pair_records: Optional[List[Tuple[str, str, str, str]]] = None,
                     id_col_override: Optional[str] = None,
                     stat_column_assignments: Optional[Dict] = None,
                     ylabel_text: str = 'Relative Abundance (%)',
                     title_wrap_width: int = 40):
    """Generate boxplots for top metabolites.

    Added per-plot size controls, font size controls, and directory de-duplication to avoid nested boxplots/boxplots.
    
    Parameters
    ----------
    annotate_comparisons : Optional[List[Tuple[str, str]]]
        If provided, only add significance stars for these specific comparisons. None means annotate all.
    using_custom_list : bool
        If True, skip generating the boxplot_significance.csv file since metabolites were manually selected.
    """

    # Directory de-duplication
    if os.path.basename(outdir.rstrip(os.sep)).lower() == 'boxplots':
        box_dir = outdir
    else:
        box_dir = os.path.join(outdir, 'boxplots')
    os.makedirs(box_dir, exist_ok=True)
    
    files_created = []
    errors = []
    
    # Get identifier column - use override if provided (must be configured via Configure Stat Columns dialog)
    id_col = id_col_override
    if not id_col or id_col not in complete_df.columns:
        # Try fallback: use first non-numeric, non-statistical column
        available_cols = [c for c in complete_df.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
        
        fallback_id = None
        for col in available_cols:
            try:
                pd.to_numeric(complete_df[col], errors='raise')
                continue  # Skip numeric columns
            except (ValueError, TypeError):
                fallback_id = col
                break
        
        if fallback_id:
            logger.warning(f"   Boxplot: Using fallback ID column '{fallback_id}' (not explicitly configured)")
            id_col = fallback_id
        else:
            logger.error(f"   Boxplot: ID column not configured. Available non-statistical columns: {available_cols[:20]}")
            return [], [f'CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use "Configure Stat Columns" dialog to specify the ID column.']
    
    sample_cols_all = [c for c in sample_cols if c in complete_df.columns]
    if not sample_cols_all:
        return [], ['No sample columns for boxplots']
    
    df_numeric = complete_df[[id_col] + sample_cols_all].copy()

    # Optional statistical filtering if adjusted p-values present (use any available pair threshold fields)
    # Heuristic: if any *_adj_p column exists and below threshold across at least one comparison, keep metabolite.
    p_cols = [c for c in complete_df.columns if c.endswith('_adj_p')]
    # Use thresholds from BoxplotParams if present on caller (passed via closure is not here). We can look for global defaults.
    # For now leave filtering to caller update in run_boxplot_analysis if needed.
    
    # Optional inclusion filter
    if include:
        wanted_lower = {w.strip().lower() for w in include if w and isinstance(w, str)}
        pre_ct = len(df_numeric)
        df_numeric = df_numeric[df_numeric[id_col].astype(str).str.lower().isin(wanted_lower)]
        if df_numeric.empty:
            return [], [f"No metabolites matched custom list ({len(include)} provided)"]

    def _normalize_id_value(value: Any) -> Optional[str]:
        """Normalize identifiers so numeric and string forms compare consistently."""
        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass

        text = str(value).strip()
        if not text:
            return None

        try:
            numeric_value = float(text)
            if numeric_value.is_integer():
                return str(int(numeric_value))
        except Exception:
            pass

        return text

    complete_id_norm = complete_df[id_col].apply(_normalize_id_value)
    
    # Select metabolites by variance
    # Treat 0 as below-detection (missing) for boxplots.
    mat_for_var = df_numeric[sample_cols_all].apply(pd.to_numeric, errors='coerce').replace(0, np.nan)
    df_numeric['__var__'] = mat_for_var.var(axis=1, skipna=True)
    if not no_limit and top_n > 0:
        df_numeric = df_numeric.sort_values('__var__', ascending=False).head(top_n)
    
    # Find pairwise p-value columns
    import itertools
    pair_pval_cols = {}
    prefer_adj = getattr(complete_df, '_viz_prefer_adj_p', True)
    
    # DEBUG: Log what we're looking for
    logger.info(f"   Boxplot: Looking for p-value columns (prefer_adj={prefer_adj})")
    logger.debug(f"   Boxplot: Available columns in complete_df: {list(complete_df.columns)}")
    
    # If pair_records provided from verified assignments, use those directly
    if pair_records:
        logger.info(f"   Boxplot: Using pre-computed pair_records ({len(pair_records)} pairs)")
        for g1, g2, fc_col, p_col in pair_records:
            pair_pval_cols[(g1, g2)] = p_col
            logger.debug(f"   Boxplot: Using pair {g1} vs {g2}: FC={fc_col}, p-value={p_col}")
    else:
        # Fallback: search for p-value columns by pattern
        logger.info(f"   Boxplot: No pair_records provided, searching for p-value columns by pattern")
        for g1, g2 in itertools.combinations(groups, 2):
            base1 = f"{g2}_vs_{g1}"
            base2 = f"{g1}_vs_{g2}"
            
            # FIXED: Use _find_p_value_column for consistent column detection
            found_col = None
            
            # Try base2 first (g1_vs_g2, which is the standard order)
            found_col = _find_p_value_column(complete_df, base2, prefer_adj=prefer_adj)
            
            # If not found, try base1 (reversed order)
            if not found_col:
                found_col = _find_p_value_column(complete_df, base1, prefer_adj=prefer_adj)
            
            if found_col:
                pair_pval_cols[(g1, g2)] = found_col
                logger.debug(f"   Boxplot: Found p-value column for {g1} vs {g2}: {found_col}")
            else:
                logger.warning(f"   Boxplot: No p-value column found for {g1} vs {g2} (tried {base2}_*, {base1}_*)")
                # List what p-value columns ARE available for this comparison
                matching = [c for c in complete_df.columns if (base1 in c or base2 in c) and ('p_' in c.lower() or 'pvalue' in c.lower())]
                if matching:
                    logger.warning(f"   Boxplot: Available p-value columns for this comparison: {matching}")
    
    # DEBUG: Show final mapping
    logger.info(f"   Boxplot: Found {len(pair_pval_cols)} p-value column mappings: {pair_pval_cols}")

    
    sig_rows = []
    saved_files = []
    
    # Helper to sanitize filenames
    import re
    def _sanitize_filename(name: str, max_len: int = 80) -> str:
        cleaned = re.sub(r'[\\/:*?"<>|]', '_', str(name))
        cleaned = cleaned.replace(' ', '_')
        cleaned = re.sub(r'_+', '_', cleaned)
        cleaned = cleaned.strip('._')
        if len(cleaned) > max_len:
            cleaned = cleaned[:max_len]
        if not cleaned:
            cleaned = 'metabolite'
        return cleaned
    
    used_filenames = set()
    
    for _, row in df_numeric.iterrows():
        metab = row[id_col]
        vals = []
        gs = []
        for c in sample_cols_all:
            g = sample_to_group[c]
            v = row[c]
            if pd.isna(v):
                continue
            try:
                v_num = float(v)
            except Exception:
                continue
            if v_num == 0:
                continue
            vals.append(v_num)
            gs.append(g)
        
        if not vals:
            continue
        
        df_long = pd.DataFrame({'group': gs, 'value': vals})
        metab_label = _normalize_id_value(metab) or str(metab).strip()
        
        plt.figure(figsize=(fig_width, fig_height))
        ax = plt.gca()
        
        palette = {g: group_color_map.get(g,'#999999') for g in groups}
        sns.boxplot(data=df_long, x='group', y='value', hue='group', dodge=False, 
                   order=groups, palette=palette, width=0.45, linewidth=1.1, 
                   showfliers=False, ax=ax, legend=False)
        
        # Style box plots
        for line in ax.lines:
            line.set_color('black')
            line.set_linewidth(1.0)
        for patch in ax.patches:
            patch.set_edgecolor('#222222')
            patch.set_linewidth(1.1)
        
        sns.stripplot(data=df_long, x='group', y='value', hue='group', order=groups,
                     ax=ax, dodge=False, palette=palette, edgecolor='black',
                     linewidth=0.6, size=4, jitter=0.15, legend=False)
        
        legend_obj = ax.get_legend()
        if legend_obj is not None:
            legend_obj.remove()
        
        # Use normalized metabolite name for title so numeric IDs stay textual.
        title_txt = str(metab_label)
        # Smart wrapping: break at spaces when exceeding character limit
        # Goal: create the LONGEST possible lines that are still under the limit
        wrap_width = title_wrap_width
        if len(title_txt) > wrap_width:
            lines = []
            remaining = title_txt
            while len(remaining) > wrap_width:
                # Look for the last space within the wrap_width limit
                chunk = remaining[:wrap_width]
                space_pos = chunk.rfind(' ')
                
                if space_pos > 0:
                    # Found a space - break there (longest possible line under limit)
                    lines.append(remaining[:space_pos])
                    remaining = remaining[space_pos + 1:]  # Skip the space
                else:
                    # No space found - force break at wrap_width
                    lines.append(remaining[:wrap_width])
                    remaining = remaining[wrap_width:]
            
            # Add remaining text
            if remaining:
                lines.append(remaining)
            
            title_txt = '\n'.join(lines)

        # Prefer placing the title as a Figure suptitle to ensure centering and avoid axes clipping.
        title_lines = title_txt.split('\n') if isinstance(title_txt, str) else [str(title_txt)]
        n_lines = len(title_lines)
        # Adaptive fontsize: single-line a bit larger, multi-line reduced but not too small
        if n_lines <= 1:
            title_fs = title_fontsize
        else:
            title_fs = max(10, title_fontsize - (n_lines - 1))  # Increased min from 9 to 10

        try:
            fig = plt.gcf()
            # Use verticalalignment='top' for better positioning
            fig.suptitle(title_txt, fontsize=title_fs, fontweight='bold', x=0.5, y=0.98, 
                        verticalalignment='top', ha='center')
            # Remove axes title duplication
            try:
                ax.set_title('')
            except Exception:
                pass
            # Reserve top margin proportional to number of title lines - more space for multi-line titles
            top_margin = max(0.70, 0.92 - 0.05 * n_lines)  # Increased spacing factor from 0.04 to 0.05
            plt.subplots_adjust(top=top_margin)
        except Exception:
            try:
                ax.set_title(title_txt, fontweight='bold', fontsize=title_fontsize, pad=12, loc='center')
            except Exception:
                pass
        ax.set_xlabel('', fontsize=xlabel_fontsize)
        ax.set_ylabel(ylabel_text, fontweight='bold', fontsize=ylabel_fontsize)
        # Explicitly ensure axis label font sizes are applied
        ax.xaxis.label.set_fontsize(xlabel_fontsize)
        ax.yaxis.label.set_fontsize(ylabel_fontsize)
        
        ax.tick_params(axis='both', labelsize=tick_fontsize, width=1.5)
        # Allow caller to control x tick rotation via DataFrame-attached hints or params fallback
        rotation = getattr(complete_df, '_viz_xtick_rotation', None)
        rotate_on = getattr(complete_df, '_viz_rotate_xticks', None)
        if rotate_on is None:
            # Fallback to 45 deg default if not provided; generate_boxplots is usually called via run_boxplot_analysis
            rotate_on = True
        if rotation is None:
            rotation = 45
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(rotation if rotate_on else 0)
            lbl.set_horizontalalignment('right' if (rotate_on and rotation != 0) else 'center')
            lbl.set_fontweight('bold')
            lbl.set_fontsize(xlabel_fontsize)
        for lbl in ax.get_yticklabels():
            lbl.set_fontweight('bold')
            lbl.set_fontsize(tick_fontsize)
        # Force redraw to ensure all tick labels update
        plt.draw()
        
        for spine in ['top','right']:
            if spine in ax.spines:
                ax.spines[spine].set_visible(False)
        for spine in ['left','bottom']:
            if spine in ax.spines:
                ax.spines[spine].set_linewidth(2.2)
        
        # Add significance annotations if requested
        if annotate and pair_pval_cols:
            y_min = df_long['value'].min()
            y_max = df_long['value'].max()
            y_range = (y_max - y_min) if y_max > y_min else 1.0
            
            # Add padding at the bottom (10% of range) so points don't sit on x-axis
            y_min_padded = y_min - (0.10 * y_range)
            
            def star(p):
                if p is None: return ''
                if p < 1e-4: return '****'
                if p < 0.001: return '***'
                if p < 0.01: return '**'
                if p < alpha: return '*'
                return ''
            
            sig_comparisons = []
            for g1, g2 in itertools.combinations(groups, 2):
                # Skip if not in annotate_comparisons list
                if annotate_comparisons is not None:
                    if not ((g1, g2) in annotate_comparisons or (g2, g1) in annotate_comparisons):
                        continue
                
                cname = pair_pval_cols.get((g1, g2)) or pair_pval_cols.get((g2, g1))
                p_val = None
                if cname:
                    metab_norm = _normalize_id_value(metab)
                    sel = complete_df[complete_id_norm == metab_norm]
                    if not sel.empty and cname in sel.columns:
                        try:
                            p_val = float(sel.iloc[0][cname])
                            # DEBUG: Log p-value retrieval
                            logger.debug(f"   Boxplot [{metab}]: {g1} vs {g2} -> column={cname}, p_val={p_val}")
                        except Exception as e:
                            logger.warning(f"   Boxplot [{metab}]: Failed to parse p-value for {g1} vs {g2} from column {cname}: {e}")
                            p_val = None
                    else:
                        logger.warning(f"   Boxplot [{metab}]: Column {cname} not found or metabolite not in data for {g1} vs {g2}")
                else:
                    logger.warning(f"   Boxplot [{metab}]: No p-value column mapping found for {g1} vs {g2}")
                s = star(p_val)
                if s:
                    sig_comparisons.append((g1, g2, p_val, s))
                    sig_rows.append({'metabolite': metab_label, 'group1': g1, 'group2': g2, 'p_value': p_val, 'stars': s})
            
            if sig_comparisons:
                # Increase spacing between stacked significance annotations to avoid overlap
                base = y_max + 0.25 * y_range
                # Step controls vertical separation between stacked comparisons (larger -> more separation)
                step = 0.25 * y_range  # larger gap to reduce collisions when many comparisons exist
                # Auto font sizing for stars when many comparisons are stacked
                star_fs = 12
                if len(sig_comparisons) >= 3:
                    star_fs = 11
                if len(sig_comparisons) >= 5:
                    star_fs = 9
                for i, (g1, g2, p_val, stars) in enumerate(sig_comparisons):
                    h_level = base + i * step
                    x1 = groups.index(g1)
                    x2 = groups.index(g2)
                    # DEBUG: Log which groups are being connected with which p-value
                    logger.info(f"   Boxplot [{metab}]: Drawing bracket from {g1} (pos {x1}) to {g2} (pos {x2}) with p={p_val:.6f} ({stars})")
                    # Keep horizontal arms short
                    arm = step * 0.08
                    # Keep stars tight to their own bracket to avoid drifting toward the next stacked line.
                    text_y = h_level + arm + (0.004 * y_range)
                    ax.plot([x1, x1, x2, x2], [h_level, h_level + arm, h_level + arm, h_level],
                        color='black', linewidth=1.2)
                    ax.text((x1 + x2) / 2, text_y, stars, ha='center', va='bottom',
                        fontsize=star_fs, fontweight='bold')

                # Dynamic top spacing: less padding for fewer annotations, more for many
                # Scale top padding based on number of comparisons (0.10 to 0.25 of range)
                if len(sig_comparisons) == 1:
                    top_padding = 0.08 * y_range  # Minimal padding for single annotation
                elif len(sig_comparisons) == 2:
                    top_padding = 0.10 * y_range  # Slightly more for 2 annotations
                else:
                    top_padding = 0.15 * y_range  # Standard padding for 3+ annotations
                
                # Expand y-limits with padding at bottom and dynamic top spacing
                ax.set_ylim(y_min_padded, base + len(sig_comparisons) * step + top_padding)
            else:
                # No significance annotations, minimal top padding
                ax.set_ylim(y_min_padded, y_max + 0.05 * y_range)
        
        plt.tight_layout()
        # If title wrapped to multiple lines, increase top margin so title doesn't collide with annotations
        # Also add right padding
        try:
            if len(title_lines) > 1:
                # Increase top spacing proportional to number of lines
                extra_top = 0.06 * len(title_lines)
                # Use a conservative top value (smaller means more room at top) and add right padding
                plt.subplots_adjust(top=max(0.75, 0.92 - extra_top), right=0.95)
            else:
                # Add right padding for single-line titles too
                plt.subplots_adjust(right=0.95)
        except Exception:
            pass
        
        # Save plot
        base_safe = _sanitize_filename(metab_label)
        fname = f"{filename_prefix}{base_safe}.png"
        if fname in used_filenames:
            i = 1
            while True:
                cand = f"{filename_prefix}{base_safe}_{i}.png"
                if cand not in used_filenames:
                    fname = cand
                    break
                i += 1
        used_filenames.add(fname)
        
        outp = os.path.join(box_dir, fname)
        try:
            plt.savefig(outp, dpi=dpi)
            saved_files.append(outp)
        except Exception as e:
            errors.append(f"Failed to save boxplot for {metab}: {e}")
        finally:
            plt.close()
    
    # Save significance table (controlled by save_excel parameter)
    # Skip if using custom list since the metabolites were manually selected, not statistically filtered
    if sig_rows and save_excel and not using_custom_list:
        sig_path = os.path.join(box_dir, 'boxplot_significance.csv')
        pd.DataFrame(sig_rows).to_csv(sig_path, index=False)
        saved_files.append(sig_path)
    
    return saved_files, errors

def run_boxplot_analysis(ctx: CommonVizContext, params: BoxplotParams) -> VizResults:
    """Boxplot analysis."""
    files_created = []
    errors = []
    
    try:
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining

        df_all = ctx.complete_df.copy()
        logger.info(f"   Boxplot: Starting with {len(df_all)} metabolites in dataset")
        prefer_adj = True
        if hasattr(ctx, 'use_adj_p'):
            prefer_adj = bool(ctx.use_adj_p)

        # Identify ID column from configuration (must be configured via Configure Stat Columns dialog)
        id_col = None
        if ctx.stat_column_assignments and ctx.stat_column_assignments.get('id_column'):
            id_col = ctx.stat_column_assignments.get('id_column')
        elif ctx.id_column:
            id_col = ctx.id_column
        
        if not id_col or id_col not in df_all.columns:
            # Try fallback: use first non-numeric, non-statistical column
            available_cols = [c for c in df_all.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
            
            fallback_id = None
            for col in available_cols:
                try:
                    pd.to_numeric(df_all[col], errors='raise')
                    continue  # Skip numeric columns
                except (ValueError, TypeError):
                    fallback_id = col
                    break
            
            if fallback_id:
                logger.warning(f"   Boxplot: Using fallback ID column '{fallback_id}' (not explicitly configured)")
                id_col = fallback_id
            else:
                logger.error(f"   Boxplot: ID column not configured. Available non-statistical columns: {available_cols[:20]}")
                return VizResults([], [f'CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use "Configure Stat Columns" dialog to specify the ID column.'], 'Boxplot analysis failed')
        
        logger.info(f"   Boxplot: Using ID column: '{id_col}'")
        logger.info(f"   Boxplot: Available ID columns in data: {[c for c in ['PubChem_CID', 'HMDB_ID', 'CAS', 'metabolite_id', 'Name'] if c in df_all.columns]}")
        
        # Debug custom list structure
        if params.include_metabolites:
            if isinstance(params.include_metabolites, dict):
                logger.info(f"   Boxplot: Custom list structure: dict with {len(params.include_metabolites.get('names', []))} names, "
                          f"{len(params.include_metabolites.get('pubchem_ids', []))} PubChem IDs, "
                          f"{len(params.include_metabolites.get('hmdb_ids', []))} HMDB IDs, "
                          f"{len(params.include_metabolites.get('cas_ids', []))} CAS IDs")
            else:
                logger.info(f"   Boxplot: Custom list structure: simple list with {len(params.include_metabolites)} items")


        # Locate all available pair columns (log2FC + adj_p)
        pair_records: List[Tuple[str,str,str,str]] = []  # (g1,g2,fc_col,p_col)
        import itertools as _it
        
        # Debug: Check if verified_assignments are being passed
        if hasattr(ctx, 'verified_assignments') and ctx.verified_assignments:
            logger.info(f"   Boxplot: Using verified column assignments: {ctx.verified_assignments}")
        else:
            logger.info(f"   Boxplot: No verified assignments provided (ctx.verified_assignments={getattr(ctx, 'verified_assignments', 'NOT SET')})")
        
        for g1, g2 in _it.combinations(ctx.groups, 2):
            found = _locate_pair_columns(df_all, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
            if found:
                fc_col, p_col, _ = found
                pair_records.append((g1, g2, fc_col, p_col))
        
        # AUTO-FIX: If custom list provided but no statistical columns found, auto-enable custom_only mode
        if params.include_metabolites and not pair_records:
            logger.info("   Boxplot: No statistical columns found. Auto-enabling custom-list-only mode.")
            params.use_custom_only = True

        # Interpret params.fc_threshold <= 0 as "no fold-change cutoff" (p-value only)
        use_fc = bool(params.fc_threshold and params.fc_threshold > 0)
        if use_fc:
            # Only compute a log2 threshold when FC cutoff is enabled and > 1
            log2fc_thresh = np.log2(params.fc_threshold) if params.fc_threshold > 1 else 0.0
        else:
            log2fc_thresh = None

        def mask_for_pair(fc_col: str, p_col: str):
            if fc_col not in df_all.columns or p_col not in df_all.columns:
                return np.zeros(len(df_all), dtype=bool)
            if not use_fc or log2fc_thresh is None:
                return (df_all[p_col] < params.p_threshold)
            return (df_all[p_col] < params.p_threshold) & (df_all[fc_col].abs() >= log2fc_thresh)

        if params.use_custom_only and params.include_metabolites:
            logger.info("   Boxplot: Using custom list ONLY mode - bypassing p-value/FC filters")
            match_mask = match_metabolites_multi_column(df_all, params.include_metabolites, id_col)
            df_filtered = df_all[match_mask].copy()
            logger.info(f"   Boxplot: After custom-only filter: {len(df_filtered)} metabolites")
        else:
            # Build list of masks according to selected pairs
            selected_pairs = pair_records
            
            # FALLBACK: If no statistical columns found at all, fall back to showing all metabolites or custom list only
            if not selected_pairs and not pair_records:
                logger.warning("   Boxplot: No statistical columns found in data. Reverting to all metabolites or custom list only.")
                if params.include_metabolites:
                    logger.info("   Boxplot: Custom list provided - using custom-only mode")
                    match_mask = match_metabolites_multi_column(df_all, params.include_metabolites, id_col)
                    df_filtered = df_all[match_mask].copy()
                    logger.info(f"   Boxplot: After custom-only filter: {len(df_filtered)} metabolites")
                    params.use_custom_only = True  # Mark for later (to skip significance file)
                else:
                    logger.info("   Boxplot: No custom list and no statistics - showing all metabolites")
                    df_filtered = df_all.copy()
            else:
                # When groups are selected, rebuild pair_records to only include pairs from selected groups
                if params.selected_groups:
                    import itertools as _it2
                    selected_pairs = []
                    # Only look at pairs between selected groups
                    for g1, g2 in _it2.combinations(params.selected_groups, 2):
                        found = _locate_pair_columns(df_all, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                        if found:
                            fc_col, p_col, _ = found
                            selected_pairs.append((g1, g2, fc_col, p_col))
                    logger.info(f"   Boxplot: Selected groups mode - using {len(selected_pairs)} pairs from selected groups")
                elif params.filter_mode == 'specific' and params.filter_pairs:
                    wanted = {tuple(sorted(p)) for p in params.filter_pairs}
                    selected_pairs = [pr for pr in pair_records if tuple(sorted(pr[:2])) in wanted]
                
                # Apply p-value and fold-change filters - ALWAYS apply these unless using custom_only mode
                masks = [pd.Series(mask_for_pair(fc, pc)).fillna(False).to_numpy(dtype=bool) for _,_,fc,pc in selected_pairs]
                if masks:
                    if params.filter_mode == 'all':
                        agg = np.logical_and.reduce(np.asarray(masks, dtype=bool))
                    else:  # any or specific default to OR
                        agg = np.logical_or.reduce(np.asarray(masks, dtype=bool))
                    df_filtered = df_all[agg].copy()
                    logger.info(f"   Boxplot: After p-value/FC filter: {len(df_filtered)} metabolites")
                else:
                    logger.info(f"   Boxplot: No pair statistics available - using all metabolites")
                    df_filtered = df_all.copy()

        # Debug / transparency summary
        box_dir = ctx.output_dir if os.path.basename(ctx.output_dir.rstrip(os.sep)).lower() == 'boxplots' else os.path.join(ctx.output_dir, 'boxplots')
        try:
            os.makedirs(box_dir, exist_ok=True)
            summary_lines = [
                f"Total metabolites: {len(df_all)}",
                f"Pair stats available: {len(pair_records)}",
                f"Filter mode: {params.filter_mode}",
                f"Custom list provided: {bool(params.include_metabolites)}",
                f"Use custom only: {params.use_custom_only}",
                f"p_threshold: {params.p_threshold}",
                f"fc_threshold: {params.fc_threshold}",
                f"log2FC threshold used: {log2fc_thresh if log2fc_thresh is not None else 'SKIPPED'}",
                f"Metabolites after filtering: {len(df_filtered)}"
            ]
            with open(os.path.join(box_dir, 'boxplot_filter_summary.txt'), 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(summary_lines))
        except Exception:
            pass

        # Attach tick rotation hints for plotting
        try:
            setattr(df_filtered, '_viz_rotate_xticks', bool(getattr(params, 'rotate_xticks', True)))
            setattr(df_filtered, '_viz_xtick_rotation', int(getattr(params, 'xtick_rotation', 45)))
        except Exception:
            pass

        # Respect Top N unless user explicitly requests unlimited plots or custom-only mode.
        # A plain custom list used as an intersection should not force unlimited plotting.
        force_no_limit = params.no_limit or bool(params.use_custom_only and params.include_metabolites)
        
        # Determine if we're using a custom list (to skip significance file generation)
        using_custom_list = bool(params.include_metabolites)
        
        # Set preference for adjusted p-values on the dataframe for generate_boxplots to use
        # FIXED: Pass prefer_adj setting to boxplot generation for correct p-value column selection
        setattr(df_filtered, '_viz_prefer_adj_p', prefer_adj)
        
        # Filter groups and samples if selected_groups is specified
        groups_to_plot = ctx.groups
        sample_cols_to_plot = ctx.sample_cols
        if params.selected_groups:
            # Filter to only include selected groups that exist in ctx.groups
            groups_to_plot = [g for g in ctx.groups if g in params.selected_groups]
            if not groups_to_plot:
                return VizResults([], ['No valid groups selected for boxplot'], 'Boxplot analysis failed: no valid groups')
            
            # Filter sample columns to only include samples from selected groups
            sample_cols_to_plot = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) in groups_to_plot]
            if not sample_cols_to_plot:
                return VizResults([], ['No samples found for selected groups'], 'Boxplot analysis failed: no samples')
            
            logger.info(f"   Boxplot: Filtered groups from {len(ctx.groups)} to {len(groups_to_plot)}: {groups_to_plot}")
            logger.info(f"   Boxplot: Filtered samples from {len(ctx.sample_cols)} to {len(sample_cols_to_plot)}")
        
        files, errs = generate_boxplots(
            df_filtered,
            groups=groups_to_plot,
            sample_cols=sample_cols_to_plot,
            sample_to_group=ctx.sample_to_group,
            outdir=ctx.output_dir,
            top_n=params.top_n,
            no_limit=force_no_limit,  # Prevent variance filtering when using custom list
            annotate=params.annotate,
            alpha=0.05,
            group_color_map=ctx.color_map,
            include=None,  # Don't re-filter; df_filtered is already filtered
            fig_width=params.fig_width,
            fig_height=params.fig_height,
            dpi=params.fig_dpi,
            xlabel_fontsize=params.xlabel_fontsize,
            ylabel_fontsize=params.ylabel_fontsize,
            title_fontsize=params.title_fontsize,
            tick_fontsize=params.tick_fontsize,
            save_excel=params.save_excel,
            annotate_comparisons=params.annotate_comparisons,
            using_custom_list=using_custom_list,
            pair_records=pair_records,
            id_col_override=id_col,
            stat_column_assignments=ctx.stat_column_assignments,
            ylabel_text=params.ylabel_text,
            title_wrap_width=params.title_wrap_width
        )
        files_created.extend(files)
        errors.extend(errs)
        
        if files_created and not errors:
            summary = f"Boxplot analysis complete: {len(files_created)} files generated"
        elif files_created and errors:
            summary = f"Boxplot analysis partial: {len(files_created)} files, {len(errors)} errors"
        else:
            summary = f"Boxplot analysis failed: {len(errors)} errors"
        
    except Exception as e:
        errors.append(f"Boxplot analysis failed: {str(e)}")
        summary = "Boxplot analysis failed"
    
    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def generate_bargraphs(complete_df: pd.DataFrame, *, groups: List[str], sample_cols: List[str],
                      sample_to_group: Dict[str, str], outdir: str, top_n: int, no_limit: bool,
                      annotate: bool, alpha: float, group_color_map: Dict[str, str],
                      include: Optional[List[str]] = None,
                      fig_width: float = 3.0, fig_height: float = 3.0, dpi: int = 240,
                      filename_prefix: str = 'bargraph_',
                      xlabel_fontsize: int = 11, ylabel_fontsize: int = 8, title_fontsize: int = 12,
                      tick_fontsize: int = 11, legend_fontsize: int = 10,
                      save_excel: bool = True, annotate_comparisons: Optional[List[Tuple[str, str]]] = None,
                      using_custom_list: bool = False,
                      pair_records: Optional[List[Tuple[str, str, str, str]]] = None,
                      id_col_override: Optional[str] = None,
                      stat_column_assignments: Optional[Dict] = None,
                      ylabel_text: str = 'Relative Abundance (%)',
                      title_wrap_width: int = 40,
                      display_mode: str = 'separate',
                      grouped_title: str = '',
                      low_value_boost_enabled: bool = False,
                      low_value_boost_threshold: float = 0.25,
                      low_value_boost_factor: float = 2.0):
    """Generate bar graphs for top metabolites with SD bars and significance brackets."""

    if os.path.basename(outdir.rstrip(os.sep)).lower() == 'bargraphs':
        bar_dir = outdir
    else:
        bar_dir = os.path.join(outdir, 'bargraphs')
    os.makedirs(bar_dir, exist_ok=True)

    files_created = []
    errors = []

    id_col = id_col_override
    if not id_col or id_col not in complete_df.columns:
        available_cols = [c for c in complete_df.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
        fallback_id = None
        for col in available_cols:
            try:
                pd.to_numeric(complete_df[col], errors='raise')
                continue
            except (ValueError, TypeError):
                fallback_id = col
                break
        if fallback_id:
            logger.warning(f"   Bargraph: Using fallback ID column '{fallback_id}' (not explicitly configured)")
            id_col = fallback_id
        else:
            logger.error(f"   Bargraph: ID column not configured. Available non-statistical columns: {available_cols[:20]}")
            return [], [f'CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use "Configure Stat Columns" dialog to specify the ID column.']

    sample_cols_all = [c for c in sample_cols if c in complete_df.columns]
    if not sample_cols_all:
        return [], ['No sample columns for bar graphs']

    df_numeric = complete_df[[id_col] + sample_cols_all].copy()

    if include:
        wanted_lower = {w.strip().lower() for w in include if w and isinstance(w, str)}
        df_numeric = df_numeric[df_numeric[id_col].astype(str).str.lower().isin(wanted_lower)]
        if df_numeric.empty:
            return [], [f"No metabolites matched custom list ({len(include)} provided)"]

    def _normalize_id_value(value: Any) -> Optional[str]:
        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass

        text = str(value).strip()
        if not text:
            return None

        try:
            numeric_value = float(text)
            if numeric_value.is_integer():
                return str(int(numeric_value))
        except Exception:
            pass

        return text

    complete_id_norm = complete_df[id_col].apply(_normalize_id_value)

    mat_for_var = df_numeric[sample_cols_all].apply(pd.to_numeric, errors='coerce').replace(0, np.nan)
    df_numeric['__var__'] = mat_for_var.var(axis=1, skipna=True)
    if not no_limit and top_n > 0:
        df_numeric = df_numeric.sort_values('__var__', ascending=False).head(top_n)

    pair_pval_cols = {}
    prefer_adj = getattr(complete_df, '_viz_prefer_adj_p', True)
    if pair_records:
        for g1, g2, _, p_col in pair_records:
            pair_pval_cols[(g1, g2)] = p_col
    else:
        for g1, g2 in itertools.combinations(groups, 2):
            base1 = f"{g2}_vs_{g1}"
            base2 = f"{g1}_vs_{g2}"
            found_col = _find_p_value_column(complete_df, base2, prefer_adj=prefer_adj)
            if not found_col:
                found_col = _find_p_value_column(complete_df, base1, prefer_adj=prefer_adj)
            if found_col:
                pair_pval_cols[(g1, g2)] = found_col

    def _star(p):
        if p is None:
            return ''
        if p < 1e-4:
            return '****'
        if p < 0.001:
            return '***'
        if p < 0.01:
            return '**'
        if p < alpha:
            return '*'
        return ''

    def _sanitize_filename(name: str, max_len: int = 80) -> str:
        import re
        cleaned = re.sub(r'[\\/:*?"<>|]', '_', str(name))
        cleaned = cleaned.replace(' ', '_')
        cleaned = re.sub(r'_+', '_', cleaned)
        cleaned = cleaned.strip('._')
        if len(cleaned) > max_len:
            cleaned = cleaned[:max_len]
        return cleaned if cleaned else 'metabolite'

    def _metabolite_group_values(row) -> Dict[str, List[float]]:
        data = {g: [] for g in groups}
        for c in sample_cols_all:
            g = sample_to_group.get(c)
            if g not in data:
                continue
            v = row[c]
            if pd.isna(v):
                continue
            try:
                v_num = float(v)
            except Exception:
                continue
            if v_num == 0:
                continue
            data[g].append(v_num)
        return data

    sig_rows = []

    if str(display_mode).lower() == 'grouped':
        selected_rows = []
        for _, row in df_numeric.iterrows():
            grouped_vals = _metabolite_group_values(row)
            if any(len(v) > 0 for v in grouped_vals.values()):
                selected_rows.append((row, grouped_vals))

        if not selected_rows:
            return [], ['No metabolites with plottable values for grouped bar graph']

        n_metabs = len(selected_rows)
        n_groups = max(1, len(groups))
        bar_width = min(0.22, 0.82 / n_groups)
        cluster_span = n_groups * bar_width
        cluster_gap = max(0.45, bar_width * 2.0)
        centers = np.arange(n_metabs) * (cluster_span + cluster_gap)

        plt.figure(figsize=(max(fig_width, 4.5), max(fig_height, 3.2)))
        ax = plt.gca()
        palette = {g: group_color_map.get(g, '#999999') for g in groups}

        mean_matrix = {}
        sd_matrix = {}
        for g in groups:
            means = []
            sds = []
            for _, grouped_vals in selected_rows:
                arr = grouped_vals.get(g, [])
                if arr:
                    means.append(float(np.mean(arr)))
                    sds.append(float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0)
                else:
                    means.append(0.0)
                    sds.append(0.0)
            mean_matrix[g] = np.asarray(means, dtype=float)
            sd_matrix[g] = np.asarray(sds, dtype=float)

        cluster_scale = np.ones(n_metabs, dtype=float)
        boost_applied_count = 0
        boost_factor = 1.0
        boost_threshold = 0.0
        if low_value_boost_enabled:
            try:
                boost_factor = max(1.0, float(low_value_boost_factor))
            except Exception:
                boost_factor = 1.0
            try:
                boost_threshold = max(0.0, float(low_value_boost_threshold))
            except Exception:
                boost_threshold = 0.0

            if boost_factor > 1.0 and boost_threshold > 0.0:
                for mi in range(n_metabs):
                    cluster_peak = max(float(mean_matrix[g][mi]) for g in groups)
                    # Scale all bars in the metabolite cluster equally when the cluster is low.
                    if cluster_peak < boost_threshold:
                        cluster_scale[mi] = boost_factor
                        boost_applied_count += 1

        bar_positions = {}
        for gi, g in enumerate(groups):
            offset = (gi - (n_groups - 1) / 2.0) * bar_width
            xpos = centers + offset
            bar_positions[g] = xpos
            means_scaled = mean_matrix[g] * cluster_scale
            sds_scaled = sd_matrix[g] * cluster_scale
            ax.bar(
                xpos,
                means_scaled,
                width=bar_width,
                yerr=sds_scaled,
                capsize=3,
                label=g,
                color=palette[g],
                edgecolor='#222222',
                linewidth=0.9
            )

        metab_labels = []
        cluster_tops = []
        for row, grouped_vals in selected_rows:
            metab_raw = row[id_col]
            metab_label = _normalize_id_value(metab_raw) or str(metab_raw).strip()
            metab_labels.append(metab_label)
            mvals = []
            for g in groups:
                arr = grouped_vals.get(g, [])
                if arr:
                    mvals.append(float(np.mean(arr)) + (float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0))
            cluster_tops.append(max(mvals) if mvals else 0.0)

        cluster_tops = [float(v) * float(cluster_scale[i]) for i, v in enumerate(cluster_tops)]

        ax.set_xticks(centers)
        ax.set_xticklabels(metab_labels)
        rotation = getattr(complete_df, '_viz_xtick_rotation', 45)
        rotate_on = getattr(complete_df, '_viz_rotate_xticks', True)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(rotation if rotate_on else 0)
            lbl.set_horizontalalignment('right' if (rotate_on and rotation != 0) else 'center')
            lbl.set_fontweight('bold')
            lbl.set_fontsize(tick_fontsize)
        for lbl in ax.get_yticklabels():
            lbl.set_fontweight('bold')
            lbl.set_fontsize(tick_fontsize)

        ax.set_ylabel(ylabel_text, fontweight='bold', fontsize=ylabel_fontsize)
        ax.set_xlabel('', fontsize=xlabel_fontsize)
        ax.tick_params(axis='both', labelsize=tick_fontsize, width=1.5)
        legend_obj = ax.legend(
            frameon=False,
            fontsize=legend_fontsize,
            ncol=1,
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),
            borderaxespad=0.0
        )

        title_txt = str(grouped_title).strip() if grouped_title is not None else ''
        if not title_txt:
            title_txt = 'Grouped Bar Graphs'
        ax.set_title(title_txt, fontweight='bold', fontsize=title_fontsize)

        if low_value_boost_enabled and boost_factor > 1.0 and boost_threshold > 0.0:
            try:
                factor_txt = str(int(boost_factor)) if float(boost_factor).is_integer() else f"{boost_factor:g}"
            except Exception:
                factor_txt = f"{boost_factor}"

            report_lines = [
                "Grouped Bar Graph - Low-Value Boost Report",
                f"threshold: < {boost_threshold:g}",
                f"scale_factor: x{factor_txt}",
                f"applied_clusters: {boost_applied_count}/{n_metabs}",
                "",
                "metabolite\tmax_mean_before\tscale_applied"
            ]

            for i, metab_label in enumerate(metab_labels):
                cluster_peak = max(float(mean_matrix[g][i]) for g in groups)
                scale_value = float(cluster_scale[i])
                scale_txt = f"x{scale_value:g}" if scale_value > 1.0 else "x1"
                report_lines.append(f"{metab_label}\t{cluster_peak:.6g}\t{scale_txt}")

            report_path = os.path.join(bar_dir, f"{filename_prefix}grouped_boost_report.txt")
            try:
                with open(report_path, 'w', encoding='utf-8') as fh:
                    fh.write('\n'.join(report_lines))
                files_created.append(report_path)
            except Exception as e:
                errors.append(f"Failed to save grouped boost report: {e}")

        if annotate and pair_pval_cols:
            max_y = max(cluster_tops) if cluster_tops else 1.0
            y_range = max(1.0, max_y)
            global_max_needed = max_y

            for mi, (row, _) in enumerate(selected_rows):
                metab_raw = row[id_col]
                metab_label = _normalize_id_value(metab_raw) or str(metab_raw).strip()
                metab_norm = _normalize_id_value(metab_raw)
                sel = complete_df[complete_id_norm == metab_norm]
                sig_comparisons = []
                for g1, g2 in itertools.combinations(groups, 2):
                    if annotate_comparisons is not None and not ((g1, g2) in annotate_comparisons or (g2, g1) in annotate_comparisons):
                        continue
                    cname = pair_pval_cols.get((g1, g2)) or pair_pval_cols.get((g2, g1))
                    p_val = None
                    if cname and not sel.empty and cname in sel.columns:
                        try:
                            p_val = float(sel.iloc[0][cname])
                        except Exception:
                            p_val = None
                    stars = _star(p_val)
                    if stars:
                        sig_comparisons.append((g1, g2, p_val, stars))
                        sig_rows.append({'metabolite': metab_label, 'group1': g1, 'group2': g2, 'p_value': p_val, 'stars': stars})

                if sig_comparisons:
                    local_top = cluster_tops[mi]
                    base = local_top + 0.12 * y_range
                    step = 0.12 * y_range
                    arm = step * 0.10
                    for si, (g1, g2, _, stars) in enumerate(sig_comparisons):
                        h = base + si * step
                        x1 = float(bar_positions[g1][mi])
                        x2 = float(bar_positions[g2][mi])
                        ax.plot([x1, x1, x2, x2], [h, h + arm, h + arm, h], color='black', linewidth=1.1)
                        ax.text((x1 + x2) / 2.0, h + arm + (0.01 * y_range), stars,
                                ha='center', va='bottom', fontsize=max(9, tick_fontsize), fontweight='bold')
                    global_max_needed = max(global_max_needed, base + len(sig_comparisons) * step + 0.08 * y_range)

            ax.set_ylim(bottom=0, top=global_max_needed)

        for spine in ['top', 'right']:
            if spine in ax.spines:
                ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            if spine in ax.spines:
                ax.spines[spine].set_linewidth(2.0)

        # Keep room on the right for the outside legend so it does not overlap the plot.
        plt.tight_layout()
        plt.subplots_adjust(right=0.82)
        outp = os.path.join(bar_dir, f"{filename_prefix}grouped.png")
        try:
            plt.savefig(outp, dpi=dpi, bbox_inches='tight', bbox_extra_artists=(legend_obj,))
            files_created.append(outp)
        except Exception as e:
            errors.append(f"Failed to save grouped bar graph: {e}")
        finally:
            plt.close()
    else:
        used_filenames = set()
        for _, row in df_numeric.iterrows():
            metab = row[id_col]
            grouped_vals = _metabolite_group_values(row)
            if not any(len(v) > 0 for v in grouped_vals.values()):
                continue

            means = []
            sds = []
            for g in groups:
                arr = grouped_vals.get(g, [])
                if arr:
                    means.append(float(np.mean(arr)))
                    sds.append(float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0)
                else:
                    means.append(0.0)
                    sds.append(0.0)

            metab_label = _normalize_id_value(metab) or str(metab).strip()
            palette = [group_color_map.get(g, '#999999') for g in groups]

            plt.figure(figsize=(fig_width, fig_height))
            ax = plt.gca()
            xpos = np.arange(len(groups), dtype=float)
            ax.bar(
                xpos,
                np.asarray(means, dtype=float),
                yerr=np.asarray(sds, dtype=float),
                capsize=4,
                color=palette,
                edgecolor='#222222',
                linewidth=0.9,
                width=0.65
            )

            ax.set_xticks(xpos)
            ax.set_xticklabels(groups)
            rotation = getattr(complete_df, '_viz_xtick_rotation', 45)
            rotate_on = getattr(complete_df, '_viz_rotate_xticks', True)
            for lbl in ax.get_xticklabels():
                lbl.set_rotation(rotation if rotate_on else 0)
                lbl.set_horizontalalignment('right' if (rotate_on and rotation != 0) else 'center')
                lbl.set_fontweight('bold')
                lbl.set_fontsize(xlabel_fontsize)
            for lbl in ax.get_yticklabels():
                lbl.set_fontweight('bold')
                lbl.set_fontsize(tick_fontsize)

            title_txt = str(metab_label)
            if len(title_txt) > title_wrap_width:
                lines = []
                remaining = title_txt
                while len(remaining) > title_wrap_width:
                    chunk = remaining[:title_wrap_width]
                    space_pos = chunk.rfind(' ')
                    if space_pos > 0:
                        lines.append(remaining[:space_pos])
                        remaining = remaining[space_pos + 1:]
                    else:
                        lines.append(remaining[:title_wrap_width])
                        remaining = remaining[title_wrap_width:]
                if remaining:
                    lines.append(remaining)
                title_txt = '\n'.join(lines)
            ax.set_title(title_txt, fontweight='bold', fontsize=title_fontsize, loc='center')

            ax.set_xlabel('', fontsize=xlabel_fontsize)
            ax.set_ylabel(ylabel_text, fontweight='bold', fontsize=ylabel_fontsize)
            ax.tick_params(axis='both', labelsize=tick_fontsize, width=1.5)

            if annotate and pair_pval_cols:
                y_max = float(np.max(np.asarray(means) + np.asarray(sds))) if means else 1.0
                y_range = max(1.0, y_max)
                metab_norm = _normalize_id_value(metab)
                sel = complete_df[complete_id_norm == metab_norm]
                sig_comparisons = []
                for g1, g2 in itertools.combinations(groups, 2):
                    if annotate_comparisons is not None and not ((g1, g2) in annotate_comparisons or (g2, g1) in annotate_comparisons):
                        continue
                    cname = pair_pval_cols.get((g1, g2)) or pair_pval_cols.get((g2, g1))
                    p_val = None
                    if cname and not sel.empty and cname in sel.columns:
                        try:
                            p_val = float(sel.iloc[0][cname])
                        except Exception:
                            p_val = None
                    stars = _star(p_val)
                    if stars:
                        sig_comparisons.append((g1, g2, p_val, stars))
                        sig_rows.append({'metabolite': metab_label, 'group1': g1, 'group2': g2, 'p_value': p_val, 'stars': stars})

                if sig_comparisons:
                    base = y_max + 0.18 * y_range
                    step = 0.18 * y_range
                    arm = step * 0.10
                    for i, (g1, g2, _, stars) in enumerate(sig_comparisons):
                        h = base + i * step
                        x1 = groups.index(g1)
                        x2 = groups.index(g2)
                        ax.plot([x1, x1, x2, x2], [h, h + arm, h + arm, h], color='black', linewidth=1.2)
                        ax.text((x1 + x2) / 2.0, h + arm + (0.01 * y_range), stars,
                                ha='center', va='bottom', fontsize=max(9, tick_fontsize), fontweight='bold')
                    ax.set_ylim(bottom=0, top=base + len(sig_comparisons) * step + 0.08 * y_range)
                else:
                    ax.set_ylim(bottom=0, top=y_max + 0.08 * y_range)

            for spine in ['top', 'right']:
                if spine in ax.spines:
                    ax.spines[spine].set_visible(False)
            for spine in ['left', 'bottom']:
                if spine in ax.spines:
                    ax.spines[spine].set_linewidth(2.0)

            plt.tight_layout()
            base_safe = _sanitize_filename(metab_label)
            fname = f"{filename_prefix}{base_safe}.png"
            if fname in used_filenames:
                i = 1
                while True:
                    cand = f"{filename_prefix}{base_safe}_{i}.png"
                    if cand not in used_filenames:
                        fname = cand
                        break
                    i += 1
            used_filenames.add(fname)

            outp = os.path.join(bar_dir, fname)
            try:
                plt.savefig(outp, dpi=dpi)
                files_created.append(outp)
            except Exception as e:
                errors.append(f"Failed to save bar graph for {metab}: {e}")
            finally:
                plt.close()

    if sig_rows and save_excel and not using_custom_list:
        sig_path = os.path.join(bar_dir, 'bargraph_significance.csv')
        pd.DataFrame(sig_rows).to_csv(sig_path, index=False)
        files_created.append(sig_path)

    return files_created, errors

def run_bargraph_analysis(ctx: CommonVizContext, params: BargraphParams) -> VizResults:
    """Bar graph analysis with the same filtering behavior as boxplots."""
    files_created = []
    errors = []

    try:
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining

        df_all = ctx.complete_df.copy()
        prefer_adj = bool(getattr(ctx, 'use_adj_p', True))

        id_col = None
        if ctx.stat_column_assignments and ctx.stat_column_assignments.get('id_column'):
            id_col = ctx.stat_column_assignments.get('id_column')
        elif ctx.id_column:
            id_col = ctx.id_column

        if not id_col or id_col not in df_all.columns:
            available_cols = [c for c in df_all.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
            fallback_id = None
            for col in available_cols:
                try:
                    pd.to_numeric(df_all[col], errors='raise')
                    continue
                except (ValueError, TypeError):
                    fallback_id = col
                    break
            if fallback_id:
                id_col = fallback_id
            else:
                return VizResults([], [f'CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use "Configure Stat Columns" dialog to specify the ID column.'], 'Bar graph analysis failed')

        pair_records: List[Tuple[str, str, str, str]] = []
        for g1, g2 in itertools.combinations(ctx.groups, 2):
            found = _locate_pair_columns(df_all, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
            if found:
                fc_col, p_col, _ = found
                pair_records.append((g1, g2, fc_col, p_col))

        use_fc = bool(params.fc_threshold and params.fc_threshold > 0)
        log2fc_thresh = np.log2(params.fc_threshold) if (use_fc and params.fc_threshold > 1) else (0.0 if use_fc else None)

        def mask_for_pair(fc_col: str, p_col: str):
            if fc_col not in df_all.columns or p_col not in df_all.columns:
                return np.zeros(len(df_all), dtype=bool)
            if not use_fc or log2fc_thresh is None:
                return (df_all[p_col] < params.p_threshold)
            return (df_all[p_col] < params.p_threshold) & (df_all[fc_col].abs() >= log2fc_thresh)

        if params.use_custom_only and params.include_metabolites:
            match_mask = match_metabolites_multi_column(df_all, params.include_metabolites, id_col)
            df_filtered = df_all[match_mask].copy()
        else:
            selected_pairs = pair_records
            if params.selected_groups:
                selected_pairs = []
                for g1, g2 in itertools.combinations(params.selected_groups, 2):
                    found = _locate_pair_columns(df_all, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                    if found:
                        fc_col, p_col, _ = found
                        selected_pairs.append((g1, g2, fc_col, p_col))
            elif params.filter_mode == 'specific' and params.filter_pairs:
                wanted = {tuple(sorted(p)) for p in params.filter_pairs}
                selected_pairs = [pr for pr in pair_records if tuple(sorted(pr[:2])) in wanted]

            masks = [pd.Series(mask_for_pair(fc, pc)).fillna(False).to_numpy(dtype=bool) for _, _, fc, pc in selected_pairs]
            if masks:
                if params.filter_mode == 'all':
                    agg = np.logical_and.reduce(np.asarray(masks, dtype=bool))
                else:
                    agg = np.logical_or.reduce(np.asarray(masks, dtype=bool))
                df_filtered = df_all[agg].copy()
            else:
                df_filtered = df_all.copy()

        try:
            setattr(df_filtered, '_viz_rotate_xticks', bool(getattr(params, 'rotate_xticks', True)))
            setattr(df_filtered, '_viz_xtick_rotation', int(getattr(params, 'xtick_rotation', 45)))
            setattr(df_filtered, '_viz_prefer_adj_p', prefer_adj)
        except Exception:
            pass

        groups_to_plot = ctx.groups
        sample_cols_to_plot = ctx.sample_cols
        if params.selected_groups:
            groups_to_plot = [g for g in ctx.groups if g in params.selected_groups]
            if not groups_to_plot:
                return VizResults([], ['No valid groups selected for bar graph'], 'Bar graph analysis failed: no valid groups')
            sample_cols_to_plot = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) in groups_to_plot]
            if not sample_cols_to_plot:
                return VizResults([], ['No samples found for selected groups'], 'Bar graph analysis failed: no samples')

        files, errs = generate_bargraphs(
            df_filtered,
            groups=groups_to_plot,
            sample_cols=sample_cols_to_plot,
            sample_to_group=ctx.sample_to_group,
            outdir=ctx.output_dir,
            top_n=params.top_n,
            no_limit=(params.no_limit or bool(params.use_custom_only and params.include_metabolites)),
            annotate=params.annotate,
            alpha=0.05,
            group_color_map=ctx.color_map,
            include=None,
            fig_width=params.fig_width,
            fig_height=params.fig_height,
            dpi=params.fig_dpi,
            xlabel_fontsize=params.xlabel_fontsize,
            ylabel_fontsize=params.ylabel_fontsize,
            title_fontsize=params.title_fontsize,
            tick_fontsize=params.tick_fontsize,
            legend_fontsize=params.legend_fontsize,
            save_excel=params.save_excel,
            annotate_comparisons=params.annotate_comparisons,
            using_custom_list=bool(params.include_metabolites),
            pair_records=pair_records,
            id_col_override=id_col,
            stat_column_assignments=ctx.stat_column_assignments,
            ylabel_text=params.ylabel_text,
            title_wrap_width=params.title_wrap_width,
            display_mode=params.display_mode,
            grouped_title=params.grouped_title,
            low_value_boost_enabled=params.low_value_boost_enabled,
            low_value_boost_threshold=params.low_value_boost_threshold,
            low_value_boost_factor=params.low_value_boost_factor
        )
        files_created.extend(files)
        errors.extend(errs)

        if files_created and not errors:
            summary = f"Bar graph analysis complete: {len(files_created)} files generated"
        elif files_created and errors:
            summary = f"Bar graph analysis partial: {len(files_created)} files, {len(errors)} errors"
        else:
            summary = f"Bar graph analysis failed: {len(errors)} errors"
    except Exception as e:
        errors.append(f"Bar graph analysis failed: {str(e)}")
        summary = "Bar graph analysis failed"

    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def _add_heatmap_splits(ax, data_matrix, sample_cols, sample_to_group, groups, 
                        log2fc_data=None, add_row_split=True, add_col_split=True):
    """Add visual splits to heatmap for up/down regulation and group boundaries.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The heatmap axes
    data_matrix : pd.DataFrame
        The data matrix being plotted
    sample_cols : list
        List of sample column names in display order
    sample_to_group : dict
        Mapping of sample columns to groups
    groups : list
        List of groups in order
    log2fc_data : pd.Series, optional
        Log2 fold change values for each row (for row split)
    add_row_split : bool
        Whether to add horizontal line separating up/down regulated
    add_col_split : bool
        Whether to add vertical lines separating groups
    """
    try:
        n_rows, n_cols = data_matrix.shape
        
        # Add fine gridlines between individual cells (subtle grid for clarity)
        # Draw thin black lines between each row and column
        for i in range(1, n_rows):
            ax.axhline(y=i - 0.5, color='black', linewidth=0.3, linestyle='-', alpha=0.3, zorder=1)
        for j in range(1, n_cols):
            ax.axvline(x=j - 0.5, color='black', linewidth=0.3, linestyle='-', alpha=0.3, zorder=1)
        
        # Add column splits (vertical lines between groups)
        if add_col_split and len(groups) > 1:
            group_boundaries = []
            current_pos = 0
            for g in groups[:-1]:  # Don't need line after last group
                # Count samples in this group
                group_size = sum(1 for c in sample_cols if sample_to_group.get(c) == g)
                current_pos += group_size
                group_boundaries.append(current_pos - 0.5)  # Position between columns
            
            # Draw vertical lines
            for boundary in group_boundaries:
                ax.axvline(x=boundary, color='white', linewidth=2.5, linestyle='-', zorder=10)
        
        # Add row split (horizontal line between up/down regulated)
        if add_row_split and log2fc_data is not None and len(log2fc_data) > 1:
            # Draw line exactly where log2FC transitions from positive to non-positive
            try:
                # Handle both Series and array-like data, preserve order
                if hasattr(log2fc_data, 'values'):
                    vals = pd.to_numeric(log2fc_data.values, errors='coerce')
                else:
                    vals = pd.to_numeric(log2fc_data, errors='coerce')
                
                # Replace NaN with 0
                vals = np.nan_to_num(vals, nan=0.0)
                
                # Find the EXACT position where positive values end and non-positive begin
                # The data should already be sorted by log2FC (descending) from _order_and_cluster_rows
                split_idx = None
                for i in range(len(vals) - 1):
                    # Split line goes between last positive and first non-positive
                    if vals[i] > 0 and vals[i + 1] <= 0:
                        split_idx = i + 0.5  # Position between row i and i+1
                        break
                
                # Draw the line if we found a transition point
                if split_idx is not None:
                    ax.axhline(y=split_idx, color='white', linewidth=3.5, linestyle='-', zorder=10)
                    
                    # Count for reporting
                    pos_count = sum(1 for v in vals if v > 0)
                    neg_count = sum(1 for v in vals if v < 0)
                    zero_count = sum(1 for v in vals if v == 0)
                    
                    logger.info(
                        f"   ✓ Split line drawn at row {split_idx:.1f} "
                        f"[Above: {int(split_idx + 0.5)} up-regulated | Below: {len(vals) - int(split_idx + 0.5)} down-regulated/unchanged]"
                        f" (Total: {pos_count} up, {neg_count} down, {zero_count} unchanged)"
                    )
                        
                else:
                    # No clear transition found; summarize counts and explain
                    pos_count = int(np.sum(vals > 0))
                    neg_count = int(np.sum(vals < 0))
                    zero_count = int(np.sum(vals == 0))

                    if pos_count == 0 and neg_count == 0:
                        logger.info(f"   No split line: All {zero_count} metabolites are unchanged (no up or down)")
                    elif pos_count == 0:
                        logger.info(f"   No split line: All {neg_count + zero_count} metabolites are downregulated/unchanged (no upregulated to split from)")
                    elif neg_count == 0:
                        logger.info(f"   No split line: All {pos_count + zero_count} metabolites are upregulated/unchanged (no downregulated to split from)")
                    else:
                        logger.debug(
                            f"   No split line: {pos_count} upregulated, {neg_count} downregulated, {zero_count} unchanged (need adjacent transition)"
                        )
            except Exception as e:
                logger.warning(f"Could not draw up/down split line: {e}")
                
                # Optional: Add text labels
                # ax.text(-0.5, upregulated_count/2, ha='right', va='center',
                #        fontweight='bold', fontsize=9, color='red', rotation=90)
                # ax.text(-0.5, upregulated_count + downregulated_count/2, 
                #        ha='right', va='center', fontweight='bold', fontsize=9, 
                #        color='green', rotation=90)
    except Exception as e:
        # Fail silently - don't break heatmap generation
        pass

def _setup_heatmap_figure(n_samples: int, n_features: int, params: HeatmapParams, 
                         row_linkage=None):
    """Setup figure with dendrogram, heatmap, and colorbar axes.
    
    Returns: (fig, ax_heatmap, ax_dendro, ax_cbar)
    """
    # Calculate dimensions - check auto_size first
    use_auto_size = getattr(params, 'auto_size', True)
    
    if use_auto_size:
        dyn_width = min(3.5 + (0.15 * n_samples) + (0.02 * n_features), 25.0)
        dyn_height = min(2.0 + (0.35 * n_features) + (0.05 * n_samples), 30.0)
    else:
        dyn_width = params.fig_width
        dyn_height = params.fig_height
    
    # Colorbar sizing
    try:
        cbar_height_inches = float(getattr(params, 'colorbar_height_inches', 0.6))
    except Exception:
        cbar_height_inches = 0.6
    cbar_height_inches = max(0.3, min(cbar_height_inches, 3.0))
    
    # FIXED: Use absolute spacing in inches (0.15 inches fixed gap)
    spacing_inches = 0.15
    show_cbar = getattr(params, 'show_colorbar', True)
    
    # Create figure with EXTRA height for colorbar in margin
    if show_cbar:
        # Add colorbar height + spacing to figure height
        total_height = dyn_height + cbar_height_inches + spacing_inches
    else:
        total_height = dyn_height
    
    fig = plt.figure(figsize=(dyn_width, total_height), constrained_layout=False)
    
    # Calculate normalized positions based on TOTAL height
    if show_cbar:
        # Heatmap occupies the bottom portion (dyn_height / total_height)
        heatmap_fraction = dyn_height / total_height
        cbar_height_norm = cbar_height_inches / total_height
        spacing_norm = spacing_inches / total_height
    else:
        heatmap_fraction = 0.98
        cbar_height_norm = 0.0
        spacing_norm = 0.0
    
    # Dendrogram width
    have_dendro = bool(params.cluster and row_linkage is not None and params.show_dendrogram)
    try:
        dendro_ratio = float(params.dendrogram_width_ratio)
    except Exception:
        dendro_ratio = 0.18
    dendro_ratio = max(0.0, min(dendro_ratio, 5.0))
    dendro_width_norm = dendro_ratio / (1.0 + dendro_ratio) if have_dendro else 0.01
    
    # Position colorbar above heatmap with fixed spacing
    if show_cbar:
        cbar_bottom = heatmap_fraction + spacing_norm
        ax_cbar = fig.add_axes((dendro_width_norm + 0.005, cbar_bottom,
                       1.0 - dendro_width_norm - 0.01, cbar_height_norm * 0.8))
        heatmap_top = heatmap_fraction - 0.02
    else:
        ax_cbar = None
        heatmap_top = 0.98
    
    # Main heatmap axis
    ax_hm = fig.add_axes((dendro_width_norm + 0.005, 0.05,
                         1.0 - dendro_width_norm - 0.01, heatmap_top - 0.05))
    
    # Dendrogram axis
    ax_dendro = None
    if have_dendro:
        ax_dendro = fig.add_axes((0.005, 0.05, dendro_width_norm - 0.01, heatmap_top - 0.05))
        try:
            from scipy.cluster.hierarchy import dendrogram
            dendrogram(row_linkage, ax=ax_dendro, orientation='left', no_labels=True,
                       color_threshold=0, above_threshold_color='black', link_color_func=lambda k: 'black')
            ax_dendro.set_xticks([])
            ax_dendro.set_yticks([])
            for sp in ax_dendro.spines.values():
                sp.set_visible(False)
        except Exception as e:
            logger.warning(f"Could not draw dendrogram: {e}")
    
    return fig, ax_hm, ax_dendro, ax_cbar

def _draw_heatmap_content(ax, mat_ord: pd.DataFrame, row_names: pd.Series,
                         cols_pair: List[str], sample_to_group: dict, groups: List[str],
                         params: HeatmapParams, vmin: float, vmax: float,
                         log2fc_ordered=None):
    """Draw the heatmap image, labels, and splits. Returns the image object."""
    from matplotlib.colors import LinearSegmentedColormap
    
    cmap = LinearSegmentedColormap.from_list('gbr', ['green','black','red'])
    im = ax.imshow(mat_ord.values, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    
    # Labels
    ax.set_yticks(range(len(row_names)))
    ax.set_yticklabels(row_names, fontsize=params.feature_fontsize, fontweight='bold')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.set_xticks(range(len(mat_ord.columns)))
    ax.set_xticklabels(mat_ord.columns, rotation=90, fontsize=params.sample_fontsize, fontweight='bold')
    ax.set_title('')
    
    # Add visual splits
    _add_heatmap_splits(
        ax, mat_ord, cols_pair, sample_to_group, groups,
        log2fc_data=log2fc_ordered,
        add_row_split=params.show_fc_divider,  # FIXED: Respect the parameter setting
        add_col_split=not params.no_col_split
    )
    
    # DEBUG: Log if split was disabled via parameter
    if not params.show_fc_divider and log2fc_ordered is not None:
        logger.info(f"   Note: FC divider line is DISABLED in settings (show_fc_divider=False)")
    
    return im

def _order_and_cluster_rows(mat_z: pd.DataFrame, sub_df: pd.DataFrame, params: HeatmapParams,
                           log2fc_series: Optional[pd.Series], id_col: str
                           ):
    """Apply row ordering: cluster or FC sort. Returns: (mat_ord, sub_df_sorted, row_linkage, log2fc_ordered)."""
    row_linkage = None
    log2fc_ordered = None
    
    if params.cluster:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            row_linkage = linkage(mat_z.values, method='average', metric='euclidean')
            row_order = leaves_list(row_linkage)
            mat_ord = mat_z.iloc[row_order]
            sub_df = sub_df.iloc[row_order]
            # FIXED: Reorder log2fc data to match clustered order for split line
            if log2fc_series is not None:
                log2fc_ordered = log2fc_series.iloc[row_order].reset_index(drop=True)
        except Exception:
            mat_ord = mat_z
            log2fc_ordered = log2fc_series
    else:
        # Sort by log2FC descending (positive first, then negative) for clear split
        if log2fc_series is not None:
            sort_idx = log2fc_series.sort_values(ascending=False).index
            sub_df = sub_df.loc[sort_idx].reset_index(drop=True)
            mat_z = mat_z.loc[sort_idx].reset_index(drop=True)
            log2fc_series = log2fc_series.loc[sort_idx].reset_index(drop=True)
        mat_ord = mat_z
        log2fc_ordered = log2fc_series
    
    return mat_ord, sub_df, row_linkage, log2fc_ordered

def run_heatmap_analysis(ctx: CommonVizContext, params: HeatmapParams) -> VizResults:
    """Heatmap analysis - basic implementation."""
    files_created = []
    errors = []
    
    try:
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining
        # Directory de-duplication
        if os.path.basename(ctx.output_dir.rstrip(os.sep)).lower() == 'heatmaps':
            hm_dir = ctx.output_dir
        else:
            hm_dir = os.path.join(ctx.output_dir, 'heatmaps')
        os.makedirs(hm_dir, exist_ok=True)
        
        # Get identifier column from configuration (must be configured via Configure Stat Columns dialog)
        id_col = None
        if ctx.stat_column_assignments and ctx.stat_column_assignments.get('id_column'):
            id_col = ctx.stat_column_assignments.get('id_column')
        elif ctx.id_column:
            id_col = ctx.id_column
        
        if not id_col or id_col not in ctx.complete_df.columns:
            # Try fallback: use first non-numeric, non-statistical column
            available_cols = [c for c in ctx.complete_df.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
            
            # Try to find a suitable fallback ID column
            fallback_id = None
            for col in available_cols:
                try:
                    pd.to_numeric(ctx.complete_df[col], errors='raise')
                    continue  # Skip numeric columns
                except (ValueError, TypeError):
                    # This is a non-numeric column - good candidate
                    fallback_id = col
                    break
            
            if fallback_id:
                logger.warning(f"   Heatmap: Using fallback ID column '{fallback_id}' (not explicitly configured)")
                id_col = fallback_id
            else:
                logger.error(f"   Heatmap: ID column not configured. Available non-statistical columns: {available_cols[:20]}")
                errors.append(f"CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use \"Configure Stat Columns\" dialog to specify the ID column.")
                return VizResults(files_created, errors, "Heatmap analysis failed: missing ID column configuration")
        
        # Simple implementation: filter significant metabolites and create basic heatmap
        from itertools import combinations
        import matplotlib as mpl
        from matplotlib.colors import LinearSegmentedColormap
        
        # Interpret params.fc_threshold <= 0 as 'no fold-change cutoff' (p-value only)
        heatmap_use_fc = bool(params.fc_threshold and params.fc_threshold > 0)
        if heatmap_use_fc:
            log2fc_threshold = float(np.log2(params.fc_threshold)) if params.fc_threshold > 1 else 0.0
        else:
            log2fc_threshold = None
        
        # Preserve explicit group order by index pairing
        ordered_groups = ctx.groups
        pair_iter = [(ordered_groups[i], ordered_groups[j]) for i in range(len(ordered_groups)) for j in range(i+1, len(ordered_groups))]

        significant_union_ids: Set[str] = set()
        pair_sig_maps = {}  # (g1,g2) -> ordered list of metabolite ids

        # Select which pairs to consider based on filter_mode and selected_comparisons
        selected_pairs = pair_iter
        if params.filter_mode == 'specific' and params.filter_pairs:
            wanted = {tuple(sorted(p)) for p in params.filter_pairs}
            selected_pairs = [(g1,g2) for (g1,g2) in pair_iter if tuple(sorted((g1,g2))) in wanted]
        
        # Apply selected_comparisons filter
        if params.selected_comparisons is not None:
            selected_pairs = [
                (g1, g2) for (g1, g2) in selected_pairs
                if (g1, g2) in params.selected_comparisons or (g2, g1) in params.selected_comparisons
            ]

        # Helper to compute significant dataframe for a pair
        def compute_sig_for_pair(g1, g2):
            # Check for per-comparison metabolite list
            comp_metab_list = _get_metabolites_for_comparison(params, g1, g2)
            
            # If skip_unlisted_comparisons is True and no list found, skip this comparison
            if params.skip_unlisted_comparisons and comp_metab_list is None:
                return None, None, "No metabolite list for this comparison (skipped)"
            
            # Determine if we're using a custom list (no statistical filtering)
            using_custom_list = (comp_metab_list is not None) or (params.use_custom_only and params.include_metabolites)
            
            # Determine preference for adjusted p-values from context
            prefer_adj = True
            if hasattr(ctx, 'use_adj_p') and ctx.use_adj_p is not None:
                prefer_adj = bool(ctx.use_adj_p)
            
            # Locate stats columns for this comparison - uses Configure Stat Columns assignments first
            pair_cols = _locate_pair_columns(ctx.complete_df, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=getattr(ctx, 'stat_column_assignments', None))
            if pair_cols:
                log2fc_col, p_col, _ = pair_cols
            else:
                # Attempt synthetic generation of pairwise stats if missing
                log2fc_col, p_col = None, None
                try:
                    g1_samples = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) == g1 and c in ctx.complete_df.columns]
                    g2_samples = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) == g2 and c in ctx.complete_df.columns]
                    if g1_samples and g2_samples:
                        # Only build if not already present
                        synth_fc = f"{g1}_vs_{g2}_log2FC"
                        synth_p = f"{g1}_vs_{g2}_p_value"  # raw p-values (Welch t-test if scipy available)
                        if synth_fc not in ctx.complete_df.columns or synth_p not in ctx.complete_df.columns:
                            import numpy as _np
                            g1_mat = ctx.complete_df[g1_samples].astype(float).values
                            g2_mat = ctx.complete_df[g2_samples].astype(float).values
                            # Means with small pseudocount to avoid division by zero
                            g1_mean = _np.nanmean(g1_mat, axis=1)
                            g2_mean = _np.nanmean(g2_mat, axis=1)
                            fc_vals = _np.log2((g2_mean + 1e-9) / (g1_mean + 1e-9))
                            p_vals = _np.ones_like(fc_vals)
                            try:
                                from scipy import stats as _stats
                                # Welch t-test per metabolite
                                # Vectorized not directly available; loop with early break if large
                                n_rows = len(fc_vals)
                                if n_rows <= 5000:
                                    p_list = []
                                    for i in range(n_rows):
                                        row1 = g1_mat[i, :]
                                        row2 = g2_mat[i, :]
                                        res = _stats.ttest_ind(row2, row1, nan_policy='omit', equal_var=False)
                                        p_val = getattr(res, 'pvalue', getattr(res, 'p', 1.0))
                                        p_list.append(p_val if p_val is not None else 1.0)
                                    p_vals = _np.array(p_list)
                                else:
                                    # For very large matrices, approximate using variance formula
                                    var1 = _np.nanvar(g1_mat, axis=1)
                                    var2 = _np.nanvar(g2_mat, axis=1)
                                    n1 = _np.sum(_np.isfinite(g1_mat), axis=1)
                                    n2 = _np.sum(_np.isfinite(g2_mat), axis=1)
                                    # Welch t statistic
                                    t_stat = (g2_mean - g1_mean) / _np.sqrt(var1 / n1 + var2 / n2 + 1e-12)
                                    # Approximate df (Welch–Satterthwaite)
                                    df_num = (var1 / n1 + var2 / n2) ** 2
                                    df_den = (var1**2 / (n1**2 * (n1 - 1 + 1e-9))) + (var2**2 / (n2**2 * (n2 - 1 + 1e-9)))
                                    df = df_num / (df_den + 1e-12)
                                    try:
                                        from scipy.stats import t as _t_dist
                                        p_vals = 2 * _t_dist.sf(_np.abs(t_stat), df)
                                    except Exception:
                                        p_vals = _np.ones_like(t_stat)
                            except Exception:
                                # scipy not available; keep p_vals as ones so filtering only by FC if requested
                                pass
                            ctx.complete_df[synth_fc] = fc_vals
                            ctx.complete_df[synth_p] = p_vals
                            logger.warning(
                                f"   Synthetic stats created for {g1} vs {g2}: columns '{synth_fc}', '{synth_p}' (raw p-values; consider running formal stats)."
                            )
                        log2fc_col, p_col = synth_fc, synth_p
                    else:
                        logger.warning(f"   Cannot synthesize stats for {g1} vs {g2} (missing sample columns for one or both groups)")
                except Exception as _e:
                    logger.warning(f"   Failed synthetic stat generation for {g1} vs {g2}: {_e}")
            
            if using_custom_list:
                # Using custom/comparison-specific list - NO statistical filtering applied
                if comp_metab_list is not None:
                    # Use comparison-specific list
                    match_mask = match_metabolites_multi_column(ctx.complete_df, comp_metab_list, id_col)
                    df_sub = ctx.complete_df[match_mask].copy()
                    logger.info(f"   Heatmap {g1} vs {g2}: Using comparison-specific list, {len(df_sub)} metabolites matched (NO p-value/FC filtering)")
                else:
                    # Use custom list
                    match_mask = match_metabolites_multi_column(ctx.complete_df, params.include_metabolites, id_col)
                    df_sub = ctx.complete_df[match_mask].copy()
                    logger.info(f"   Heatmap {g1} vs {g2}: Using custom list only, {len(df_sub)} metabolites matched (NO p-value/FC filtering)")
                
                if df_sub.empty:
                    logger.warning(f"   Heatmap {g1} vs {g2}: No metabolites from custom list found in data!")
                    return None, None, "No metabolites from custom list found in data"
                
                # Sort by p-value if available (for display order, NOT for filtering)
                if p_col is not None and p_col in df_sub.columns:
                    df_sub = df_sub.sort_values(p_col, ascending=True)
                
                # Apply max_metabolites limit if specified
                if params.max_metabolites > 0:
                    df_sub = df_sub.head(params.max_metabolites)
                
                return df_sub, (log2fc_col, p_col), None
            
            else:
                # NOT using custom list - apply statistical filtering
                if not pair_cols:
                    # Provide specific error about what's missing
                    expected_log2fc = f"{g1}_vs_{g2}_log2FC"
                    expected_pval = f"{g1}_vs_{g2}_pvalue or {g1}_vs_{g2}_adj_p"
                    actual_cols = [c for c in ctx.complete_df.columns if '_vs_' in c and (g1 in c or g2 in c)]
                    if actual_cols:
                        hint = f"Found similar columns: {actual_cols[:5]}"
                    else:
                        hint = "No pairwise columns found for this comparison"
                    return None, None, f"Missing stat columns for {g1}_vs_{g2}. Expected: {expected_log2fc} and {expected_pval}. {hint}. Use 'Configure Stat Columns' to map columns."
                
                # Use p-value/FC threshold filtering
                if log2fc_threshold is None:
                    mask = (ctx.complete_df[p_col] < params.p_threshold)
                else:
                    mask = (ctx.complete_df[p_col] < params.p_threshold) & (ctx.complete_df[log2fc_col].abs() >= log2fc_threshold)
                
                df_sub = ctx.complete_df[mask].copy()
                logger.info(f"   Heatmap {g1} vs {g2}: After p-value/FC filter: {len(df_sub)} metabolites")
                
                if df_sub.empty:
                    logger.warning(f"   Heatmap {g1} vs {g2}: No metabolites pass statistical filters!")
                    return None, None, "No metabolites pass filter"
                
                p_col_name = p_col if isinstance(p_col, str) else None
                if p_col_name is not None:
                    df_sub = df_sub.sort_values(by=p_col_name, ascending=True)
                
                if params.max_metabolites > 0:
                    df_sub = df_sub.head(params.max_metabolites)
                
                return df_sub, (log2fc_col, p_col), None

        per_pair_records = []  # for summary

        for g1, g2 in selected_pairs:
            df_sig, pair_cols, err_msg = compute_sig_for_pair(g1, g2)
            if df_sig is None:
                errors.append(f"Skipping heatmap {g1} vs {g2}: {err_msg}")
                continue
            row_ids = df_sig[id_col].astype(str).tolist()
            pair_sig_maps[(g1,g2)] = row_ids
            per_pair_records.append({'pair': f"{g1} vs {g2}", 'count': len(row_ids)})

        # Aggregate set based on filter_mode (for combined heatmap logic later)
        if pair_sig_maps:
            if params.filter_mode == 'all' and selected_pairs:
                sets = [set(pair_sig_maps[p]) for p in pair_sig_maps]
                aggregate_ids = set.intersection(*sets) if sets else set()
            else:  # any or specific
                aggregate_ids = set().union(*[set(v) for v in pair_sig_maps.values()])
        else:
            aggregate_ids = set()
        significant_union_ids.update(aggregate_ids)

        # Helper to order sample columns by group while preserving the user's
        # original per-sample ordering within each group
        def _order_cols_by_groups(cols: List[str], groups_seq: List[str]) -> List[str]:
            ordered: List[str] = []
            for g in groups_seq:
                ordered.extend([c for c in ctx.sample_cols if c in cols and ctx.sample_to_group.get(c) == g])
            # Include any leftovers (safety)
            leftovers = [c for c in cols if c not in ordered]
            return ordered + leftovers

        # Now generate per-pair heatmaps using filtered lists (pair_sig_maps)
        for (g1,g2), id_list in pair_sig_maps.items():
            # Respect group ordering: all samples from g1 followed by all from g2
            cols_pair = [c for c in ctx.sample_cols if ctx.sample_to_group.get(c) in (g1, g2) and c in ctx.complete_df.columns]
            cols_pair = _order_cols_by_groups(cols_pair, [g1, g2])
            if not cols_pair:
                errors.append(f"No sample columns for {g1} vs {g2}")
                continue
            sub_df = ctx.complete_df[ctx.complete_df[id_col].astype(str).isin(id_list)].copy()
            if sub_df.empty:
                errors.append(f"No data matrix for {g1} vs {g2}")
                continue

            # Compute an effect z-score per row to use for splitting (mean difference / per-row std)
            # This is more robust when rows are z-scored and fold-change columns are unreliable.
            cols_g1 = [c for c in cols_pair if ctx.sample_to_group.get(c) == g1]
            cols_g2 = [c for c in cols_pair if ctx.sample_to_group.get(c) == g2]
            effect_z = None
            try:
                if cols_g1 and cols_g2:
                    mean_g1 = sub_df[cols_g1].apply(pd.to_numeric, errors='coerce').mean(axis=1)
                    mean_g2 = sub_df[cols_g2].apply(pd.to_numeric, errors='coerce').mean(axis=1)
                    pooled_std = sub_df[cols_pair].apply(pd.to_numeric, errors='coerce').std(axis=1)
                    # Avoid zero division
                    pooled_std = pooled_std.replace(0, np.nan).fillna(1.0)
                    effect_z = (mean_g2 - mean_g1) / pooled_std
                    
                    # Debug logging for split line positioning
                    pos_count = (effect_z > 0).sum()
                    neg_count = (effect_z < 0).sum()
                    zero_count = (effect_z == 0).sum()
                    logger.debug(f"   Heatmap {g1} vs {g2} effect_z: "
                               f"{pos_count} upregulated, {neg_count} downregulated, {zero_count} no change "
                               f"(range: {effect_z.min():.3f} to {effect_z.max():.3f})")
                else:
                    effect_z = None
                    logger.warning(f"   Heatmap {g1} vs {g2}: Could not compute effect_z (missing group columns)")
            except Exception as e:
                effect_z = None
                logger.warning(f"   Heatmap {g1} vs {g2}: Error computing effect_z: {e}")

            # Prepare matrix
            mat = sub_df[cols_pair].apply(pd.to_numeric, errors='coerce').fillna(0)
            if params.scale == 'row':
                mat_z = (mat.T - mat.mean(axis=1)).T
                stds = mat.std(axis=1).replace(0, np.nan)
                mat_z = (mat_z.T / stds).T.fillna(0)
            else:
                mat_z = mat
            
            # Order and cluster rows using helper function (pass effect_z so splitting can use it)
            mat_ord, sub_df, row_linkage, log2fc_ordered = _order_and_cluster_rows(
                mat_z, sub_df, params, effect_z, id_col
            )
            row_names = sub_df[id_col].astype(str)
            n_samples = len(mat_ord.columns)
            n_features = len(mat_ord)
            
            # Setup figure and axes using helper function
            fig, ax, dendro_ax, ax_cbar = _setup_heatmap_figure(n_samples, n_features, params, row_linkage)
            
            # Prepare color scale - respect user's choice hierarchy
            # Priority: 1) Fixed scale (-3 to 3), 2) Auto scale, 3) Manual scale
            if getattr(params, 'use_fixed_scale', True):
                # Default: use fixed -3 to 3 scale
                vmin, vmax = -3.0, 3.0
                logger.debug(f"   Using fixed scale: {vmin} to {vmax}")
            elif params.auto_scale:
                # Auto scale to 5th-95th percentile
                data_vals = mat_ord.values.flatten()
                finite_vals = data_vals[np.isfinite(data_vals)]
                if finite_vals.size > 0:
                    vmin = np.percentile(finite_vals, 5)
                    vmax = np.percentile(finite_vals, 95)
                    if vmin >= vmax:
                        vmin, vmax = finite_vals.min(), finite_vals.max()
                    logger.debug(f"   Auto scale: {vmin:.2f} to {vmax:.2f}")
                else:
                    vmin, vmax = -3, 3
                    logger.debug(f"   Auto scale failed, using default: {vmin} to {vmax}")
            else:
                # Manual scale from user input
                vmin, vmax = params.vmin, params.vmax
                logger.debug(f"   Manual scale: {vmin} to {vmax}")
            
            # Draw heatmap content using helper function
            im = _draw_heatmap_content(
                ax, mat_ord, row_names, cols_pair, ctx.sample_to_group, [g1, g2],
                params, vmin, vmax, log2fc_ordered
            )
            
            # Move colorbar to the top, horizontal orientation using dedicated axis
            show_cbar = getattr(params, 'show_colorbar', True)
            if show_cbar and ax_cbar is not None:
                # Colorbar: use fixed ticks; when using fixed z-range show -3/0/3
                cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
                cbar.set_ticks([vmin, 0.0, vmax])
                cbar.set_ticklabels([f'{vmin:.0f}', '0', f'{vmax:.0f}'])
                cbar.ax.tick_params(labelsize=max(14, params.sample_fontsize+2))  # Increased font size for tick labels
                cbar.ax.xaxis.set_ticks_position('top')  # Put ticks on top
                cbar.ax.xaxis.set_label_position('top')
                
                try:
                    for lbl in cbar.ax.get_xticklabels():
                        lbl.set_fontweight('bold')
                except Exception:
                    pass
            # Title positioned above colorbar (which is now above y=1.0)
            # Recalculate spacing for title positioning
            try:
                cbar_height_inches = float(getattr(params, 'colorbar_height_inches', 0.6))
            except Exception:
                cbar_height_inches = 0.6
            dyn_height = fig.get_figheight()
            cbar_height_norm = cbar_height_inches / dyn_height
            spacing_norm = 0.15 / dyn_height
            title_y = 1.0 + spacing_norm + cbar_height_norm + 0.02 if show_cbar else 1.01
            fig.suptitle(f'Heatmap: {g1} vs {g2}', fontweight='bold', fontsize=params.title_fontsize, y=title_y)
            tag = f"{g1}_vs_{g2}".replace(' ', '_')
            out_png = os.path.join(hm_dir, f'HM_{tag}.png')
            plt.savefig(out_png, dpi=params.fig_dpi, bbox_inches='tight')
            files_created.append(out_png)
            
            # Save metabolite list (controlled by save_excel parameter)
            if params.save_excel:
                list_path = os.path.join(hm_dir, f'HM_{tag}_metabolites.csv')
                # Include log2FC and p-value for verification of split and significance
                # Get pair columns info for this comparison
                # Determine preference for adjusted p-values from context
                prefer_adj = True
                if hasattr(ctx, 'use_adj_p') and ctx.use_adj_p is not None:
                    prefer_adj = bool(ctx.use_adj_p)
                pair_cols_info = _locate_pair_columns(ctx.complete_df, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                p_col = pair_cols_info[1] if pair_cols_info else None
                try:
                    df_out = pd.DataFrame({
                        'metabolite_id': sub_df[id_col].astype(str).values,
                        'effect_z': (log2fc_ordered.values if log2fc_ordered is not None else np.full(len(sub_df), np.nan)),
                        'p_value': (sub_df[p_col].values if (p_col is not None and p_col in sub_df.columns) else np.full(len(sub_df), np.nan)),
                    })
                    df_out['Regulation'] = np.where(df_out['effect_z'] > 0, 'Upregulated', np.where(df_out['effect_z'] < 0, 'Downregulated', 'No change'))
                    df_out['Rank'] = np.arange(1, len(df_out) + 1)
                    df_out.to_csv(list_path, index=False)
                except Exception:
                    # Fallback to simple export if any issue
                    pd.DataFrame({'metabolite_id': row_names}).to_csv(list_path, index=False)
                files_created.append(list_path)
            plt.close()
        
        # Combined union heatmap if requested
        if params.combined and significant_union_ids:
            mode = getattr(params, 'combined_mode', 'union').lower()
            # Order all sample columns by explicit group order
            sample_cols_all = _order_cols_by_groups([c for c in ctx.sample_cols if c in ctx.complete_df.columns], ctx.groups)

            if mode == 'intersection':
                # Only metabolites present (significant) in all pair comparisons that had at least one sig metabolite
                involved_pairs = [p for p, ids in pair_sig_maps.items() if ids]
                if involved_pairs:
                    intersect_ids = set(pair_sig_maps[involved_pairs[0]])
                    for p in involved_pairs[1:]:
                        intersect_ids &= set(pair_sig_maps[p])
                else:
                    intersect_ids = set()
                candidate_ids = list(intersect_ids)
            elif mode == 'concatenate':
                # Concatenate metabolites per pair in group order, de-duplicating while preserving first appearance
                seen = set()
                ordered_concat = []
                for (g1, g2) in pair_iter:
                    ids_list = pair_sig_maps.get((g1, g2), [])
                    for m in ids_list:
                        if m not in seen:
                            seen.add(m)
                            ordered_concat.append(m)
                candidate_ids = ordered_concat
            else:  # union (default)
                # Order by minimum p-value across pairs
                # Determine preference for adjusted p-values from context
                prefer_adj = True
                if hasattr(ctx, 'use_adj_p') and ctx.use_adj_p is not None:
                    prefer_adj = bool(ctx.use_adj_p)
                
                pval_records = []
                for metab in significant_union_ids:
                    min_p = np.inf
                    for (g1, g2), ids_list in pair_sig_maps.items():
                        if metab in ids_list:
                            pair_cols = _locate_pair_columns(ctx.complete_df, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                            if pair_cols:
                                _, p_col, _ = pair_cols
                                row = ctx.complete_df[ctx.complete_df[id_col].astype(str) == metab]
                                if not row.empty and p_col in row.columns:
                                    try:
                                        pval = float(row.iloc[0][p_col])
                                        if pval < min_p:
                                            min_p = pval
                                    except Exception:
                                        pass
                    if not np.isfinite(min_p):
                        min_p = 1.0
                    pval_records.append((metab, min_p))
                pval_records.sort(key=lambda x: x[1])
                candidate_ids = [m for m, _ in pval_records]

            # Apply max limit
            if params.max_metabolites > 0:
                candidate_ids = candidate_ids[:params.max_metabolites]

            if candidate_ids:
                sub = ctx.complete_df[ctx.complete_df[id_col].astype(str).isin(candidate_ids)].copy()

                # Use the FIRST comparison's log2FC for sorting (to maintain consistent direction)
                # Get the first pair from selected pairs
                first_pair = selected_pairs[0] if selected_pairs else None
                log2fc_for_sorting = {}
                
                # Determine preference for adjusted p-values from context
                prefer_adj = True
                if hasattr(ctx, 'use_adj_p') and ctx.use_adj_p is not None:
                    prefer_adj = bool(ctx.use_adj_p)
                
                if first_pair:
                    g1, g2 = first_pair
                    pair_cols_info = _locate_pair_columns(ctx.complete_df, g1, g2, prefer_adj=prefer_adj, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                    if pair_cols_info:
                        fc_col, _, _ = pair_cols_info
                        if fc_col in ctx.complete_df.columns:
                            for metab in candidate_ids:
                                row = ctx.complete_df[ctx.complete_df[id_col].astype(str) == metab]
                                if not row.empty:
                                    try:
                                        fc_val = float(row.iloc[0][fc_col])
                                        if np.isfinite(fc_val):
                                            log2fc_for_sorting[metab] = fc_val
                                        else:
                                            log2fc_for_sorting[metab] = 0.0
                                    except Exception:
                                        log2fc_for_sorting[metab] = 0.0
                                else:
                                    log2fc_for_sorting[metab] = 0.0
                else:
                    # Fallback to average if no pairs available
                    for metab in candidate_ids:
                        fc_values = []
                        for (g1, g2) in pair_sig_maps.keys():
                            pair_cols_info = _locate_pair_columns(ctx.complete_df, g1, g2, verified_assignments=ctx.verified_assignments, stat_column_assignments=ctx.stat_column_assignments)
                            if pair_cols_info:
                                fc_col, _, _ = pair_cols_info
                                if fc_col in ctx.complete_df.columns:
                                    row = ctx.complete_df[ctx.complete_df[id_col].astype(str) == metab]
                                    if not row.empty:
                                        try:
                                            fc_val = float(row.iloc[0][fc_col])
                                            if np.isfinite(fc_val):
                                                fc_values.append(fc_val)
                                        except Exception:
                                            pass
                        log2fc_for_sorting[metab] = np.mean(fc_values) if fc_values else 0.0

                mat = sub[sample_cols_all].apply(pd.to_numeric, errors='coerce')
                mat = mat.fillna(mat.median())
                
                # Add log2fc to sub for sorting
                sub['__log2fc__'] = [log2fc_for_sorting.get(m, 0.0) for m in sub[id_col].astype(str)]
                # Sort by log2fc descending (positive first, then negative)
                sub = sub.sort_values('__log2fc__', ascending=False).reset_index(drop=True)
                log2fc_sorted = sub['__log2fc__']
                
                # Re-align mat with sorted sub
                mat = sub[sample_cols_all].apply(pd.to_numeric, errors='coerce')
                mat = mat.fillna(mat.median())

                # Convert average log2FC to an effect z-score by dividing by per-row std
                try:
                    per_row_std = mat.std(axis=1).replace(0, np.nan).fillna(1.0)
                    effect_z_combined = sub['__log2fc__'] / per_row_std
                    
                    # Debug logging for split line positioning
                    pos_count = (effect_z_combined > 0).sum()
                    neg_count = (effect_z_combined < 0).sum()
                    zero_count = (effect_z_combined == 0).sum()
                    logger.debug(f"   Combined heatmap ({mode}) effect_z: "
                               f"{pos_count} upregulated, {neg_count} downregulated, {zero_count} no change "
                               f"(range: {effect_z_combined.min():.3f} to {effect_z_combined.max():.3f})")
                except Exception as e:
                    effect_z_combined = sub['__log2fc__']
                    logger.warning(f"   Combined heatmap ({mode}): Could not normalize effect_z: {e}")
                
                if params.scale == 'row':
                    mat_z = (mat.T - mat.mean(axis=1)).T
                    stds = mat.std(axis=1).replace(0, np.nan)
                    mat_z = (mat_z.T / stds).T.fillna(0)
                else:
                    mat_z = mat
                
                # FIXED: Define title_mode BEFORE clustering/drawing section
                title_mode = mode.capitalize()
                
                # Order and cluster rows using helper function 
                # IMPORTANT: Pass raw log2fc_sorted for split line, not effect_z which can have sign issues
                mat_ord, sub, row_linkage, log2fc_combined = _order_and_cluster_rows(
                    mat_z, sub, params, log2fc_sorted, id_col
                )
                row_names = sub[id_col].astype(str)
                n_samples = len(mat_ord.columns)
                n_features = len(mat_ord)
                
                # Setup figure and axes using helper function
                fig, ax, dendro_ax, ax_cbar = _setup_heatmap_figure(n_samples, n_features, params, row_linkage)
                
                # Prepare color scale
                use_fixed_zrange = (getattr(params, 'scale', None) == 'row') or getattr(params, 'force_z_range', False)
                if use_fixed_zrange:
                    vmin, vmax = -3.0, 3.0
                else:
                    if params.auto_scale:
                        data_vals = mat_ord.values.flatten()
                        finite_vals = data_vals[np.isfinite(data_vals)]
                        if finite_vals.size > 0:
                            vmin = np.percentile(finite_vals, 5)
                            vmax = np.percentile(finite_vals, 95)
                            if vmin >= vmax:
                                vmin, vmax = finite_vals.min(), finite_vals.max()
                        else:
                            vmin, vmax = -3, 3
                    else:
                        vmin, vmax = params.vmin, params.vmax
                
                # Draw heatmap content using helper function
                im = _draw_heatmap_content(
                    ax, mat_ord, row_names, sample_cols_all, ctx.sample_to_group, ctx.groups,
                    params, vmin, vmax, log2fc_combined
                )
                
                # Move colorbar to the top (dedicated axis)
                if getattr(params, 'show_colorbar', True) and ax_cbar is not None:
                    cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
                    # If fixed z-range is used, show -3/0/3 ticks; otherwise show computed extremes
                    if (getattr(params, 'scale', None) == 'row') or getattr(params, 'force_z_range', False):
                        cbar.set_ticks([vmin, 0.0, vmax])
                        cbar.set_ticklabels([f'{vmin:.0f}', '0', f'{vmax:.0f}'])
                    else:
                        cbar.set_ticks([vmin, 0.0, vmax])
                        cbar.set_ticklabels([f'{vmin:.2f}', '0', f'{vmax:.2f}'])
                    cbar.ax.tick_params(labelsize=max(14, params.sample_fontsize+2))  # Increased font size for tick labels
                    cbar.ax.xaxis.set_ticks_position('top')
                    cbar.ax.xaxis.set_label_position('top')
                    # Removed log2(FC) label as requested
                    try:
                        for lbl in cbar.ax.get_xticklabels():
                            lbl.set_fontweight('bold')
                    except Exception:
                        pass
                # Title positioned above colorbar
                title_y = 1.0 + spacing_norm + cbar_height_norm + 0.02 if show_cbar else 1.01
                fig.suptitle(f'Heatmap: Combined ({title_mode})', fontweight='bold', fontsize=params.title_fontsize, y=title_y)
                out_png = os.path.join(hm_dir, f'HM_combined_{mode}.png')
                plt.savefig(out_png, dpi=params.fig_dpi, bbox_inches='tight')
                files_created.append(out_png)

                # Save combined metabolite list with the log2FC used for sorting
                if params.save_excel:
                    try:
                        list_path = os.path.join(hm_dir, f'HM_combined_{mode}_metabolites.csv')
                        df_out = pd.DataFrame({
                            'metabolite_id': sub[id_col].astype(str).values,
                            'log2FC': sub['__log2fc__'].values,
                        })
                        df_out['Regulation'] = np.where(df_out['log2FC'] > 0, 'Upregulated', np.where(df_out['log2FC'] < 0, 'Downregulated', 'No change'))
                        df_out['Rank'] = np.arange(1, len(df_out) + 1)
                        df_out.to_csv(list_path, index=False)
                        files_created.append(list_path)
                    except Exception:
                        pass
                
                # Save metabolite list (controlled by save_excel parameter)
                if params.save_excel:
                    list_path = os.path.join(hm_dir, f'HM_combined_{mode}_metabolites.csv')
                    pd.DataFrame({'metabolite_id': row_names}).to_csv(list_path, index=False)
                    files_created.append(list_path)
                plt.close()

        # Write summary file for transparency
        try:
            summary_lines = [
                f"Total metabolites: {len(ctx.complete_df)}",
                f"Pairs considered: {len(selected_pairs)}",
                f"Filter mode: {params.filter_mode}",
                f"Use custom only: {params.use_custom_only}",
                f"p_threshold: {params.p_threshold}",
                f"fc_threshold: {params.fc_threshold}",
                f"Metabolites in aggregate set: {len(significant_union_ids)}"
            ]
            if per_pair_records:
                summary_lines.append("Per-pair counts:")
                for rec in per_pair_records:
                    summary_lines.append(f"  {rec['pair']}: {rec['count']}")
            with open(os.path.join(hm_dir, 'heatmap_filter_summary.txt'), 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(summary_lines))
        except Exception:
            pass
        
        if files_created and not errors:
            summary = f"Heatmap analysis complete: {len(files_created)} files generated"
        elif files_created and errors:
            summary = f"Heatmap analysis partial: {len(files_created)} files, {len(errors)} errors"
        else:
            summary = f"Heatmap analysis failed: {len(errors)} errors"
        
    except Exception as e:
        errors.append(f"Heatmap analysis failed: {str(e)}")
        summary = "Heatmap analysis failed"
    
    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def run_roc_analysis(ctx: CommonVizContext, params: ROCParams) -> VizResults:
    """ROC analysis."""
    files_created = []
    errors = []
    
    try:
        if ctx.preferred_group_order:
            ordered = [g for g in ctx.preferred_group_order if g in ctx.groups]
            remaining = [g for g in ctx.groups if g not in ordered]
            ctx.groups = ordered + remaining
        from itertools import combinations
        from sklearn.linear_model import LogisticRegression
        
        # Directory de-duplication
        if os.path.basename(ctx.output_dir.rstrip(os.sep)).lower() == 'roc':
            roc_dir = ctx.output_dir
        else:
            roc_dir = os.path.join(ctx.output_dir, 'roc')
        os.makedirs(roc_dir, exist_ok=True)
        
        # Get identifier column from configuration (must be configured via Configure Stat Columns dialog)
        id_col = None
        if ctx.stat_column_assignments and ctx.stat_column_assignments.get('id_column'):
            id_col = ctx.stat_column_assignments.get('id_column')
        elif ctx.id_column:
            id_col = ctx.id_column
        
        if not id_col or id_col not in ctx.complete_df.columns:
            # Try fallback: use first non-numeric, non-statistical column
            available_cols = [c for c in ctx.complete_df.columns if not any(g in str(c) for g in ['_vs_', '_log2FC', '_adj_p', '_p_value', '_pvalue'])]
            
            fallback_id = None
            for col in available_cols:
                try:
                    pd.to_numeric(ctx.complete_df[col], errors='raise')
                    continue  # Skip numeric columns
                except (ValueError, TypeError):
                    fallback_id = col
                    break
            
            if fallback_id:
                logger.warning(f"   ROC: Using fallback ID column '{fallback_id}' (not explicitly configured)")
                id_col = fallback_id
            else:
                logger.error(f"   ROC: ID column not configured. Available non-statistical columns: {available_cols[:20]}")
                errors.append(f"CONFIGURATION ERROR: ID column not configured. Available columns: {", ".join(available_cols[:10])}... Use \"Configure Stat Columns\" dialog to specify the ID column.")
                return VizResults(files_created, errors, "ROC analysis failed: missing ID column configuration")
        
        if params.all_pairs:
            import itertools
            # Select pairs
            all_pairs = list(combinations(ctx.groups, 2))
            selected_pairs = all_pairs
            if params.filter_mode == 'specific' and params.filter_pairs:
                wanted = {tuple(sorted(p)) for p in params.filter_pairs}
                selected_pairs = [p for p in all_pairs if tuple(sorted(p)) in wanted]
            
            # Apply selected_comparisons filter
            if params.selected_comparisons is not None:
                selected_pairs = [
                    (g1, g2) for (g1, g2) in selected_pairs
                    if (g1, g2) in params.selected_comparisons or (g2, g1) in params.selected_comparisons
                ]
            
            # Filtering logic: treat fc_threshold <= 0 as 'no fold-change cutoff' (p-value only)
            roc_use_fc = bool(params.fc_threshold and params.fc_threshold > 0)
            if roc_use_fc:
                log2fc_thresh = np.log2(params.fc_threshold) if params.fc_threshold > 1 else 0.0
            else:
                log2fc_thresh = None
            summary_lines = [f"Total metabolites: {len(ctx.complete_df)}",
                             f"Pairs considered: {len(selected_pairs)}",
                             f"Filter mode: {params.filter_mode}",
                             f"Use custom only: {params.use_custom_only}",
                             f"Max metabolites: {params.max_metabolites}",
                             f"Excel only: {params.excel_only}",
                             f"p_threshold: {params.p_threshold}",
                             f"fc_threshold: {params.fc_threshold}"]
            logger.info(f"   ROC: max_metabolites={params.max_metabolites}, use_custom_only={params.use_custom_only}, excel_only={params.excel_only}")
            per_pair_counts = []
            for g1, g2 in selected_pairs:
                # Check for per-comparison metabolite list
                comp_metab_list = _get_metabolites_for_comparison(params, g1, g2)
                
                # If skip_unlisted_comparisons is True and no list found, skip this comparison
                if params.skip_unlisted_comparisons and comp_metab_list is None:
                    continue
                
                cols_pos = [c for c in ctx.sample_cols if ctx.sample_to_group[c] == g2]
                cols_neg = [c for c in ctx.sample_cols if ctx.sample_to_group[c] == g1]
                if len(cols_pos) == 0 or len(cols_neg) == 0:
                    errors.append(f"Insufficient samples for ROC {g1} vs {g2}")
                    continue
                
                # Metabolite selection - prioritize per-comparison list
                if comp_metab_list is not None:
                    metabolites = comp_metab_list
                    # Custom list: no max_metabolites limit
                elif params.use_custom_only and params.metabolites:
                    metabolites = params.metabolites
                    # Custom list: no max_metabolites limit
                else:
                    # Find pairwise stat columns (and ignore custom list when use_custom_only is False)
                    base = f"{g2}_vs_{g1}"
                    p_col = f"{base}_adj_p"
                    fc_col = f"{base}_log2FC"
                    if p_col in ctx.complete_df.columns and fc_col in ctx.complete_df.columns:
                        if log2fc_thresh is None:
                            mask = (ctx.complete_df[p_col] < params.p_threshold)
                        else:
                            mask = (ctx.complete_df[p_col] < params.p_threshold) & (ctx.complete_df[fc_col].abs() >= log2fc_thresh)
                        df_sel = ctx.complete_df[mask].copy()
                        top_df = df_sel.nsmallest(params.max_metabolites, p_col)
                        metabolites = top_df[id_col].astype(str).tolist()
                    else:
                        # Fallback: variance
                        sample_cols_pair = cols_pos + cols_neg
                        mat = ctx.complete_df[sample_cols_pair].apply(pd.to_numeric, errors='coerce')
                        ctx.complete_df['__var__'] = mat.var(axis=1, skipna=True)
                        top_df = ctx.complete_df.nlargest(params.max_metabolites, '__var__')
                        metabolites = top_df[id_col].astype(str).tolist()
                per_pair_counts.append({'pair': f"{g1} vs {g2}", 'count': len(metabolites)})
                if not metabolites:
                    errors.append(f"No metabolites selected for ROC {g1} vs {g2}")
                    continue
                y_true = np.array([1]*len(cols_pos) + [0]*len(cols_neg))
                plt.figure(figsize=(params.fig_width, params.fig_height))
                ax = plt.gca()
                
                # First pass: calculate all AUCs and store data
                metabolite_data = []  # Store (metabolite, auc_score, fpr, tpr)
                feature_matrix = []  # Store metabolite values for combined ROC (n_metabolites x n_samples)
                
                # When using custom list, process ALL metabolites (no limit)
                # When using statistical filtering, respect max_metabolites limit
                using_custom = (comp_metab_list is not None) or (params.use_custom_only and params.metabolites)
                if using_custom:
                    metabolites_to_process = metabolites  # No limit for custom lists
                else:
                    metabolites_to_process = metabolites[:params.max_metabolites] if params.max_metabolites > 0 else metabolites
                
                for metab in metabolites_to_process:
                    met_rows = ctx.complete_df[ctx.complete_df[id_col].astype(str) == str(metab)]
                    if met_rows.empty:
                        continue
                    met_row = met_rows.iloc[0]
                    values_pos = np.asarray(pd.to_numeric(met_row[cols_pos], errors='coerce').to_numpy(), dtype=float)
                    values_neg = np.asarray(pd.to_numeric(met_row[cols_neg], errors='coerce').to_numpy(), dtype=float)
                    scores = np.array(list(values_pos) + list(values_neg), dtype=float)

                    # Treat 0 as below-detection (missing) for ROC.
                    scores[scores == 0] = np.nan

                    # If all missing after treating zeros as NaN, skip.
                    finite = scores[np.isfinite(scores)]
                    if finite.size == 0:
                        continue

                    # Impute remaining NaNs with within-metabolite median.
                    if np.isnan(scores).any():
                        fill_val = float(np.nanmedian(scores))
                        if not np.isfinite(fill_val):
                            continue
                        scores = np.where(np.isnan(scores), fill_val, scores)
                    if np.allclose(scores, scores[0]):
                        continue
                    fpr, tpr, _ = roc_curve(y_true, scores)
                    auc_score = auc(fpr, tpr)
                    fpr2, tpr2, _ = roc_curve(y_true, -scores)
                    auc2 = auc(fpr2, tpr2)
                    if auc2 > auc_score:
                        fpr, tpr, auc_score = fpr2, tpr2, auc2
                    if auc_score < params.min_auc:
                        continue
                    metabolite_data.append((metab, auc_score, fpr, tpr))
                    feature_matrix.append(scores)  # Store original scores for logistic regression
                
                if len(metabolite_data) == 0:
                    errors.append(f"No valid ROC curves for {g1} vs {g2}")
                    plt.close()
                    continue
                
                # Sort by AUC (ascending order - lowest to highest)
                metabolite_data.sort(key=lambda x: x[1])
                
                # Excel-only mode: save AUC values and skip plot generation
                if params.excel_only:
                    auc_entries = [{'metabolite': metab, 'AUC': auc_score} for metab, auc_score, fpr, tpr in metabolite_data]
                    tag = f"{g2}_vs_{g1}".replace(' ', '_')
                    auc_path = os.path.join(roc_dir, f'roc_{tag}_AUC.csv')
                    pd.DataFrame(auc_entries).to_csv(auc_path, index=False)
                    files_created.append(auc_path)
                    per_pair_counts.append({'pair': f"{g1} vs {g2}", 'count': len(metabolite_data)})
                    plt.close()
                    continue
                
                # Second pass: plot in sorted order
                auc_entries = []
                for i, (metab, auc_score, fpr, tpr) in enumerate(metabolite_data):
                    color = sns.color_palette('tab20')[i % 20]
                    ax.plot(fpr, tpr, label=f"{metab} (AUC={auc_score:.2f})", color=color, linewidth=2.5, alpha=0.9)
                    auc_entries.append({'metabolite': metab, 'AUC': auc_score})
                
                kept_curves = len(metabolite_data)
                
                # Add combined ROC if requested and we have multiple metabolites (plot last)
                if params.include_combined and len(feature_matrix) > 1:
                    try:
                        # Build feature matrix: transpose so shape is (n_samples, n_metabolites)
                        X = np.array(feature_matrix).T  # Shape: (n_samples, n_metabolites)
                        
                        # Standardize features for better logistic regression performance
                        from sklearn.preprocessing import StandardScaler
                        scaler = StandardScaler()
                        X_scaled = scaler.fit_transform(X)
                        
                        # Fit logistic regression to learn optimal metabolite weights
                        lr_model = LogisticRegression(max_iter=1000, solver='lbfgs', random_state=42)
                        lr_model.fit(X_scaled, y_true)
                        
                        # Get predicted probabilities for combined ROC
                        predicted_probs = lr_model.predict_proba(X_scaled)[:, 1]
                        
                        # Calculate combined ROC curve
                        fpr_combined, tpr_combined, _ = roc_curve(y_true, predicted_probs)
                        auc_combined = auc(fpr_combined, tpr_combined)
                        
                        # Plot combined as bold black line (plotted last so it appears last in legend)
                        ax.plot(fpr_combined, tpr_combined, label=f"Combined (AUC={auc_combined:.2f})", 
                            color='black', linewidth=3.5, alpha=0.95, linestyle='-')
                        auc_entries.append({'metabolite': 'Combined', 'AUC': auc_combined})
                        kept_curves += 1
                    except Exception as e:
                        logger.warning(f"Failed to calculate combined ROC: {e}")
                
                # Add random line (no label so it doesn't appear in legend)
                ax.plot([0,1], [0,1], color='lightgray', linestyle=':', linewidth=1.6)
                
                ax.set_xlabel('1 - Specificity', fontweight='bold', fontsize=params.xlabel_fontsize)
                ax.set_ylabel('Sensitivity', fontweight='bold', fontsize=params.ylabel_fontsize)
                ax.set_title(f'ROC Curves: {g2} vs {g1}', fontweight='bold', fontsize=params.title_fontsize)
                # Explicitly ensure axis label font sizes are applied
                ax.xaxis.label.set_fontsize(params.xlabel_fontsize)
                ax.yaxis.label.set_fontsize(params.ylabel_fontsize)
                ax.set_xlim(0,1)
                ax.set_ylim(0,1)
                ax.tick_params(axis='both', labelsize=params.tick_fontsize, width=1.5)
                for spine in ['top','right']:
                    if spine in ax.spines:
                        ax.spines[spine].set_visible(False)
                for spine in ['left','bottom']:
                    if spine in ax.spines:
                        ax.spines[spine].set_linewidth(2.2)
                for lbl in ax.get_xticklabels()+ax.get_yticklabels():
                    lbl.set_fontweight('bold')
                
                # Move legend outside if more than 11 curves
                if kept_curves > 11:
                    leg = ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=params.legend_fontsize)
                else:
                    leg = ax.legend(frameon=False, loc='lower right', fontsize=params.legend_fontsize)
                for t in leg.get_texts():
                    t.set_fontweight('bold')
                plt.tight_layout()
                tag = f"{g2}_vs_{g1}".replace(' ', '_')
                out_png = os.path.join(roc_dir, f'roc_{tag}.png')
                plt.savefig(out_png, dpi=params.fig_dpi)
                files_created.append(out_png)
                
                # Save AUC values (controlled by save_excel parameter)
                if auc_entries and params.save_excel:
                    auc_df = pd.DataFrame(auc_entries)
                    auc_path = os.path.join(roc_dir, f'roc_{tag}_AUC.csv')
                    auc_df.to_csv(auc_path, index=False)
                    files_created.append(auc_path)
                plt.close()
            # Write summary file
            try:
                summary_lines.append(f"Metabolites per pair:")
                for rec in per_pair_counts:
                    summary_lines.append(f"  {rec['pair']}: {rec['count']}")
                with open(os.path.join(roc_dir, 'roc_filter_summary.txt'), 'w', encoding='utf-8') as fh:
                    fh.write('\n'.join(summary_lines))
            except Exception:
                pass
        
        if files_created and not errors:
            summary = f"ROC analysis complete: {len(files_created)} files generated"
        elif files_created and errors:
            summary = f"ROC analysis partial: {len(files_created)} files, {len(errors)} errors"
        else:
            summary = f"ROC analysis failed: {len(errors)} errors"
        
    except Exception as e:
        errors.append(f"ROC analysis failed: {str(e)}")
        summary = "ROC analysis failed"
    
    return VizResults(
        files_created=files_created,
        errors=errors,
        summary=summary
    )

def validate_visualization_ready(complete_df: pd.DataFrame, groups: List[str], sample_cols: List[str]) -> Tuple[bool, str]:
    """Validate that data is ready for visualization."""
    if len(groups) < 2:
        return False, "At least 2 groups required for visualization"
    
    if len(sample_cols) < 3:
        return False, "Insufficient sample columns for analysis"
    
    if complete_df.empty:
        return False, "No data available for visualization"
    
    return True, "Data ready for visualization"
