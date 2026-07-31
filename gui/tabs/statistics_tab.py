import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import _tkinter
import logging
import os
import time
import pandas as pd
import json
import threading
import traceback
import numpy as np
from collections import Counter


# Import shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.shared.utils import resolve_runtime_config_path, is_statistics_metadata_col
from gui.shared.column_assignment import show_column_assignment_dialog

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class StatisticsTab(BaseTab):
    """Statistics Tab - Placeholder for statistical analysis
    
    This tab will handle statistical analysis of metabolite data.
    """
    
    def __init__(self, parent, data_manager):
        """Initialize Statistics tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles (runs only once)
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Setup memory_store as reference to data_manager's memory store
        self.memory_store = self.data_manager.memory_store
        
        # Create UI
        self.setup_ui()
        print("[OK] Statistics Tab initialized")
    
    def setup_ui(self):
        """Create the tab interface"""
        # Call the full statistics setup
        self.setup_statistics_tab()
        
        logger.info("Statistics tab UI created")
    
    @staticmethod
    def _remove_readonly_attribute(filepath):
        """Remove read-only attribute from a file (Windows-specific fix)"""
        try:
            import os
            import stat
            os.chmod(filepath, stat.S_IWRITE)
        except Exception as e:
            logger.debug(f"Could not remove read-only attribute from {filepath}: {e}")
    
    def ordered_groups(self, groups):
        """Order groups according to group_definitions or group_order setting"""
        if not groups:
            return []
        
        # Try to use custom group order if specified
        try:
            if hasattr(self, 'statistics_group_order_var'):
                order_str = self.statistics_group_order_var.get()
                if order_str and order_str.strip():
                    ordered = order_str.split(',')
                    ordered = [g.strip() for g in ordered if g.strip()]
                    # Return ordered list, filtering for groups that exist in data
                    return [g for g in ordered if g in groups]
        except Exception:
            pass
        
        # Fallback: use group_definitions order if available
        if hasattr(self, 'group_definitions'):
            defined_order = list(self.group_definitions.keys())
            return [g for g in defined_order if g in groups]
        
        # Fallback: just return as-is
        return list(groups)
    
    def on_data_ready(self, data_key: str = None):
        """Called when ID Annotation tab has finished and data is ready"""
        print(f"✅ Statistics tab notified: Data ready (key: {data_key})")
        # Trigger auto-load of metabolite data
        try:
            self.root.after(100, self._auto_load_metabolite_data)
        except Exception as e:
            print(f"Warning: Could not auto-load metabolite data: {e}")
    
    def auto_load_annotated_file(self, file_path: str):
        """Auto-load annotated file from ID Annotation tab"""
        try:
            if file_path and os.path.exists(file_path):
                # Store in memory for _auto_load_metabolite_data to find
                self.annotated_metabolites_excel_path = file_path
                self.id_annotated_excel_path = file_path
                self.annotated_ids_excel_path = file_path
                
                print(f"✅ Statistics tab auto-loading annotated file: {file_path}")
                
                # Actually load the data
                self._auto_load_metabolite_data()
            else:
                print(f"⚠️ Annotated file not found: {file_path}")
        except Exception as e:
            print(f"Error auto-loading annotated file: {e}")
    
    def _auto_load_to_visualization(self):
        """Auto-load statistical results to Visualization tab without switching.

        Per user request: autoload immediately when stats finish, but only auto-switch
        to the Visualization tab after the user exports statistical results.
        """
        try:
            # Notify Visualization tab that results are ready
            self.notify_data_ready("📊 Visualization", "statistical_results")
            
            # Get Visualization tab and trigger its data loading
            viz_tab = self.get_tab_by_name("📊 Visualization")
            if viz_tab:
                # Call visualization tab's data update method if available
                if hasattr(viz_tab, 'update_viz_data_status'):
                    viz_tab.update_viz_data_status()
                    print(f"✅ Updated Visualization tab data status")
                # Do not auto-switch here; switch occurs upon export
                self._thread_safe_log("✅ Visualization data autoloaded (no auto-switch).\n")
            # else:
            #     print(f"⚠️ Visualization tab not found")
        except Exception as e:
            print(f"Warning: Could not auto-load to visualization: {e}")



    def setup_statistics_tab(self):
        """Create the enhanced Statistics tab with scrollbar, action buttons, and right-side log"""
        import importlib

        # Layout defaults (can be overridden before setup or at runtime
        # via set_statistics_layout)
        # Use a more reasonable default for the statistics middle height to
        # ensure configuration content is visible without excessive scrolling
        self.stats_middle_height = getattr(self, 'stats_middle_height', 800)
        self.stats_middle_ratio = getattr(self, 'stats_middle_ratio', (1, 4))
        self.stats_column_mins = getattr(self, 'stats_column_mins', (300, 500, 280))
        # Number of text lines for the right-hand stats log. Set to None to
        # allow the ScrolledText to size by grid instead of a fixed line height.
        self.stats_log_lines = getattr(self, 'stats_log_lines', 35)
        # Backend toggle: keep lipid class sheets hidden by default.
        # Set this to True in code if you want to load workbook class sheets
        # instead of deriving class outputs after Step 4.
        self.use_lipid_class_sheet_backend = getattr(self, 'use_lipid_class_sheet_backend', False)

        # Define config persistence methods early so widgets can use them
        if not hasattr(self, '_stats_config_file'):
            def _stats_config_file():
                try:
                    base_dir = os.path.dirname(os.path.abspath(__file__))
                except Exception:
                    base_dir = os.getcwd()
                return resolve_runtime_config_path('statistics_config.json', {base_dir})
            self._stats_config_file = _stats_config_file
        
        if not hasattr(self, '_log_stats'):
            def _log_stats(msg: str):
                if hasattr(self, 'stats_log'):
                    ts = time.strftime('%H:%M:%S')
                    self.stats_log.insert(tk.END, f"[{ts}] {msg}\n")
                    self.stats_log.see(tk.END)
            self._log_stats = _log_stats
        
        if not hasattr(self, '_remove_readonly_attribute'):
            def _remove_readonly_attribute(file_path: str):
                try:
                    import stat
                    if os.path.exists(file_path):
                        os.chmod(file_path, stat.S_IWRITE | stat.S_IREAD)
                except Exception:
                    pass
            self._remove_readonly_attribute = _remove_readonly_attribute

        # Helper: apply custom pairwise p-value adjustment methods beyond BH
        if not hasattr(self, '_apply_custom_pairwise_adjustment'):
            def _apply_custom_pairwise_adjustment(df, method: str, scope: str):
                """Apply alternative multiple-testing corrections to pairwise results.

                Parameters:
                    df (pd.DataFrame): Pairwise results with at least 'p_value' or 'p_value_adj'.
                    method (str): One of {'Bonferroni','Holm','Hochberg','BY','BH','None'}.
                    scope (str): FDR scope ('per-comparison','per-feature','global').
                Notes:
                    - We re-compute adjusted p-values from raw 'p_value' (BH removed) to keep semantics.
                    - For Tukey/Dunn outputs having 'p_value_adj' only, we derive raw vector from available columns if possible.
                """
                import numpy as _np
                import pandas as _pd
                if method in ('BH','None'):
                    return  # Already handled or user requested none
                if df is None or df.empty:
                    return
                raw_col = 'p_value'
                if raw_col not in df.columns:
                    # Fallback: attempt to reconstruct from existing adjusted p-values (not strictly correct)
                    raw_col = 'p_value_adj' if 'p_value_adj' in df.columns else None
                if raw_col is None:
                    return
                # Decide grouping based on scope
                def _adjust_vector(pvals: list[float]):
                    m = len(pvals)
                    if m == 0:
                        return []
                    try:
                        from statsmodels.stats.multitest import multipletests
                        method_map = {
                            'Bonferroni': 'bonferroni',
                            'Holm': 'holm',
                            'Hochberg': 'simes-hochberg',  # statsmodels uses 'simes-hochberg'
                            'BY': 'fdr_by'
                        }
                        if method in method_map:
                            adj = multipletests(pvals, method=method_map[method])[1]
                            return adj.tolist()
                    except Exception:
                        pass  # Fallback manual implementations below
                    p = _np.array(pvals, dtype=float)
                    if method == 'Bonferroni':
                        return (_np.minimum(p * m, 1.0)).tolist()
                    if method == 'Holm':
                        order = _np.argsort(p)
                        adj = _np.empty_like(p)
                        for i, idx in enumerate(order):
                            adj[idx] = min((m - i) * p[idx], 1.0)
                        # enforce monotonicity
                        rev_order = _np.argsort(order)
                        for i in range(1, m):
                            prev = order[i-1]; cur = order[i]
                            adj[cur] = max(adj[cur], adj[prev])
                        return adj.tolist()
                    if method == 'Hochberg':
                        order = _np.argsort(p)[::-1]  # descending
                        adj = _np.empty_like(p)
                        cumulative = 1.0
                        for i, idx in enumerate(order):
                            rank = i + 1
                            val = p[idx] * rank
                            cumulative = min(cumulative, val)
                            adj[idx] = min(val, 1.0)
                        # monotonic from smallest p to largest
                        for i in range(m-2, -1, -1):
                            cur = order[i]; nxt = order[i+1]
                            adj[cur] = min(adj[cur], adj[nxt])
                        return adj.tolist()
                    if method == 'BY':
                        # Benjamini-Yekutieli: BH * harmonic factor
                        harmonic = sum(1.0 / k for k in range(1, m+1))
                        order = _np.argsort(p)
                        adj = _np.empty_like(p)
                        for i, idx in enumerate(order):
                            rank = i + 1
                            adj[idx] = min(p[idx] * m * harmonic / rank, 1.0)
                        # enforce monotonicity (from largest rank to smallest)
                        for i in range(m-2, -1, -1):
                            cur = order[i]; nxt = order[i+1]
                            adj[cur] = min(adj[cur], adj[nxt])
                        return adj.tolist()
                    return p.tolist()
                # Grouping logic
                if scope == 'per-comparison' and {'group1','group2'} <= set(df.columns):
                    for (g1, g2), sub_idx in df.groupby(['group1','group2']).groups.items():
                        sub = df.loc[sub_idx]
                        pvals = [float(x) if _pd.notna(x) else 1.0 for x in sub[raw_col].tolist()]
                        adj = _adjust_vector(pvals)
                        df.loc[sub_idx, 'p_value_adj'] = adj
                elif scope == 'per-metabolite' and 'metabolite' in df.columns:
                    for metab, sub_idx in df.groupby('metabolite').groups.items():
                        sub = df.loc[sub_idx]
                        pvals = [float(x) if _pd.notna(x) else 1.0 for x in sub[raw_col].tolist()]
                        adj = _adjust_vector(pvals)
                        df.loc[sub_idx, 'p_value_adj'] = adj
                else:  # global
                    pvals = [float(x) if _pd.notna(x) else 1.0 for x in df[raw_col].tolist()]
                    adj = _adjust_vector(pvals)
                    df['p_value_adj'] = adj
                # Recompute neg_log10 and Expression if log2_fold_change present
                if 'p_value_adj' in df.columns:
                    df['neg_log10_p_adj'] = -_np.log10(df['p_value_adj'].replace(0, _np.finfo(float).eps))
                    if 'log2_fold_change' in df.columns and 'p_value_adj' in df.columns:
                        mask_sig = (df['p_value_adj'] <= 0.05) & df['log2_fold_change'].notna()
                        df.loc[mask_sig & (df['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
                        df.loc[mask_sig & (df['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
            self._apply_custom_pairwise_adjustment = _apply_custom_pairwise_adjustment
        
        if not hasattr(self, '_gather_statistics_config'):
            def _gather_statistics_config():
                return {
                    'group_definitions': self.group_definitions,
                    'group_count': self.group_count,
                    'base_group': self.stat_base_group.get() if hasattr(self, 'stat_base_group') else '',
                    'stat_norm_method': self.stat_norm_method.get() if hasattr(self, 'stat_norm_method') else 'none',
                    'stat_test_type': self.stat_test_type.get() if hasattr(self, 'stat_test_type') else 'overall',
                    'overall_test': self.stat_overall_test.get() if hasattr(self, 'stat_overall_test') else 'anova',
                    'pairwise_test': self.stat_pairwise_test.get() if hasattr(self, 'stat_pairwise_test') else 'welch',
                    'fdr_scope': self.fdr_scope_var.get() if hasattr(self, 'fdr_scope_var') else 'per-comparison',
                    'alpha': '0.05',  # Hardcoded; user cannot change

                    'pairwise_p_adjust_method': self.pairwise_p_adjust_method.get() if hasattr(self, 'pairwise_p_adjust_method') else 'BH',
                    'data_mode': self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite',
                    'custom_comparisons': self.custom_comparisons_var.get() if hasattr(self, 'custom_comparisons_var') else '',
                    'filter_timing': self.filter_timing_var.get() if hasattr(self, 'filter_timing_var') else 'before',
                    'enable_imputation_before': bool(self.enable_imputation_before_var.get()) if hasattr(self, 'enable_imputation_before_var') else False,
                    'imputation_before_method': self.imputation_before_method_var.get() if hasattr(self, 'imputation_before_method_var') else 'half_min',
                    'imputation_before_knn_neighbors': self.imputation_before_knn_neighbors_var.get() if hasattr(self, 'imputation_before_knn_neighbors_var') else '5',
                    'enable_variability_filter': bool(self.enable_variability_filter_var.get()) if hasattr(self, 'enable_variability_filter_var') else False,
                    'variability_percentile': self.variability_percent_var.get() if hasattr(self, 'variability_percent_var') else '10',
                    'enable_imputation': bool(self.enable_imputation_var.get()) if hasattr(self, 'enable_imputation_var') else False,
                    'imputation_method': self.imputation_method_var.get() if hasattr(self, 'imputation_method_var') else 'half_min',
                    'knn_neighbors': self.knn_neighbors_var.get() if hasattr(self, 'knn_neighbors_var') else '5',
                    'imputation_min_group_percent': self.imputation_min_group_percent_var.get() if hasattr(self, 'imputation_min_group_percent_var') else '50.0',
                    'imputation_prefilter_scope': self.imputation_prefilter_scope_var.get() if hasattr(self, 'imputation_prefilter_scope_var') else 'per_group',
                    'enable_pca_outlier': bool(self.enable_pca_outlier_var.get()) if hasattr(self, 'enable_pca_outlier_var') else False,
                    'sample_group_assignments': {c: v.get() for c, v in getattr(self, 'sample_group_vars', {}).items()},
                    'n_jobs': self._get_workers_count(getattr(self, 'stats_workers', None), default=3)
                }
            self._gather_statistics_config = _gather_statistics_config
        
        if not hasattr(self, '_save_statistics_config'):
            def _save_statistics_config():
                # Don't save during initialization before config is loaded
                if not hasattr(self, '_stats_config_loaded') or not self._stats_config_loaded:
                    return
                    
                path = self._stats_config_file()
                try:
                    data = self._gather_statistics_config()
                    # Include auto-assign patterns in the config (ALWAYS save, even if empty)
                    if hasattr(self, '_auto_assign_patterns'):
                        data['auto_assign_patterns'] = self._auto_assign_patterns
                    else:
                        data['auto_assign_patterns'] = {}  # Save empty dict if not set
                    with open(path, 'w', encoding='utf-8') as f:
                        json.dump(data, f, indent=2)
                    # Remove read-only attribute
                    self._remove_readonly_attribute(path)
                    print(f"💾 Saved statistics config with auto-assign patterns: {list(data.get('auto_assign_patterns', {}).keys())}")
                except Exception as e:
                    self._log_stats(f"Failed to save stats config: {e}")
            self._save_statistics_config = _save_statistics_config
        
        if not hasattr(self, '_apply_statistics_config'):
            def _apply_statistics_config(data: dict):
                gd = data.get('group_definitions') or {}
                if gd:
                    self.group_definitions = gd
                    self.group_count = data.get('group_count', len(gd))
                    self.group_id_vars = {gid: tk.StringVar(value=label) for gid, label in self.group_definitions.items()}
                    if hasattr(self, 'refresh_group_ui'):
                        self.refresh_group_ui()
                mapping = [
                    ('stat_norm_method','stat_norm_method','none'),
                    ('stat_test_type','stat_test_type','overall'),
                    ('stat_overall_test','overall_test','anova'),
                    ('stat_pairwise_test','pairwise_test','welch'),
                    ('fdr_scope_var','fdr_scope','per-comparison'),
                    # alpha removed (hardcoded 0.05)
                    ('pairwise_p_adjust_method','pairwise_p_adjust_method','BH'),
                    ('statistics_data_mode','data_mode','metabolite'),
                    ('custom_comparisons_var','custom_comparisons',''),
                    ('filter_timing_var','filter_timing','before')
                ]
                for attr, key, default in mapping:
                    if hasattr(self, attr) and key in data:
                        try:
                            getattr(self, attr).set(data.get(key, default))
                        except Exception:
                            pass
                # Workers (n_jobs)
                try:
                    if 'n_jobs' in data:
                        if not hasattr(self, 'stats_workers'):
                            self.stats_workers = tk.StringVar(value=str(int(data.get('n_jobs', 3) or 3)))
                        else:
                            self.stats_workers.set(str(int(data.get('n_jobs', 3) or 3)))
                except Exception:
                    pass
                # Boolean settings
                if 'use_adj_p' in data and hasattr(self, 'use_adj_p_var'):
                    try:
                        pass  # use_adj_p removed - not used
                    except Exception:
                        pass
                if 'base_group' in data and hasattr(self, 'stat_base_group') and data['base_group'] in self.group_definitions:
                    self.stat_base_group.set(data['base_group'])
                if hasattr(self, 'enable_variability_filter_var') and 'enable_variability_filter' in data:
                    self.enable_variability_filter_var.set(bool(data.get('enable_variability_filter', False)))
                if hasattr(self, 'variability_percent_var') and 'variability_percentile' in data:
                    self.variability_percent_var.set(str(data.get('variability_percentile', '10')))
                if hasattr(self, 'enable_imputation_before_var') and 'enable_imputation_before' in data:
                    self.enable_imputation_before_var.set(bool(data.get('enable_imputation_before', False)))
                if hasattr(self, 'imputation_before_method_var') and 'imputation_before_method' in data:
                    self.imputation_before_method_var.set(str(data.get('imputation_before_method', 'half_min')))
                if hasattr(self, 'imputation_before_knn_neighbors_var') and 'imputation_before_knn_neighbors' in data:
                    self.imputation_before_knn_neighbors_var.set(str(data.get('imputation_before_knn_neighbors', '5')))
                if hasattr(self, 'enable_imputation_var') and 'enable_imputation' in data:
                    self.enable_imputation_var.set(bool(data.get('enable_imputation', False)))
                if hasattr(self, 'imputation_method_var') and 'imputation_method' in data:
                    self.imputation_method_var.set(str(data.get('imputation_method', 'half_min')))
                if hasattr(self, 'knn_neighbors_var') and 'knn_neighbors' in data:
                    self.knn_neighbors_var.set(str(data.get('knn_neighbors', '5')))
                if hasattr(self, 'imputation_min_group_percent_var') and 'imputation_min_group_percent' in data:
                    self.imputation_min_group_percent_var.set(str(data.get('imputation_min_group_percent', '50.0')))
                if hasattr(self, 'imputation_prefilter_scope_var') and 'imputation_prefilter_scope' in data:
                    scope = str(data.get('imputation_prefilter_scope', 'per_group')).strip().lower()
                    if scope not in ('per_group', 'all_groups'):
                        scope = 'per_group'
                    self.imputation_prefilter_scope_var.set(scope)
                if hasattr(self, 'enable_pca_outlier_var') and 'enable_pca_outlier' in data:
                    self.enable_pca_outlier_var.set(bool(data.get('enable_pca_outlier', False)))
                sga = data.get('sample_group_assignments') or {}
                if sga and hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                    # Only apply saved assignments to columns that exist in the CURRENT data
                    # This prevents old group assignments from persisting when new data is loaded
                    # with different sample column names
                    current_cols = set(self.sample_group_vars.keys())
                    for col, grp in sga.items():
                        if col in current_cols and col in self.sample_group_vars:
                            # Handle both internal IDs (Group1) and display labels (Control)
                            if grp in self.group_definitions:
                                # It's an internal ID, convert to display label
                                self.sample_group_vars[col].set(self.group_definitions[grp])
                            elif grp in self.group_definitions.values():
                                # It's already a display label, use as-is
                                self.sample_group_vars[col].set(grp)
                            # If neither, skip (invalid saved value)
                # Two-way ANOVA factor assignments no longer persisted (console-based)
                # Load auto-assign patterns (even if empty, overwrite existing)
                patterns = data.get('auto_assign_patterns', {})
                self._auto_assign_patterns = patterns  # Always set, even if empty dict
                if patterns:
                    print(f"✅ Loaded auto-assign patterns: {list(patterns.keys())}")
                    self._log_stats(f'Loaded auto-assign patterns for {len(patterns)} groups.')
                else:
                    print(f"ℹ️ No auto-assign patterns in saved config")
                # Update filter timing explanation text
                if hasattr(self, '_update_filter_timing_explanation'):
                    self._update_filter_timing_explanation()
                self._log_stats('Previous statistics configuration loaded.')
            self._apply_statistics_config = _apply_statistics_config
        
        if not hasattr(self, '_load_statistics_config'):
            def _load_statistics_config():
                path = self._stats_config_file()
                if os.path.exists(path):
                    try:
                        with open(path, 'r', encoding='utf-8') as f:
                            data = json.load(f)
                        self._apply_statistics_config(data)
                    except Exception as e:
                        self._log_stats(f"Failed to load stats config: {e}")
                else:
                    self._log_stats('No previous statistics configuration found.')
            self._load_statistics_config = _load_statistics_config
        
        # Ensure stats config change handler calls save
        if not hasattr(self, '_stats_config_changed'):
            def _stats_config_changed(log: str | None = None):
                try:
                    self._save_statistics_config()
                    print(f"💾 Config saved with group_definitions: {self.group_definitions}")
                except Exception as e:
                    print(f"❌ Failed to save config: {e}")
                if log:
                    self._log_stats(log)
            self._stats_config_changed = _stats_config_changed

        # FDR scope change handler with warning logic
        def _on_fdr_scope_changed():
            scope = self.fdr_scope_var.get()
            self._stats_config_changed(log=f"FDR scope: {scope}")
            
            # Update warning visibility based on selection
            if hasattr(self, 'fdr_scope_warning'):
                if scope == 'per-metabolite':
                    self.fdr_scope_warning.pack(fill='x', pady=(2, 0))
                else:
                    self.fdr_scope_warning.pack_forget()
        
        self._on_fdr_scope_changed = _on_fdr_scope_changed

        # Use the frame from BaseTab as the main container
        # (this tab is already inside the main GUI notebook)
        main_container = tk.Frame(self.frame, bg='#f0f0f0')
        main_container.pack(fill='both', expand=True, padx=5, pady=5)

        # Let the canvas expand to the container size instead of fixed dimensions
        self.stats_canvas = tk.Canvas(main_container, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(main_container, orient='vertical', command=self.stats_canvas.yview)
        scrollable_frame = tk.Frame(self.stats_canvas, bg='#f0f0f0')
        scrollable_frame.bind("<Configure>", lambda e: self.stats_canvas.configure(scrollregion=self.stats_canvas.bbox("all")))
        canvas_frame = self.stats_canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        self.stats_canvas.configure(yscrollcommand=scrollbar.set)

        def configure_canvas(event):
            # Update both width and height of the canvas window so the
            # scrollable_frame fills the available notebook area when resized
            try:
                # Only set the width of the inner window. Do NOT force its
                # height to the canvas height — that prevents vertical
                # scrolling because the inner window would always match the
                # visible canvas area.
                self.stats_canvas.itemconfig(canvas_frame, width=event.width)
            except Exception:
                try:
                    self.stats_canvas.itemconfig(canvas_frame, width=event.width)
                except Exception:
                    pass
            # Update scrollregion to match the content
            try:
                self.stats_canvas.configure(scrollregion=self.stats_canvas.bbox('all'))
            except Exception:
                pass
        self.stats_canvas.bind('<Configure>', configure_canvas)

        # Avoid forcing a fixed pixel height on the canvas. Let the
        # notebook / container manage sizing so vertical scrolling works
        # normally and the tab does not create excessive blank space.
        # If callers want to enforce a minimum, they can call
        # set_statistics_layout() which will adjust grid minsize instead.

        self.stats_canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        def _on_mousewheel(event):
            # Normalize wheel delta for Windows (event.delta is multiple of 120)
            try:
                step = int(-1 * (event.delta / 120))
            except Exception:
                # Fallback if event.delta not available
                step = -1 if getattr(event, 'num', None) == 5 else 1
            try:
                self.stats_canvas.yview_scroll(step, 'units')
            except Exception:
                pass

        # Bind mousewheel locally to the canvas and its inner frame so that
        # scrolling only affects the pane under the cursor (avoids global
        # conflicts with other canvases in the app).
        def _bind_mousewheel_recursive(widget):
            """Recursively bind mousewheel to all child widgets"""
            widget.bind('<MouseWheel>', _on_mousewheel)
            for child in widget.winfo_children():
                _bind_mousewheel_recursive(child)
        
        self.stats_canvas.bind('<MouseWheel>', _on_mousewheel)
        # Schedule recursive binding after all widgets are created
        scrollable_frame.after(100, lambda: _bind_mousewheel_recursive(scrollable_frame))

        container = tk.Frame(scrollable_frame, bg='#f0f0f0')
        container.pack(fill='both', expand=True, padx=10, pady=10)
        #tk.Label(container, text='📊 Statistics & Normalization', font=('Arial', 16, 'bold'), bg='#f0f0f0').pack(pady=(0, 10))

        body = tk.Frame(container, bg='#f0f0f0')
        body.pack(fill='both', expand=True)
        # keep a reference for runtime layout updates
        self.stats_body = body

        # NEW 3-COLUMN LAYOUT: Prioritize left and middle columns
        # Left: Steps 1-3 (Mode, Import, Verify, Groups) - 40%
        # Middle: Steps 4-5 (Normalization, Statistics) - 40%
        # Right: Progress & Results - 20%
        body.grid_columnconfigure(0, weight=2)  # Left column - Steps 1-3 (higher priority)
        body.grid_columnconfigure(1, weight=2)  # Middle column - Steps 4-5 (higher priority)
        body.grid_columnconfigure(2, weight=1)  # Right column - Progress (lower priority)
        body.grid_rowconfigure(0, weight=1, minsize=800)

        # ========================================================================
        # LEFT COLUMN: Steps 1-3 (Mode Selection, Import, Verify, Group Config)
        # ========================================================================
        left_col = tk.LabelFrame(body, text='📋 Steps 1-3: Data & Groups', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        left_col.grid(row=0, column=0, sticky='nsew', padx=(0, 3))
        
        # Canvas with scrollbars for left column
        left_canvas = tk.Canvas(left_col, bg='#f0f0f0', highlightthickness=0, height=700)
        left_scrollbar_y = ttk.Scrollbar(left_col, orient="vertical", command=left_canvas.yview)
        left_scrollbar_x = ttk.Scrollbar(left_col, orient="horizontal", command=left_canvas.xview)
        left_scrollable = tk.Frame(left_canvas, bg='#f0f0f0')
        
        left_scrollable.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )
        
        left_canvas_window = left_canvas.create_window((0, 0), window=left_scrollable, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar_y.set, xscrollcommand=left_scrollbar_x.set)
        
        def configure_left_scroll(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
            left_canvas.itemconfig(left_canvas_window, width=event.width)
        
        left_canvas.bind('<Configure>', configure_left_scroll)
        
        left_scrollbar_y.pack(side="right", fill="y", padx=(2, 0))
        left_scrollbar_x.pack(side="bottom", fill="x", pady=(2, 0))
        left_canvas.pack(side="left", fill="both", expand=True)
        
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_scrollable.bind("<MouseWheel>", _on_left_mousewheel)
        
        # ========== STEP 1: Data Mode Selection ==========
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        
        step1_frame = tk.LabelFrame(left_scrollable, text='Step 1: Select Data Mode', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        step1_frame.pack(fill='x', padx=5, pady=(5, 10))
        
        # Data Mode Selection (Metabolite vs Lipid vs Custom)
        mode_frame = tk.Frame(step1_frame, bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Label(mode_frame, text='Select mode:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(side='left', padx=(0, 10))
        self.statistics_data_mode = tk.StringVar(value='metabolite')
        self.statistics_data_mode.trace_add('write', lambda *a: self._stats_config_changed(log=f"Data mode: {self.statistics_data_mode.get()}"))
        tk.Radiobutton(mode_frame, text='Metabolite', variable=self.statistics_data_mode,
                       value='metabolite', bg='#f0f0f0', command=self.on_statistics_mode_change).pack(side='left', padx=5)
        tk.Radiobutton(mode_frame, text='Lipid', variable=self.statistics_data_mode,
                       value='lipid', bg='#f0f0f0', command=self.on_statistics_mode_change).pack(side='left', padx=5)
        tk.Radiobutton(mode_frame, text='Custom', variable=self.statistics_data_mode,
                       value='custom', bg='#f0f0f0', command=self.on_statistics_mode_change).pack(side='left', padx=5)
        
        # Mode description label
        self.stats_mode_desc_label = tk.Label(
            step1_frame,
            text='Metabolite mode expects Pos_id/Neg_id sheets',
            bg='#e3f2fd',
            fg='#1565c0',
            font=('Arial', 8, 'italic'),
            wraplength=280,
            justify='left',
            padx=5,
            pady=5
        )
        self.stats_mode_desc_label.pack(fill='x', padx=5, pady=(0, 5))
        
        # ========== STEP 2: Import & Verify Columns ==========
        step2_frame = tk.LabelFrame(left_scrollable, text='Step 2: Import & Verify Columns', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        step2_frame.pack(fill='x', padx=5, pady=(0, 10))
        
        tk.Button(step2_frame, text='📂 Import Excel Data (Pos/Neg)', command=self.import_statistics_excel,
                  bg='#8e44ad', fg='white', **btn_style).pack(fill='x', padx=5, pady=5)
        
        # Important note about sample name consistency
        note_frame = tk.Frame(step2_frame, bg='#fff3cd', relief='solid', bd=1)
        note_frame.pack(fill='x', padx=5, pady=(0, 5))
        tk.Label(note_frame, text='⚠️', bg='#fff3cd', font=('Arial', 12), fg='#856404').pack(side='left', padx=(5, 0))
        tk.Label(note_frame, text='Important: Ensure sample names in Positive sheet match exactly\nwith sample names in Negative sheet for accurate analysis.',
                bg='#fff3cd', font=('Arial', 8), fg='#856404', justify='left', wraplength=380).pack(side='left', padx=(5, 5), pady=5)
        
        tk.Button(step2_frame, text='🔍 Verify Columns', command=self.verify_statistics_columns,
                  bg='#2980b9', fg='white', **btn_style).pack(fill='x', padx=5, pady=5)
        
        # ========== STEP 3: Configure Groups ==========
        step3_frame = tk.LabelFrame(left_scrollable, text='Step 3: Configure Groups', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        step3_frame.pack(fill='x', padx=5, pady=(0, 10))
        
        # Configure Groups button (disabled until verification complete)
        self.configure_groups_btn = tk.Button(step3_frame, text='⚙️ Configure Groups', command=self.auto_assign_groups,
                  bg='#9b59b6', fg='white', state='disabled', **btn_style)
        self.configure_groups_btn.pack(fill='x', padx=5, pady=5)
        
        # Group IDs & Labels section (compact)
        group_ids_subsection = tk.LabelFrame(step3_frame, text='Group IDs & Labels', bg='#f0f0f0')
        group_ids_subsection.pack(fill='x', padx=5, pady=(5, 5))
        
        # Initialize group data
        if not hasattr(self, 'group_definitions') or not self.group_definitions:
            self.group_definitions = {'Group1': 'Control', 'Group2': 'Disease', 'Group3': 'Treatment', 'Group4': 'Other'}
            self.group_count = 4
        if hasattr(self, 'base_group_combo'):
            self.base_group_combo['values'] = [''] + list(self.group_definitions.keys())
        
        # Add/Remove group buttons at the top
        buttons_frame = tk.Frame(group_ids_subsection, bg='#f0f0f0')
        buttons_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Button(buttons_frame, text='+ Add Group', command=self.add_group,
                  bg='#27ae60', fg='white', font=('Arial', 8, 'bold'), pady=2, padx=8).pack(side='left', padx=2)
        tk.Button(buttons_frame, text='- Remove', command=self.remove_group,
                  bg='#e74c3c', fg='white', font=('Arial', 8, 'bold'), pady=2, padx=8).pack(side='left', padx=2)
        
        # Canvas for group IDs (compact scrollable list)
        groups_canvas = tk.Canvas(group_ids_subsection, bg='#f0f0f0', highlightthickness=0, height=120)
        groups_scrollbar = ttk.Scrollbar(group_ids_subsection, orient='vertical', command=groups_canvas.yview)
        self.groups_scrollable_frame = tk.Frame(groups_canvas, bg='#f0f0f0')
        
        self.groups_scrollable_frame.bind(
            "<Configure>",
            lambda e: groups_canvas.configure(scrollregion=groups_canvas.bbox("all"))
        )
        
        groups_canvas_window = groups_canvas.create_window((0, 0), window=self.groups_scrollable_frame, anchor='nw')
        groups_canvas.configure(yscrollcommand=groups_scrollbar.set)
        
        def configure_groups_scroll(event):
            groups_canvas.configure(scrollregion=groups_canvas.bbox('all'))
            groups_canvas.itemconfig(groups_canvas_window, width=event.width)
        
        groups_canvas.bind('<Configure>', configure_groups_scroll)
        groups_scrollbar.pack(side="right", fill="y")
        groups_canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        self.groups_canvas = groups_canvas
        
        def _on_groups_mousewheel(event):
            groups_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        groups_canvas.bind('<MouseWheel>', _on_groups_mousewheel)
        self.groups_scrollable_frame.bind('<MouseWheel>', _on_groups_mousewheel)
        
        # Refresh group UI
        self.group_id_vars = {}
        self.refresh_group_ui()
        
        # ========== Replicate Filtering Configuration ==========
        filtering_subsection = tk.LabelFrame(step3_frame, text='Replicate Filtering', bg='#f0f0f0')
        filtering_subsection.pack(fill='x', padx=5, pady=(5, 5))
        
        tk.Label(filtering_subsection, text='Apply filtering:', bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 2))
        
        self.filter_timing_var = tk.StringVar(value='before')
        self.filter_timing_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Filter timing: {self.filter_timing_var.get()}"))
        
        filter_radio_frame = tk.Frame(filtering_subsection, bg='#f0f0f0')
        filter_radio_frame.pack(fill='x', padx=10, pady=(0, 5))
        
        tk.Radiobutton(filter_radio_frame, text='Before Normalization', variable=self.filter_timing_var,
                      value='before', bg='#f0f0f0', command=self._on_filter_timing_change).pack(anchor='w')
        tk.Radiobutton(filter_radio_frame, text='After Normalization', variable=self.filter_timing_var,
                      value='after', bg='#f0f0f0', command=self._on_filter_timing_change).pack(anchor='w')
        
        # Explanation label
        self.filter_timing_explanation = tk.Label(filtering_subsection, text='', bg='#fff3cd', fg='#856404',
                                                 font=('Arial', 8), wraplength=280, justify='left', padx=5, pady=5)
        self.filter_timing_explanation.pack(fill='x', padx=5, pady=(0, 5))
        self._update_filter_timing_explanation()
        
        # Minimum samples per group
        min_samples_frame = tk.Frame(filtering_subsection, bg='#f0f0f0')
        min_samples_frame.pack(fill='x', padx=5, pady=(0, 5))
        self.min_samples_controls_container = min_samples_frame
        
        self.min_samples_label = tk.Label(min_samples_frame, text='Min samples per group:', bg='#f0f0f0', font=('Arial', 9))
        self.min_samples_label.pack(anchor='w')
        
        self.min_samples_type_var = tk.StringVar(value='absolute')
        self.min_samples_type_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Min samples type: {self.min_samples_type_var.get()}"))
        
        type_frame = tk.Frame(min_samples_frame, bg='#f0f0f0')
        type_frame.pack(fill='x')
        
        self.min_samples_absolute_rb = tk.Radiobutton(type_frame, text='Absolute', variable=self.min_samples_type_var,
                  value='absolute', bg='#f0f0f0')
        self.min_samples_absolute_rb.pack(side='left')
        self.min_samples_percentage_rb = tk.Radiobutton(type_frame, text='Percentage', variable=self.min_samples_type_var,
                  value='percentage', bg='#f0f0f0')
        self.min_samples_percentage_rb.pack(side='left')
        
        value_frame = tk.Frame(min_samples_frame, bg='#f0f0f0')
        value_frame.pack(fill='x', pady=(2, 0))
        
        self.min_samples_count_label = tk.Label(value_frame, text='Count:', bg='#f0f0f0', font=('Arial', 8))
        self.min_samples_count_label.pack(side='left')
        self.min_samples_per_group_var = tk.StringVar(value='2')
        self.min_samples_per_group_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Min samples count: {self.min_samples_per_group_var.get()}"))
        self.min_samples_count_entry = tk.Entry(value_frame, textvariable=self.min_samples_per_group_var, width=8)
        self.min_samples_count_entry.pack(side='left', padx=5)
        
        self.min_samples_percent_label = tk.Label(value_frame, text='%:', bg='#f0f0f0', font=('Arial', 8))
        self.min_samples_percent_label.pack(side='left', padx=(10, 0))
        self.min_samples_percent_var = tk.StringVar(value='50.0')
        self.min_samples_percent_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Min samples percent: {self.min_samples_percent_var.get()}"))
        self.min_samples_percent_entry = tk.Entry(value_frame, textvariable=self.min_samples_percent_var, width=8)
        self.min_samples_percent_entry.pack(side='left', padx=5)
        
        # Group assignment instructions (no UI, handled in Verify Columns)
        info_label = tk.Label(step3_frame, 
                             text="ℹ️ Use 'Configure Groups' to assign samples to groups for analysis.",
                             bg='#e3f2fd', fg='#1565c0', font=('Arial', 8), 
                             wraplength=280, justify='left', padx=5, pady=5)
        info_label.pack(fill='x', padx=5, pady=5)
        
        # ========== Base Group ==========
        base_frame = tk.LabelFrame(step3_frame, text='Base Group (Optional)', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        base_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Label(base_frame, text='Compare all groups vs this base only', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2,2))
        self.stat_base_group = tk.StringVar(value='')
        self.base_group_combo = ttk.Combobox(base_frame, values=[''], textvariable=self.stat_base_group, state='readonly')
        self.base_group_combo.pack(fill='x', padx=5, pady=(2, 5))
        self.base_group_combo.bind('<<ComboboxSelected>>', lambda e: self._stats_config_changed(log=f"Base group: {self.stat_base_group.get() or '[None]'}"))
        
        # ========== Custom Comparisons (moved from Step 5) ==========
        custom_comp_frame = tk.LabelFrame(step3_frame, text='Custom Comparisons (Optional)', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        custom_comp_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Label(custom_comp_frame, text='E.g., "Group1-Group2,Group3-Group4"', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2,2))
        self.custom_comparisons_var = tk.StringVar(value='')
        self.custom_comparisons_var.trace_add('write', lambda *a: self._stats_config_changed())
        tk.Entry(custom_comp_frame, textvariable=self.custom_comparisons_var, font=('Arial', 9)).pack(fill='x', padx=5, pady=(2, 2))
        tk.Label(custom_comp_frame, text='Leave empty for all pairwise comparisons', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(0, 5))
        
        # ========== Group Order (moved from Step 5) ==========
        group_order_frame = tk.LabelFrame(step3_frame, text='Group Order (Optional)', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        group_order_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Label(group_order_frame, text='Comma separated labels e.g. Control,Disease,Treatment', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2,2))
        self.statistics_group_order_var = tk.StringVar(value='')
        self.statistics_group_order_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Group order updated"))
        tk.Entry(group_order_frame, textvariable=self.statistics_group_order_var, font=('Arial', 9)).pack(fill='x', padx=5, pady=(2, 2))
        tk.Label(group_order_frame, text='If empty, order follows group definitions listing.', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(0, 5))
        
        # ========================================================================
        # MIDDLE COLUMN: Steps 4-5 (Normalization & Statistical Tests)
        # ========================================================================
        middle_col = tk.LabelFrame(body, text='🔬 Steps 4-5: Analysis', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        middle_col.grid(row=0, column=1, sticky='nsew', padx=3)
        
        # Canvas with scrollbars for middle column
        middle_canvas = tk.Canvas(middle_col, bg='#f0f0f0', highlightthickness=0, height=700)
        middle_scrollbar_y = ttk.Scrollbar(middle_col, orient="vertical", command=middle_canvas.yview)
        middle_scrollbar_x = ttk.Scrollbar(middle_col, orient="horizontal", command=middle_canvas.xview)
        middle_scrollable = tk.Frame(middle_canvas, bg='#f0f0f0')
        
        middle_scrollable.bind(
            "<Configure>",
            lambda e: middle_canvas.configure(scrollregion=middle_canvas.bbox("all"))
        )
        
        middle_canvas_window = middle_canvas.create_window((0, 0), window=middle_scrollable, anchor="nw")
        middle_canvas.configure(yscrollcommand=middle_scrollbar_y.set, xscrollcommand=middle_scrollbar_x.set)
        
        def configure_middle_scroll(event):
            middle_canvas.configure(scrollregion=middle_canvas.bbox("all"))
            middle_canvas.itemconfig(middle_canvas_window, width=event.width)
        
        middle_canvas.bind('<Configure>', configure_middle_scroll)
        middle_scrollbar_y.pack(side="right", fill="y", padx=(2, 0))
        middle_scrollbar_x.pack(side="bottom", fill="x", pady=(2, 0))
        middle_canvas.pack(side="left", fill="both", expand=True)
        
        def _on_middle_mousewheel(event):
            middle_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        middle_canvas.bind("<MouseWheel>", _on_middle_mousewheel)
        middle_scrollable.bind("<MouseWheel>", _on_middle_mousewheel)
        
        # ========== STEP 4: Normalization ==========
        step4_frame = tk.LabelFrame(middle_scrollable, text='Step 4: Normalization', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        step4_frame.pack(fill='x', padx=5, pady=(0, 10))
        
        tk.Label(step4_frame, text='Normalization Methods:', bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 2))
        
        # Normalization methods configuration
        self.norm_method_vars = {}
        self.norm_method_order = {}
        self.norm_selection_counter = [0]
        self.stat_norm_method = tk.StringVar(value='none')
        
        self.norm_methods_list = [
            ('Median', 'median'),
            ('TIC', 'TIC'),
            ('PQN', 'PQN'),
            ('IS', 'IS'),
            ('LOESS_QC', 'LOESS_QC'),
            ('Rel_Abundance(%)', 'Rel_Abundance(%)'),
            ('Quantile', 'quantile'),
            ('VSN', 'VSN'),
            ('CLR', 'CLR'),
            ('Z-Score', 'zscore'),
            ('Log2', 'log2'),
        ]
        
        # Initialize vars
        for display_name, method_name in self.norm_methods_list:
            self.norm_method_vars[method_name] = tk.BooleanVar(value=False)
            self.norm_method_order[method_name] = None
        
        # Dropdown button with current selection display
        dropdown_frame = tk.Frame(step4_frame, bg='#f0f0f0')
        dropdown_frame.pack(fill='x', padx=5, pady=(0, 5))
        
        self.norm_dropdown_btn = tk.Button(dropdown_frame, text='Select Methods ▼', 
                                          command=self._show_norm_dropdown,
                                          bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                                          relief='raised', bd=2, pady=5, anchor='w')
        self.norm_dropdown_btn.pack(fill='x')
        
        # Display current chain below button
        chain_display_frame = tk.Frame(step4_frame, bg='#e8f5e9', relief='solid', bd=1)
        chain_display_frame.pack(fill='x', padx=5, pady=(0, 5))
        tk.Label(chain_display_frame, text='Current Chain:', bg='#e8f5e9', font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.norm_chain_display = tk.Label(chain_display_frame, text='none', bg='#e8f5e9', font=('Arial', 9, 'italic'), fg='#27ae60')
        self.norm_chain_display.pack(side='left', padx=5)
        
        tk.Label(
            step4_frame,
            text='💡 Methods applied in order selected',
            bg='#f0f0f0',
            font=('Arial', 8, 'italic'),
            fg='#555'
        ).pack(anchor='w', padx=5, pady=(0, 5))

        # Step 4a: Imputation before normalization (optional)
        imputation_before_frame = tk.LabelFrame(
            step4_frame,
            text='Step 4a: Optional Pre-Normalization Imputation',
            bg='#f0f0f0'
        )
        imputation_before_frame.pack(fill='x', padx=5, pady=(6, 4))

        self.enable_imputation_before_var = tk.BooleanVar(value=False)
        self.imputation_before_method_var = tk.StringVar(value='half_min')
        self.imputation_before_knn_neighbors_var = tk.StringVar(value='5')

        self.enable_imputation_before_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Pre-normalization imputation: {'on' if self.enable_imputation_before_var.get() else 'off'}"))
        self.imputation_before_method_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.imputation_before_knn_neighbors_var.trace_add('write', lambda *a: self._stats_config_changed())

        imp_before_frame = tk.Frame(imputation_before_frame, bg='#f0f0f0')
        imp_before_frame.pack(fill='x', padx=5, pady=(4, 2))
        tk.Checkbutton(imp_before_frame, text='Imputation Before Normalization', variable=self.enable_imputation_before_var, bg='#f0f0f0').pack(side='left')
        imp_before_method_label = tk.Label(imp_before_frame, text='Method:', bg='#f0f0f0', font=('Arial', 8))
        imp_before_method_label.pack(side='left', padx=(8, 3))
        self.imputation_before_method_combo = ttk.Combobox(
            imp_before_frame,
            textvariable=self.imputation_before_method_var,
            values=['half_min', 'median_per_group', 'median_global', 'knn'],
            width=16,
            state='readonly'
        )
        self.imputation_before_method_combo.pack(side='left')
        self.knn_k_before_label = tk.Label(imp_before_frame, text='KNN k:', bg='#f0f0f0', font=('Arial', 8))
        self.knn_k_before_label.pack(side='left', padx=(8, 3))
        self.knn_k_before_entry = tk.Entry(imp_before_frame, textvariable=self.imputation_before_knn_neighbors_var, width=4)
        self.knn_k_before_entry.pack(side='left')

        self._create_tooltip(
            imp_before_method_label,
            "Pre-normalization imputation method:\n"
            "- half_min: fills missing with half of row minimum positive value\n"
            "- median_per_group: fills using row-wise median positive value\n"
            "- median_global: fills using one global median positive value\n"
            "- knn: multivariate KNN imputation\n\n"
            "This runs BEFORE normalization to complete the dataset first."
        )

        optional_proc_frame = tk.LabelFrame(
            step4_frame,
            text='Step 4b: Optional Post-Normalization Steps',
            bg='#f0f0f0'
        )
        optional_proc_frame.pack(fill='x', padx=5, pady=(6, 4))

        self.enable_variability_filter_var = tk.BooleanVar(value=False)
        self.variability_percent_var = tk.StringVar(value='10')
        self.enable_imputation_var = tk.BooleanVar(value=False)
        self.imputation_method_var = tk.StringVar(value='half_min')
        self.knn_neighbors_var = tk.StringVar(value='5')
        self.imputation_min_group_percent_var = tk.StringVar(value='50.0')
        self.imputation_prefilter_scope_var = tk.StringVar(value='per_group')
        self.enable_pca_outlier_var = tk.BooleanVar(value=False)

        self.enable_variability_filter_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Variability filter: {'on' if self.enable_variability_filter_var.get() else 'off'}"))
        self.variability_percent_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.enable_imputation_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Imputation: {'on' if self.enable_imputation_var.get() else 'off'}"))
        self.imputation_method_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.knn_neighbors_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.imputation_min_group_percent_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.imputation_prefilter_scope_var.trace_add('write', lambda *a: self._stats_config_changed())
        self.enable_pca_outlier_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"PCA outlier filter: {'on' if self.enable_pca_outlier_var.get() else 'off'}"))

        var_frame = tk.Frame(optional_proc_frame, bg='#f0f0f0')
        var_frame.pack(fill='x', padx=5, pady=(4, 2))
        tk.Checkbutton(var_frame, text='Variability filter', variable=self.enable_variability_filter_var, bg='#f0f0f0').pack(side='left')
        tk.Label(var_frame, text='Variance percentile cutoff (%):', bg='#f0f0f0', font=('Arial', 8)).pack(side='left', padx=(8, 3))
        ttk.Combobox(var_frame, textvariable=self.variability_percent_var, values=['5', '10', '15', '20'], width=5, state='readonly').pack(side='left')
        # tk.Label(optional_proc_frame,
        #      text='Note: This is a percentile threshold, not exact count removal; ties at the threshold may remove slightly more/fewer rows.',
        #      bg='#f0f0f0', fg='#666666', font=('Arial', 8, 'italic'), wraplength=460, justify='left').pack(fill='x', padx=5, pady=(0, 2))

        imp_frame = tk.Frame(optional_proc_frame, bg='#f0f0f0')
        imp_frame.pack(fill='x', padx=5, pady=2)
        tk.Checkbutton(imp_frame, text='Imputation', variable=self.enable_imputation_var, bg='#f0f0f0').pack(side='left')
        imp_method_label = tk.Label(imp_frame, text='Method:', bg='#f0f0f0', font=('Arial', 8))
        imp_method_label.pack(side='left', padx=(8, 3))
        self.imputation_method_combo = ttk.Combobox(
            imp_frame,
            textvariable=self.imputation_method_var,
            values=['half_min', 'median_per_group', 'median_global', 'knn'],
            width=16,
            state='readonly'
        )
        self.imputation_method_combo.pack(side='left')
        self.knn_k_label = tk.Label(imp_frame, text='KNN k:', bg='#f0f0f0', font=('Arial', 8))
        self.knn_k_label.pack(side='left', padx=(8, 3))
        self.knn_k_entry = tk.Entry(imp_frame, textvariable=self.knn_neighbors_var, width=4)
        self.knn_k_entry.pack(side='left')

        self._create_tooltip(
            imp_method_label,
            "Imputation method:\n"
            "- half_min: fills missing with half of row minimum positive value\n"
            "- median_per_group: fills using row-wise median positive value\n"
            "- median_global: fills using one global median positive value\n"
            "- knn: multivariate KNN imputation"
        )
        self._create_tooltip(
            self.imputation_method_combo,
            "Choose the missing-value imputation strategy for Step 4b.\n"
            "Use KNN for MAR-like missingness; half_min/median methods are simple deterministic options."
        )
        self._create_tooltip(
            self.knn_k_label,
            "Number of nearest neighbors used for KNN imputation.\n"
            "Typical values: 3-10. This setting is only used when method=knn."
        )
        self._create_tooltip(
            self.knn_k_entry,
            "KNN neighborhood size (k). Higher k smooths more; lower k keeps local structure."
        )

        imp_prefilter_frame = tk.Frame(optional_proc_frame, bg='#f0f0f0')
        imp_prefilter_frame.pack(fill='x', padx=5, pady=(0, 2))
        imp_prefilter_label = tk.Label(imp_prefilter_frame, text='Min %/group before impute:', bg='#f0f0f0', font=('Arial', 8))
        imp_prefilter_label.pack(side='left', padx=(24, 3))
        self.imp_prefilter_percent_entry = tk.Entry(imp_prefilter_frame, textvariable=self.imputation_min_group_percent_var, width=6)
        self.imp_prefilter_percent_entry.pack(side='left')
        imp_scope_label = tk.Label(imp_prefilter_frame, text='Scope:', bg='#f0f0f0', font=('Arial', 8))
        imp_scope_label.pack(side='left', padx=(8, 3))
        self.imputation_prefilter_scope_combo = ttk.Combobox(
            imp_prefilter_frame,
            textvariable=self.imputation_prefilter_scope_var,
            values=['per_group', 'all_groups'],
            width=11,
            state='readonly'
        )
        self.imputation_prefilter_scope_combo.pack(side='left')

        self._create_tooltip(
            imp_prefilter_label,
            "Pre-filter threshold before imputation.\n"
            "A value is counted as valid only if it is non-zero and non-missing."
        )
        self._create_tooltip(
            self.imp_prefilter_percent_entry,
            "Minimum valid percentage required before imputation (0-100)."
        )
        self._create_tooltip(
            imp_scope_label,
            "Scope for applying the threshold:\n"
            "- per_group: pass if threshold met in at least one group\n"
            "- all_groups: require threshold in every group"
        )
        self._create_tooltip(
            self.imputation_prefilter_scope_combo,
            "Switch between per_group and all_groups filtering rules before imputation."
        )

        def _set_widget_state_recursive(widget, state: str):
            for child in widget.winfo_children():
                _set_widget_state_recursive(child, state)
            try:
                widget.configure(state=state)
            except Exception:
                pass

        def _update_min_samples_controls_state(*_args):
            imputation_on = bool(self.enable_imputation_var.get())
            new_state = 'disabled' if imputation_on else 'normal'
            if hasattr(self, 'min_samples_controls_container'):
                _set_widget_state_recursive(self.min_samples_controls_container, new_state)
            if hasattr(self, 'min_samples_label'):
                self.min_samples_label.configure(fg=('#9e9e9e' if imputation_on else 'black'))

            # KNN parameter is only relevant when imputation is enabled AND method is KNN.
            knn_active = imputation_on and (self.imputation_method_var.get().strip().lower() == 'knn')
            if hasattr(self, 'knn_k_entry'):
                try:
                    self.knn_k_entry.configure(state=('normal' if knn_active else 'disabled'))
                except Exception:
                    pass
            if hasattr(self, 'knn_k_label'):
                try:
                    self.knn_k_label.configure(fg=('black' if knn_active else '#9e9e9e'))
                except Exception:
                    pass

        self._update_min_samples_controls_state = _update_min_samples_controls_state
        self.enable_imputation_var.trace_add('write', self._update_min_samples_controls_state)
        self.imputation_method_var.trace_add('write', self._update_min_samples_controls_state)
        self._update_min_samples_controls_state()

        pca_frame = tk.Frame(optional_proc_frame, bg='#f0f0f0')
        pca_frame.pack(fill='x', padx=5, pady=2)
        tk.Checkbutton(pca_frame, text='PCA-based sample outlier removal', variable=self.enable_pca_outlier_var, bg='#f0f0f0').pack(side='left')

        # tk.Label(optional_proc_frame,
        #          text='Optional only: disabled by default. \nIf imputation is enabled, min-samples-per-group filtering is bypassed and imputation pre-filter is used.\nReview PCA score plots before final reporting.',
        #          bg='#fff3cd', fg='#856404', font=('Arial', 8), wraplength=460, justify='left').pack(fill='x', padx=5, pady=(2, 5))

        self.run_norm_btn = tk.Button(step4_frame, text='⚙️ Run Normalization & Test Normality', command=self.run_statistics_pipeline,
              bg='#27ae60', fg='white', state='disabled', **btn_style)
        self.run_norm_btn.pack(fill='x', padx=5, pady=2)
        tk.Button(step4_frame, text='💾 Export Normalized Data', command=self.export_normalized_results,
              bg='#8e44ad', fg='white', **btn_style).pack(fill='x', padx=5, pady=2)
        
        # ========== STEP 5: Statistical Tests ==========
        step5_frame = tk.LabelFrame(middle_scrollable, text='Step 5: Statistical Tests', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        step5_frame.pack(fill='x', padx=5, pady=(0, 10))
        
        # Run Statistics and Export buttons at top of Step 5
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        self.run_stats_btn = tk.Button(step5_frame, text='🧪 Run Statistics', command=self.run_statistical_tests,
                  bg='#27ae60', fg='white', state='disabled', **btn_style)
        self.run_stats_btn.pack(fill='x', padx=5, pady=(5, 2))
        tk.Button(step5_frame, text='📈 Export Stat Results', command=self.export_statistical_results,
                  bg='#8e44ad', fg='white', **btn_style).pack(fill='x', padx=5, pady=2)
        
        # Test type selection
        self.stat_test_type = tk.StringVar(value='overall')
        
        overall_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        overall_frame.pack(fill='x', padx=5, pady=(5, 2))
        tk.Radiobutton(overall_frame, text='Overall (>2 groups):', variable=self.stat_test_type,
                       value='overall', bg='#f0f0f0', command=self.on_test_type_change).pack(side='left')
        
        self.stat_overall_test = tk.StringVar(value='anova')
        self.stat_overall_test.trace_add('write', lambda *a: self._stats_config_changed(log=f"Overall test: {self.stat_overall_test.get()}"))
        self.stat_overall_test.trace_add('write', self._update_two_way_button_state)
        # Overall test options with descriptions and Configure button on same row
        test_info_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        test_info_frame.pack(fill='x', padx=15, pady=(0, 2))
        self.overall_combo = ttk.Combobox(test_info_frame, 
                                         values=['anova', 'kruskal', 'two_way_anova', 'nonparametric_two_way'], 
                                         textvariable=self.stat_overall_test, state='readonly', width=22)
        self.overall_combo.pack(side='left', padx=(0, 5))
        
        # Info label that updates based on selection (placed after button, expands to the right)
        self.overall_test_info = tk.Label(test_info_frame, text='', bg='#f0f0f0',
                                          font=('Arial', 8, 'italic'), fg='#666', wraplength=400, justify='left')
        self.overall_test_info.pack(side='left', padx=10, fill='x', expand=True)

        def _update_overall_test_info(*args):
            test = self.stat_overall_test.get()
            info_text = {
                'One-way-anova': '→ Parametric',
                'Kruskal-Wallis ': 'Non-parametric',
                'Two_way_anova': 'Parametric',
                'Nonparametric_two_way': '→ ART/Rank'
            }.get(test, '')
            self.overall_test_info.config(text=info_text)
        
        self.stat_overall_test.trace_add('write', _update_overall_test_info)
        _update_overall_test_info()  # Set initial text

        # Two-Way ANOVA now launches a console for automated factor detection
        def _open_two_way_anova_console():
            # Reuse the same prerequisite flow as Run Statistics:
            # if groups are not configured yet, open Configure Groups first,
            # then continue directly to Two-Way ANOVA configuration.
            if not self._ensure_groups_ready(
                after_config_callback=_open_two_way_anova_console,
                action_label='two-way ANOVA configuration'
            ):
                return
            
            from main_script.factor_mapping_manager import FactorMappingConfig
            
            win = tk.Toplevel(self.root)
            win.title('Two-Way ANOVA Factor Configuration')
            win.geometry('900x900')
            try:
                win.minsize(750, 700)
                win.resizable(True, True)
            except Exception:
                pass
            win.configure(bg='#f0f0f0')
            
            # Title
            title_frame = tk.Frame(win, bg='#3498db')
            title_frame.pack(fill='x')
            tk.Label(title_frame, text='Two-Way ANOVA Factor Configuration', font=('Arial',14,'bold'), 
                     bg='#3498db', fg='white').pack(pady=10)
            
            # Top action buttons (Confirm/Cancel) fixed at top
            button_frame_top = tk.Frame(win, bg='#f0f0f0')
            button_frame_top.pack(side='top', fill='x', padx=10, pady=(8,6))
            btn_confirm = tk.Button(button_frame_top, text='Confirm', bg='#27ae60', fg='white', width=15, command=lambda: None)
            btn_confirm.pack(side='left', padx=5)
            btn_cancel = tk.Button(button_frame_top, text='Cancel', bg='#e74c3c', fg='white', width=15, command=lambda: None)
            btn_cancel.pack(side='left', padx=5)

            # Scrollable content area below buttons
            container = tk.Frame(win, bg='#f0f0f0')
            container.pack(side='top', fill='both', expand=True)
            canvas = tk.Canvas(container, bg='#f0f0f0', highlightthickness=0)
            vscroll = tk.Scrollbar(container, orient='vertical', command=canvas.yview)
            vscroll.pack(side='right', fill='y')
            canvas.pack(side='left', fill='both', expand=True)
            canvas.configure(yscrollcommand=vscroll.set)
            scroll_content = tk.Frame(canvas, bg='#f0f0f0')
            canvas.create_window((0,0), window=scroll_content, anchor='nw')
            def _on_sc_content_configure(event):
                try:
                    canvas.configure(scrollregion=canvas.bbox('all'))
                except Exception:
                    pass
            scroll_content.bind('<Configure>', _on_sc_content_configure)
            
            # Get initial heuristic mappings from column/group names
            groups = sorted({v.get() for v in self.sample_group_vars.values() if v.get()})
            
            # Enhanced heuristic for Factor A (Treatment/Injury)
            def _factA(g):
                gl = g.lower()
                # Check for TBI/Sham patterns
                if 'tbi' in gl:
                    return 'TBI'
                if 'sham' in gl:
                    return 'Sham'
                # Check for treatment patterns
                if 'treat' in gl or 'drug' in gl:
                    return 'Treatment'
                # Check for control patterns
                if 'ctl' in gl or 'control' in gl or 'ctrl' in gl:
                    return 'Control'
                return 'OTHER'
            
            # Enhanced heuristic for Factor B (Diet/Condition)
            def _factB(g):
                gl = g.lower()
                # Check for HFD/NC patterns (matching R script)
                if 'hfd' in gl or 'high' in gl:
                    return 'HFD'
                if 'nc' in gl or 'normal' in gl or 'regular' in gl or 'control' in gl:
                    return 'NC'
                # Other diet patterns
                if 'chow' in gl:
                    return 'Chow'
                return 'Regular'
            
            factorA_map_group = {g:_factA(g) for g in groups}
            factorB_map_group = {g:_factB(g) for g in groups}
            
            # Build sample-level maps
            sample_factorA = {}
            sample_factorB = {}
            for col, gvar in self.sample_group_vars.items():
                g = gvar.get()
                sample_factorA[col] = factorA_map_group.get(g,'OTHER')
                sample_factorB[col] = factorB_map_group.get(g,'Regular')
            
            # Initialize config manager
            config_manager = FactorMappingConfig()
            saved_configs = config_manager.list_mappings()
            
            # ===== Load Previous Config Section =====
            load_frame = tk.LabelFrame(scroll_content, text='📂 Load Previous Configuration (Optional)', 
                                      bg='#f0f0f0', font=('Arial',10,'bold'))
            load_frame.pack(fill='x', padx=10, pady=(10,5))
            
            config_var = tk.StringVar(value='<new>')
            if saved_configs:
                def _load_saved_config(*_):
                    cfg_name = config_var.get()
                    if cfg_name and cfg_name != '<new>':
                        cfg = config_manager.load_mapping(cfg_name)
                        if cfg:
                            # Update Factor A name and mappings
                            factorA_name_entry.delete(0, tk.END)
                            factorA_name_entry.insert(0, cfg['factorA']['name'])
                            for grp in groups:
                                if grp in cfg['factorA']['groups']:
                                    factorA_edits[grp].delete(0, tk.END)
                                    factorA_edits[grp].insert(0, cfg['factorA']['groups'][grp])
                            
                            # Update Factor B name and mappings
                            factorB_name_entry.delete(0, tk.END)
                            factorB_name_entry.insert(0, cfg['factorB']['name'])
                            for grp in groups:
                                if grp in cfg['factorB']['groups']:
                                    factorB_edits[grp].delete(0, tk.END)
                                    factorB_edits[grp].insert(0, cfg['factorB']['groups'][grp])
                            
                            _update_design_summary()
                
                config_combo = ttk.Combobox(load_frame, values=['<new>'] + saved_configs, 
                                           textvariable=config_var, state='readonly', width=40)
                config_combo.pack(side='left', padx=10, pady=5)
                config_combo.bind('<<ComboboxSelected>>', _load_saved_config)
            
            # ===== Factor A Configuration =====
            factorA_frame = tk.LabelFrame(scroll_content, text='⚙️ Factor A Configuration', 
                                         bg='#f0f0f0', font=('Arial',10,'bold'))
            factorA_frame.pack(fill='x', padx=10, pady=5)
            
            # Factor A name with better defaults
            name_frame = tk.Frame(factorA_frame, bg='#f0f0f0')
            name_frame.pack(fill='x', padx=10, pady=5)
            tk.Label(name_frame, text='Factor A Name (e.g., "Treatment", "Injury"):', bg='#f0f0f0').pack(side='left')
            factorA_name_entry = tk.Entry(name_frame, width=30)
            factorA_name_entry.pack(side='left', padx=5)
            # Smart default based on detected patterns
            default_factorA = 'Treatment' if any('TBI' in _factA(g) or 'Sham' in _factA(g) for g in groups) else 'Factor A'
            factorA_name_entry.insert(0, default_factorA)
            
            # Factor A group assignments
            factorA_edits = {}
            for grp in groups:
                grp_frame = tk.Frame(factorA_frame, bg='#f0f0f0')
                grp_frame.pack(fill='x', padx=20, pady=3)
                tk.Label(grp_frame, text=f'{grp}:', width=15, bg='#f0f0f0').pack(side='left')
                edit = tk.Entry(grp_frame, width=15)
                edit.pack(side='left', padx=5)
                edit.insert(0, factorA_map_group[grp])
                factorA_edits[grp] = edit
            
            # ===== Factor B Configuration =====
            factorB_frame = tk.LabelFrame(scroll_content, text='⚙️ Factor B Configuration', 
                                         bg='#f0f0f0', font=('Arial',10,'bold'))
            factorB_frame.pack(fill='x', padx=10, pady=5)
            
            # Factor B name with better defaults
            name_frame_b = tk.Frame(factorB_frame, bg='#f0f0f0')
            name_frame_b.pack(fill='x', padx=10, pady=5)
            tk.Label(name_frame_b, text='Factor B Name (e.g., "Diet", "Condition"):', bg='#f0f0f0').pack(side='left')
            factorB_name_entry = tk.Entry(name_frame_b, width=30)
            factorB_name_entry.pack(side='left', padx=5)
            # Smart default based on detected patterns
            default_factorB = 'Diet' if any('HFD' in _factB(g) or 'NC' in _factB(g) for g in groups) else 'Factor B'
            factorB_name_entry.insert(0, default_factorB)
            
            # Factor B group assignments
            factorB_edits = {}
            for grp in groups:
                grp_frame = tk.Frame(factorB_frame, bg='#f0f0f0')
                grp_frame.pack(fill='x', padx=20, pady=3)
                tk.Label(grp_frame, text=f'{grp}:', width=15, bg='#f0f0f0').pack(side='left')
                edit = tk.Entry(grp_frame, width=15)
                edit.pack(side='left', padx=5)
                edit.insert(0, factorB_map_group[grp])
                factorB_edits[grp] = edit
            
            # ===== Design Summary Section =====
            summary_frame = tk.LabelFrame(scroll_content, text='📊 Design Summary (Two-Way ANOVA Structure)', 
                                         bg='#f0f0f0', font=('Arial',10,'bold'))
            summary_frame.pack(fill='both', expand=True, padx=10, pady=5)
            
            # Add helpful hint
            hint_label = tk.Label(summary_frame, 
                                 text="This shows your experimental design matching R's aov(value ~ Factor_A * Factor_B) approach",
                                 bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#555', wraplength=700, justify='left')
            hint_label.pack(fill='x', padx=5, pady=(5,0))
            
            # Add a vertical scrollbar to the summary text area
            summary_scroll = tk.Scrollbar(summary_frame, orient='vertical')
            summary_scroll.pack(side='right', fill='y')
            summary_text = tk.Text(summary_frame, height=12, width=80, bg='white', yscrollcommand=summary_scroll.set,
                                  font=('Consolas', 9))
            summary_text.pack(fill='both', expand=True, padx=5, pady=5)
            summary_scroll.config(command=summary_text.yview)
            
            def _update_design_summary():
                # Collect current values (group -> level mapping)
                current_factorA = {grp: factorA_edits[grp].get() for grp in groups}
                current_factorB = {grp: factorB_edits[grp].get() for grp in groups}

                # Compute sample-level cell counts using current mappings
                from collections import Counter
                cell_counts = Counter()
                group_counts = Counter()  # Track how many samples per named group
                total_samples = 0
                try:
                    for col, gvar in getattr(self, 'sample_group_vars', {}).items():
                        grp = gvar.get()
                        if not grp:
                            continue
                        a_lvl = current_factorA.get(grp, 'OTHER')
                        b_lvl = current_factorB.get(grp, 'Regular')
                        cell_counts[(a_lvl, b_lvl)] += 1
                        group_counts[grp] += 1
                        total_samples += 1
                except Exception:
                    pass

                factorA_levels = sorted(set(current_factorA.values()))
                factorB_levels = sorted(set(current_factorB.values()))
                counts = list(cell_counts.values())
                min_reps = min(counts) if counts else 0
                max_reps = max(counts) if counts else 0
                balanced = (len(set(counts)) == 1) if counts else False

                # Build combined group labels (matching R: paste(Diet, Treatment, sep="_"))
                combined_groups = {}  # (a_lvl, b_lvl) -> list of original group names
                for grp in groups:
                    a_lvl = current_factorA.get(grp, 'OTHER')
                    b_lvl = current_factorB.get(grp, 'Regular')
                    key = (a_lvl, b_lvl)
                    if key not in combined_groups:
                        combined_groups[key] = []
                    combined_groups[key].append(grp)

                summary_text.config(state='normal')
                summary_text.delete('1.0', tk.END)
                
                # Get factor names
                fa_name = factorA_name_entry.get().strip() or 'Factor A'
                fb_name = factorB_name_entry.get().strip() or 'Factor B'
                
                summary_text.insert('end', f"═══════ Two-Way ANOVA Design ═══════\n\n", 'bold')
                summary_text.insert('end', f"{fa_name} Levels: {', '.join(factorA_levels)}\n")
                summary_text.insert('end', f"{fb_name} Levels: {', '.join(factorB_levels)}\n")
                summary_text.insert('end', f"Total Samples: {total_samples}\n")
                summary_text.insert('end', f"Design: {'✓ Balanced' if balanced else '⚠ Unbalanced'}\n")
                summary_text.insert('end', f"Replicates per Cell: {min_reps}-{max_reps}\n\n")
                
                summary_text.insert('end', f"─── Cell Structure ({fa_name} × {fb_name}) ───\n")
                for (a, b), count in sorted(cell_counts.items()):
                    combined_name = f"{a}_{b}"
                    orig_groups = combined_groups.get((a, b), [])
                    summary_text.insert('end', f"  {combined_name}: {count} samples")
                    if orig_groups:
                        summary_text.insert('end', f" (from {', '.join(orig_groups)})")
                    summary_text.insert('end', '\n')
                
                summary_text.insert('end', f"\n─── ANOVA Formula ───\n")
                summary_text.insert('end', f"  value ~ {fa_name} * {fb_name}\n")
                summary_text.insert('end', f"\n─── Post-hoc Test ───\n")
                summary_text.insert('end', f"  Tukey HSD on Group factor\n")
                summary_text.insert('end', f"  ({len(combined_groups)} group comparisons)\n")
                
                summary_text.config(state='disabled')
            
            # Apply tag for bold text
            summary_text.tag_config('bold', font=('Arial', 9, 'bold'))
            
            # Initial summary update
            _update_design_summary()
            
            # Bind entries to update summary on change
            for edit in list(factorA_edits.values()) + list(factorB_edits.values()):
                edit.bind('<KeyRelease>', lambda e: _update_design_summary())
            
            # ===== Save Configuration Section =====
            save_frame = tk.Frame(scroll_content, bg='#f0f0f0')
            save_frame.pack(fill='x', padx=10, pady=5)
            tk.Label(save_frame, text='Save Configuration Name (optional):', bg='#f0f0f0').pack(side='left')
            save_name_entry = tk.Entry(save_frame, width=30)
            save_name_entry.pack(side='left', padx=5)
            
            def _apply_and_save():
                # Collect final mappings
                final_factorA_map = {grp: factorA_edits[grp].get() for grp in groups}
                final_factorB_map = {grp: factorB_edits[grp].get() for grp in groups}
                
                # Build sample-level maps
                self.sample_factorA_vars = {}
                self.sample_factorB_vars = {}
                for col, gvar in self.sample_group_vars.items():
                    g = gvar.get()
                    factorA_level = final_factorA_map.get(g, 'OTHER')
                    factorB_level = final_factorB_map.get(g, 'Regular')
                    self.sample_factorA_vars[col] = tk.StringVar(value=factorA_level)
                    self.sample_factorB_vars[col] = tk.StringVar(value=factorB_level)

                # Persist factor display names on the instance for downstream labeling
                try:
                    self.factorA_name = factorA_name_entry.get().strip() or 'Factor A'
                    self.factorB_name = factorB_name_entry.get().strip() or 'Factor B'
                except Exception:
                    self.factorA_name = 'Factor A'
                    self.factorB_name = 'Factor B'
                
                # Save configuration if name provided
                save_name = save_name_entry.get().strip()
                if save_name:
                    config_manager.save_mapping(
                        save_name,
                        factorA_name_entry.get(),
                        factorB_name_entry.get(),
                        final_factorA_map,
                        final_factorB_map,
                        notes=f"Auto-saved configuration from {len(groups)} groups"
                    )
                    self._thread_safe_log(f'✅ Configuration "{save_name}" saved.\n')
                
                self._thread_safe_log('✅ Factor assignments confirmed. Use "Run Statistics" button to execute Two-Way ANOVA.\n')
                win.destroy()
            
            def _cancel():
                win.destroy()
            
            # Wire up top buttons now that handlers exist
            try:
                btn_confirm.config(command=_apply_and_save)
                btn_cancel.config(command=_cancel)
            except Exception:
                pass
        def _on_overall_test_change(*_):
            val = self.stat_overall_test.get().lower()
            if val == 'two_way_anova':
                # Defer console until user explicitly runs stats (avoid auto-popup on config load)
                pass
        self.stat_overall_test.trace_add('write', _on_overall_test_change)
        
        # Add explicit button to launch two-way ANOVA console (works for both parametric and non-parametric)
        def _launch_console_button():
            current_test = self.stat_overall_test.get().lower()
            if current_test in ['two_way_anova', 'nonparametric_two_way']:
                _open_two_way_anova_console()
            else:
                messagebox.showinfo('Two-Way ANOVA', 'Select "two_way_anova" or "nonparametric_two_way" from Overall test dropdown first.')
        
        # Configure button sits immediately to the right of the dropdown for left alignment
        self.two_way_config_btn = tk.Button(
            test_info_frame,
            text='🧬 Configure TWO-Way ANOVA',
            command=_launch_console_button,
            bg='#3498db', fg='white', font=('Arial', 8, 'bold'),
            relief='raised', bd=2, padx=8, pady=3, state='disabled'
        )
        # Place before the info label by using pack with before parameter
        self.two_way_config_btn.pack(side='left', padx=(5, 10), before=self.overall_test_info)
        
        # ========== NON-PARAMETRIC TWO-WAY CONFIGURATION ==========
        nonparam_config_frame = tk.LabelFrame(step5_frame, text='⚗️ Non-Parametric Two-Way Settings', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        # Only show this frame when nonparametric_two_way is selected
        
        def _update_nonparam_visibility(*args):
            if self.stat_overall_test.get() == 'nonparametric_two_way':
                nonparam_config_frame.pack(fill='x', padx=5, pady=(5,5), after=step5_frame.winfo_children()[-1])
            else:
                nonparam_config_frame.pack_forget()
        
        self.stat_overall_test.trace_add('write', _update_nonparam_visibility)
        
        # Method selection
        method_frame = tk.Frame(nonparam_config_frame, bg='#f0f0f0')
        method_frame.pack(fill='x', padx=10, pady=(5,2))
        tk.Label(method_frame, text='Analysis Method:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w')
        
        self.nonparam_method = tk.StringVar(value='art')
        self.nonparam_method.trace_add('write', lambda *a: self._stats_config_changed(log=f"Non-parametric method: {self.nonparam_method.get()}"))
        
        # Radio buttons with descriptions
        tk.Radiobutton(method_frame, text='ART (Aligned Rank Transform)', variable=self.nonparam_method, 
                      value='art', bg='#f0f0f0').pack(anchor='w', padx=10)
        tk.Label(method_frame, text='   → RECOMMENDED: Tests main effects + interaction, preserves power', 
                bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#27ae60', wraplength=400, justify='left').pack(anchor='w', padx=20)
        
        tk.Radiobutton(method_frame, text='Rank-Transformed ANOVA', variable=self.nonparam_method, 
                      value='rank', bg='#f0f0f0').pack(anchor='w', padx=10, pady=(5,0))
        tk.Label(method_frame, text='   → Simpler rank approach', 
                bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#666', wraplength=400, justify='left').pack(anchor='w', padx=20)
        
        # # Use same factor configuration as parametric two-way
        # nonparam_note = tk.Label(nonparam_config_frame, 
        #                         text='📌 Factor assignments: Use "Configure Two-Way ANOVA" button above to set up factors\n'
        #                              '   (same configuration applies to both parametric and non-parametric methods)\n'
        #                              '   Post-hoc tests use Dunn method with BH correction (configure in console if needed)',
        #                         bg='#fff3cd', fg='#856404', font=('Arial', 8), wraplength=450, justify='left', padx=8, pady=5)
        # nonparam_note.pack(fill='x', padx=10, pady=(5,10))
        
        # Post-hoc is always enabled with BH correction by default (can be changed in configure console)
        self.nonparam_posthoc = tk.BooleanVar(value=True)
        
        pairwise_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        pairwise_frame.pack(fill='x', padx=5, pady=2)
        tk.Radiobutton(pairwise_frame, text='Pairwise:', variable=self.stat_test_type,
                       value='pairwise', bg='#f0f0f0', command=self.on_test_type_change).pack(side='left')
        
        # Pairwise test combo with Configure ROTS button on same row
        pairwise_combo_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        pairwise_combo_frame.pack(fill='x', padx=15, pady=(0, 5))
        
        self.stat_pairwise_test = tk.StringVar(value='welch')
        self.stat_pairwise_test.trace_add('write', lambda *a: self._stats_config_changed(log=f"Pairwise test: {self.stat_pairwise_test.get()}"))
        self.stat_pairwise_test.trace_add('write', lambda *a: self._toggle_rots_parameters())
        self.stat_pairwise_test.trace_add('write', self._update_rots_button_state)
        self.pairwise_combo = ttk.Combobox(pairwise_combo_frame, values=['welch', 'mannwhitney', 'rots', 'limma'], textvariable=self.stat_pairwise_test, state='readonly', width=22)
        self.pairwise_combo.pack(side='left', padx=(0, 5))
        # Start disabled - will be enabled when pairwise is selected
        self.pairwise_combo.config(state='disabled')
        
        # ROTS Configure button (placed beside pairwise dropdown)
        self.rots_config_button = tk.Button(
            pairwise_combo_frame,
            text='⚙️ Configure ROTS Parameters',
            command=self._open_rots_config_dialog,
            bg='#3498db', fg='white', font=('Arial', 8, 'bold'),
            relief='raised', bd=2, padx=8, pady=3, state='disabled'
        )
        self.rots_config_button.pack(side='left', padx=(5, 0))
        
        # ROTS Parameters (hidden - accessed via dialog)
        # Store as instance variables for the dialog to access
        if not hasattr(self, 'rots_B'):
            self.rots_B = tk.StringVar(value='1000')
        if not hasattr(self, 'rots_K'):
            self.rots_K = tk.StringVar(value='100')
        if not hasattr(self, 'rots_alpha'):
            self.rots_alpha = tk.StringVar(value='0.1')
        if not hasattr(self, 'rots_seed'):
            self.rots_seed = tk.StringVar(value='42')

        # Pairwise p-value adjustment method (applies ONLY to pairwise comparisons)
        adj_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        adj_frame.pack(fill='x', padx=5, pady=(4, 2))
        tk.Label(adj_frame, text='Pairwise p-value adjustment:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w')
        tk.Label(adj_frame, text='Applies only to pairwise tests (BH default)', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w')
        self.pairwise_p_adjust_method = tk.StringVar(value='BH')
        self.pairwise_p_adjust_method.trace_add('write', lambda *a: self._stats_config_changed(log=f"Pairwise p-adjust: {self.pairwise_p_adjust_method.get()}"))
        ttk.Combobox(adj_frame, values=['BH','Bonferroni','Holm','Hochberg','BY','None'], textvariable=self.pairwise_p_adjust_method, state='readonly').pack(fill='x', padx=15, pady=(2, 5))

        # Note: Processing Settings, Base Group, and Custom Comparisons moved to Step 3 (left column)
        
        # FDR Scope
        fdr_frame = tk.Frame(step5_frame, bg='#f0f0f0')
        fdr_frame.pack(fill='x', padx=5, pady=(6, 2))
        tk.Label(fdr_frame, text='FDR Scope:', bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w', pady=(0, 2))
        self.fdr_scope_var = tk.StringVar(value='per-comparison')
        self.fdr_scope_var.trace_add('write', lambda *a: self._on_fdr_scope_changed())
        fdr_radio_frame = tk.Frame(fdr_frame, bg='#f0f0f0')
        fdr_radio_frame.pack(fill='x', pady=(0, 5))
        for txt, val in [('Per-Comparison', 'per-comparison'), ('Per-Feature', 'per-metabolite')]:
            tk.Radiobutton(fdr_radio_frame, text=txt, variable=self.fdr_scope_var, value=val, bg='#f0f0f0').pack(side='left', padx=2)
        
        # Add help text for FDR scope
        fdr_help = tk.Label(fdr_frame, text='⚠️ Only use when comparing 3+ groups \n (returns identical p-values with 2 groups)', 
                           bg='#fff3cd', fg='#856404', font=('Arial', 8), wraplength=450, justify='left', padx=5, pady=3)
        # Don't pack initially - will be shown only when per-feature is selected
        self.fdr_scope_warning = fdr_help  # Store reference for dynamic updates
        
        # Note: Group Order moved to Step 3 (left column)

        # Adjusted p-values toggle removed - not used

        # ========== STEP 5b: Covariate Adjustment Button (Optional) ==========
        covariate_frame = tk.LabelFrame(middle_scrollable, text='🎯 Step 5b: Covariate Adjustment (Optional)',
                                       bg='#f0f0f0', font=('Arial', 10, 'bold'))
        covariate_frame.pack(fill='x', padx=5, pady=(10, 5))
        
        # Info label
        info_label = tk.Label(
            covariate_frame,
            text='Adjust for covariates (Age, Sex, BMI, etc.) using linear regression',
            bg='#e3f2fd', fg='#1565c0', font=('Arial', 9), wraplength=450, justify='left'
        )
        info_label.pack(fill='x', padx=5, pady=5)
        
        # Button to open covariate dialog
        self.open_covariate_btn = tk.Button(
            covariate_frame,
            text='🎯 Run Stat with Covariate Adjustment',
            command=self.open_covariate_dialog,
            bg='#8e44ad', fg='white', font=('Arial', 10, 'bold'),
            relief='raised', bd=3, pady=8
        )
        self.open_covariate_btn.pack(fill='x', padx=10, pady=10)

        # ========================================================================
        # RIGHT COLUMN: Progress & Results (with scrollbar for small screens)
        # ========================================================================
        right_col = tk.LabelFrame(body, text='📊 Statistics Log', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        right_col.grid(row=0, column=2, sticky='nsew', padx=(3, 0))
        
        # Canvas with scrollbar for right column
        right_canvas = tk.Canvas(right_col, bg='#f0f0f0', highlightthickness=0, height=700)
        right_scrollbar = ttk.Scrollbar(right_col, orient="vertical", command=right_canvas.yview)
        right_scrollable = tk.Frame(right_canvas, bg='#f0f0f0')
        
        right_scrollable.bind(
            "<Configure>",
            lambda e: right_canvas.configure(scrollregion=right_canvas.bbox("all"))
        )
        
        right_canvas_window = right_canvas.create_window((0, 0), window=right_scrollable, anchor="nw")
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        def configure_right_scroll(event):
            right_canvas.configure(scrollregion=right_canvas.bbox("all"))
            right_canvas.itemconfig(right_canvas_window, width=event.width)
        
        right_canvas.bind('<Configure>', configure_right_scroll)
        
        right_scrollbar.pack(side="right", fill="y", padx=(2, 0))
        right_canvas.pack(side="left", fill="both", expand=True)
        
        def _on_right_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        right_canvas.bind("<MouseWheel>", _on_right_mousewheel)
        right_scrollable.bind("<MouseWheel>", _on_right_mousewheel)
        
        # Use right_scrollable as parent for all right column content
        results = right_scrollable
        
        # ========== OPEN FOLDER BUTTON ==========
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        tk.Button(results, text='📂 Open Folder', command=self.open_statistics_output_folder,
                  bg='#9b59b6', fg='white', **btn_style).pack(fill='x', padx=5, pady=(5, 5))
        
        # Add progress bar and label above the log
        self.stats_progress_label = tk.Label(results, text="", bg='#f0f0f0', font=('Arial', 9))
        self.stats_progress_label.pack(fill='x', padx=5, pady=(5, 0))
        self.stats_progress_label.pack_forget()  # Hide initially
        
        self.stats_progress = ttk.Progressbar(results, mode='indeterminate', length=400)
        self.stats_progress.pack(fill='x', padx=5, pady=(2, 2))
        self.stats_progress.pack_forget()  # Hide initially
        
        # Statistics log - increased height for maximum space
        self.stats_log = scrolledtext.ScrolledText(results, font=('Courier', 9), wrap=tk.WORD, height=35)
        self.stats_log.pack(fill='both', expand=True, padx=5, pady=(2, 5))
        self.stats_log.insert(tk.END, '📊 Statistics Log Ready\n')
        self.stats_log.insert(tk.END, 'Ready for statistical analysis operations...\n\n')
        
        # ========== Processing Settings (below log) ==========
        proc_frame_right = tk.LabelFrame(results, text='⚙️ Processing Settings', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        proc_frame_right.pack(fill='x', padx=5, pady=(5, 5))
        try:
            max_workers = os.cpu_count() or 4
        except Exception:
            max_workers = 4
        tk.Label(proc_frame_right, text=f"System has {max_workers} CPU cores available", bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(4,2))
        workers_row_right = tk.Frame(proc_frame_right, bg='#f0f0f0')
        workers_row_right.pack(fill='x', padx=5, pady=(0,6))
        tk.Label(workers_row_right, text='Parallel workers:', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
        if not hasattr(self, 'stats_workers'):
            self.stats_workers = tk.StringVar(value='3')
        workers_spin_right = tk.Spinbox(workers_row_right, from_=1, to=max_workers, textvariable=self.stats_workers, width=5)
        workers_spin_right.pack(side='left', padx=(6,0))
        tk.Button(workers_row_right, text='Auto', command=lambda: self.stats_workers.set(str(max_workers)), bg='#9b59b6', fg='white', font=('Arial', 8, 'bold')).pack(side='left', padx=(8,0))
        
        try:
            try:
                importlib.import_module('main_script.metabolite_statistics_analysis')
            except Exception:
                importlib.import_module('metabolite_statistics_analysis')
        except Exception as e:
            self.stats_log.insert(tk.END, f'Import error metabolite_statistics_analysis: {e}\n')
        self.statistics_results = {}

        # Provide a lightweight listbox for sample columns if not yet created (avoid attribute errors)
        if not hasattr(self, 'sample_cols_list'):
            # Hidden container (not displayed) but ensures attribute exists for methods expecting it
            hidden_frame = tk.Frame(container, height=1, width=1)
            hidden_frame.pack_forget()
            self.sample_cols_list = tk.Listbox(hidden_frame)

        # Ensure main scrollregion is set
        self._update_main_scroll_region()
        
        # Load persisted configuration now that widgets are ready
        try:
            self._load_statistics_config()
            self._stats_config_loaded = True
        except Exception as e:
            self.stats_log.insert(tk.END, f"⚠️ Could not auto-load statistics configuration: {e}\n")
        
        # Set initial combo state based on test_type (overall/pairwise)
        self.on_test_type_change()
        
        self.stats_log.insert(tk.END, '🔁 Statistics configuration auto-loaded (if available).\n')

    def show_stats_progress(self, message="Processing..."):
        """Show (or configure) the statistics progress bar.

        If a total step count was set previously via update_stats_progress it will
        switch to determinate mode; otherwise it stays indeterminate.
        """
        if hasattr(self, 'stats_progress'):
            # If a total was defined, use determinate mode; else indeterminate spinner
            if getattr(self, '_stats_total_steps', None):
                self.stats_progress.configure(mode='determinate', maximum=self._stats_total_steps, value=getattr(self, '_stats_current_step', 0))
            else:
                self.stats_progress.configure(mode='indeterminate')
                self.stats_progress.start(10)
            self.stats_progress_label.config(text=message)
            # Show progress widgets (they were pack_forget initially)
            if not self.stats_progress_label.winfo_ismapped():
                self.stats_progress_label.pack(fill='x', padx=5, pady=(5, 0))
            if not self.stats_progress.winfo_ismapped():
                self.stats_progress.pack(fill='x', padx=5, pady=(2, 2))
            self.stats_log.insert(tk.END, f"🔄 {message}\n")
            self.stats_log.see(tk.END)
            self.root.update()

    def set_statistics_layout(self, *, middle_height: int | None = None, log_lines: int | None = None,
                              column_mins: tuple[int, int, int] | None = None,
                              top_bottom_ratio: tuple[int, int] | None = None):
        """Update statistics tab layout values at runtime.

        Parameters (all optional):
        - middle_height: minimum height in pixels for the central statistics area
        - log_lines: number of text lines for the right-hand log (set None to allow grid sizing)
        - column_mins: (min_col0, min_col1, min_col2) minimum widths for columns
        - top_bottom_ratio: (top, bottom) integer ratio for group management rows
        """
        try:
            if middle_height is not None:
                self.stats_middle_height = int(middle_height)
                # Apply minimum height constraint to ensure content visibility
                if hasattr(self, 'stats_body'):
                    self.stats_body.grid_rowconfigure(0, weight=1, minsize=self.stats_middle_height)

            if log_lines is not None:
                self.stats_log_lines = int(log_lines) if log_lines else None
                if hasattr(self, 'stats_log'):
                    try:
                        # Recreate or reconfigure the scrolled text height setting
                        if self.stats_log_lines:
                            self.stats_log.configure(height=self.stats_log_lines)
                        else:
                            # Remove explicit height by setting to a small value and rely on grid
                            self.stats_log.configure(height=1)
                    except Exception:
                        pass

            if column_mins is not None:
                try:
                    self.stats_column_mins = tuple(int(x) for x in column_mins)
                    # Don't apply minsize constraints - let columns expand freely
                    # Removed: self.stats_body.grid_columnconfigure(..., minsize=...)
                except Exception:
                    pass

            if top_bottom_ratio is not None:
                try:
                    self.stats_middle_ratio = (int(top_bottom_ratio[0]), int(top_bottom_ratio[1]))
                    if hasattr(self, 'stats_grp_mgmt'):
                        try:
                            tr, br = self.stats_middle_ratio
                            # Only set weight ratios, don't set minsize constraints
                            self.stats_grp_mgmt.grid_rowconfigure(0, weight=tr)
                            self.stats_grp_mgmt.grid_rowconfigure(1, weight=br)
                        except Exception:
                            pass
                except Exception:
                    pass

            # Refresh scroll regions/layout helpers if present
            try:
                if hasattr(self, '_update_assignment_scroll_region'):
                    self._update_assignment_scroll_region()
                if hasattr(self, '_update_main_scroll_region'):
                    self._update_main_scroll_region()
            except Exception:
                pass
        except Exception:
            pass

    def set_stats_log_button_heights(self, check_btn_pady=3, action_btn_pady=3, 
                                      progress_label_height=None, progress_bar_height=None):
        """
        Adjust the fixed heights of buttons and progress elements in the Statistics Log panel.
        
        Parameters:
            check_btn_pady (int): Vertical padding for "Check Available Data" button (default: 3)
            action_btn_pady (int): Vertical padding for action buttons (Run Tests, Export) (default: 3)
            progress_label_height (int|None): Height for progress label in pixels (None = auto)
            progress_bar_height (int|None): Height for progress bar in pixels (None = auto)
        
        Example usage:
            # Make buttons more compact
            self.set_stats_log_button_heights(check_btn_pady=2, action_btn_pady=2)
            
            # Make buttons taller
            self.set_stats_log_button_heights(check_btn_pady=6, action_btn_pady=6)
        """
        try:
            # Update button style if stats tab exists
            if hasattr(self, 'stats_tab'):
                # The buttons are already created, so we'd need to recreate them
                # or just store these values for next creation
                self._stats_check_btn_pady = check_btn_pady
                self._stats_action_btn_pady = action_btn_pady
                self._stats_progress_label_height = progress_label_height
                self._stats_progress_bar_height = progress_bar_height
                
                # Log the change
                if hasattr(self, 'stats_log'):
                    import time
                    ts = time.strftime('%H:%M:%S')
                    self.stats_log.insert('end', 
                        f"[{ts}] Layout updated: button pady={action_btn_pady}, "
                        f"progress heights adjusted\n")
                    self.stats_log.see('end')
                
                print(f"Stats log layout settings updated. Restart tab or app to see changes.")
                print(f"  Check button pady: {check_btn_pady}")
                print(f"  Action buttons pady: {action_btn_pady}")
                if progress_label_height:
                    print(f"  Progress label height: {progress_label_height}px")
                if progress_bar_height:
                    print(f"  Progress bar height: {progress_bar_height}px")
        except Exception as e:
            print(f"Error adjusting stats log layout: {e}")

    def hide_stats_progress(self):
        """Hide the progress bar"""
        if hasattr(self, 'stats_progress'):
            self.stats_progress.stop()
            self.stats_progress.pack_forget()
            self.stats_progress_label.pack_forget()
            self.root.update()

    def update_stats_progress(self, step: int, message: str = ""):
        """Update determinate statistics progress bar.

        Call once before steps with step == 0 (or 1) to initialize total via
        setting self._stats_total_steps externally, then increment.
        """
        # Default total if not explicitly set
        if not hasattr(self, '_stats_total_steps') or self._stats_total_steps is None:
            # Provide a sane default so bar renders determinately
            self._stats_total_steps = max(step, 1)
        self._stats_current_step = max(0, min(step, self._stats_total_steps))
        if hasattr(self, 'stats_progress'):
            try:
                self.stats_progress.configure(mode='determinate', maximum=self._stats_total_steps)
                self.stats_progress['value'] = self._stats_current_step
            except Exception:
                pass
        if message:
            if hasattr(self, 'stats_progress_label'):
                self.stats_progress_label.config(text=message)
            if hasattr(self, 'stats_log'):
                self.stats_log.insert(tk.END, f"[{self._stats_current_step}/{self._stats_total_steps}] {message}\n")
                self.stats_log.see(tk.END)
        self.root.update_idletasks()

    def on_statistics_mode_change(self):
        """Handle data mode change (Metabolite <-> Lipid <-> Custom) in Statistics tab."""
        mode = self.statistics_data_mode.get()
        if hasattr(self, 'stats_log'):
            self.stats_log.insert(tk.END, f"\n🔄 Data mode changed to: {mode.upper()}\n")
            self.stats_log.see(tk.END)
        
        # Update mode description label
        if hasattr(self, 'stats_mode_desc_label'):
            if mode == 'metabolite':
                desc = 'Metabolite mode expects Pos_id/Neg_id sheets'
            elif mode == 'lipid':
                desc = 'Lipid mode expects Positive_Lipids/Negative_Lipids sheets'
            else:  # custom
                desc = 'Custom mode for preprocessed/combined data (single sheet)'
            self.stats_mode_desc_label.config(text=desc)
        
        # If switching to custom mode, set preprocessed flag
        if mode == 'custom':
            self.memory_store['is_preprocessed_custom'] = True
        else:
            # Clear custom mode flag when switching to other modes
            if 'is_preprocessed_custom' in self.memory_store:
                del self.memory_store['is_preprocessed_custom']
    
    def _clear_statistics_memory(self):
        """Clear all cached statistics data from memory before importing new file."""
        # Clear memory store
        keys_to_clear = [
            'preprocessed_combined_df',
            'preprocessed_feature_cols',
            'preprocessed_sample_cols',
            'preprocessed_verified_assignments',
            'preprocessed_mode',
            'is_preprocessed_backdoor',
            'is_preprocessed_custom'
        ]
        for key in keys_to_clear:
            if key in self.memory_store:
                del self.memory_store[key]
        
        # Clear normalized dataframes
        self.normalized_combined_df = None
        self.normalized_positive_df = None
        self.normalized_negative_df = None
        
        # Clear verified column assignments
        if hasattr(self, 'verified_pos_sample_cols'):
            self.verified_pos_sample_cols = []
        if hasattr(self, 'verified_neg_sample_cols'):
            self.verified_neg_sample_cols = []
        if hasattr(self, 'verified_pos_lipid_sample_cols'):
            self.verified_pos_lipid_sample_cols = []
        if hasattr(self, 'verified_neg_lipid_sample_cols'):
            self.verified_neg_lipid_sample_cols = []
        
        # Clear group assignments
        if hasattr(self, 'sample_group_vars'):
            for var in self.sample_group_vars.values():
                var.set('')
            self.sample_group_vars.clear()
        
        # Clear normality test results
        if hasattr(self, 'normality_test_results'):
            self.normality_test_results.clear()
        if hasattr(self, 'normality_test_targets'):
            self.normality_test_targets.clear()
    
    # verify_statistics_columns method removed - automated detection handles all column mapping

    def import_statistics_excel(self):
        """Load Excel file with Pos/Neg sheets for Statistics tab based on selected data mode."""
        try:
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            
            # Set dialog title based on mode
            if mode == 'custom':
                title = 'Select Preprocessed Excel (Custom Mode)'
            elif mode == 'lipid':
                title = 'Select Lipid Data Excel'
            else:
                title = 'Select ID Annotated Excel'
            
            path = filedialog.askopenfilename(title=title, filetypes=[('Excel Files','*.xlsx *.xls')])
            if not path:
                return
            if not hasattr(self, 'stats_log'):
                messagebox.showerror('Error', 'Statistics log not initialized.')
                return
            
            # CRITICAL: Clear all previous data from memory before importing new file
            self._clear_statistics_memory()
            
            self.stats_log.insert(tk.END, f"\n===== Importing {mode.upper()} Excel =====\n{os.path.basename(path)}\n")
            self.stats_log.see(tk.END)
            
            # Show progress bar
            self.show_stats_progress(f"Importing {mode} Excel file...")
            
            xl = pd.ExcelFile(path)
            
            # Handle Custom mode directly - import preprocessed/combined data
            if mode == 'custom':
                self.stats_log.insert(tk.END, f'🔄 Importing custom preprocessed data...\n')
                self._import_custom_combined_sheet(xl, path)
                return
            
            pos_sheet = None
            neg_sheet = None
            pos_class_sheet = None
            neg_class_sheet = None
            
            # Define sheet names based on data mode
            if mode == 'lipid':
                # For lipid mode: look for lipid and class sheets
                for candidate in ['Positive_Lipids', 'Positive_Lipid', 'Pos_Lipids', 'Pos_Lipid']:
                    if candidate in xl.sheet_names:
                        pos_sheet = candidate
                        break
                for candidate in ['Negative_Lipids', 'Negative_Lipid', 'Neg_Lipids', 'Neg_Lipid']:
                    if candidate in xl.sheet_names:
                        neg_sheet = candidate
                        break
                if self.use_lipid_class_sheet_backend:
                    for candidate in ['Positive_Lipid_Class', 'Pos_Lipid_Class', 'Positive_Class']:
                        if candidate in xl.sheet_names:
                            pos_class_sheet = candidate
                            break
                    for candidate in ['Negative_Lipid_Class', 'Neg_Lipid_Class', 'Negative_Class']:
                        if candidate in xl.sheet_names:
                            neg_class_sheet = candidate
                            break
                else:
                    self.stats_log.insert(tk.END, 'ℹ️  Lipid class sheets are disabled by backend setting and will be derived after normalization.\n')
            else:
                # For metabolite mode: look for ID annotated sheets
                for candidate in ['Pos_id','Positive','Pos','POS']:
                    if candidate in xl.sheet_names:
                        pos_sheet = candidate
                        break
                for candidate in ['Neg_id','Negative','Neg','NEG']:
                    if candidate in xl.sheet_names:
                        neg_sheet = candidate
                        break

            if not pos_sheet and not neg_sheet:
                # BACKDOOR: If no expected sheets, check if this is a preprocessed combined sheet
                sheet_desc = 'lipid sheets' if mode == 'lipid' else 'Pos_id or Neg_id sheets'
                self.stats_log.insert(tk.END, f'⚠️  No expected {sheet_desc} found.\n')
                self.stats_log.insert(tk.END, f'🔍 Checking for preprocessed combined data...\n')
                
                # Try to use first sheet as preprocessed combined data
                if len(xl.sheet_names) > 0:
                    first_sheet = xl.sheet_names[0]
                    response = messagebox.askyesno(
                        'Preprocessed Data Detected?',
                        f'No standard {sheet_desc} found.\n\n'
                        f'Found sheet: "{first_sheet}"\n\n'
                        f'This appears to be preprocessed/combined data.\n'
                        f'Would you like to import it as a single combined sheet?\n\n'
                        f'You will need to:\n'
                        f'1. Select which columns are feature IDs\n'
                        f'2. Apply normalization\n'
                        f'3. Proceed to statistics (no merging needed)',
                        icon='question'
                    )
                    
                    if response:
                        # User wants to use backdoor - import as combined preprocessed data
                        self._import_preprocessed_combined_sheet(xl, first_sheet, path, mode)
                        return
                    else:
                        self.stats_log.insert(tk.END, f'❌ Import cancelled by user.\n')
                        self.hide_stats_progress()
                        return
                else:
                    messagebox.showwarning('No Sheets Found', f'Could not find expected {sheet_desc} in the selected file.')
                    self.stats_log.insert(tk.END, f'❌ No sheets available in file.\n')
                    self.hide_stats_progress()
                    return
            
            # Load main polarity sheets
            pos_df = xl.parse(pos_sheet) if pos_sheet else None
            neg_df = xl.parse(neg_sheet) if neg_sheet else None
            pos_class_df = xl.parse(pos_class_sheet) if pos_class_sheet else None
            neg_class_df = xl.parse(neg_class_sheet) if neg_class_sheet else None
            
            # Store in memory based on mode
            if mode == 'lipid':
                if pos_df is not None:
                    self.memory_store['pos_lipid_df'] = pos_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {pos_sheet}: {len(pos_df)} rows, {len(pos_df.columns)} columns.\n')
                if neg_df is not None:
                    self.memory_store['neg_lipid_df'] = neg_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {neg_sheet}: {len(neg_df)} rows, {len(neg_df.columns)} columns.\n')
                if pos_class_df is not None:
                    self.memory_store['pos_lipid_class_df'] = pos_class_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {pos_class_sheet}: {len(pos_class_df)} rows, {len(pos_class_df.columns)} columns.\n')
                if neg_class_df is not None:
                    self.memory_store['neg_lipid_class_df'] = neg_class_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {neg_class_sheet}: {len(neg_class_df)} rows, {len(neg_class_df.columns)} columns.\n')
            else:
                if pos_df is not None:
                    self.memory_store['pos_id_df'] = pos_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {pos_sheet}: {len(pos_df)} rows, {len(pos_df.columns)} columns.\n')
                if neg_df is not None:
                    self.memory_store['neg_id_df'] = neg_df
                    self.stats_log.insert(tk.END, f'✅ Loaded {neg_sheet}: {len(neg_df)} rows, {len(neg_df.columns)} columns.\n')
            
            # Attempt automatic sample column detection from union
            union_sample_cols = []
            seen = set()
            
            # Define feature columns based on mode
            if mode == 'lipid':
                # Process lipid sheets
                pos_sample_cols = []
                neg_sample_cols = []
                for label, df in [('Positive Lipids', pos_df), ('Negative Lipids', neg_df)]:
                    if df is None:
                        self.stats_log.insert(tk.END, f'{label}: DataFrame is None, skipping.\n')
                        continue
                    try:
                        # Use robust feature detection
                        sample_cols = []
                        excluded_feature_cols = []
                        excluded_metadata_cols = []
                        excluded_non_numeric = []
                        
                        for col in df.columns:
                            # Skip non-string columns
                            if not isinstance(col, str):
                                continue
                            # treat as feature if normalized name matches canonical lipid features
                            if self._is_lipid_feature_col(col):
                                excluded_feature_cols.append(col)
                                continue
                            if is_statistics_metadata_col(col):
                                excluded_metadata_cols.append(col)
                                continue
                            # numeric columns not identified as features are treated as sample intensity columns
                            if pd.api.types.is_numeric_dtype(df[col]):
                                sample_cols.append(col)
                            else:
                                excluded_non_numeric.append(col)
                        
                        # Store per-sheet sample columns
                        if label == 'Positive Lipids':
                            pos_sample_cols = sample_cols
                        else:
                            neg_sample_cols = sample_cols
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        import traceback
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
                        self.stats_log.insert(tk.END, f'Traceback: {traceback.format_exc()}\n')
                # Store detected sample columns for each sheet
                self.detected_pos_lipid_sample_cols = pos_sample_cols
                self.detected_neg_lipid_sample_cols = neg_sample_cols
                
                # Process class sheets
                for label, df in [('Positive Class', pos_class_df), ('Negative Class', neg_class_df)]:
                    if df is None:
                        continue
                    try:
                        sample_cols = []
                        for col in df.columns:
                            if col is None:
                                continue
                            # Use normalized matching for Class column
                            if self._normalize_col(col) == 'class':
                                continue
                            if is_statistics_metadata_col(col):
                                continue
                            if pd.api.types.is_numeric_dtype(df[col]):
                                sample_cols.append(col)
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
            else:
                # Use existing metabolite detection
                from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
                pos_sample_cols = []
                neg_sample_cols = []
                for label, df in [('Positive', pos_df), ('Negative', neg_df)]:
                    if df is None:
                        continue
                    try:
                        feature_cols, sample_cols = detect_feature_and_sample_columns(df)
                        # Store per-sheet sample columns
                        if label == 'Positive':
                            pos_sample_cols = sample_cols
                        else:
                            neg_sample_cols = sample_cols
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
                # Store detected sample columns for each sheet
                self.detected_pos_sample_cols = pos_sample_cols
                self.detected_neg_sample_cols = neg_sample_cols
            
            if union_sample_cols:
                # Store sample columns for later use
                self.detected_sample_cols = union_sample_cols
                # Ensure listbox exists (hidden, but needed for compatibility)
                if not hasattr(self, 'sample_cols_list'):
                    hidden_frame = tk.Frame(self.root, height=1, width=1)
                    hidden_frame.pack_forget()
                    self.sample_cols_list = tk.Listbox(hidden_frame)
                self.sample_cols_list.delete(0, tk.END)
                for c in union_sample_cols:
                    self.sample_cols_list.insert(tk.END, c)
                self.populate_sample_assignments(union_sample_cols)
                self.stats_log.insert(tk.END, f'Union sample columns loaded: {len(union_sample_cols)} columns.\n')
            else:
                self.stats_log.insert(tk.END, 'No sample columns auto-detected from imported file.\n')
            
            self.stats_log.insert(tk.END, 'Ready for normalization. Configure groups and click "Normalization & Test Normality".\n')
            self.stats_log.see(tk.END)
            
            # Hide progress bar on success
            self.hide_stats_progress()
            
        except Exception as e:
            if hasattr(self, 'stats_log'):
                self.stats_log.insert(tk.END, f'❌ Import failed: {e}\n')
                self.stats_log.see(tk.END)
            self.hide_stats_progress()
            messagebox.showerror('Import Failed', str(e))
    
    def import_statistics_id_excel(self):
        """Legacy wrapper for backward compatibility - calls new import_statistics_excel."""
        self.import_statistics_excel()
    
    def verify_statistics_columns(self):
        """Verify and assign columns for statistics analysis using unified dialog"""
        import threading
        from gui.shared.column_assignment import show_column_assignment_dialog
        
        def _load_and_verify():
            """Background thread worker for loading files and showing dialogs"""
            try:
                mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
                
                # Determine the correct tab_type based on mode
                if mode == 'custom':
                    tab_type = 'statistics_metabolite'  # Use metabolite-style dialog for custom
                else:
                    tab_type = 'statistics_metabolite' if mode == 'metabolite' else 'statistics_lipid'
                
                # 📊 CUSTOM MODE or BACKDOOR MODE: Handle preprocessed combined data
                is_custom_mode = self.memory_store.get('is_preprocessed_custom', False) or self.memory_store.get('is_preprocessed_backdoor', False)
                
                if is_custom_mode:
                    mode_label = "CUSTOM" if self.memory_store.get('is_preprocessed_custom', False) else "PREPROCESSED"
                    self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                        f"\n📊 {mode_label} MODE: Verifying preprocessed combined data...\n"))
                    self.root.after(0, lambda: self.stats_log.see(tk.END))
                    self.root.after(0, lambda: self.show_stats_progress("Verifying preprocessed columns..."))
                    
                    # Get preprocessed data
                    combined_df = self.memory_store.get('preprocessed_combined_df')
                    feature_cols = self.memory_store.get('preprocessed_feature_cols', [])
                    sample_cols = self.memory_store.get('preprocessed_sample_cols', [])
                    
                    if combined_df is None:
                        self.root.after(0, lambda: messagebox.showerror(
                            "No Preprocessed Data",
                            "Preprocessed data not found. Please import Excel data again."
                        ))
                        self.root.after(0, lambda: self.hide_stats_progress())
                        return
                    
                    # Show column assignment dialog for preprocessed data
                    result = show_column_assignment_dialog(
                        parent=self.root,
                        df=combined_df,
                        tab_type=tab_type,
                        auto_calculate=False,
                        dialog_title=f"Custom Mode - Verify Columns",
                        detected_sample_cols=sample_cols,
                    )
                    
                    if not result:
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "❌ Verification cancelled\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        self.root.after(0, lambda: self.hide_stats_progress())
                        return
                    
                    # Update preprocessed data with verified columns
                    verified_sample_cols = result.get('sample_cols', [])
                    verified_assignments = result.get('assignments', {})
                    
                    # CRITICAL: Use the original feature column selected during import
                    # The dialog might return standardized keys, but we need the actual column name
                    feature_cols = self.memory_store.get('preprocessed_feature_cols', [])
                    feature_id_key = 'Feature ID'  # Always use 'Feature ID' for custom mode
                    
                    # If user selected a feature column during import, use that as Feature ID
                    if feature_cols and len(feature_cols) > 0:
                        actual_feature_id_col = feature_cols[0]  # Use first feature column as ID
                        verified_assignments[feature_id_key] = actual_feature_id_col

                    # Critical safety: remove any feature-assigned columns from group-assignable samples
                    forbidden_cols = self._extract_feature_assigned_columns(verified_assignments)
                    verified_sample_cols = [c for c in verified_sample_cols if c not in forbidden_cols]
                    
                    # Store verified data back to memory
                    self.memory_store['preprocessed_sample_cols'] = verified_sample_cols
                    self.memory_store['preprocessed_verified_assignments'] = verified_assignments
                    
                    # Populate sample assignments for Configure Groups
                    if verified_sample_cols:
                        self.root.after(0, lambda cols=verified_sample_cols: self.populate_sample_assignments(cols))
                        if hasattr(self, 'sample_cols_list'):
                            self.root.after(0, lambda: self.sample_cols_list.delete(0, tk.END))
                            for col in verified_sample_cols:
                                self.root.after(0, lambda c=col: self.sample_cols_list.insert(tk.END, c))
                    
                    # Get Feature ID from assignments
                    feature_id_col = verified_assignments.get(feature_id_key, 'N/A')
                    
                    # Log success
                    self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                        f"\n✅ CUSTOM MODE: Columns verified!\n"
                        f"• Sample columns: {len(verified_sample_cols)}\n"
                        f"• Feature ID: {feature_id_col}\n\n"))
                    self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                        "Next steps:\n"
                        "1. Configure sample groups\n"
                        "2. Select normalization method\n"
                        "3. Click 'Run Normalization & Test Normality'\n\n"))
                    self.root.after(0, lambda: self.stats_log.see(tk.END))
                    self.root.after(0, lambda: self.hide_stats_progress())
                    
                    # Enable Configure Groups and Run Normalization buttons
                    self.root.after(0, lambda: self.configure_groups_btn.configure(state='normal'))
                    self.root.after(0, lambda: self.run_norm_btn.configure(state='normal'))
                    
                    self.root.after(0, lambda: messagebox.showinfo(
                        "Custom Mode - Columns Verified",
                        f"✅ Columns verified successfully!\n\n"
                        f"Verified sample columns: {len(verified_sample_cols)}\n"
                        f"Feature ID: {feature_id_col}\n\n"
                        f"You can now:\n"
                        f"• Configure sample groups\n"
                        f"• Run normalization and analysis"
                    ))
                    return
                
                # Normal mode: Check if data has been imported
                if mode == 'custom':
                    # Custom mode should have triggered above, but handle edge case
                    self.root.after(0, lambda: messagebox.showerror(
                        "No Custom Data",
                        "Custom mode selected but no preprocessed data found.\nPlease import Excel data first."
                    ))
                    return
                    
                if mode == 'lipid':
                    pos_df = self.memory_store.get('pos_lipid_df')
                    neg_df = self.memory_store.get('neg_lipid_df')
                else:
                    pos_df = self.memory_store.get('pos_id_df')
                    neg_df = self.memory_store.get('neg_id_df')
                
                if pos_df is None and neg_df is None:
                    self.root.after(0, lambda: messagebox.showerror(
                        "No Data Loaded",
                        "Please import Excel data first using the 'Import Excel Data' button."
                    ))
                    return
                
                # Update status
                self.root.after(0, lambda: self.show_stats_progress("Verifying columns..."))
                self.root.after(0, lambda: self.stats_log.insert(tk.END, f"\n🔍 Verifying {mode} columns...\n"))
                self.root.after(0, lambda: self.stats_log.see(tk.END))
                
                # Verify POSITIVE sheet if it exists
                if pos_df is not None:
                    self.root.after(0, lambda: self.stats_log.insert(tk.END, "📂 Verifying Positive sheet...\n"))
                    self.root.after(0, lambda: self.stats_log.see(tk.END))
                    
                    # Get detected sample columns for positive sheet
                    if mode == 'lipid':
                        detected_samples = getattr(self, 'detected_pos_lipid_sample_cols', [])
                    else:
                        detected_samples = getattr(self, 'detected_pos_sample_cols', [])
                    
                    pos_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=pos_df,
                        tab_type=tab_type,
                        auto_calculate=False,
                        dialog_title=f"Positive {mode.capitalize()} - Column Assignment",
                        detected_sample_cols=detected_samples,
                        allow_skip=True,
                    )
                    
                    # Handle Cancel vs Skip vs Confirm
                    if pos_result is None:
                        # User cancelled the entire verification
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "❌ Verification cancelled\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        self.root.after(0, lambda: self.hide_stats_progress())
                        return
                    if pos_result.get('skipped'):
                        # Skip using positive sheet entirely
                        if mode == 'lipid':
                            if 'pos_lipid_df' in self.memory_store:
                                self.memory_store['pos_lipid_df'] = None
                            self.verified_pos_lipid_assignments = {}
                            self.verified_pos_lipid_sample_cols = []
                        else:
                            if 'pos_id_df' in self.memory_store:
                                self.memory_store['pos_id_df'] = None
                            self.verified_pos_assignments = {}
                            self.verified_pos_sample_cols = []
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "⏭ Skipped Positive sheet (will not be used).\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                    else:
                        # Store positive assignments
                        if mode == 'lipid':
                            self.verified_pos_lipid_assignments = pos_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_pos_lipid_assignments)
                            self.verified_pos_lipid_sample_cols = [
                                c for c in pos_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                        else:
                            self.verified_pos_assignments = pos_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_pos_assignments)
                            self.verified_pos_sample_cols = [
                                c for c in pos_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                        
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                            f"✅ Positive: {len(pos_result.get('sample_cols', []))} sample columns verified\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                
                # Verify NEGATIVE sheet if it exists
                if neg_df is not None:
                    self.root.after(0, lambda: self.stats_log.insert(tk.END, "📂 Verifying Negative sheet...\n"))
                    self.root.after(0, lambda: self.stats_log.see(tk.END))
                    
                    # Get detected sample columns for negative sheet
                    if mode == 'lipid':
                        detected_samples = getattr(self, 'detected_neg_lipid_sample_cols', [])
                    else:
                        detected_samples = getattr(self, 'detected_neg_sample_cols', [])
                    
                    neg_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=neg_df,
                        tab_type=tab_type,
                        auto_calculate=False,
                        dialog_title=f"Negative {mode.capitalize()} - Column Assignment",
                        detected_sample_cols=detected_samples,
                        allow_skip=True,
                    )
                    
                    if neg_result is None:
                        # User cancelled the entire verification
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "❌ Verification cancelled\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        self.root.after(0, lambda: self.hide_stats_progress())
                        return
                    if neg_result.get('skipped'):
                        # Skip using negative sheet entirely
                        if mode == 'lipid':
                            if 'neg_lipid_df' in self.memory_store:
                                self.memory_store['neg_lipid_df'] = None
                            self.verified_neg_lipid_assignments = {}
                            self.verified_neg_lipid_sample_cols = []
                        else:
                            if 'neg_id_df' in self.memory_store:
                                self.memory_store['neg_id_df'] = None
                            self.verified_neg_assignments = {}
                            self.verified_neg_sample_cols = []
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "⏭ Skipped Negative sheet (will not be used).\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                    else:
                        # Store negative assignments
                        if mode == 'lipid':
                            self.verified_neg_lipid_assignments = neg_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_neg_lipid_assignments)
                            self.verified_neg_lipid_sample_cols = [
                                c for c in neg_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                        else:
                            self.verified_neg_assignments = neg_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_neg_assignments)
                            self.verified_neg_sample_cols = [
                                c for c in neg_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                        
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                            f"✅ Negative: {len(neg_result.get('sample_cols', []))} sample columns verified\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                
                # Verify CLASS sheets if they exist (lipid mode only)
                if mode == 'lipid':
                    pos_class_df = self.memory_store.get('pos_lipid_class_df')
                    neg_class_df = self.memory_store.get('neg_lipid_class_df')
                    
                    # Verify POSITIVE_CLASS sheet if it exists
                    if pos_class_df is not None and not pos_class_df.empty:
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "📂 Verifying Positive_Lipid_Class sheet...\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        
                        # Use the SAME detected sample columns from main positive sheet
                        # Class sheets have identical sample columns (aggregated from main sheet)
                        class_sample_cols = getattr(self, 'detected_pos_lipid_sample_cols', [])
                        if not class_sample_cols:
                            # Fallback: detect from verified columns if available
                            class_sample_cols = getattr(self, 'verified_pos_lipid_sample_cols', [])
                        
                        self.root.after(0, lambda cols=class_sample_cols: self.stats_log.insert(tk.END, 
                            f"   ℹ️  Inheriting {len(cols)} sample columns from Positive sheet\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        
                        pos_class_result = show_column_assignment_dialog(
                            parent=self.root,
                            df=pos_class_df,
                            tab_type='statistics_lipid',  # Same as regular lipid, but Class replaces LipidID
                            auto_calculate=False,
                            dialog_title="Positive Lipid Class - Column Verification",
                            detected_sample_cols=class_sample_cols,
                            allow_skip=True,
                        )
                        
                        if pos_class_result is None:
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, "❌ Class sheet verification cancelled\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                            self.root.after(0, lambda: self.hide_stats_progress())
                            return
                        if pos_class_result.get('skipped'):
                            # Skip class sheet
                            if 'pos_lipid_class_df' in self.memory_store:
                                self.memory_store['pos_lipid_class_df'] = None
                            self.verified_pos_lipid_class_assignments = {}
                            self.verified_pos_lipid_class_sample_cols = []
                            self.pos_class_inherits_grouping = False
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, "⏭ Skipped Positive_Lipid_Class sheet (will not be used).\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                        else:
                            # Store positive class assignments
                            self.verified_pos_lipid_class_assignments = pos_class_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_pos_lipid_class_assignments)
                            self.verified_pos_lipid_class_sample_cols = [
                                c for c in pos_class_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                            # Mark that class uses same grouping as parent sheet
                            self.pos_class_inherits_grouping = True
                            
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                                f"✅ Positive_Class: {len(pos_class_result.get('sample_cols', []))} sample columns verified\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                    
                    # Verify NEGATIVE_CLASS sheet if it exists
                    if neg_class_df is not None and not neg_class_df.empty:
                        self.root.after(0, lambda: self.stats_log.insert(tk.END, "📂 Verifying Negative_Lipid_Class sheet...\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        
                        # Use the SAME detected sample columns from main negative sheet
                        # Class sheets have identical sample columns (aggregated from main sheet)
                        class_sample_cols = getattr(self, 'detected_neg_lipid_sample_cols', [])
                        if not class_sample_cols:
                            # Fallback: detect from verified columns if available
                            class_sample_cols = getattr(self, 'verified_neg_lipid_sample_cols', [])
                        
                        self.root.after(0, lambda cols=class_sample_cols: self.stats_log.insert(tk.END, 
                            f"   ℹ️  Inheriting {len(cols)} sample columns from Negative sheet\n"))
                        self.root.after(0, lambda: self.stats_log.see(tk.END))
                        
                        neg_class_result = show_column_assignment_dialog(
                            parent=self.root,
                            df=neg_class_df,
                            tab_type='statistics_lipid',  # Same as regular lipid, but Class replaces LipidID
                            auto_calculate=False,
                            dialog_title="Negative Lipid Class - Column Verification",
                            detected_sample_cols=class_sample_cols,
                            allow_skip=True,
                        )
                        
                        if neg_class_result is None:
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, "❌ Class sheet verification cancelled\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                            self.root.after(0, lambda: self.hide_stats_progress())
                            return
                        if neg_class_result.get('skipped'):
                            # Skip class sheet
                            if 'neg_lipid_class_df' in self.memory_store:
                                self.memory_store['neg_lipid_class_df'] = None
                            self.verified_neg_lipid_class_assignments = {}
                            self.verified_neg_lipid_class_sample_cols = []
                            self.neg_class_inherits_grouping = False
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, "⏭ Skipped Negative_Lipid_Class sheet (will not be used).\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                        else:
                            # Store negative class assignments
                            self.verified_neg_lipid_class_assignments = neg_class_result['assignments']
                            forbidden_cols = self._extract_feature_assigned_columns(self.verified_neg_lipid_class_assignments)
                            self.verified_neg_lipid_class_sample_cols = [
                                c for c in neg_class_result.get('sample_cols', []) if c not in forbidden_cols
                            ]
                            # Mark that class uses same grouping as parent sheet
                            self.neg_class_inherits_grouping = True
                            
                            self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                                f"✅ Negative_Class: {len(neg_class_result.get('sample_cols', []))} sample columns verified\n"))
                            self.root.after(0, lambda: self.stats_log.see(tk.END))
                
                # Show summary
                self.root.after(0, lambda: self.stats_log.insert(tk.END, "\n✅ Column verification complete!\n"))
                self.root.after(0, lambda: self.stats_log.insert(tk.END, 
                    "Verified columns will be used during normalization.\n\n"))
                self.root.after(0, lambda: self.stats_log.see(tk.END))
                self.root.after(0, lambda: self.hide_stats_progress())
                
                # Re-populate sample assignments with ONLY verified columns
                # This ensures Configure Groups only sees the columns user verified
                verified_union_cols = []
                if mode == 'lipid':
                    if hasattr(self, 'verified_pos_lipid_sample_cols') and self.verified_pos_lipid_sample_cols:
                        verified_union_cols.extend(self.verified_pos_lipid_sample_cols)
                    if hasattr(self, 'verified_neg_lipid_sample_cols') and self.verified_neg_lipid_sample_cols:
                        verified_union_cols.extend(self.verified_neg_lipid_sample_cols)
                else:
                    if hasattr(self, 'verified_pos_sample_cols') and self.verified_pos_sample_cols:
                        verified_union_cols.extend(self.verified_pos_sample_cols)
                    if hasattr(self, 'verified_neg_sample_cols') and self.verified_neg_sample_cols:
                        verified_union_cols.extend(self.verified_neg_sample_cols)
                
                # Remove duplicates while preserving order
                seen = set()
                verified_union_cols = [col for col in verified_union_cols if not (col in seen or seen.add(col))]
                
                if verified_union_cols:
                    self.root.after(0, lambda cols=verified_union_cols: self.populate_sample_assignments(cols))
                    self.root.after(0, lambda cols=verified_union_cols: self.stats_log.insert(tk.END, 
                        f"✅ Updated group configuration with {len(cols)} verified sample columns.\n\n"))
                    self.root.after(0, lambda: self.stats_log.see(tk.END))
                
                # Enable Configure Groups and Run Normalization buttons after successful verification
                self.root.after(0, lambda: self.configure_groups_btn.configure(state='normal'))
                self.root.after(0, lambda: self.run_norm_btn.configure(state='normal'))
                
                self.root.after(0, lambda: messagebox.showinfo(
                    "Verification Complete",
                    "Column verification completed successfully.\n\n"
                    "Next steps:\n"
                    "• Configure Groups (assign samples to groups)\n"
                    "• Run Normalization & Test Normality\n"
                    "• Run Statistical Tests"
                ))
                
            except Exception as e:
                logger.error(f"Error verifying statistics columns: {e}")
                import traceback
                traceback.print_exc()
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to verify columns: {str(e)}"))
                self.root.after(0, lambda: self.stats_log.insert(tk.END, f"❌ Error: {str(e)}\n"))
                self.root.after(0, lambda: self.stats_log.see(tk.END))
                self.root.after(0, lambda: self.hide_stats_progress())
        
        # Start verification in background thread
        verification_thread = threading.Thread(target=_load_and_verify, daemon=True)
        verification_thread.start()
    
    def _import_custom_combined_sheet(self, xl: pd.ExcelFile, file_path: str):
        """
        CUSTOM MODE: Import a single preprocessed/combined sheet.
        User selects which sheet to use and which columns are feature columns.
        The rest of the numeric columns are treated as sample data.
        """
        try:
            self.stats_log.insert(tk.END, f'\n📊 CUSTOM MODE: Importing preprocessed data...\n')
            
            # If multiple sheets, let user choose
            if len(xl.sheet_names) == 0:
                messagebox.showerror("No Sheets", "No sheets found in Excel file")
                self.hide_stats_progress()
                return
            
            sheet_name = str(xl.sheet_names[0])
            if len(xl.sheet_names) > 1:
                from tkinter import simpledialog
                sheets_list = "\n".join([str(s) for s in xl.sheet_names[:10]])
                if len(xl.sheet_names) > 10:
                    sheets_list += f"\n... and {len(xl.sheet_names) - 10} more"
                sheet_name = simpledialog.askstring(
                    "Select Sheet",
                    f"Multiple sheets found:\n{sheets_list}\n\nEnter sheet name:",
                    initialvalue=sheet_name
                )
                if not sheet_name or sheet_name not in xl.sheet_names:
                    self.stats_log.insert(tk.END, "❌ Invalid sheet selected\n")
                    self.hide_stats_progress()
                    return
            
            self.stats_log.insert(tk.END, f'📋 Sheet: "{sheet_name}"\n')
            
            # Load the sheet
            combined_df = xl.parse(sheet_name)
            self.stats_log.insert(tk.END, f'✅ Loaded {len(combined_df)} rows, {len(combined_df.columns)} columns.\n')
            
            # Show column selector dialog for user to pick feature columns
            feature_cols = self._select_feature_columns_dialog(combined_df, 'custom')
            
            if feature_cols is None:
                self.stats_log.insert(tk.END, f'❌ Feature selection cancelled.\n')
                self.hide_stats_progress()
                return
            
            self.stats_log.insert(tk.END, f'📋 Feature columns selected: {len(feature_cols)}\n')
            for fc in feature_cols:
                self.stats_log.insert(tk.END, f'   - {fc}\n')
            
            # Identify sample columns (numeric columns that aren't features)
            sample_cols = []
            for col in combined_df.columns:
                if col in feature_cols:
                    continue
                # Skip known metadata columns
                if is_statistics_metadata_col(col):
                    continue
                # Skip lipid feature columns
                try:
                    if self._is_lipid_feature_col(col):
                        continue
                except Exception:
                    pass
                # Include if numeric
                if pd.api.types.is_numeric_dtype(combined_df[col]):
                    sample_cols.append(col)
            
            self.stats_log.insert(tk.END, f'📊 Sample columns detected: {len(sample_cols)}\n')
            
            if not sample_cols:
                messagebox.showerror('No Sample Columns', 
                    'No numeric sample columns found after removing feature columns.\n'
                    'Please check your data.')
                self.stats_log.insert(tk.END, f'❌ No sample columns available.\n')
                self.hide_stats_progress()
                return
            
            # Store main feature ID column (first selected feature column)
            main_feature_col = feature_cols[0] if feature_cols else None
            
            # Store the combined dataframe in memory with custom mode flags
            self.memory_store['preprocessed_combined_df'] = combined_df
            self.memory_store['preprocessed_feature_cols'] = feature_cols
            self.memory_store['preprocessed_sample_cols'] = sample_cols
            self.memory_store['preprocessed_mode'] = 'custom'
            self.memory_store['is_preprocessed_custom'] = True  # Custom mode flag
            self.memory_store['is_preprocessed_backdoor'] = True  # Keep for backward compatibility
            
            # CRITICAL: Pre-store verified assignments with the feature ID
            # This ensures the Feature ID is available when normalization runs
            verified_assignments = {'Feature ID': main_feature_col}
            self.memory_store['preprocessed_verified_assignments'] = verified_assignments
            
            self.stats_log.insert(tk.END, f'✅ Preprocessed data loaded into memory.\n')
            self.stats_log.insert(tk.END, f'✅ Feature ID column: {main_feature_col}\n')
            
            # Populate sample columns list
            if hasattr(self, 'sample_cols_list'):
                self.sample_cols_list.delete(0, tk.END)
                for c in sample_cols:
                    self.sample_cols_list.insert(tk.END, c)
                self.populate_sample_assignments(sample_cols)
            
            # Show instructions
            self.stats_log.insert(tk.END, f'\n✨ CUSTOM MODE READY ✨\n')
            self.stats_log.insert(tk.END, f'✅ Detected {len(sample_cols)} sample columns. Use "Auto-Assign by Pattern" button to configure groups.\n')
            self.stats_log.insert(tk.END, f'\nNext steps:\n')
            self.stats_log.insert(tk.END, f'1. Click "🔍 Verify Columns" to confirm column assignments\n')
            self.stats_log.insert(tk.END, f'2. Configure sample groups using "⚙️ Configure Groups"\n')
            self.stats_log.insert(tk.END, f'3. Select normalization method\n')
            self.stats_log.insert(tk.END, f'4. Click "Run Normalization & Test Normality"\n')
            self.stats_log.see(tk.END)
            
            self.hide_stats_progress()
            
            # Enable Configure Groups button since we have sample columns
            if hasattr(self, 'configure_groups_btn'):
                self.configure_groups_btn.configure(state='normal')
            
            messagebox.showinfo('Custom Data Loaded',
                f'Successfully loaded preprocessed data!\n\n'
                f'Rows: {len(combined_df)}\n'
                f'Feature ID: {main_feature_col}\n'
                f'Sample columns: {len(sample_cols)}\n\n'
                f'Next: Click "🔍 Verify Columns" then configure groups.')
            
        except Exception as e:
            self.stats_log.insert(tk.END, f'❌ Custom import failed: {e}\n')
            self.stats_log.see(tk.END)
            self.hide_stats_progress()
            messagebox.showerror('Import Failed', f'Failed to import custom data:\n{e}')

    def _import_preprocessed_combined_sheet(self, xl, sheet_name, file_path, mode):
        """
        BACKDOOR: Import a single preprocessed/combined sheet.
        Allows users to import data that's already combined (not split into Pos/Neg).
        User selects feature columns, rest are treated as numeric sample columns.
        Data is normalized and ready for statistics without merging.
        """
        try:
            self.stats_log.insert(tk.END, f'\n🔓 BACKDOOR MODE: Importing preprocessed combined data...\n')
            self.stats_log.insert(tk.END, f'📊 Sheet: "{sheet_name}"\n')
            
            # Load the sheet
            combined_df = xl.parse(sheet_name)
            self.stats_log.insert(tk.END, f'✅ Loaded {len(combined_df)} rows, {len(combined_df.columns)} columns.\n')
            
            # Show column selector dialog for user to pick feature columns
            feature_cols = self._select_feature_columns_dialog(combined_df, mode)
            
            if feature_cols is None:
                self.stats_log.insert(tk.END, f'❌ Feature selection cancelled.\n')
                self.hide_stats_progress()
                return
            
            self.stats_log.insert(tk.END, f'📋 Feature columns selected: {len(feature_cols)}\n')
            for fc in feature_cols:
                self.stats_log.insert(tk.END, f'   - {fc}\n')
            
            # Identify sample columns (numeric columns that aren't features)
            sample_cols = []
            for col in combined_df.columns:
                # Exclude lipid feature canonical columns from being considered numeric sample columns
                if col not in feature_cols:
                    try:
                        if self._is_lipid_feature_col(col):
                            continue
                    except Exception:
                        pass
                if is_statistics_metadata_col(col):
                    continue
                if col not in feature_cols and pd.api.types.is_numeric_dtype(combined_df[col]):
                    sample_cols.append(col)
            
            self.stats_log.insert(tk.END, f'📊 Sample columns detected: {len(sample_cols)}\n')
            
            if not sample_cols:
                messagebox.showerror('No Sample Columns', 
                    'No numeric sample columns found after removing feature columns.\n'
                    'Please check your data.')
                self.stats_log.insert(tk.END, f'❌ No sample columns available.\n')
                self.hide_stats_progress()
                return
            
            # Store the combined dataframe in memory with a special flag
            self.memory_store['preprocessed_combined_df'] = combined_df
            self.memory_store['preprocessed_feature_cols'] = feature_cols
            self.memory_store['preprocessed_sample_cols'] = sample_cols
            self.memory_store['preprocessed_mode'] = mode
            self.memory_store['is_preprocessed_backdoor'] = True  # Flag for later processing
            
            self.stats_log.insert(tk.END, f'✅ Preprocessed data loaded into memory.\n')
            
            # Populate sample columns list
            if hasattr(self, 'sample_cols_list'):
                self.sample_cols_list.delete(0, tk.END)
                for c in sample_cols:
                    self.sample_cols_list.insert(tk.END, c)
                self.populate_sample_assignments(sample_cols)
            
            # Show instructions
            self.stats_log.insert(tk.END, f'\n✨ BACKDOOR MODE READY ✨\n')
            self.stats_log.insert(tk.END, f'Next steps:\n')
            self.stats_log.insert(tk.END, f'1. Configure sample groups in "Assign Samples to Groups"\n')
            self.stats_log.insert(tk.END, f'2. Select normalization method\n')
            self.stats_log.insert(tk.END, f'3. Click "Normalization & Test Normality" (pipeline + normality)\n')
            self.stats_log.insert(tk.END, f'4. Data will be normalized and ready for statistics!\n')
            self.stats_log.insert(tk.END, f'\n⚠️  Note: Pos/Neg splitting is bypassed in backdoor mode.\n')
            self.stats_log.see(tk.END)
            
            self.hide_stats_progress()
            
            messagebox.showinfo('Preprocessed Data Loaded',
                f'Successfully loaded preprocessed combined data!\n\n'
                f'Rows: {len(combined_df)}\n'
                f'Feature columns: {len(feature_cols)}\n'
                f'Sample columns: {len(sample_cols)}\n\n'
                f'Configure groups and run normalization.')
            
        except Exception as e:
            self.stats_log.insert(tk.END, f'❌ Backdoor import failed: {e}\n')
            self.stats_log.see(tk.END)
            self.hide_stats_progress()
            messagebox.showerror('Import Failed', f'Failed to import preprocessed data:\n{e}')
    
    def _select_feature_columns_dialog(self, df, mode):
        """
        Show dialog for user to select which columns are feature IDs.
        Returns list of selected column names, or None if cancelled.
        """
        dialog = tk.Toplevel(self.root)
        dialog.title(f"Select Feature Columns - {mode.capitalize()} Mode")
        dialog.geometry("600x500")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Header
        header_text = (
            f"Select which columns contain feature identifiers.\n"
            f"All other numeric columns will be treated as sample data.\n\n"
            f"Common feature columns:\n"
            f"  Metabolite: Name, Compound, metabolite_id, Feature_ID\n"
            f"  Lipid: lipidid, LipidName, class, subclass"
        )
        header = ttk.Label(dialog, text=header_text, justify='left', 
                          font=('TkDefaultFont', 9), wraplength=580)
        header.pack(pady=10, padx=10, anchor='w')
        
        # Listbox with scrollbar for column selection
        list_frame = ttk.LabelFrame(dialog, text=f"Columns in '{df.columns[0] if len(df.columns) > 0 else 'data'}'...", padding=10)
        list_frame.pack(fill='both', expand=True, padx=10, pady=(0,10))
        
        scrollbar = ttk.Scrollbar(list_frame, orient='vertical')
        listbox = tk.Listbox(
            list_frame,
            selectmode='multiple',
            exportselection=False,
            yscrollcommand=scrollbar.set,
            height=15,
            font=('TkDefaultFont', 9)
        )
        scrollbar.config(command=listbox.yview)
        scrollbar.pack(side='right', fill='y')
        listbox.pack(side='left', fill='both', expand=True)
        
        # Populate with columns
        for col in df.columns:
            listbox.insert(tk.END, col)
        
        # Auto-select likely feature columns
        likely_features_metabolite = ['name', 'compound', 'metabolite_id', 'metaboliteid', 'feature_id', 
                                      'featureid', 'id', 'identifier', 'metabolite', 'feature']
        likely_features_lipid = ['lipidid', 'lipid_id', 'lipidname', 'lipid_name', 'class', 
                                'subclass', 'lipid', 'lipidgroup', 'charge', 'calcmz', 'basert', 
                                'obsrt', 'ppmdiff', 'mz', 'retention', 'rt', 'mass']
        likely_features = likely_features_lipid if mode == 'lipid' else likely_features_metabolite
        
        for idx, col in enumerate(df.columns):
            col_lower = str(col).lower().replace('_', '').replace(' ', '')
            if any(lf in col_lower for lf in likely_features):
                listbox.selection_set(idx)
        
        # Info label
        info_label = ttk.Label(dialog, text="Tip: Use Ctrl+Click to select multiple columns", 
                              foreground='gray', font=('TkDefaultFont', 8))
        info_label.pack(pady=(0,10))
        
        # Result holder
        selected_cols = []
        
        def on_ok():
            nonlocal selected_cols
            indices = listbox.curselection()
            if not indices:
                messagebox.showwarning("No Selection", 
                    "Please select at least one feature column.\n"
                    "If all columns are samples, click Cancel and report this use case.")
                return
            selected_cols = [listbox.get(i) for i in indices]
            dialog.destroy()
        
        def on_cancel():
            nonlocal selected_cols
            selected_cols = None
            dialog.destroy()
        
        # Buttons
        button_frame = ttk.Frame(dialog)
        button_frame.pack(fill='x', padx=10, pady=10)
        
        ttk.Button(button_frame, text="Cancel", command=on_cancel).pack(side='right', padx=5)
        ttk.Button(button_frame, text="OK - Use Selected", command=on_ok).pack(side='right', padx=5)
        
        # Select All / Deselect All helpers
        helper_frame = ttk.Frame(button_frame)
        helper_frame.pack(side='left')
        ttk.Button(helper_frame, text="Select All", 
                  command=lambda: listbox.selection_set(0, tk.END)).pack(side='left', padx=2)
        ttk.Button(helper_frame, text="Clear All", 
                  command=lambda: listbox.selection_clear(0, tk.END)).pack(side='left', padx=2)
        
        dialog.wait_window()
        return selected_cols if selected_cols else None

    def _update_assignment_scroll_region(self):
        """Update the scroll region for the sample assignment panel - now a no-op since UI removed."""
        # Sample assignment UI was removed, this method is kept for compatibility
        pass

    def _update_main_scroll_region(self):
        """Update the scroll region for the main statistics canvas."""
        try:
            self.stats_canvas.configure(scrollregion=self.stats_canvas.bbox('all'))
        except Exception:
            pass

    def _refresh_sample_columns_view(self):
        """Refresh listbox and assignment panel based on selected polarity."""
        if not hasattr(self, 'detected_sample_cols'):
            return
        pol = getattr(self, 'sample_polarity_var', None)
        if pol:
            pol = self.sample_polarity_var.get()
        else:
            # Fallback: choose available
            pol = 'positive' if self.detected_sample_cols.get('positive') else 'negative'
        sample_cols = self.detected_sample_cols.get(pol, [])
        # Update listbox
        self.sample_cols_list.delete(0, tk.END)
        for c in sample_cols:
            self.sample_cols_list.insert(tk.END, c)
        # Rebuild assignment rows
        self.populate_sample_assignments(sample_cols)
        self.stats_log.insert(tk.END, f'Active polarity view: {pol} ({len(sample_cols)} sample columns).\n')
        self.stats_log.see(tk.END)

    def detect_statistics_sample_columns(self):
        """Detect sample columns - now requires completed normalization."""
        # Check if normalization has been completed
        if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
            messagebox.showwarning('Normalization Required', 
                                 'Please complete data normalization first before detecting sample columns.\n\n'
                                 'Use "Normalization & Test Normality" button to process your data.')
            return
        
        # Use the normalized data for sample column detection
        self.populate_sample_assignments_from_normalized_data()

    def run(self):
        """Run the GUI application"""
        self.root.mainloop()

    def populate_sample_assignments_from_normalized_data(self):
        """Populate sample assignments using sample columns from normalized data"""
        if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
            self.stats_log.insert(tk.END, 'No normalized data available for sample assignment.\n')
            return
        
        # Use verified sample columns from column verification step
        from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
        try:
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            
            # PRIORITY 1: Use verified sample columns if available
            sample_cols = []
            if mode == 'lipid':
                # Combine verified sample columns from both polarities
                if hasattr(self, 'verified_pos_lipid_sample_cols') and self.verified_pos_lipid_sample_cols:
                    sample_cols.extend([col for col in self.verified_pos_lipid_sample_cols if col in self.normalized_combined_df.columns])
                if hasattr(self, 'verified_neg_lipid_sample_cols') and self.verified_neg_lipid_sample_cols:
                    sample_cols.extend([col for col in self.verified_neg_lipid_sample_cols if col in self.normalized_combined_df.columns])
            else:
                # Combine verified sample columns from both polarities
                if hasattr(self, 'verified_pos_sample_cols') and self.verified_pos_sample_cols:
                    sample_cols.extend([col for col in self.verified_pos_sample_cols if col in self.normalized_combined_df.columns])
                if hasattr(self, 'verified_neg_sample_cols') and self.verified_neg_sample_cols:
                    sample_cols.extend([col for col in self.verified_neg_sample_cols if col in self.normalized_combined_df.columns])
            
            # Remove duplicates while preserving order
            seen = set()
            sample_cols = [col for col in sample_cols if not (col in seen or seen.add(col))]
            
            # FALLBACK: If no verified columns, auto-detect from normalized data
            if not sample_cols:
                self.stats_log.insert(tk.END, '⚠️ No verified sample columns found, auto-detecting from normalized data.\n')
                if mode == 'lipid':
                    # For lipid data, use robust feature detection
                    for col in self.normalized_combined_df.columns:
                        # Use normalized matching for lipid features
                        if not self._is_lipid_feature_col(col) and pd.api.types.is_numeric_dtype(self.normalized_combined_df[col]):
                            sample_cols.append(col)
                else:
                    # For metabolite data, use standard detection
                    feature_cols, sample_cols = detect_feature_and_sample_columns(self.normalized_combined_df)
            
            self.stats_log.insert(tk.END, f'✅ Using {len(sample_cols)} verified sample columns for group assignment.\n')
            
            # Populate assignment dropdowns
            self.populate_sample_assignments(sample_cols)
            
            self.stats_log.insert(tk.END, 'Sample assignments populated from verified columns. Ready for group assignment.\n')
            self.stats_log.see(tk.END)
            
        except Exception as e:
            self.stats_log.insert(tk.END, f'Error detecting sample columns: {str(e)}\n')
            self.stats_log.see(tk.END)

    def _extract_feature_assigned_columns(self, assignments) -> set[str]:
        """Extract columns mapped to Feature ID/Feature Column from assignment dicts.

        Supports both mapping shapes:
        - {"Feature ID": "Molecule"}
        - {"Molecule": "Feature ID"}
        """
        feature_cols: set[str] = set()
        if not isinstance(assignments, dict):
            return feature_cols

        def _is_feature_token(text: str) -> bool:
            t = str(text).strip().lower().replace('_', ' ')
            return ('feature' in t) or (t == 'id')

        for key, value in assignments.items():
            key_is_feature = _is_feature_token(key)
            val_is_feature = isinstance(value, str) and _is_feature_token(value)

            if key_is_feature:
                if isinstance(value, str) and value.strip():
                    feature_cols.add(value)
                elif isinstance(value, (list, tuple, set)):
                    for v in value:
                        if isinstance(v, str) and v.strip():
                            feature_cols.add(v)

            if val_is_feature and isinstance(key, str) and key.strip():
                feature_cols.add(key)

        return feature_cols

    def _get_all_assigned_feature_columns(self) -> set[str]:
        """Collect all known feature-assigned columns across modes/sheets."""
        assigned: set[str] = set()

        # Custom/preprocessed persisted assignments
        pre = self.memory_store.get('preprocessed_verified_assignments', {}) if hasattr(self, 'memory_store') else {}
        assigned.update(self._extract_feature_assigned_columns(pre))

        # Explicitly selected feature columns during custom import/backdoor import
        if hasattr(self, 'memory_store'):
            for col in self.memory_store.get('preprocessed_feature_cols', []) or []:
                if isinstance(col, str) and col.strip():
                    assigned.add(col)

        # Sheet-level verified assignments (metabolite/lipid/class)
        assignment_attrs = [
            'verified_pos_assignments',
            'verified_neg_assignments',
            'verified_pos_lipid_assignments',
            'verified_neg_lipid_assignments',
            'verified_pos_lipid_class_assignments',
            'verified_neg_lipid_class_assignments',
        ]
        for attr in assignment_attrs:
            assigned.update(self._extract_feature_assigned_columns(getattr(self, attr, {})))

        return assigned

    def _sanitize_group_sample_columns(self, sample_cols) -> list[str]:
        """Remove any columns assigned as Feature ID/Feature Column from group configuration."""
        cols = [c for c in (sample_cols or []) if isinstance(c, str) and c.strip()]
        forbidden = self._get_all_assigned_feature_columns()
        if not forbidden:
            return cols

        filtered = [c for c in cols if c not in forbidden]
        removed = [c for c in cols if c in forbidden]
        if removed and hasattr(self, 'stats_log'):
            self.stats_log.insert(
                tk.END,
                f"⚠️ Removed {len(removed)} feature-assigned column(s) from group configuration: {', '.join(removed)}\n"
            )
            self.stats_log.see(tk.END)
        return filtered

    def populate_sample_assignments(self, sample_cols):
        """Populate sample assignments using verified columns from Verify Columns step.
        
        All feature/sample detection is done in the Verify Columns dialog.
        This method trusts the user's verified column assignments and does NOT 
        apply any additional filtering.
        """
        # Final safety guard: never allow Feature ID/Feature Column assignments
        # to appear in group configuration, even if verification output is noisy.
        sample_cols = self._sanitize_group_sample_columns(sample_cols)
        self.detected_sample_cols = sample_cols
        self.sample_group_vars = {}
        
        # Initialize excluded samples store if missing
        if not hasattr(self, 'excluded_samples'):
            self.excluded_samples: set[str] = set()
        
        # Initialize sample_group_vars as empty StringVars for each sample
        for col in sample_cols:
            if col not in getattr(self, 'excluded_samples', set()):
                self.sample_group_vars[col] = tk.StringVar(value='')
        
        self.stats_log.insert(tk.END, f'✅ Using {len(sample_cols)} verified sample columns. Use "Auto-Assign by Pattern" button to configure groups.\n')
        self.stats_log.see(tk.END)

    def _on_filter_timing_change(self):
        """Handle filter timing selection change"""
        self._update_filter_timing_explanation()
        # Save to config
        if hasattr(self, '_stats_config_changed'):
            timing = self.filter_timing_var.get()
            self._stats_config_changed(log=f"Filter timing: {timing} normalization")
    
    def _update_filter_timing_explanation(self):
        """Update the explanation text based on selected filter timing"""
        if not hasattr(self, 'filter_timing_explanation'):
            return
        
        timing = self.filter_timing_var.get() if hasattr(self, 'filter_timing_var') else 'before'
        
        if timing == 'before':
            text = ("Before Normalization: Features with insufficient samples per group are zeroed out,\n"
                   "then normalized. Rows with no valid data are removed. This is more conservative.")
        else:
            text = ("After Normalization: All data is normalized first. During statistical tests, metabolites\n"
                   "that don't meet the threshold are skipped. No data modification occurs.")
        
        self.filter_timing_explanation.config(text=text)

    def _detect_sample_columns_for_optional_processing(self, df: pd.DataFrame, mode: str) -> list[str]:
        """Detect sample columns in a merged dataframe for optional processing steps.
        
        NOTE: This is a fallback only. Verified sample columns should be passed directly
        to _apply_optional_post_normalization_processing via the verified_sample_cols parameter.
        """
        if mode == 'lipid':
            return [
                col for col in df.columns
                if (not self._is_lipid_feature_col(col))
                and not is_statistics_metadata_col(col)
                and pd.api.types.is_numeric_dtype(df[col])
            ]

        from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
        _, sample_cols = detect_feature_and_sample_columns(df)
        return sample_cols

    def _get_verified_sample_columns_for_current_mode(self, mode: str, polarity: str = None, class_level: bool = False) -> list[str]:
        """Return verified sample columns for the current statistics mode.

        Verified columns are the source of truth. Class-level columns are preferred
        when available, otherwise the corresponding verified per-polarity sample
        columns are used as a fallback.
        """
        verified_cols: list[str] = []

        def extend_unique(values):
            for value in values or []:
                if value not in verified_cols:
                    verified_cols.append(value)

        if mode == 'lipid':
            if polarity == 'positive':
                if class_level:
                    extend_unique(getattr(self, 'verified_pos_lipid_class_sample_cols', []))
                extend_unique(getattr(self, 'verified_pos_lipid_sample_cols', []))
            elif polarity == 'negative':
                if class_level:
                    extend_unique(getattr(self, 'verified_neg_lipid_class_sample_cols', []))
                extend_unique(getattr(self, 'verified_neg_lipid_sample_cols', []))
            else:
                if class_level:
                    extend_unique(getattr(self, 'verified_pos_lipid_class_sample_cols', []))
                    extend_unique(getattr(self, 'verified_neg_lipid_class_sample_cols', []))
                extend_unique(getattr(self, 'verified_pos_lipid_sample_cols', []))
                extend_unique(getattr(self, 'verified_neg_lipid_sample_cols', []))
        else:
            if polarity == 'positive':
                extend_unique(getattr(self, 'verified_pos_sample_cols', []))
            elif polarity == 'negative':
                extend_unique(getattr(self, 'verified_neg_sample_cols', []))
            else:
                extend_unique(getattr(self, 'verified_pos_sample_cols', []))
                extend_unique(getattr(self, 'verified_neg_sample_cols', []))

        return verified_cols

    def _build_lipid_class_dataframe(self, df: pd.DataFrame, original_data: pd.DataFrame = None) -> tuple[pd.DataFrame | None, dict]:
        """Aggregate a lipid dataframe into one row per class using sample means.
        
        Parameters
        ----------
        df : pd.DataFrame
            Processed lipid dataframe (may contain imputed values or filtered rows)
        original_data : pd.DataFrame, optional
            Original lipid data ONLY for negative-polarity classes in raw input phase.
            Should be None when imputation is enabled to prevent imputed-only lipids
            from being excluded during class aggregation.
            
            When provided (only for raw input), filters out imputed values from class
            aggregation to prevent classes from appearing in polarities where they 
            don't naturally exist.
            
            IMPORTANT: Do NOT pass original_data when imputation is enabled.
            Imputation already handles preserving valid lipids for aggregation.
            Passing original_data after imputation causes ~49 lipids to disappear,
            leading to inconsistent row counts between imputation-disabled and 
            imputation-enabled workflows (480 vs 431 lipids).
        """
        if df is None or getattr(df, 'empty', True):
            return None, {'applied': False, 'reason': 'no_data'}

        class_col = None
        for candidate in ['Class', 'Lipid_Class', 'Class_name']:
            if candidate in df.columns:
                class_col = candidate
                break

        if not class_col:
            return None, {'applied': False, 'reason': 'no_class_column'}

        sample_cols = [
            col for col in df.columns
            if col != class_col and pd.api.types.is_numeric_dtype(df[col])
        ]
        if not sample_cols:
            return None, {'applied': False, 'reason': 'no_sample_columns', 'class_column': class_col}

        # Diagnostics: how many source rows/classes feed the class aggregation.
        class_series = df[class_col]
        non_null_class_mask = class_series.notna() & class_series.astype(str).str.strip().ne('')
        source_rows = int(len(df))
        valid_class_rows = int(non_null_class_mask.sum())
        unique_classes_in_source = int(class_series[non_null_class_mask].astype(str).nunique())

        try:
            from main_script.metabolites_visualization import aggregate_by_lipid_class
            class_df = aggregate_by_lipid_class(df, sample_cols, class_col=class_col, original_data=original_data)
        except Exception as e:
            return None, {'applied': False, 'reason': f'aggregation_failed: {e}', 'class_column': class_col}

        if class_col != 'Lipid_Class' and class_col in class_df.columns:
            class_df = class_df.rename(columns={class_col: 'Lipid_Class'})
        if 'lipid_class_id' in class_df.columns:
            if 'Lipid_Class' not in class_df.columns:
                class_df.insert(0, 'Lipid_Class', class_df['lipid_class_id'])
            class_df = class_df.drop(columns=['lipid_class_id'], errors='ignore')

        if 'Lipid_Class' in class_df.columns:
            ordered_cols = ['Lipid_Class'] + [c for c in class_df.columns if c != 'Lipid_Class']
            class_df = class_df[ordered_cols]

        class_count = int(len(class_df))
        collapsed_rows = max(0, valid_class_rows - class_count)

        return class_df.reset_index(drop=True), {
            'applied': True,
            'class_column': class_col,
            'class_count': class_count,
            'sample_column_count': int(len(sample_cols)),
            'source_rows': source_rows,
            'valid_class_rows': valid_class_rows,
            'unique_classes_in_source': unique_classes_in_source,
            'collapsed_rows': int(collapsed_rows),
        }

    def _export_pre_normalization_debug_workbook(self, raw_exports: dict) -> str | None:
        """Write a temporary workbook with pre-normalization lipid-ion and lipid-class sheets."""
        if not raw_exports:
            return None

        try:
            from datetime import datetime
            import shutil

            backup_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'logs'))
            os.makedirs(backup_dir, exist_ok=True)

            export_dir = None
            if hasattr(self, 'stats_results_folder') and self.stats_results_folder.get():
                candidate_dir = os.path.abspath(self.stats_results_folder.get())
                if os.path.isdir(candidate_dir):
                    export_dir = candidate_dir

            if export_dir is None:
                input_path = None
                for attr_name in ['annotated_metabolites_excel_path', 'id_annotated_excel_path', 'annotated_ids_excel_path']:
                    candidate = getattr(self, attr_name, None)
                    if candidate:
                        input_path = candidate
                        break
                if input_path:
                    export_dir = os.path.dirname(os.path.abspath(input_path))

            if export_dir is None:
                export_dir = backup_dir

            os.makedirs(export_dir, exist_ok=True)
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            export_path = os.path.join(export_dir, f'pre_normalization_raw_debug_{timestamp}.xlsx')
            backup_path = os.path.join(backup_dir, f'pre_normalization_raw_debug_{timestamp}.xlsx')

            summary_rows = []
            with pd.ExcelWriter(export_path, engine='openpyxl') as writer:
                for polarity in ('positive', 'negative'):
                    payload = raw_exports.get(polarity) or {}
                    ions_df = payload.get('ions')
                    ions_input_df = payload.get('ions_input')
                    class_df = payload.get('class')
                    class_input_df = payload.get('class_input')
                    imputed_df = payload.get('imputed')

                    if ions_df is not None and not getattr(ions_df, 'empty', True):
                        sheet_name = f'{polarity.capitalize()}_LipidIons'
                        ions_df.to_excel(writer, sheet_name=sheet_name, index=False)
                        summary_rows.append({
                            'Polarity': polarity.capitalize(),
                            'Sheet': sheet_name,
                            'Type': 'LipidIons',
                            'Rows': int(len(ions_df)),
                            'Stage': 'Filtered before imputation/normalization'
                        })

                    if ions_input_df is not None and not getattr(ions_input_df, 'empty', True):
                        summary_rows.append({
                            'Polarity': polarity.capitalize(),
                            'Sheet': f'{polarity.capitalize()}_LipidIons_Input',
                            'Type': 'LipidIons',
                            'Rows': int(len(ions_input_df)),
                            'Stage': 'Before pre-normalization filtering'
                        })

                    if class_df is not None and not getattr(class_df, 'empty', True):
                        sheet_name = f'{polarity.capitalize()}_LipidClass'
                        class_df.to_excel(writer, sheet_name=sheet_name, index=False)
                        summary_rows.append({
                            'Polarity': polarity.capitalize(),
                            'Sheet': sheet_name,
                            'Type': 'LipidClass',
                            'Rows': int(len(class_df)),
                            'Stage': 'Created before normalization'
                        })

                    if class_input_df is not None and not getattr(class_input_df, 'empty', True):
                        summary_rows.append({
                            'Polarity': polarity.capitalize(),
                            'Sheet': f'{polarity.capitalize()}_LipidClass_Input',
                            'Type': 'LipidClass',
                            'Rows': int(len(class_input_df)),
                            'Stage': 'Class built from unfiltered input'
                        })

                    if imputed_df is not None and not getattr(imputed_df, 'empty', True):
                        summary_rows.append({
                            'Polarity': polarity.capitalize(),
                            'Sheet': f'{polarity.capitalize()}_Imputed',
                            'Type': 'LipidIons',
                            'Rows': int(len(imputed_df)),
                            'Stage': 'After pre-normalization imputation'
                        })

                if summary_rows:
                    pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)

            if os.path.abspath(export_path) != os.path.abspath(backup_path):
                try:
                    shutil.copy2(export_path, backup_path)
                except Exception as copy_error:
                    self._thread_safe_log(f'⚠️ Could not create backup raw workbook copy: {copy_error}\n')

            return export_path
        except Exception as e:
            self._thread_safe_log(f'⚠️ Temporary raw workbook export failed: {e}\n')
            return None

    def _apply_optional_post_normalization_processing(self, combined_df: pd.DataFrame, mode: str, apply_imputation_prefilter: bool = False, verified_sample_cols: list[str] = None) -> tuple[pd.DataFrame, dict]:
        """Apply optional variability filtering, imputation, and PCA outlier removal.
        
        Parameters
        ----------
        verified_sample_cols : list[str], optional
            Pre-verified sample columns to use. If provided, only these columns are used for processing.
            If None, falls back to auto-detection (not recommended).
        """
        from main_script.metabolite_statistics_analysis import (
            apply_variability_filter,
            apply_imputation,
            apply_pca_outlier_filter,
        )

        processed_df = combined_df.copy()
        # Use verified sample columns if provided, otherwise auto-detect (fallback only)
        if verified_sample_cols:
            sample_cols = [col for col in verified_sample_cols if col in processed_df.columns]
        else:
            sample_cols = self._detect_sample_columns_for_optional_processing(processed_df, mode)
        group_map_full = self._parse_group_assignments() if hasattr(self, 'sample_group_vars') else {}
        alias_map = getattr(self, '_current_sample_column_aliases', {}) or {}
        group_map = {}
        for col_name, group_name in group_map_full.items():
            if not group_name:
                continue
            if col_name in sample_cols:
                group_map[col_name] = group_name
                continue
            aliased_name = alias_map.get(col_name)
            if aliased_name in sample_cols:
                group_map[aliased_name] = group_name

        enable_var = bool(self.enable_variability_filter_var.get()) if hasattr(self, 'enable_variability_filter_var') else False
        enable_imp = bool(self.enable_imputation_var.get()) if hasattr(self, 'enable_imputation_var') else False
        enable_pca = bool(self.enable_pca_outlier_var.get()) if hasattr(self, 'enable_pca_outlier_var') else False

        report = {
            'applied': bool(enable_var or enable_imp or enable_pca),
            'sample_column_count': len(sample_cols),
            'variability_filter': None,
            'imputation': None,
            'pca_outlier': None,
        }

        if not report['applied']:
            return processed_df, report

        self._thread_safe_log('\n🧪 Running optional post-normalization processing...\n')

        if enable_var:
            try:
                var_pct = float(self.variability_percent_var.get()) if hasattr(self, 'variability_percent_var') else 10.0
            except Exception:
                var_pct = 10.0
            processed_df, var_report = apply_variability_filter(
                processed_df,
                sample_cols,
                group_map=group_map,
                variance_percentile=var_pct,
                compute_anova=True,
                require_testable_rows=False
            )
            report['variability_filter'] = var_report
            removed = var_report.get('removed_features', 0)
            removed_low_var = var_report.get('removed_low_variance_features', removed)
            anova_untested = ((var_report.get('anova_summary') or {}).get('rows_untested', 0))
            self._thread_safe_log(
                f"  • Variability filter: removed {removed} total features "
                f"({removed_low_var} low-variance only; percentile cutoff={var_pct:.0f}% of row variance, ties may affect exact count). "
                f"ANOVA diagnostics: {anova_untested} rows were not ANOVA-testable and were kept.\n"
            )
            sample_cols = self._detect_sample_columns_for_optional_processing(processed_df, mode)

        if enable_imp:
            if apply_imputation_prefilter and group_map and len(group_map) == len(sample_cols):
                imp_min_pct = self._get_imputation_min_group_percent()
                imp_scope = self._get_imputation_prefilter_scope()
                scope_text = 'at least one group' if imp_scope == 'per_group' else 'all groups'
                processed_df, pre_imp_report = self._apply_imputation_prefilter(
                    processed_df,
                    sample_cols,
                    group_map,
                    imp_min_pct,
                    scope=imp_scope,
                )
                report['imputation_prefilter'] = pre_imp_report
                self._thread_safe_log(
                    f"  • Imputation pre-filter: min {imp_min_pct:.1f}% valid in {scope_text}; "
                    f"removed {pre_imp_report.get('removed', 0)} rows, kept {pre_imp_report.get('kept', 0)}.\n"
                )
            elif apply_imputation_prefilter:
                self._thread_safe_log(
                    "  • Imputation pre-filter skipped: incomplete/absent group assignments for current sample columns.\n"
                )
            else:
                if self._is_imputation_before_normalization_enabled():
                    self._thread_safe_log(
                        "  • Imputation pre-filter: already applied before normalization.\n"
                    )
                else:
                    self._thread_safe_log(
                        "  • Imputation pre-filter: not applied before normalization; using verified sample columns only.\n"
                    )

            imp_method = self.imputation_method_var.get().strip().lower() if hasattr(self, 'imputation_method_var') else 'half_min'
            try:
                knn_k = int(self.knn_neighbors_var.get()) if hasattr(self, 'knn_neighbors_var') else 5
            except Exception:
                knn_k = 5
            processed_df, imp_report = apply_imputation(
                processed_df,
                sample_cols,
                method=imp_method,
                knn_neighbors=knn_k,
                debug=True,
                log_fn=self._thread_safe_log,
            )
            report['imputation'] = imp_report
            self._thread_safe_log(
                f"  • Imputation ({imp_method}): filled {imp_report.get('imputed_cells', 0)} cells; "
                f"missing_after={imp_report.get('missing_after', 'NA')}; "
                f"zeros_after={imp_report.get('zeros_after', 'NA')}"
                + (
                    f"; unresolved_before_fallback={imp_report.get('unresolved_before_fallback', 0)}; "
                    f"unresolved_after_fallback={imp_report.get('unresolved_after_fallback', 0)}"
                    if imp_method == 'knn' else ""
                )
                + "\n"
            )

        if enable_pca:
            matrix = processed_df[sample_cols].apply(pd.to_numeric, errors='coerce') if sample_cols else pd.DataFrame()
            has_missing = bool(matrix.isna().any().any() or (matrix == 0).any().any()) if not matrix.empty else True
            if has_missing and not enable_imp:
                self._thread_safe_log('  • PCA outlier step skipped: missing values detected and imputation is disabled.\n')
                report['pca_outlier'] = {'applied': False, 'reason': 'missing_values_without_imputation'}
            else:
                pca_df, pca_report = apply_pca_outlier_filter(
                    processed_df,
                    sample_cols,
                    group_map=group_map if group_map else None,
                    threshold_sd=3.0,
                )
                processed_df = pca_df
                report['pca_outlier'] = {
                    k: v for k, v in pca_report.items() if k != 'scores'
                }
                self.last_pca_outlier_scores = pca_report.get('scores')
                removed_samples = pca_report.get('removed_samples', [])
                if removed_samples:
                    self._thread_safe_log(f"  • PCA outlier removal: removed {len(removed_samples)} samples: {', '.join(removed_samples)}\n")
                else:
                    self._thread_safe_log('  • PCA outlier removal: no samples exceeded threshold.\n')

        self._thread_safe_log('✅ Optional post-normalization processing complete.\n')
        return processed_df, report

    def on_test_type_change(self):
        """Handle radio button change for test type selection"""
        test_type = self.stat_test_type.get()
        if test_type == 'overall':
            # Enable overall combo, disable pairwise
            self.overall_combo.config(state='readonly')
            self.pairwise_combo.config(state='disabled')
        else:  # pairwise
            # Enable pairwise combo, disable overall
            self.overall_combo.config(state='disabled')
            self.pairwise_combo.config(state='readonly')
        # Check if ROTS parameters should be shown
        self._toggle_rots_parameters()
    
    def _toggle_rots_parameters(self):
        """Update ROTS configure button state based on current pairwise test."""
        # This function is kept for compatibility but ROTS button state
        # is now managed by _update_rots_button_state
        pass
    
    def _update_two_way_button_state(self, *args):
        """Enable Two-Way ANOVA config button only when two_way_anova or nonparametric_two_way is selected."""
        try:
            if hasattr(self, 'two_way_config_btn'):
                test = self.stat_overall_test.get().lower()
                if test in ['two_way_anova', 'nonparametric_two_way']:
                    self.two_way_config_btn.config(state='normal')
                else:
                    self.two_way_config_btn.config(state='disabled')
        except Exception:
            pass
    
    def _update_rots_button_state(self, *args):
        """Enable ROTS config button only when rots is selected as pairwise test."""
        try:
            if hasattr(self, 'rots_config_button'):
                test = self.stat_pairwise_test.get().lower()
                if test == 'rots':
                    self.rots_config_button.config(state='normal')
                else:
                    self.rots_config_button.config(state='disabled')
        except Exception:
            pass
    
    def _open_rots_config_dialog(self):
        """Open a dialog to configure ROTS parameters"""
        import tkinter as tk
        from tkinter import ttk, messagebox
        
        # Create dialog window
        dialog = tk.Toplevel(self.root)
        dialog.title("⚙️ Configure ROTS Parameters")
        dialog.geometry("600x550")  # Increased height to show all content
        dialog.resizable(True, True)  # Allow resizing
        dialog.minsize(550, 500)  # Set minimum size
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Center the dialog
        dialog.update_idletasks()
        x = (dialog.winfo_screenwidth() // 2) - (600 // 2)
        y = (dialog.winfo_screenheight() // 2) - (550 // 2)
        dialog.geometry(f"+{x}+{y}")
        
        # Main frame
        main_frame = tk.Frame(dialog, bg='white', padx=20, pady=20)
        main_frame.pack(fill='both', expand=True)
        
        # Title
        title_label = tk.Label(main_frame, text="ROTS (Reproducibility-Optimized Test Statistic)", 
                               bg='white', font=('Arial', 12, 'bold'), fg='#2c3e50')
        title_label.pack(pady=(0, 10))
        
        # Description
        desc_text = ("ROTS uses bootstrap resampling and reproducibility optimization\n"
                    "to find the optimal variance stabilization parameter.\n"
                    "Adjust parameters below for your analysis needs.")
        desc_label = tk.Label(main_frame, text=desc_text, bg='white', 
                             font=('Arial', 9), fg='#7f8c8d', justify='left')
        desc_label.pack(pady=(0, 15))
        
        # Parameters frame
        params_frame = tk.LabelFrame(main_frame, text="Parameters", bg='white', 
                                     font=('Arial', 10, 'bold'), padx=15, pady=15)
        params_frame.pack(fill='both', expand=True, pady=(0, 15))
        
        # Bootstrap iterations
        b_frame = tk.Frame(params_frame, bg='white')
        b_frame.pack(fill='x', pady=8)
        tk.Label(b_frame, text='Bootstrap iterations (B):', bg='white', 
                font=('Arial', 10), width=22, anchor='w').pack(side='left')
        b_spin = tk.Spinbox(b_frame, from_=100, to=10000, increment=100, 
                           textvariable=self.rots_B, width=10, font=('Arial', 10))
        b_spin.pack(side='left', padx=(10, 10))
        tk.Label(b_frame, text='(Higher = more accurate)', bg='white', 
                font=('Arial', 8, 'italic'), fg='#95a5a6').pack(side='left')
        
        # K parameter
        k_frame = tk.Frame(params_frame, bg='white')
        k_frame.pack(fill='x', pady=8)
        tk.Label(k_frame, text='Top K features:', bg='white', 
                font=('Arial', 10), width=22, anchor='w').pack(side='left')
        k_spin = tk.Spinbox(k_frame, from_=10, to=1000, increment=10, 
                           textvariable=self.rots_K, width=10, font=('Arial', 10))
        k_spin.pack(side='left', padx=(10, 10))
        tk.Label(k_frame, text='(For reproducibility check)', bg='white', 
                font=('Arial', 8, 'italic'), fg='#95a5a6').pack(side='left')
        
        # Alpha parameter
        alpha_frame = tk.Frame(params_frame, bg='white')
        alpha_frame.pack(fill='x', pady=8)
        tk.Label(alpha_frame, text='Alpha (top proportion):', bg='white', 
                font=('Arial', 10), width=22, anchor='w').pack(side='left')
        alpha_spin = tk.Spinbox(alpha_frame, from_=0.01, to=0.5, increment=0.05, 
                               textvariable=self.rots_alpha, width=10, 
                               font=('Arial', 10), format='%.2f')
        alpha_spin.pack(side='left', padx=(10, 10))
        tk.Label(alpha_frame, text='(Typically 0.1-0.3)', bg='white', 
                font=('Arial', 8, 'italic'), fg='#95a5a6').pack(side='left')
        
        # Random seed
        seed_frame = tk.Frame(params_frame, bg='white')
        seed_frame.pack(fill='x', pady=8)
        tk.Label(seed_frame, text='Random seed:', bg='white', 
                font=('Arial', 10), width=22, anchor='w').pack(side='left')
        seed_entry = tk.Entry(seed_frame, textvariable=self.rots_seed, 
                             width=12, font=('Arial', 10))
        seed_entry.pack(side='left', padx=(10, 10))
        tk.Label(seed_frame, text='(Empty = random)', bg='white', 
                font=('Arial', 8, 'italic'), fg='#95a5a6').pack(side='left')
        
        # Info box
        info_frame = tk.Frame(main_frame, bg='#e8f4f8', relief='solid', bd=1)
        info_frame.pack(fill='x', pady=(0, 15))
        info_icon = tk.Label(info_frame, text='ℹ️', bg='#e8f4f8', font=('Arial', 12))
        info_icon.pack(side='left', padx=(10, 5), pady=8)
        info_text = tk.Label(info_frame, 
                            text="Setting a seed ensures reproducible results across runs.\n"
                                 "Higher bootstrap iterations increase accuracy but take longer.",
                            bg='#e8f4f8', font=('Arial', 8), fg='#2980b9', 
                            justify='left', wraplength=400)
        info_text.pack(side='left', padx=(5, 10), pady=8)
        
        # Buttons
        button_frame = tk.Frame(main_frame, bg='white')
        button_frame.pack(fill='x')
        
        def on_ok():
            dialog.destroy()
            self._stats_config_changed(log=f"ROTS params: B={self.rots_B.get()}, K={self.rots_K.get()}, alpha={self.rots_alpha.get()}, seed={self.rots_seed.get()}")
        
        def on_reset():
            self.rots_B.set('1000')
            self.rots_K.set('100')
            self.rots_alpha.set('0.1')
            self.rots_seed.set('42')
        
        tk.Button(button_frame, text='Reset to Defaults', command=on_reset,
                 bg='#95a5a6', fg='white', font=('Arial', 9, 'bold'),
                 padx=15, pady=8).pack(side='left')
        tk.Button(button_frame, text='Cancel', command=dialog.destroy,
                 bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'),
                 padx=20, pady=8).pack(side='right', padx=(5, 0))
        tk.Button(button_frame, text='OK', command=on_ok,
                 bg='#27ae60', fg='white', font=('Arial', 9, 'bold'),
                 padx=30, pady=8).pack(side='right')
    
    def _get_rots_parameters(self):
        """Get ROTS parameters from GUI controls"""
        rots_params = {}
        try:
            # Bootstrap iterations
            rots_params['rots_B'] = int(self.rots_B.get()) if hasattr(self, 'rots_B') else 1000
        except Exception:
            rots_params['rots_B'] = 1000
        
        try:
            # Top K features
            rots_params['rots_K'] = int(self.rots_K.get()) if hasattr(self, 'rots_K') else 100
        except Exception:
            rots_params['rots_K'] = 100
        
        try:
            # Alpha parameter
            rots_params['rots_alpha'] = float(self.rots_alpha.get()) if hasattr(self, 'rots_alpha') else 0.1
        except Exception:
            rots_params['rots_alpha'] = 0.1
        
        try:
            # Random seed
            seed_str = self.rots_seed.get() if hasattr(self, 'rots_seed') else '42'
            if seed_str and seed_str.strip():
                rots_params['rots_seed'] = int(seed_str)
            else:
                rots_params['rots_seed'] = None
        except Exception:
            rots_params['rots_seed'] = None
        
        return rots_params
    
    def _update_groups_scroll_region(self):
        """Update the scroll region for the groups panel."""
        if not hasattr(self, 'groups_scrollable_frame') or not hasattr(self, 'groups_canvas'):
            return
        try:
            self.groups_canvas.update_idletasks()
            self.groups_scrollable_frame.update_idletasks()
            width = self.groups_scrollable_frame.winfo_reqwidth()
            height = self.groups_scrollable_frame.winfo_reqheight()
            self.groups_canvas.configure(scrollregion=(0, 0, width, height))
        except Exception:
            # Best-effort
            pass

    def refresh_group_ui(self):
        """Refresh the group UI with current group definitions"""
        # Clear existing group entries
        for widget in self.groups_scrollable_frame.winfo_children():
            widget.destroy()
        
        # Recreate group entries
        for i, (group_id, default_label) in enumerate(self.group_definitions.items()):
            id_frame = tk.Frame(self.groups_scrollable_frame, bg='#f0f0f0')
            id_frame.pack(fill='x', padx=3, pady=2)
            tk.Label(id_frame, text=f'{group_id}:', bg='#f0f0f0', width=8).pack(side='left')
            
            if group_id not in self.group_id_vars:
                self.group_id_vars[group_id] = tk.StringVar(value=default_label)
            
            entry_var = self.group_id_vars[group_id]
            # Limit entry width to half by setting explicit width instead of expand
            tk.Entry(id_frame, textvariable=entry_var, font=('Arial', 9), width=20).pack(side='left', padx=(5,5))
            
            # Immediate trace to auto-refresh assignments when user edits a label
            if not hasattr(self, '_group_label_trace_ids'):
                self._group_label_trace_ids = {}
            if group_id not in self._group_label_trace_ids:
                def _on_label_change(*args, gid=group_id):
                    if getattr(self, '_suppress_group_label_trace', False) or getattr(self, '_shutting_down', False):
                        return
                    try:
                        if getattr(self, '_group_label_change_after', None):
                            self.root.after_cancel(self._group_label_change_after)
                    except Exception:
                        pass
                    # Immediate update: update group_definitions and refresh all UI elements
                    try:
                        # Update the group definition immediately
                        new_label = self.group_id_vars[gid].get().strip()
                        if new_label:
                            old_label = self.group_definitions.get(gid, gid)
                            self.group_definitions[gid] = new_label
                            print(f"🔄 Group label changed: {gid} '{old_label}' → '{new_label}'")
                            
                            # Update base group combo if it exists
                            display_labels = [self.group_definitions[g] for g in self.group_definitions.keys()]
                            if hasattr(self, 'base_group_combo'):
                                try:
                                    current = self.base_group_combo.get()
                                    self.base_group_combo['values'] = [''] + display_labels
                                    if current == old_label or current == gid:
                                        self.base_group_combo.set(new_label)
                                    self.base_group_combo.update_idletasks()
                                except Exception:
                                    pass
                            
                            # Update all sample assignment dropdowns with new group labels
                            if hasattr(self, 'sample_group_combos'):
                                try:
                                    for sample_col, combo in self.sample_group_combos.items():
                                        if combo and combo.winfo_exists():
                                            combo['values'] = [''] + display_labels
                                            # If current selection was the old label, update to new label
                                            current_selection = self.sample_group_vars.get(sample_col)
                                            if current_selection and current_selection.get() == old_label:
                                                current_selection.set(new_label)
                                            combo.update_idletasks()
                                except Exception:
                                    pass
                    except Exception:
                        pass
                    # Short debounce for config persistence to avoid too frequent saves
                    try:
                        self._group_label_change_after = self.root.after(200, self._stats_config_changed)
                    except Exception:
                        pass
                trace_id = entry_var.trace_add('write', _on_label_change)
                self._group_label_trace_ids[group_id] = trace_id
        
        # Update scroll region after refreshing groups
        self._update_groups_scroll_region()
    
    def _show_norm_dropdown(self):
        """Show popup window with checkboxes for normalization method selection"""
        # Create popup window
        popup = tk.Toplevel(self.frame)
        popup.title("Select Normalization Methods")
        popup.geometry("550x500")
        popup.transient(self.frame)
        popup.grab_set()
        
        # Center the popup
        popup.update_idletasks()
        x = self.frame.winfo_rootx() + (self.frame.winfo_width() // 2) - (popup.winfo_width() // 2)
        y = self.frame.winfo_rooty() + (self.frame.winfo_height() // 2) - (popup.winfo_height() // 2)
        popup.geometry(f"+{x}+{y}")
        
        # Header
        header = tk.Frame(popup, bg='#3498db', pady=10)
        header.pack(fill='x')
        tk.Label(header, text="Select methods in the order you want them applied",
                bg='#3498db', fg='white', font=('Arial', 10, 'bold')).pack()
        
        # Scrollable frame for checkboxes
        canvas = tk.Canvas(popup, bg='white', highlightthickness=0)
        scrollbar = ttk.Scrollbar(popup, orient='vertical', command=canvas.yview)
        scrollable = tk.Frame(canvas, bg='white')
        
        scrollable.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side='left', fill='both', expand=True, padx=5, pady=5)
        scrollbar.pack(side='right', fill='y')
        
        # Tooltips for methods
        tooltips = {
            'median': 'Robust to outliers',
            'TIC': 'Total Ion Current',
            'PQN': 'Probabilistic Quotient Normalization',
            'IS': 'Internal Standard (requires IS identification)',
            'LOESS_QC': 'QC drift correction',
            'Rel_Abundance(%)': 'Relative abundance percentage',
            'quantile': 'Quantile normalization',
            'VSN': 'Variance Stabilizing Normalization (arcsinh)',
            'CLR': 'Centered Log-Ratio transformation',
            'zscore': 'Z-score standardization',
            'log2': 'Log2 transformation'
        }
        
        # Create checkboxes
        for display_name, method_name in self.norm_methods_list:
            frame = tk.Frame(scrollable, bg='white', pady=3)
            frame.pack(fill='x', padx=10)
            
            var = self.norm_method_vars[method_name]
            
            def make_callback(m=method_name):
                def callback(*args):
                    v = self.norm_method_vars[m]
                    if v.get():
                        self.norm_selection_counter[0] += 1
                        self.norm_method_order[m] = self.norm_selection_counter[0]
                    else:
                        self.norm_method_order[m] = None
                    self._update_norm_chain_display()
                return callback
            
            # Temporarily remove trace, update, then re-add
            cb = tk.Checkbutton(frame, text=display_name, variable=var, 
                               bg='white', font=('Arial', 10), anchor='w',
                               command=make_callback())
            cb.pack(side='left', fill='x', expand=True)
            
            # Order indicator
            order = self.norm_method_order.get(method_name)
            order_text = f"[{order}]" if order else ""
            order_label = tk.Label(frame, text=order_text, bg='white', 
                                  fg='#2980b9', font=('Arial', 9, 'bold'), width=4)
            order_label.pack(side='left')
            
            # Store reference to update later
            setattr(self, f'_popup_order_{method_name}', order_label)
            
            # Tooltip
            tooltip_text = tooltips.get(method_name, '')
            tk.Label(frame, text=f"({tooltip_text})", bg='white',
                    fg='#7f8c8d', font=('Arial', 8)).pack(side='left', padx=(5, 0))
        
        # Buttons at bottom
        button_frame = tk.Frame(popup, bg='#ecf0f1', pady=10)
        button_frame.pack(fill='x', side='bottom')
        
        btn_style = {'font': ('Arial', 9, 'bold'), 'pady': 5}
        
        tk.Button(button_frame, text='✓ Done', command=popup.destroy,
                 bg='#27ae60', fg='white', **btn_style).pack(fill='x', padx=20, pady=(5, 2))
        
        tk.Button(button_frame, text='🔄 Reset All', command=lambda: [self._reset_norm_selection(), self._refresh_popup_orders()],
                 bg='#95a5a6', fg='white', **btn_style).pack(fill='x', padx=20, pady=(2, 5))
        
        def _refresh_popup_orders():
            """Refresh order labels in popup"""
            for method_name in self.norm_method_order.keys():
                label = getattr(self, f'_popup_order_{method_name}', None)
                if label:
                    order = self.norm_method_order.get(method_name)
                    label.config(text=f"[{order}]" if order else "")
        
        self._refresh_popup_orders = _refresh_popup_orders
    
    def _update_norm_chain_display(self):
        """Update the normalization chain display based on selected methods"""
        # Get methods with their selection order
        ordered_methods = []
        for method, order in self.norm_method_order.items():
            if order is not None:
                ordered_methods.append((order, method))
        
        # Sort by order
        ordered_methods.sort()
        
        # Update order labels
        for method in self.norm_method_order.keys():
            label = getattr(self, f'norm_order_label_{method}', None)
            if label:
                order = self.norm_method_order[method]
                if order is not None:
                    label.config(text=f'[{order}]')
                else:
                    label.config(text='')
        
        # Build chain string
        if ordered_methods:
            chain = '+'.join([m[1] for m in ordered_methods])
            self.norm_chain_display.config(text=chain)
            self.stat_norm_method.set(chain)
        else:
            self.norm_chain_display.config(text='none')
            self.stat_norm_method.set('none')
    
    def _reset_norm_selection(self):
        """Reset all normalization selections"""
        for method, var in self.norm_method_vars.items():
            var.set(False)
        
        self.norm_selection_counter[0] = 0
        for method in self.norm_method_order.keys():
            self.norm_method_order[method] = None
        
        self._update_norm_chain_display()
        self._stats_config_changed(log="Normalization selection reset")
    
    def add_group(self):
        """Add a new group to the group definitions"""
        self.group_count += 1
        new_group_id = f'Group{self.group_count}'
        self.group_definitions[new_group_id] = f'Group{self.group_count}'
        self.group_id_vars[new_group_id] = tk.StringVar(value=f'Group{self.group_count}')
        self.refresh_group_ui()
        self.refresh_group_assignments()
    
    def remove_group(self):
        """Remove the last group from group definitions"""
        if self.group_count <= 2:  # Don't allow less than 2 groups
            messagebox.showwarning("Minimum Groups", "At least 2 groups are required for statistical analysis.")
            return
        
        last_group_id = f'Group{self.group_count}'
        if last_group_id in self.group_definitions:
            del self.group_definitions[last_group_id]
            if last_group_id in self.group_id_vars:
                del self.group_id_vars[last_group_id]
            self.group_count -= 1
            self.refresh_group_ui()
            self.refresh_group_assignments()

    def refresh_group_assignments(self):
        """Refresh group assignments - now handled via auto-assign pattern"""
        # Update group definitions from entry fields
        prev_suppress = getattr(self, '_suppress_group_label_trace', False)
        self._suppress_group_label_trace = True
        print(f"🔄 Refreshing group assignments from UI...")
        for group_id, var in self.group_id_vars.items():
            label = var.get().strip()
            if label:
                old_label = self.group_definitions.get(group_id, '')
                self.group_definitions[group_id] = label
                print(f"   {group_id}: '{old_label}' → '{label}'")
        
        # Update base group combo if it exists
        if hasattr(self, 'base_group_combo'):
            display_labels = [self.group_definitions[g] for g in self.group_definitions.keys()]
            self.base_group_combo['values'] = [''] + display_labels
        
        self._suppress_group_label_trace = prev_suppress
        
        # Log the update
        self.stats_log.insert(tk.END, f'✅ Updated group definitions: {", ".join(self.group_definitions.values())}\n')
        self.stats_log.see(tk.END)

    def auto_assign_groups(self, on_done=None):
        """Auto-assign groups based on common naming patterns."""
        if not hasattr(self, 'sample_group_vars') or not self.sample_group_vars:
            messagebox.showwarning('No Columns', 'No sample columns detected yet. Complete normalization and run "Assign Groups" first.')
            return

        # Track pending callback if invoked from a guarded action (Run Statistics / Covariate)
        self._pending_group_config_callback = on_done
        if getattr(self, '_group_config_window_open', False):
            try:
                if hasattr(self, '_group_config_window') and self._group_config_window:
                    self._group_config_window.lift()
                    self._group_config_window.focus_force()
            except Exception:
                pass
            return
        
        # Force update group_definitions from the UI entry fields before opening dialog
        self.refresh_group_assignments()
        
        # DEBUG: Print current group definitions
        print(f"🔍 DEBUG: Opening auto-assign dialog with group_definitions: {self.group_definitions}")
        
        callback_after_close = on_done

        pattern_window = tk.Toplevel(self.root)
        pattern_window.title('Auto-Assign Groups by Pattern')
        pattern_window.geometry('550x500')
        pattern_window.configure(bg='#f0f0f0')
        self._group_config_window = pattern_window
        self._group_config_window_open = True
        
        # Main frame with scrollbar
        main_frame = tk.Frame(pattern_window, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        # Create canvas and scrollbar
        canvas = tk.Canvas(main_frame, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#f0f0f0')
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        def configure_scroll(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(canvas_window, width=event.width)
        
        canvas.bind('<Configure>', configure_scroll)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind("<MouseWheel>", _on_mousewheel)
        scrollable_frame.bind("<MouseWheel>", _on_mousewheel)
        
        # Pattern definition frame
        # Top button frame for Apply and Cancel
        top_btn_frame = tk.Frame(pattern_window, bg='#f0f0f0')
        top_btn_frame.pack(fill='x', padx=10, pady=(5, 0))
        
        tk.Label(top_btn_frame, text='Auto-Assignment Patterns', font=('Arial', 12, 'bold'), bg='#f0f0f0').pack(side='left')
        
        # Pattern definition frame
        pattern_frame = tk.LabelFrame(scrollable_frame, text='Define Patterns', bg='#f0f0f0')
        pattern_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        tk.Label(pattern_frame, text='Enter keywords/patterns for each group (one per line):', bg='#f0f0f0').pack(anchor='w', padx=5)
        
        # Load saved patterns from config or instance variable
        saved_patterns = None  # None means never saved before, {} means user cleared patterns
        
        # First try to load from config file
        try:
            config_path = self._stats_config_file()
            print(f"📁 Config path: {config_path}")
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                    print(f"📄 Config file loaded successfully")
                    print(f"   Full config keys: {list(config_data.keys())}")
                    print(f"   Has 'auto_assign_patterns': {'auto_assign_patterns' in config_data}")
                    
                    if 'auto_assign_patterns' in config_data:
                        saved_patterns = config_data['auto_assign_patterns']
                        print(f"✅ Raw auto_assign_patterns from file: {saved_patterns}")
                        print(f"   Type: {type(saved_patterns)}")
                        print(f"   Keys: {list(saved_patterns.keys()) if saved_patterns else 'None'}")
                        print(f"   Values: {list(saved_patterns.values()) if saved_patterns else 'None'}")
                    else:
                        print(f"❌ 'auto_assign_patterns' key not found in config")
            else:
                print(f"⚠️ Config file does not exist: {config_path}")
        except Exception as e:
            print(f"⚠️ Error loading auto_assign_patterns from config: {e}")
            import traceback
            traceback.print_exc()
            pass
        
        # Fallback to instance variable if config not found
        if saved_patterns is None and hasattr(self, '_auto_assign_patterns'):
            saved_patterns = self._auto_assign_patterns
            print(f"✅ Using instance variable _auto_assign_patterns: {list(saved_patterns.keys()) if saved_patterns else 'empty'}")
        
        # Track if this is first time use (no saved patterns at all)
        first_time_use = (saved_patterns is None)
        if first_time_use:
            saved_patterns = {}
            print("📝 First time use - will show defaults and auto-discover")
        else:
            print(f"♻️ Using saved patterns from previous session")
        
        pattern_vars = {}
        for group_id in self.group_definitions.keys():
            # Always use the current label from group_definitions
            current_label = self.group_definitions[group_id]
            group_frame = tk.LabelFrame(pattern_frame, text=f'{group_id}: {current_label}', bg='#f0f0f0')
            group_frame.pack(fill='x', padx=5, pady=5)
            
            pattern_text = tk.Text(group_frame, height=3, font=('Arial', 9))
            pattern_text.pack(fill='x', padx=5, pady=5)
            pattern_vars[group_id] = pattern_text
            
            # Load saved patterns if available (this includes empty strings if user cleared them)
            if group_id in saved_patterns:
                pattern_content = saved_patterns[group_id]
                print(f"🔵 Inserting pattern for {group_id}: '{pattern_content}' (length={len(pattern_content)})")
                pattern_text.insert(tk.END, pattern_content)
                # Verify insertion
                actual_content = pattern_text.get('1.0', tk.END).strip()
                print(f"✅ Verification - Text widget now contains: '{actual_content}' (length={len(actual_content)})")
            # No defaults - leave empty for user to define

        # Auto-save patterns when text is modified
        def save_patterns_from_widgets():
            """Save current patterns from all text widgets to instance variable and config."""
            try:
                self._auto_assign_patterns = {}
                for gid, ptext in pattern_vars.items():
                    pattern_content = ptext.get('1.0', tk.END).strip()
                    self._auto_assign_patterns[gid] = pattern_content
                # Save to config file
                if hasattr(self, '_save_statistics_config'):
                    self._save_statistics_config()
                    print(f"💾 Auto-saved patterns: {list(self._auto_assign_patterns.keys())}")
            except Exception as e:
                print(f"⚠️ Auto-save failed: {e}")
        
        # Bind all pattern text widgets to auto-save on modification
        # Only bind KeyRelease (user typing), not <<Modified>> (fires during programmatic insert)
        for group_id, pattern_text in pattern_vars.items():
            def make_save_callback(gid):
                def callback(event):
                    print(f"🔔 Event triggered for {gid}: {event.type}")
                    save_patterns_from_widgets()
                return callback
            
            pattern_text.bind('<KeyRelease>', make_save_callback(group_id))
            # Don't bind <<Modified>> - it fires during .insert() and causes premature saves
        
        # ========== CURRENT GROUP ASSIGNMENTS DISPLAY ==========
        assignments_frame = tk.LabelFrame(scrollable_frame, text='📋 Current Group Assignments', 
                                         bg='#f0f0f0', font=('Arial',10,'bold'))
        assignments_frame.pack(fill='both', expand=True, padx=10, pady=(10,5))
        
        # Create canvas for scrollable assignment list
        assign_canvas = tk.Canvas(assignments_frame, bg='white', highlightthickness=1, height=300)
        assign_scrollbar = tk.Scrollbar(assignments_frame, orient='vertical', command=assign_canvas.yview)
        assign_scrollable = tk.Frame(assign_canvas, bg='white')
        
        assign_scrollable.bind(
            "<Configure>",
            lambda e: assign_canvas.configure(scrollregion=assign_canvas.bbox("all"))
        )
        
        assign_canvas_window = assign_canvas.create_window((0, 0), window=assign_scrollable, anchor="nw")
        assign_canvas.configure(yscrollcommand=assign_scrollbar.set)
        
        def configure_assign_scroll(event):
            assign_canvas.configure(scrollregion=assign_canvas.bbox("all"))
            assign_canvas.itemconfig(assign_canvas_window, width=event.width)
        
        assign_canvas.bind('<Configure>', configure_assign_scroll)
        assign_scrollbar.pack(side="right", fill="y")
        assign_canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        
        # Header
        header_frame = tk.Frame(assign_scrollable, bg='#e8f4f8')
        header_frame.pack(fill='x', pady=(0, 2))
        tk.Label(header_frame, text='Sample Column', bg='#e8f4f8', font=('Arial', 9, 'bold'), 
                width=25, anchor='w').pack(side='left', padx=5)
        tk.Label(header_frame, text='Assigned Group', bg='#e8f4f8', font=('Arial', 9, 'bold'), 
                width=20, anchor='w').pack(side='left', padx=5)
        
        # Display current assignments with dropdowns for manual editing
        assignment_combos = {}  # Store comboboxes for updating
        
        # Get group options for dropdowns (group labels, not IDs)
        group_options = [''] + [self.group_definitions.get(gid, gid) for gid in sorted(self.group_definitions.keys())]
        
        for col_name in sorted(self.sample_group_vars.keys()):
            row_frame = tk.Frame(assign_scrollable, bg='white')
            row_frame.pack(fill='x', pady=1)
            
            # Sample column name
            tk.Label(row_frame, text=col_name, bg='white', font=('Arial', 8), 
                    width=25, anchor='w').pack(side='left', padx=5)
            
            # Dropdown for group assignment (manual override)
            current_group = self.sample_group_vars[col_name].get() or ''
            group_combo = ttk.Combobox(row_frame, values=group_options, 
                                      textvariable=self.sample_group_vars[col_name],
                                      state='readonly', width=18, font=('Arial', 8))
            group_combo.pack(side='left', padx=5)
            assignment_combos[col_name] = group_combo
        
        # --- AUTO DISCOVER PREFIX TOKENS FROM SAMPLE COLUMN NAMES ---
        # Only run auto-discover on first time use (when no config exists at all)
        if first_time_use:
            try:
                sample_cols_existing = list(self.sample_group_vars.keys())
                prefix_counts: dict[str,int] = {}
                for col_name in sample_cols_existing:
                    # Token = up to first underscore OR first digit block
                    token = col_name.split('_')[0]
                    if not token:
                        continue
                    if len(token) > 30:  # avoid accidentally huge tokens
                        continue
                    prefix_counts[token] = prefix_counts.get(token, 0) + 1
                # Sort by frequency (desc)
                sorted_tokens = sorted(prefix_counts.items(), key=lambda x: (-x[1], x[0]))
                # Map tokens heuristically to groups where not already present in defaults
                for token, freq in sorted_tokens:
                    token_upper = token.upper()
                    # Skip purely numeric tokens
                    if token_upper.isdigit():
                        continue
                    # Check if token already present in any pattern box
                    already_present = False
                    for tv in pattern_vars.values():
                        if token in tv.get('1.0', tk.END):
                            already_present = True
                            break
                    if already_present:
                        continue
                    # Heuristic assignment rules
                    target_group_id = None
                    if token_upper.startswith('CTR') or token_upper.startswith('CTL') or token_upper == 'C':
                        target_group_id = 'Group1'
                    elif token_upper.startswith('ORTH') or token_upper.startswith('ORT') or token_upper.startswith('OR'): 
                        target_group_id = 'Group2'
                    elif token_upper.startswith('TBI') or token_upper.startswith('TBI') or token_upper.startswith('TB'):
                        target_group_id = 'Group3'
                    elif token_upper in {'TRT','TX','TREAT','DRUG'}:
                        target_group_id = 'Group3'
                    # If we found a target group, append token to its pattern box
                    if target_group_id and target_group_id in pattern_vars:
                        tv = pattern_vars[target_group_id]
                        existing_text = tv.get('1.0', tk.END).strip()
                        if existing_text:
                            tv.insert(tk.END, f'\n{token}')
                        else:
                            tv.insert(tk.END, token)
                # Log what was auto-detected
                detected_tokens_msg = ', '.join(f"{t}({c})" for t,c in sorted_tokens[:10])
                self.stats_log.insert(tk.END, f"🔍 Auto-detected sample name prefixes: {detected_tokens_msg}\n")
                self.stats_log.see(tk.END)
            except Exception as e:
                try:
                    self.stats_log.insert(tk.END, f"⚠️ Prefix auto-detect failed: {e}\n")
                    self.stats_log.see(tk.END)
                except Exception:
                    pass
        else:
            # Saved patterns loaded - skip auto-discovery
            try:
                self.stats_log.insert(tk.END, f"✅ Loaded saved auto-assign patterns from configuration.\n")
                self.stats_log.see(tk.END)
            except Exception:
                pass

        # Config persistence methods already defined at tab initialization
        # Buttons
        def apply_patterns():
            # Save patterns for persistence
            self._auto_assign_patterns = {}
            for group_id, pattern_text in pattern_vars.items():
                self._auto_assign_patterns[group_id] = pattern_text.get('1.0', tk.END).strip()
            
            # Save to config immediately
            try:
                self._save_statistics_config()
                self._log_stats('Auto-assign patterns saved to configuration.')
            except Exception as e:
                self._log_stats(f'Failed to save auto-assign patterns: {e}')
            
            assignments_made = 0
            # Determine if combobox values are labels (most recent refresh) or internal IDs
            internal_ids = set(self.group_definitions.keys())
            labels_set = set(self.group_definitions.values())
            use_labels = True  # default to labels (safer)
            try:
                first_val = next(iter(self.sample_group_vars.values())).get()
                if first_val in internal_ids and first_val not in labels_set:
                    use_labels = False
            except Exception:
                pass

            for col_name, group_var in self.sample_group_vars.items():
                original_value = group_var.get()
                matched = False
                #New
                all_pattern_candidates = []
                for group_id, pattern_text in pattern_vars.items():
                    patterns = [p.strip() for p in pattern_text.get('1.0', tk.END).splitlines() if p.strip()]
                    # if not patterns:
                    #     continue
                    # col_lower = col_name.lower()
                    for pattern in patterns:
                    #     if pattern.lower() and pattern.lower() in col_lower:
                    #         # Assign using label (preferred) unless UI still using IDs
                    #         assigned_value = self.group_definitions.get(group_id, group_id)
                    #         if not use_labels:
                    #             # Fallback: UI values still internal IDs
                    #             assigned_value = group_id
                    #         group_var.set(assigned_value)
                    #         assignments_made += 1
                    #         matched = True
                    #         break
                    # if matched:
                    #     break  # stop checking other groups for this column
                         if pattern.lower():
                            all_pattern_candidates.append((group_id, pattern, len(pattern)))
                
                # Sort by length descending (longest patterns first for more specific matching)
                all_pattern_candidates.sort(key=lambda x: -x[2])
                
                col_lower = col_name.lower()
                for group_id, pattern, _ in all_pattern_candidates:
                    if pattern.lower() in col_lower:
                        # Assign using label (preferred) unless UI still using IDs
                        assigned_value = self.group_definitions.get(group_id, group_id)
                        if not use_labels:
                            # Fallback: UI values still internal IDs
                            assigned_value = group_id
                        group_var.set(assigned_value)
                        assignments_made += 1
                        matched = True
                        break  # stop after first match
                
                # If nothing matched leave user's original manual assignment intact
                if not matched:
                    group_var.set(original_value)

            # Update assignment display dropdowns (comboboxes auto-update via StringVar binding)
            # Just ensure the current value is in the dropdown options
            for col_name, combo_widget in assignment_combos.items():
                current_group = self.sample_group_vars[col_name].get()
                if current_group and current_group not in group_options:
                    # Add it temporarily if it's a custom value
                    combo_widget['values'] = group_options + [current_group]
            
            # Post-assignment diagnostics: count groups and separate by polarity
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            
            # Determine which columns are positive vs negative
            pos_cols = []
            neg_cols = []
            
            if mode == 'lipid':
                # For lipids, use stored detection results
                pos_cols = getattr(self, 'detected_pos_lipid_sample_cols', [])
                neg_cols = getattr(self, 'detected_neg_lipid_sample_cols', [])
            else:
                # For metabolites, use stored detection results
                pos_cols = getattr(self, 'detected_pos_sample_cols', [])
                neg_cols = getattr(self, 'detected_neg_sample_cols', [])
            
            # Count groups separately for positive and negative
            pos_group_counts: dict[str,int] = {}
            neg_group_counts: dict[str,int] = {}
            
            for col_name, gv in self.sample_group_vars.items():
                grp = gv.get()
                if col_name in pos_cols:
                    pos_group_counts[grp] = pos_group_counts.get(grp, 0) + 1
                elif col_name in neg_cols:
                    neg_group_counts[grp] = neg_group_counts.get(grp, 0) + 1
            
            # Log results with polarity breakdown
            self.stats_log.insert(tk.END, f"\n📊 Auto-assignment completed: {assignments_made} columns assigned\n")
            
            if pos_group_counts:
                pos_summary = ', '.join(f"{gid}:{cnt}" for gid,cnt in pos_group_counts.items())
                self.stats_log.insert(tk.END, f"   ✅ Positive: {pos_summary} (Total: {sum(pos_group_counts.values())})\n")
            
            if neg_group_counts:
                neg_summary = ', '.join(f"{gid}:{cnt}" for gid,cnt in neg_group_counts.items())
                self.stats_log.insert(tk.END, f"   ✅ Negative: {neg_summary} (Total: {sum(neg_group_counts.values())})\n")
            
            # Inform user that class sheets will use same grouping (lipid mode only)
            if mode == 'lipid':
                has_class_sheets = (
                    (hasattr(self, 'pos_class_inherits_grouping') and self.pos_class_inherits_grouping) or
                    (hasattr(self, 'neg_class_inherits_grouping') and self.neg_class_inherits_grouping)
                )
                if has_class_sheets:
                    self.stats_log.insert(tk.END, f"   ℹ️  Class sheets will use the same group assignments\n")
            
            # Warn if any expected tokens (e.g., TBI) not assigned
            if 'Group3' in self.group_definitions:
                # Check if any column assigned to Group3 by either ID or label
                group3_label = self.group_definitions['Group3']
                any_group3 = any(gv.get() in {'Group3', group3_label} for gv in self.sample_group_vars.values())
                if not any_group3:
                    self.stats_log.insert(tk.END, f"⚠️ No columns matched Group3 (e.g., {group3_label}). Add a pattern like 'TBI' and re-apply if needed.\n")
            
            self.stats_log.see(tk.END)
            
            # Save the updated configuration after assignment
            try:
                self._stats_config_changed(log='Auto-assignment complete. Configuration saved.')
            except Exception:
                pass
            
            # Don't close the window - let user review and manually adjust if needed
        
        def _close_without_action():
            """Close the window without triggering post-actions."""
            self._pending_group_config_callback = None
            self._group_config_window_open = False
            self._group_config_window = None
            try:
                pattern_window.destroy()
            except Exception:
                pass

        def done_and_close():
            """Save patterns and close the window."""
            save_patterns_from_widgets()
            
            # Save final group assignments to config
            try:
                self._stats_config_changed(log='Group assignments saved.')
            except Exception:
                pass
            
            # Enable Run Statistics button after groups are configured
            if hasattr(self, 'run_stats_btn'):
                self.run_stats_btn.config(state='normal')
                self._log_stats('✅ Groups configured successfully. Run Statistics button enabled.')
            
            self._groups_configured_once = True
            self._group_config_window_open = False
            self._group_config_window = None
            try:
                pattern_window.destroy()
            except Exception:
                pass

            cb = getattr(self, '_pending_group_config_callback', None) or callback_after_close
            self._pending_group_config_callback = None
            if callable(cb):
                self.root.after(0, cb)
        
        # Bottom button frame
        bottom_btn_frame = tk.Frame(pattern_window, bg='#f0f0f0')
        bottom_btn_frame.pack(side='bottom', fill='x', padx=10, pady=(5, 10))
        
        tk.Button(bottom_btn_frame, text='Apply Patterns', command=apply_patterns, 
            bg='#27ae60', fg='white', font=('Arial', 10, 'bold'), width=15).pack(side='left', padx=5)
        tk.Button(bottom_btn_frame, text='Done', command=done_and_close, 
            bg='#3498db', fg='white', font=('Arial', 10, 'bold'), width=15).pack(side='left', padx=5)
        tk.Button(bottom_btn_frame, text='Cancel', command=_close_without_action, 
            bg='#e74c3c', fg='white', font=('Arial', 10, 'bold'), width=15).pack(side='right', padx=5)

        pattern_window.protocol("WM_DELETE_WINDOW", _close_without_action)

    def _parse_group_assignments(self):
        """Parse group assignments from the dropdown selections"""
        mapping = {}
        if hasattr(self, 'sample_group_vars'):
            # Build reverse lookup: label -> label to keep exported semantics stable
            current_labels = {gid: self.group_definitions[gid] for gid in self.group_definitions}
            label_set = set(current_labels.values())
            for col_name, group_var in self.sample_group_vars.items():
                sel = group_var.get()
                # If selection is an internal ID, convert to its label
                if sel in current_labels:
                    mapping[col_name] = current_labels[sel]
                elif sel in label_set:
                    mapping[col_name] = sel
                else:
                    # Fallback to first defined label if something unexpected
                    mapping[col_name] = next(iter(label_set)) if label_set else sel
        return mapping

    def _get_min_group_size(self) -> int:
        """Return the minimum samples per group threshold configured by the user (absolute count)."""
        try:
            value = int(self.min_samples_per_group_var.get()) if hasattr(self, 'min_samples_per_group_var') else 2
        except Exception:
            value = 2
        return max(1, value)
    
    def _get_min_group_size_type(self) -> str:
        """Return threshold type: 'absolute' or 'percentage'."""
        try:
            return self.min_samples_type_var.get() if hasattr(self, 'min_samples_type_var') else 'absolute'
        except Exception:
            return 'absolute'
    
    def _get_threshold_type(self) -> str:
        """Return threshold type: 'percentage' or 'count'."""
        try:
            val = self.min_samples_type_var.get() if hasattr(self, 'min_samples_type_var') else 'percentage'
            # Normalize to 'percentage' or 'count' (handles both 'percentage'/'absolute' and 'percentage'/'count')
            if val in ('percentage', '%'):
                return 'percentage'
            else:
                return 'count'
        except Exception:
            return 'percentage'
    
    def _get_min_group_size_percent(self) -> float:
        """Return minimum percentage threshold."""
        try:
            value = float(self.min_samples_percent_var.get()) if hasattr(self, 'min_samples_percent_var') else 50.0
        except Exception:
            value = 50.0
        return max(0.0, min(100.0, value))

    def _is_imputation_enabled(self) -> bool:
        """Return whether Step 4b imputation is enabled."""
        try:
            return bool(self.enable_imputation_var.get()) if hasattr(self, 'enable_imputation_var') else False
        except Exception:
            return False

    def _is_imputation_before_normalization_enabled(self) -> bool:
        """Return whether pre-normalization imputation (Step 4a) is enabled."""
        try:
            return bool(self.enable_imputation_before_var.get()) if hasattr(self, 'enable_imputation_before_var') else False
        except Exception:
            return False

    def _get_imputation_min_group_percent(self) -> float:
        """Return imputation pre-filter threshold (% non-missing per group)."""
        try:
            value = float(self.imputation_min_group_percent_var.get()) if hasattr(self, 'imputation_min_group_percent_var') else 50.0
        except Exception:
            value = 50.0
        return max(0.0, min(100.0, value))

    def _get_imputation_prefilter_scope(self) -> str:
        """Return imputation pre-filter scope: per_group or all_groups."""
        try:
            value = str(self.imputation_prefilter_scope_var.get()) if hasattr(self, 'imputation_prefilter_scope_var') else 'per_group'
        except Exception:
            value = 'per_group'
        value = value.strip().lower()
        return value if value in ('per_group', 'all_groups') else 'per_group'

    def _apply_imputation_prefilter(
        self,
        df: pd.DataFrame,
        sample_cols: list[str],
        group_map: dict[str, str],
        min_percent: float,
        scope: str = 'per_group'
    ) -> tuple[pd.DataFrame, dict]:
        """Exclude features that fail minimum valid-sample percentage in all groups before imputation.

        A value is valid if it is non-missing and non-zero.
        """
        if df is None or df.empty or not sample_cols or not group_map:
            return df, {'applied': False, 'reason': 'missing_inputs', 'kept': len(df) if df is not None else 0, 'removed': 0}

        group_cols: dict[str, list[str]] = {}
        for c in sample_cols:
            g = group_map.get(c)
            if g:
                group_cols.setdefault(g, []).append(c)

        if not group_cols:
            return df, {'applied': False, 'reason': 'no_group_columns', 'kept': len(df), 'removed': 0}

        vals = df[sample_cols].apply(pd.to_numeric, errors='coerce')
        valid_mask = vals.notna() & (vals != 0)

        scope = str(scope or 'per_group').strip().lower()
        if scope not in ('per_group', 'all_groups'):
            scope = 'per_group'

        keep_mask = pd.Series(False, index=df.index) if scope == 'per_group' else pd.Series(True, index=df.index)
        per_group_thresholds = {}
        for grp, cols in group_cols.items():
            n = len(cols)
            required = int(np.ceil(n * (min_percent / 100.0)))
            required = max(1, required)
            per_group_thresholds[grp] = required
            grp_valid_counts = valid_mask[cols].sum(axis=1)
            grp_pass = (grp_valid_counts >= required)
            if scope == 'per_group':
                keep_mask |= grp_pass
            else:
                keep_mask &= grp_pass

        filtered = df.loc[keep_mask].copy()
        return filtered, {
            'applied': True,
            'min_percent': float(min_percent),
            'scope': scope,
            'group_thresholds': per_group_thresholds,
            'total_before': int(len(df)),
            'kept': int(len(filtered)),
            'removed': int(len(df) - len(filtered)),
        }

    def _filter_groups_by_min_samples(self, group_map: dict[str, str], min_required: int) -> tuple[dict[str, str], dict[str, int], dict[str, int]]:
        """Filter out groups whose sample counts fall below the required threshold.

        Returns a tuple of (filtered_map, counts_by_group, excluded_groups).
        """
        if not group_map:
            return {}, {}, {}
        try:
            threshold = max(1, int(min_required))
        except Exception:
            threshold = 1
        counts = Counter(group_map.values())
        excluded = {grp: cnt for grp, cnt in counts.items() if cnt < threshold}
        filtered = {sample: grp for sample, grp in group_map.items() if counts.get(grp, 0) >= threshold}
        return filtered, dict(counts), excluded

    def _get_pos_neg_for_stats(self):
        """Get positive and negative DataFrames based on current data mode (metabolite/lipid)."""
        pos_df = None
        neg_df = None
        pos_class_df = None
        neg_class_df = None
        
        # Determine current mode
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        
        # Debug: Log current memory store status
        if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
            self.stats_log.insert(tk.END, f"🔍 DEBUG: Mode={mode}, Memory store keys: {list(self.memory_store.keys())}\n")
            
            if mode == 'lipid':
                # Lipid mode: look for lipid and class DataFrames (support multiple naming conventions)
                lipid_pos_keys = ['pos_lipid_df', 'lipid_pos_df']
                for key in lipid_pos_keys:
                    if key in self.memory_store and self.memory_store.get(key) is not None:
                        pos_df = self.memory_store.get(key)
                        self.stats_log.insert(tk.END, f"✅ Found positive lipid data in '{key}': {len(pos_df)} rows, {len(pos_df.columns)} columns\n")
                        break
                
                lipid_neg_keys = ['neg_lipid_df', 'lipid_neg_df']
                for key in lipid_neg_keys:
                    if key in self.memory_store and self.memory_store.get(key) is not None:
                        neg_df = self.memory_store.get(key)
                        self.stats_log.insert(tk.END, f"✅ Found negative lipid data in '{key}': {len(neg_df)} rows, {len(neg_df.columns)} columns\n")
                        break
                
                # Load class DataFrames if available (support multiple naming conventions)
                class_pos_keys = ['pos_lipid_class_df', 'lipid_pos_class_df']
                for key in class_pos_keys:
                    if key in self.memory_store and self.memory_store[key] is not None:
                        pos_class_df = self.memory_store[key]
                        self.stats_log.insert(tk.END, f"✅ Found positive lipid class data in '{key}': {len(pos_class_df)} rows\n")
                        break
                
                class_neg_keys = ['neg_lipid_class_df', 'lipid_neg_class_df']
                for key in class_neg_keys:
                    if key in self.memory_store and self.memory_store[key] is not None:
                        neg_class_df = self.memory_store[key]
                        self.stats_log.insert(tk.END, f"✅ Found negative lipid class data in '{key}': {len(neg_class_df)} rows\n")
                        break
                
                if pos_df is None and neg_df is None:
                    self.stats_log.insert(tk.END, "❌ No lipid data found in memory store\n")
                    self.stats_log.insert(tk.END, "💡 REQUIRED: Import Excel with Positive_Lipids and/or Negative_Lipids sheets (at least one required)\n")
            else:
                # Metabolite mode: look for ID DataFrames
                id_pos_keys = ['pos_id_df', 'clean_pos_id_df']
                for key in id_pos_keys:
                    if key in self.memory_store and self.memory_store.get(key) is not None:
                        pos_df = self.memory_store.get(key)
                        self.stats_log.insert(tk.END, f"✅ Found positive ID data in '{key}': {len(pos_df)} rows, {len(pos_df.columns)} columns\n")
                        break
                
                if pos_df is None:
                    self.stats_log.insert(tk.END, "❌ No positive ID data found in memory store\n")
                    self.stats_log.insert(tk.END, "💡 TIP: To include positive data, complete ID Annotation or import Excel with a Pos_id sheet (optional).\n")
                    
                # Check for ID DataFrames (from ID annotation workflow OR Excel import)
                id_neg_keys = ['neg_id_df', 'clean_neg_id_df']
                for key in id_neg_keys:
                    if key in self.memory_store and self.memory_store.get(key) is not None:
                        neg_df = self.memory_store.get(key)
                        self.stats_log.insert(tk.END, f"✅ Found negative ID data in '{key}': {len(neg_df)} rows, {len(neg_df.columns)} columns\n")
                        break
                        
                if neg_df is None:
                    self.stats_log.insert(tk.END, "❌ No negative ID data found in memory store\n")
                    self.stats_log.insert(tk.END, "💡 TIP: To include negative data, complete ID Annotation or import Excel with a Neg_id sheet (optional).\n")
                        
        else:
            self.stats_log.insert(tk.END, "❌ Memory store not found or not a dictionary\n")
            
        self.stats_log.see(tk.END)
        
        # Store class DataFrames for later use if in lipid mode
        if mode == 'lipid':
            self._current_pos_class_df = pos_class_df
            self._current_neg_class_df = neg_class_df
        
        return pos_df, neg_df

    def check_memory_store_status(self):
        """Check and display the current status of the memory store for debugging."""
        # Determine current mode
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        
        self.stats_log.insert(tk.END, f"\n🔍 === MEMORY STORE STATUS CHECK ({mode.upper()} MODE) ===\n")
        
        if not hasattr(self, 'memory_store'):
            self.stats_log.insert(tk.END, "❌ ERROR: memory_store attribute does not exist\n")
            self.stats_log.see(tk.END)
            return
            
        if not isinstance(self.memory_store, dict):
            self.stats_log.insert(tk.END, f"❌ ERROR: memory_store is not a dictionary (type: {type(self.memory_store)})\n")
            self.stats_log.see(tk.END)
            return
            
        if not self.memory_store:
            self.stats_log.insert(tk.END, "❌ WARNING: memory_store is empty\n")
            if mode == 'lipid':
                self.stats_log.insert(tk.END, "💡 SOLUTION: Please run Lipid Data Cleaning first\n")
            else:
                self.stats_log.insert(tk.END, "💡 SOLUTION: Please run Data Cleaning and ID Annotation first\n")
            self.stats_log.see(tk.END)
            return
            
        self.stats_log.insert(tk.END, f"✅ Memory store contains {len(self.memory_store)} items:\n")
        
        for key, value in self.memory_store.items():
            if value is None:
                self.stats_log.insert(tk.END, f"  • {key}: None\n")
            elif hasattr(value, '__len__') and hasattr(value, 'columns'):  # DataFrame-like
                self.stats_log.insert(tk.END, f"  • {key}: DataFrame with {len(value)} rows, {len(value.columns)} columns\n")
            else:
                self.stats_log.insert(tk.END, f"  • {key}: {type(value).__name__}\n")
        
        # Check for required DataFrames based on mode
        if mode == 'lipid':
            # Check for lipid data
            required_lipid_keys = ['lipid_pos_df', 'pos_lipid_df']
            required_lipid_neg_keys = ['lipid_neg_df', 'neg_lipid_df']
            optional_class_keys = ['lipid_pos_class_df', 'pos_lipid_class_df', 'lipid_neg_class_df', 'neg_lipid_class_df']
            
            found_pos_lipid = []
            found_neg_lipid = []
            found_class_data = []
            
            # Check for positive lipid data
            for key in required_lipid_keys:
                if key in self.memory_store and self.memory_store[key] is not None:
                    found_pos_lipid.append(key)
            
            # Check for negative lipid data
            for key in required_lipid_neg_keys:
                if key in self.memory_store and self.memory_store[key] is not None:
                    found_neg_lipid.append(key)
            
            # Check for class data
            for key in optional_class_keys:
                if key in self.memory_store and self.memory_store[key] is not None:
                    found_class_data.append(key)
            
            if found_pos_lipid or found_neg_lipid:
                self.stats_log.insert(tk.END, f"✅ Lipid data ready:\n")
                if found_pos_lipid:
                    self.stats_log.insert(tk.END, f"   • Positive: {', '.join(found_pos_lipid)}\n")
                if found_neg_lipid:
                    self.stats_log.insert(tk.END, f"   • Negative: {', '.join(found_neg_lipid)}\n")
                if found_class_data:
                    self.stats_log.insert(tk.END, f"   • Class data: {', '.join(found_class_data)}\n")
                self.stats_log.insert(tk.END, "✅ Statistics tab ready to use\n")
            else:
                self.stats_log.insert(tk.END, "❌ NO LIPID DATA FOUND\n")
                self.stats_log.insert(tk.END, "💡 REQUIRED: Complete Lipid Data Cleaning OR import Excel with Positive_Lipids and/or Negative_Lipids sheets\n")
        else:
            # Check for metabolite ID data
            required_id_keys = ['pos_id_df', 'neg_id_df']
            legacy_keys = ['pos_enhanced_df', 'neg_enhanced_df', 'pos_result', 'neg_result']
            found_id_data = []
            found_legacy_data = []
            
            # Check for required ID DataFrames
            for key in required_id_keys:
                if key in self.memory_store and self.memory_store[key] is not None:
                    found_id_data.append(key)
            
            # Check for legacy DataFrames (informational only)
            for key in legacy_keys:
                if key in self.memory_store and self.memory_store[key] is not None:
                    found_legacy_data.append(key)
            
            if found_id_data:
                self.stats_log.insert(tk.END, f"✅ ID-annotated data ready: {', '.join(found_id_data)}\n")
                self.stats_log.insert(tk.END, "✅ Statistics tab ready to use\n")
            else:
                self.stats_log.insert(tk.END, "❌ NO ID-ANNOTATED DATA FOUND\n")
                self.stats_log.insert(tk.END, "💡 REQUIRED: Complete ID Annotation OR import Excel with Pos_id/Neg_id sheets\n")
                if found_legacy_data:
                    self.stats_log.insert(tk.END, f"⚠️ Legacy data found but not supported: {', '.join(found_legacy_data)}\n")
                    self.stats_log.insert(tk.END, "⚠️ Statistics now requires ID-annotated data only\n")
        
        self.stats_log.insert(tk.END, "===================================\n\n")
        self.stats_log.see(tk.END)

    def run_statistics_pipeline(self):
        """Launch statistics pipeline in a separate thread to prevent GUI freezing."""

        # Workflow notice: when Step 4b imputation is enabled, skip original min-samples filtering.
        if self._is_imputation_enabled():
            imp_pct = self._get_imputation_min_group_percent()
            imp_scope = self._get_imputation_prefilter_scope()
            scope_text = 'at least one group' if imp_scope == 'per_group' else 'all groups'
            messagebox.showinfo(
                "Imputation Workflow Active",
                "Step 4b Imputation is enabled.\n\n"
                "• Original 'Min samples per group' filtering will be skipped.\n"
                f"• Imputation pre-filter will be used instead ({imp_pct:.1f}% valid values in {scope_text}).\n"
                "• Valid means non-zero and non-missing.\n\n"
                "This avoids double filtering before statistical tests."
            )
        
        # Check if group assignments exist
        has_assignments = False
        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
            assigned_count = sum(1 for var in self.sample_group_vars.values() if var.get() and var.get() != 'Unassigned')
            has_assignments = assigned_count > 0
        
        # Prompt user if no group assignments
        if not has_assignments:
            messagebox.showerror(
                "Configure Groups Required",
                "❌ No group assignments detected!\n\n"
                "REQUIRED STEPS:\n"
                "1. Click 'Configure Groups' button\n"
                "2. Use 'Apply Patterns' to auto-assign groups\n"
                "   OR manually assign each sample to a group\n"
                "3. Click 'Done' to save assignments\n"
                "4. Return here and run normalization\n\n"
                "⚠️ Pre-normalization filtering requires group assignments\n"
                "to remove features with insufficient replicates per group.\n"
                "This is critical for accurate statistical results.\n\n"
                "Please configure groups before proceeding."
            )
            self.stats_log.insert(tk.END, "\n❌ Normalization cancelled. Configure groups first using 'Configure Groups' button.\n\n")
            self.stats_log.see(tk.END)
            return
        
        # Disable the button to prevent multiple simultaneous runs
        for widget in self.frame.winfo_children():
            if isinstance(widget, tk.Frame):
                for child in widget.winfo_children():
                    if isinstance(child, tk.Button) and "Run Normalization" in child.cget('text'):
                        child.config(state='disabled')
        
        # Start the actual pipeline in a separate thread
        threading.Thread(target=self._run_statistics_pipeline_threaded, daemon=True).start()
    
    def _normalize_col(self, col):
        """Normalize a column name for comparison (lowercase, no spaces/underscores)."""
        if col is None:
            return ''
        return str(col).lower().replace('_', '').replace(' ', '')
    
    def _get_verified_id_column(self, dataframe=None):
        """
        Get the verified ID column name from column verification.
        Returns the actual column name that should be used for metabolite/lipid IDs.
        
        Args:
            dataframe: Optional dataframe to verify column exists in
            
        Returns:
            str: The verified ID column name ('Name' for metabolites, 'LipidID' for lipids)
        """
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        
        # Default column names
        default_col = 'LipidID' if mode == 'lipid' else 'Name'
        
        # Try to get from verified assignments
        if mode == 'lipid':
            # Check positive lipid assignments
            if hasattr(self, 'verified_pos_lipid_assignments'):
                if 'LipidID' in self.verified_pos_lipid_assignments:
                    verified_col = self.verified_pos_lipid_assignments['LipidID']
                    if verified_col and verified_col != 'Ignore':
                        return verified_col
            # Check negative lipid assignments
            if hasattr(self, 'verified_neg_lipid_assignments'):
                if 'LipidID' in self.verified_neg_lipid_assignments:
                    verified_col = self.verified_neg_lipid_assignments['LipidID']
                    if verified_col and verified_col != 'Ignore':
                        return verified_col
        else:
            # Check positive metabolite assignments
            if hasattr(self, 'verified_pos_assignments'):
                if 'Name' in self.verified_pos_assignments:
                    verified_col = self.verified_pos_assignments['Name']
                    if verified_col and verified_col != 'Ignore':
                        return verified_col
            # Check negative metabolite assignments
            if hasattr(self, 'verified_neg_assignments'):
                if 'Name' in self.verified_neg_assignments:
                    verified_col = self.verified_neg_assignments['Name']
                    if verified_col and verified_col != 'Ignore':
                        return verified_col
        
        # If dataframe provided, verify the column exists
        if dataframe is not None:
            if default_col in dataframe.columns:
                return default_col
            # Try to find any similar column
            possible_cols = ['Name', 'LipidID', 'Lipid_ID', 'Metabolite', 'Compound', 'Feature ID']
            for col in possible_cols:
                if col in dataframe.columns:
                    return col
        
        return default_col
    
    def _is_lipid_feature_col(self, col):
        """Check if a column is a canonical lipid feature column."""
        if col is None:
            return False

        normalized = self._normalize_col(col)
        if not normalized:
            return False

        # Strict metadata matching only.
        # IMPORTANT: Do not match on generic "lipid" substring because
        # sample names commonly include it (e.g., "..._lipid_Neg").
        exact_features = {
            'lipidid', 'lipidname',
            'class', 'classname', 'subclass', 'superclass',
            'adduct', 'adduction', 'adductionion', 'adductiontype', 'adductionmode',
            'calcmz', 'obsmz', 'basert', 'obsrt', 'ppmdiff',
            'polarity', 'mz', 'retention', 'retentiontime', 'rt',
            'charge', 'mass', 'neutralmass'
        }

        prefix_features = (
            'lipidid', 'lipidname', 'class', 'subclass', 'superclass',
            'adduct', 'calcmz', 'obsmz', 'basert', 'obsrt',
            'ppmdiff', 'retention', 'polarity'
        )

        return normalized in exact_features or normalized.startswith(prefix_features)
    
    def _thread_safe_log(self, message):
        """Thread-safe logging to stats_log"""
        if getattr(self, '_shutting_down', False):
            return
        try:
            self.root.after(0, lambda: (not getattr(self, '_shutting_down', False)) and self.stats_log.insert(tk.END, message))
            self.root.after(0, lambda: (not getattr(self, '_shutting_down', False)) and self.stats_log.see(tk.END))
        except Exception:
            pass
    
    def _thread_safe_progress(self, message=None):
        """Thread-safe progress bar update"""
        if getattr(self, '_shutting_down', False):
            return
        if message:
            try:
                self.root.after(0, lambda: (not getattr(self, '_shutting_down', False)) and self.show_stats_progress(message))
            except Exception:
                pass
        else:
            try:
                self.root.after(0, lambda: (not getattr(self, '_shutting_down', False)) and self.hide_stats_progress())
            except Exception:
                pass

    def _thread_safe_progress_step(self, step, total, message):
        """Thread-safe wrapper to update determinate statistics progress."""
        # Cache total so show_stats_progress knows to use determinate mode
        if not hasattr(self, '_stats_total_steps') or self._stats_total_steps != total:
            self._stats_total_steps = total
        if getattr(self, '_shutting_down', False):
            return
        try:
            self.root.after(0, lambda: (not getattr(self, '_shutting_down', False)) and self.update_stats_progress(step, message))
        except Exception:
            pass
    
    def _clear_downstream_data(self):
        """Clear all downstream data (statistics, visualization, pathway) when new data is uploaded.
        This prevents old data from previous runs contaminating new analyses."""
        # Clear normalized data
        if hasattr(self, 'normalized_combined_df'):
            self.normalized_combined_df = None
        if hasattr(self, 'raw_normalized_combined_df'):
            self.raw_normalized_combined_df = None
        if hasattr(self, 'optional_processing_report'):
            self.optional_processing_report = {}
        if hasattr(self, 'optional_processing_applied'):
            self.optional_processing_applied = False
        if hasattr(self, 'last_pca_outlier_scores'):
            self.last_pca_outlier_scores = None
        
        # Clear statistical results
        self._clear_statistical_results()
        
        # Clear pathway data
        if hasattr(self, 'pathway_filtered_metabolites_data'):
            self.pathway_filtered_metabolites_data = None
        if hasattr(self, 'pathway_filtered_pathways_data'):
            self.pathway_filtered_pathways_data = None
        if hasattr(self, 'pathway_original_metabolites_data'):
            self.pathway_original_metabolites_data = None
        if hasattr(self, 'pathway_original_pathways_data'):
            self.pathway_original_pathways_data = None
        
        # Clear visualization imported data
        if hasattr(self, 'imported_complete_df'):
            self.imported_complete_df = None
        if hasattr(self, 'imported_lipid_class_df'):
            self.imported_lipid_class_df = None
    
    def _clear_statistical_results(self):
        """Clear statistical test results to prevent visualization from using stale data."""
        if hasattr(self, 'statistical_test_results'):
            self.statistical_test_results = None
        if hasattr(self, 'statistical_test_results_class'):
            self.statistical_test_results_class = None
    
    def _re_enable_stats_button(self):
        """Re-enable the statistics button after completion"""
        def enable_button():
            for widget in self.frame.winfo_children():
                if isinstance(widget, tk.Frame):
                    for child in widget.winfo_children():
                        if isinstance(child, tk.Button) and "Run Normalization" in child.cget('text'):
                            child.config(state='normal')
        self.root.after(0, enable_button)
    
    def _run_statistics_pipeline_threaded(self):
        """Enhanced normalization pipeline with proper ion ordering and column name cleaning."""
        try:
            # Call the original pipeline logic
            self._run_statistics_pipeline_core()
        except Exception as e:
            self._thread_safe_log(f"❌ CRITICAL ERROR: {str(e)}\n")
            self._thread_safe_log(f"Stack trace: {traceback.format_exc()}\n")
        finally:
            self._thread_safe_progress()  # Hide progress bar
            self._re_enable_stats_button()  # Re-enable button
    
    def _run_statistics_pipeline_core(self):
        """Core statistics pipeline logic (runs in thread)"""
        from main_script.metabolite_statistics_analysis import normalize_dataframe, detect_feature_and_sample_columns, clean_sample_column_names
        
        # Show progress bar (thread-safe UI update)
        total_steps = 10
        self._thread_safe_progress_step(1, total_steps, "Initializing normalization pipeline...")
        
        # Determine current mode
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        self._thread_safe_log(f"\n🔄 Starting normalization pipeline in {mode.upper()} mode...\n")
        
        # STRICT VALIDATION: Only proceed with required data
        if not hasattr(self, 'memory_store') or not isinstance(self.memory_store, dict):
            error_msg = 'Memory store not available. Please complete Data Cleaning and ID Annotation first.'
            self._thread_safe_log(f"❌ ERROR: {error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('No Memory Store', error_msg))
            return
        
        # � CUSTOM MODE: Check if preprocessed combined data was loaded (custom or backdoor)
        is_custom_mode = self.memory_store.get('is_preprocessed_custom', False) or self.memory_store.get('is_preprocessed_backdoor', False)
        if is_custom_mode:
            mode_label = "CUSTOM" if self.memory_store.get('is_preprocessed_custom', False) else "PREPROCESSED"
            self._thread_safe_log(f"\n📊 {mode_label} MODE ACTIVE: Using preprocessed combined data...\n")
            self._run_preprocessed_pipeline()
            return
        
        # Check for required DataFrames based on mode
        if mode == 'lipid':
            has_pos = any(key in self.memory_store and self.memory_store[key] is not None 
                         for key in ['pos_lipid_df', 'lipid_pos_df'])
            has_neg = any(key in self.memory_store and self.memory_store[key] is not None 
                         for key in ['neg_lipid_df', 'lipid_neg_df'])
            data_type = 'lipid data'
            import_instruction = 'Import Excel file with Positive_Lipids/Negative_Lipids sheets'
        else:
            has_pos = any(key in self.memory_store and self.memory_store[key] is not None 
                            for key in ['pos_id_df', 'clean_pos_id_df'])
            has_neg = any(key in self.memory_store and self.memory_store[key] is not None 
                            for key in ['neg_id_df', 'clean_neg_id_df'])
            data_type = 'metabolite data'
            import_instruction = 'Complete ID Annotation in Tab 2 first, OR\n2. Import Excel file with Pos_id/Neg_id sheets'
        
        if not has_pos and not has_neg:
            error_msg = (f'❌ REQUIREMENT: {data_type} required for statistics.\n\n'
                        f'Found in memory: No data for {mode} mode\n\n'
                        f'SOLUTIONS:\n'
                        f'1. {import_instruction} using Import button\n\n')
            self._thread_safe_log(f"{error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('Data Required', error_msg))
            return
        
        self._thread_safe_progress_step(2, total_steps, "Loading polarity datasets...")
        pos_df, neg_df = self._get_pos_neg_for_stats()
        if pos_df is None and neg_df is None:
            error_msg = (f'Unexpected error: {data_type} detected but could not be loaded.\n\n'
                        'Please try importing your Excel file again.')
            self._thread_safe_log(f"❌ ERROR: {error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('Data Loading Error', error_msg))
            return
        
        # Process each polarity separately with proper column cleaning
        normalized_frames = []
        raw_pre_normalization_exports = {}
        norm_method = self.stat_norm_method.get()
        optional_4b_enabled = any([
            bool(self.enable_variability_filter_var.get()) if hasattr(self, 'enable_variability_filter_var') else False,
            bool(self.enable_imputation_var.get()) if hasattr(self, 'enable_imputation_var') else False,
            bool(self.enable_pca_outlier_var.get()) if hasattr(self, 'enable_pca_outlier_var') else False,
        ])
        polarity_processing_reports = {}
        
        # Initialize individual polarity storage
        normalized_positive_df = None
        normalized_negative_df = None
        raw_positive_class_df = None
        raw_negative_class_df = None
        class_column_mappings = {}
        
        step_index = 2
        for polarity, df in [('positive', pos_df), ('negative', neg_df)]:
            if df is None:
                continue
                
            label = polarity.capitalize()
            step_index += 1
            self._thread_safe_progress_step(step_index, total_steps, f"Detecting sample columns ({polarity})...")
            
            # Check if we have verified assignments for this polarity
            verified_sample_cols = None
            verified_feature_id = None
            
            if mode == 'lipid':
                if polarity == 'positive' and hasattr(self, 'verified_pos_lipid_sample_cols'):
                    verified_sample_cols = self.verified_pos_lipid_sample_cols
                    verified_feature_id = getattr(self, 'verified_pos_lipid_assignments', {}).get('Feature ID')
                elif polarity == 'negative' and hasattr(self, 'verified_neg_lipid_sample_cols'):
                    verified_sample_cols = self.verified_neg_lipid_sample_cols
                    verified_feature_id = getattr(self, 'verified_neg_lipid_assignments', {}).get('Feature ID')
            else:
                if polarity == 'positive' and hasattr(self, 'verified_pos_sample_cols'):
                    verified_sample_cols = self.verified_pos_sample_cols
                    verified_feature_id = getattr(self, 'verified_pos_assignments', {}).get('Feature ID')
                elif polarity == 'negative' and hasattr(self, 'verified_neg_sample_cols'):
                    verified_sample_cols = self.verified_neg_sample_cols
                    verified_feature_id = getattr(self, 'verified_neg_assignments', {}).get('Feature ID')
            
            # Use verified sample columns if available, otherwise auto-detect
            if verified_sample_cols:
                self._thread_safe_log(f'{label}: Using {len(verified_sample_cols)} verified sample columns.\n')
                sample_cols = [col for col in verified_sample_cols if col in df.columns]
                
                # For feature cols, if we have verified Feature ID, use it
                if verified_feature_id and verified_feature_id in df.columns:
                    # Feature columns = all non-sample numeric columns + Feature ID
                    feature_cols = [verified_feature_id]
                    for col in df.columns:
                        if col != verified_feature_id and col not in sample_cols:
                            if not pd.api.types.is_numeric_dtype(df[col]) or self._is_lipid_feature_col(col) if mode == 'lipid' else True:
                                feature_cols.append(col)
                else:
                    # No verified feature ID, use all non-sample columns as features
                    feature_cols = [col for col in df.columns if col not in sample_cols]
            else:
                # Auto-detect columns based on mode
                if mode == 'lipid':
                    # For lipid data, use robust feature detection
                    feature_cols = []
                    sample_cols = []
                    for col in df.columns:
                        # Use normalized matching for lipid features
                        if self._is_lipid_feature_col(col):
                            feature_cols.append(col)
                        elif pd.api.types.is_numeric_dtype(df[col]):
                            sample_cols.append(col)
                else:
                    # For metabolite data, use existing detection
                    feature_cols, sample_cols = detect_feature_and_sample_columns(df)
                
                self._thread_safe_log(f'{label}: Auto-detected {len(feature_cols)} feature cols, {len(sample_cols)} sample cols.\n')
            
            self._thread_safe_log(f'{label}: {len(feature_cols)} feature cols, {len(sample_cols)} sample cols confirmed.\n')
            
            if not sample_cols:
                self._thread_safe_log(f'{label}: No sample columns detected, skipping.\n')
                continue
            
            # CRITICAL: Ensure sample_cols only contains columns from THIS polarity's DataFrame
            # Filter out any columns that don't actually exist in the current df
            df_cols_set = set(df.columns)
            sample_cols = [col for col in sample_cols if col in df_cols_set]
            
            if not sample_cols:
                self._thread_safe_log(f'{label}: No valid sample columns in this polarity, skipping.\n')
                continue
            
            # Clean sample column names before normalization
            self._thread_safe_progress(f"Processing {label} data...")
            self._thread_safe_progress_step(step_index, total_steps, f"Cleaning column names ({polarity})...")
            column_mapping = clean_sample_column_names(sample_cols, polarity)
            cleaned_samples = [column_mapping[col] for col in sample_cols]
            self._current_sample_column_aliases = column_mapping.copy()
            
            # Log cleaning details
            renamed_count = sum(1 for orig, clean in column_mapping.items() if orig != clean)
            self._thread_safe_log(f'{label}: Cleaned {renamed_count}/{len(sample_cols)} sample column names.\n')
            
            # Show sample of name changes for debugging
            if renamed_count > 0:
                sample_changes = [(orig, clean) for orig, clean in column_mapping.items() if orig != clean][:5]
                for orig, clean in sample_changes:
                    self._thread_safe_log(f'  📝 {orig} → {clean}\n')
                if renamed_count > 5:
                    self._thread_safe_log(f'  ... and {renamed_count - 5} more\n')
            
            # CRITICAL: Check if pre-normalization filtering should be applied
            # Get filter timing preference (before or after normalization)
            filter_timing = self.filter_timing_var.get() if hasattr(self, 'filter_timing_var') else 'before'
            
            # Parse current group assignments from sample_group_vars
            # IMPORTANT: Only include samples from THIS polarity for per-polarity filtering
            temp_group_map = {}
            if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                # Get current group assignments (using labels not IDs)
                full_assignments = self._parse_group_assignments()
                
                # CRITICAL FIX: Only map samples that exist in THIS polarity's DataFrame
                # We must check against the actual column names in the current df, not cleaned names
                # This prevents samples from other polarities being included
                df_cols_set = set(df.columns)
                
                for orig_col in sample_cols:
                    # Only include if this column exists in the current polarity's DataFrame
                    if orig_col in df_cols_set:
                        cleaned_col = column_mapping[orig_col]
                        # Look up group assignment using both cleaned name AND original name
                        # (group assignments may be stored with either name depending on when they were set)
                        if cleaned_col in full_assignments:
                            temp_group_map[orig_col] = full_assignments[cleaned_col]
                        elif orig_col in full_assignments:
                            temp_group_map[orig_col] = full_assignments[orig_col]
            
            # Snapshot the input before any pre-normalization filtering or imputation.
            raw_input_ions_df = df.reset_index(drop=True).copy()

            # Workflow order when imputation is enabled:
            # 1) filtering, 2) normalization, 3) variance filter, 4) imputation, 5) PCA.
            if self._is_imputation_enabled():
                if temp_group_map and len(temp_group_map) == len(sample_cols):
                    imp_min_pct = self._get_imputation_min_group_percent()
                    imp_scope = self._get_imputation_prefilter_scope()
                    scope_text = 'at least one group' if imp_scope == 'per_group' else 'all groups'
                    self._thread_safe_log(f'\n🔍 Applying imputation pre-filter BEFORE normalization ({label})...\n')
                    self._thread_safe_log(f'   Rule: keep rows with at least {imp_min_pct:.1f}% valid values in {scope_text}\n')
                    filtered_df, pre_imp_report = self._apply_imputation_prefilter(
                        df,
                        sample_cols,
                        temp_group_map,
                        imp_min_pct,
                        scope=imp_scope,
                    )
                    self._thread_safe_log(
                        f"   ✅ Imputation pre-filter (pre-normalization): removed {pre_imp_report.get('removed', 0)} rows, "
                        f"kept {pre_imp_report.get('kept', 0)}.\n"
                    )
                    df = filtered_df
                else:
                    self._thread_safe_log(
                        f'⚠️ {label}: Imputation pre-filter before normalization skipped due to incomplete group assignments '
                        f'({len(temp_group_map)}/{len(sample_cols)}).\n'
                    )
            # Only apply original pre-normalization filtering if:
            # 1. User selected "before" timing
            # 2. All samples have group assignments
            # 3. Imputation workflow is NOT enabled
            elif filter_timing == 'before' and temp_group_map and len(temp_group_map) == len(sample_cols):
                # All samples have group assignments - apply pre-normalization filtering
                from main_script.metabolite_statistics_analysis import apply_min_group_size_filter
                
                # Get threshold settings with safe fallbacks
                min_type = self.min_samples_type_var.get() if hasattr(self, 'min_samples_type_var') else 'absolute'
                
                # Safe getter for min_count - handle empty string
                try:
                    min_count_str = self.min_samples_per_group_var.get() if hasattr(self, 'min_samples_per_group_var') else '2'
                    min_count = int(min_count_str) if min_count_str else 2
                except (ValueError, _tkinter.TclError):
                    min_count = 2
                
                # Safe getter for min_percent - handle empty string
                try:
                    min_percent_str = self.min_samples_percent_var.get() if hasattr(self, 'min_samples_percent_var') else '50.0'
                    min_percent = float(min_percent_str) if min_percent_str else 50.0
                except (ValueError, _tkinter.TclError):
                    min_percent = 50.0
                
                self._thread_safe_log(f'\n🔍 Applying pre-normalization filter ({label})...\n')
                self._thread_safe_log(f'   ⚡ Per-polarity filtering: Only {label} samples used for threshold calculation\n')
                self._thread_safe_log(f'   🔢 Sample count for {label}: {len(sample_cols)} samples in total\n')
                self._thread_safe_log(f'   📋 Group mapping size: {len(temp_group_map)} sample assignments\n')
                
                # Log group distribution
                from collections import Counter
                group_counts = Counter(temp_group_map.values())
                self._thread_safe_log(f'   👥 Groups detected: {", ".join([f"{g}={c}" for g, c in sorted(group_counts.items())])}\n')
                
                if min_type == 'percentage':
                    self._thread_safe_log(f'   Threshold: {min_percent}% of group size (per {label} samples only)\n')
                else:
                    self._thread_safe_log(f'   Threshold: {min_count} samples per group (per {label} samples only)\n')
                
                filtered_df, group_stats = apply_min_group_size_filter(
                    df, sample_cols, temp_group_map, 
                    min_count, min_type, min_percent
                )
                
                # Log results per group
                self._thread_safe_log(f'\n📊 Pre-normalization filtering results:\n')
                for grp in sorted(group_stats.keys()):
                    stats = group_stats[grp]
                    total_samples = stats['total_samples']
                    threshold = stats['threshold']
                    before = stats['before']
                    after = stats['after']
                    removed = stats['removed']
                    self._thread_safe_log(f'   • {grp} (n={total_samples}, threshold={threshold}):\n')
                    self._thread_safe_log(f'     - {after} features retained')
                    if removed > 0:
                        self._thread_safe_log(f' ({removed} removed)')
                    self._thread_safe_log(f'\n')
                
                df = filtered_df
                
                # Remove rows that are completely zero across all sample columns
                rows_before = len(df)
                sample_data = df[sample_cols].apply(pd.to_numeric, errors='coerce')
                # A row is kept if it has at least one non-zero, non-NaN value across all samples
                rows_with_data = (sample_data != 0).any(axis=1) & (~sample_data.isna().all(axis=1))
                df = df[rows_with_data].copy()
                rows_after = len(df)
                rows_removed = rows_before - rows_after
                
                if rows_removed > 0:
                    self._thread_safe_log(f'🗑️  Removed {rows_removed} features with no valid data across all groups.\n')
                    self._thread_safe_log(f'   Remaining: {rows_after} features\n')
                
                self._thread_safe_log(f'{label}: Pre-normalization filtering complete.\n')
            elif filter_timing == 'after':
                # User chose to apply filtering during statistics (after normalization)
                self._thread_safe_log(f'ℹ️  {label}: Filtering will be applied during statistical tests (after normalization).\n')
                self._thread_safe_log(f'   All data will be normalized without pre-filtering.\n')
            else:
                # Log more detailed information about missing assignments (only for "before" mode)
                if not hasattr(self, 'sample_group_vars') or not self.sample_group_vars:
                    self._thread_safe_log(f'⚠️ {label}: No group assignment controls found.\n')
                    self._thread_safe_log(f'   This should not happen - please report this issue.\n')
                elif not temp_group_map:
                    self._thread_safe_log(f'⚠️ {label}: No samples have group assignments - skipping pre-normalization filtering.\n')
                    self._thread_safe_log(f'   Found {len(sample_cols)} samples but 0 with group assignments.\n')
                    self._thread_safe_log(f'   Hint: Use "Auto-Assign by Pattern" or manually assign groups.\n')
                else:
                    assigned = len(temp_group_map)
                    total = len(sample_cols)
                    self._thread_safe_log(f'⚠️ {label}: Incomplete group assignments ({assigned}/{total}) - skipping pre-normalization filtering.\n')
                    self._thread_safe_log(f'   All samples must have group assignments for filtering to be applied.\n')
                    # Show which samples are missing assignments
                    missing = [column_mapping[col] for col in sample_cols if col not in temp_group_map]
                    if missing:
                        self._thread_safe_log(f'   Missing assignments for: {", ".join(missing[:5])}')
                        if len(missing) > 5:
                            self._thread_safe_log(f' ... and {len(missing)-5} more')
                        self._thread_safe_log(f'\n')

            # Capture the filtered pre-imputation state for debugging.
            raw_filtered_ions_df = df.reset_index(drop=True).copy()

            raw_pre_normalization_exports[polarity] = {
                'ions_input': raw_input_ions_df,
                'ions': raw_filtered_ions_df,
                'class_input': None,
                'class': None,
                'imputed': None,
            }
            
            # Store raw filtered ion tables for Step 7 (class grouping will happen AFTER dedup and split)
            if mode == 'lipid':
                if polarity == 'positive':
                    raw_positive_class_df = raw_filtered_ions_df.reset_index(drop=True).copy()
                elif polarity == 'negative':
                    raw_negative_class_df = raw_filtered_ions_df.reset_index(drop=True).copy()
            
            # Apply pre-normalization imputation if enabled
            if self._is_imputation_before_normalization_enabled():
                from main_script.metabolite_statistics_analysis import apply_imputation
                try:
                    imp_method = self.imputation_before_method_var.get().strip().lower() if hasattr(self, 'imputation_before_method_var') else 'half_min'
                    try:
                        knn_k = int(self.imputation_before_knn_neighbors_var.get()) if hasattr(self, 'imputation_before_knn_neighbors_var') else 5
                    except Exception:
                        knn_k = 5
                    
                    self._thread_safe_progress_step(step_index, total_steps, f"Imputing missing values BEFORE normalization ({label})...")
                    self._thread_safe_log(f'\n🔧 Applying pre-normalization imputation ({label})...\n')
                    self._thread_safe_log(f'   Method: {imp_method}\n')
                    
                    # Apply imputation to individual data
                    df, imp_report = apply_imputation(
                        df,
                        sample_cols,
                        method=imp_method,
                        knn_neighbors=knn_k,
                        debug=True,
                        log_fn=self._thread_safe_log,
                    )
                    raw_pre_normalization_exports[polarity]['imputed'] = df.reset_index(drop=True).copy()
                    self._thread_safe_log(
                        f"   ✅ Individual data imputation: filled {imp_report.get('imputed_cells', 0)} cells; "
                        f"missing_after={imp_report.get('missing_after', 'NA')}; "
                        f"zeros_after={imp_report.get('zeros_after', 'NA')}\n"
                    )
                except Exception as e:
                    self._thread_safe_log(f'⚠️ Pre-normalization imputation failed: {e}\n')

            # Ion normalization happens next (class already created from raw filtered data above)
            
            # Normalize the data
            try:
                self._thread_safe_progress_step(step_index, total_steps, f"Normalizing {polarity} data ({norm_method})...")
                
                # Validate normalization requirements
                is_feature_name = None
                qc_sample_cols = None
                
                # Check if Internal Standard normalization is selected
                if norm_method.lower() in ('is', 'internal_standard', 'istd'):
                    # Get Internal Standard from verified assignments
                    verified_is = None
                    if mode == 'lipid':
                        if polarity == 'positive' and hasattr(self, 'verified_pos_lipid_assignments'):
                            verified_is = self.verified_pos_lipid_assignments.get('Internal Standard')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_lipid_assignments'):
                            verified_is = self.verified_neg_lipid_assignments.get('Internal Standard')
                    else:
                        if polarity == 'positive' and hasattr(self, 'verified_pos_assignments'):
                            verified_is = self.verified_pos_assignments.get('Internal Standard')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_assignments'):
                            verified_is = self.verified_neg_assignments.get('Internal Standard')
                    
                    if not verified_is:
                        # IS normalization requires Internal Standard to be identified
                        error_msg = (f'❌ Internal Standard (IS) Normalization Requires IS Feature\n\n'
                                   f'{label} polarity: No Internal Standard identified.\n\n'
                                   f'SOLUTION:\n'
                                   f'1. Click "Verify Columns" button\n'
                                   f'2. In the column assignment dialog, identify your Internal Standard feature\n'
                                   f'   (Look for compounds like ISTD, Internal_Standard, or reference standards)\n'
                                   f'3. Assign it as "Internal Standard" type\n'
                                   f'4. Re-run normalization\n\n'
                                   f'NOTE: Internal Standard normalization divides all features by the IS intensity.')
                        self._thread_safe_log(f"\n{error_msg}\n")
                        self.root.after(0, lambda msg=error_msg: messagebox.showerror(
                            'Internal Standard Required', msg))
                        return
                    else:
                        is_feature_name = verified_is
                        self._thread_safe_log(f'{label}: Using Internal Standard: {is_feature_name}\n')
                
                # Check if LOESS QC normalization is selected
                elif norm_method.lower() in ('loess_qc', 'loess', 'qc_correction'):
                    # Get QC columns from verified assignments
                    verified_qc = None
                    if mode == 'lipid':
                        if polarity == 'positive' and hasattr(self, 'verified_pos_lipid_assignments'):
                            verified_qc = self.verified_pos_lipid_assignments.get('_qc_columns')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_lipid_assignments'):
                            verified_qc = self.verified_neg_lipid_assignments.get('_qc_columns')
                    else:
                        if polarity == 'positive' and hasattr(self, 'verified_pos_assignments'):
                            verified_qc = self.verified_pos_assignments.get('_qc_columns')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_assignments'):
                            verified_qc = self.verified_neg_assignments.get('_qc_columns')
                    
                    if verified_qc and len(verified_qc) >= 3:
                        qc_sample_cols = verified_qc
                        self._thread_safe_log(f'{label}: Using {len(qc_sample_cols)} QC samples for drift correction\n')
                    else:
                        # Warn but allow fallback to auto-detection
                        self._thread_safe_log(f'⚠️ {label}: No QC samples identified in verified columns.\n')
                        self._thread_safe_log(f'   LOESS correction will attempt auto-detection or use fallback.\n')
                        self._thread_safe_log(f'   For best results, verify columns and mark QC samples.\n')
                
                # Perform normalization with appropriate parameters
                norm_df = normalize_dataframe(
                    df, sample_cols, norm_method,
                    is_feature_name=is_feature_name,
                    qc_sample_cols=qc_sample_cols
                )
                
                # Apply same normalization to class data if it exists
                # DEFERRED: Class normalization will happen after merge/split workflow
                # Only normalize class data if NOT in lipid mode (metabolites normalize here)
            except Exception as e:
                messagebox.showerror('Normalization Error', f'{label} failed: {e}')
                return
            self._thread_safe_log(f'{label}: Applied normalization ({norm_method}).\n')
            
            # Rename sample columns to cleaned names BEFORE normality testing
            rename_dict = {orig: clean for orig, clean in column_mapping.items()}
            class_column_mappings[polarity] = rename_dict.copy()
            norm_df = norm_df.rename(columns=rename_dict)
            cleaned_sample_cols = [rename_dict[col] for col in sample_cols]
            self._current_sample_column_aliases = rename_dict.copy()
            
            # Update temp_group_map to use cleaned column names
            if temp_group_map:
                cleaned_group_map = {}
                for orig_col, group in temp_group_map.items():
                    cleaned_col = rename_dict.get(orig_col, orig_col)
                    cleaned_group_map[cleaned_col] = group
                temp_group_map = cleaned_group_map

            # Step 4b optional processing must run before normality testing.
            if optional_4b_enabled:
                norm_df, pol_report = self._apply_optional_post_normalization_processing(
                    norm_df.reset_index(drop=True),
                    mode,
                    verified_sample_cols=cleaned_sample_cols
                )
                polarity_processing_reports[polarity] = pol_report
                cleaned_sample_cols = [c for c in cleaned_sample_cols if c in norm_df.columns]
            
            # Perform normality tests (Shapiro-Wilk only) after normalization
            # Tests PER FEATURE across all samples, not per sample
            # QQ plots will be generated during export only (max 8 plots)
            if norm_method != 'none':
                try:
                    from scipy.stats import shapiro
                    import numpy as _np
                    import pandas as _pd
                except Exception:
                    self._thread_safe_log(f'Normality test dependencies (scipy) missing; skipping tests.\n')
                    normality_results = None
                else:
                    # Prepare results list
                    records = []
                    def _clean_vals(series):
                        vals = _pd.to_numeric(series, errors='coerce')
                        vals = vals.dropna()
                        return vals[vals != 0]
                    
                    # Test each feature across all sample columns
                    # Find feature ID column
                    metabolite_id_col = None
                    for col_name in ['Name', 'Metabolite', 'Feature ID', 'metabolite_id']:
                        if col_name in norm_df.columns:
                            metabolite_id_col = col_name
                            break
                    
                    # Determine sample columns to use for testing
                    test_sample_cols = cleaned_sample_cols
                    
                    # Store feature indices and results for QQ plot filtering
                    normal_metabolites = []  # Will store (idx, feature_id, values)
                    not_normal_metabolites = []  # Will store (idx, feature_id, values)
                    
                    total_metabolites = len(norm_df)
                    callback_freq = max(1, total_metabolites // 100)
                    
                    for seq_idx, (idx, row) in enumerate(norm_df.iterrows(), start=1):
                        if seq_idx % callback_freq == 0 or seq_idx == 1 or seq_idx == total_metabolites:
                            self._thread_safe_log(f'Testing normality: {seq_idx}/{total_metabolites}\n')
                        
                        # Get metabolite identifier
                        if metabolite_id_col and metabolite_id_col in norm_df.columns:
                            metabolite_id = row[metabolite_id_col]
                        else:
                            metabolite_id = f'Feature_{idx}'
                        
                        # Get values for this feature across all samples
                        vals_list = []
                        for col in test_sample_cols:
                            if col in norm_df.columns:
                                vals_list.append(_clean_vals(row[[col]]))
                        
                        if not vals_list:
                            continue
                        
                        all_vals = _pd.concat(vals_list, ignore_index=True)
                        n = len(all_vals)
                        
                        if n < 3:
                            records.append({
                                'Metabolite': metabolite_id,
                                'N': n,
                                'Shapiro_p': _np.nan,
                                'Is_Normal': 'TooFew'
                            })
                            continue
                        
                        # Shapiro-Wilk test (cap at 5000 per SciPy recommendation)
                        shapiro_sample = all_vals.sample(min(n, 5000), random_state=42) if n > 5000 else all_vals
                        try:
                            sh_stat, sh_p = shapiro(shapiro_sample)
                        except Exception:
                            sh_p = _np.nan
                        
                        # Determine normality: p > 0.05 is normal
                        is_norm = 'Yes' if (not _np.isnan(sh_p) and sh_p >= 0.05) else 'No' if not _np.isnan(sh_p) else 'Unknown'
                        
                        records.append({
                            'Metabolite': metabolite_id,
                            'N': n,
                            'Shapiro_p': sh_p,
                            'Is_Normal': is_norm
                        })
                        
                        # Store for QQ plot selection (keep track of actual values for later)
                        if is_norm == 'Yes':
                            normal_metabolites.append((idx, metabolite_id, all_vals.values))
                        elif is_norm == 'No':
                            not_normal_metabolites.append((idx, metabolite_id, all_vals.values))
                    
                    normality_results = _pd.DataFrame(records) if records else _pd.DataFrame(columns=['Metabolite','N','Shapiro_p','Is_Normal'])
                    if not hasattr(self, 'normality_test_results'):
                        self.normality_test_results = {}
                    self.normality_test_results[polarity] = normality_results
                    
                    # Store selected metabolites for QQ plot generation (max 2 normal + 2 not-normal per polarity)
                    if not hasattr(self, 'normality_test_targets'):
                        self.normality_test_targets = {}
                    
                    selected_plots = []
                    # Select up to 2 normal metabolites
                    if len(normal_metabolites) > 0:
                        import random
                        selected_normal = random.sample(normal_metabolites, min(2, len(normal_metabolites)))
                        selected_plots.extend([('NORMAL', m_id, vals) for _, m_id, vals in selected_normal])
                    # Select up to 2 not-normal metabolites
                    if len(not_normal_metabolites) > 0:
                        import random
                        selected_not_normal = random.sample(not_normal_metabolites, min(2, len(not_normal_metabolites)))
                        selected_plots.extend([('NOT_NORMAL', m_id, vals) for _, m_id, vals in selected_not_normal])
                    
                    self.normality_test_targets[polarity] = {
                        'plot_data': selected_plots,
                        'normal_count': len(normal_metabolites),
                        'not_normal_count': len(not_normal_metabolites)
                    }
                    
                    # Log summary of normality test results
                    if not normality_results.empty:
                        normal_count = len(normality_results[normality_results['Is_Normal'] == 'Yes'])
                        not_normal_count = len(normality_results[normality_results['Is_Normal'] == 'No'])
                        too_few_count = len(normality_results[normality_results['Is_Normal'] == 'TooFew'])
                        total = len(normality_results)
                        self._thread_safe_log(f'📊 Normality Test Summary ({polarity.capitalize()}):\n')
                        # Percentages
                        pct_norm = (normal_count / total * 100.0) if total else 0.0
                        pct_not = (not_normal_count / total * 100.0) if total else 0.0
                        self._thread_safe_log(f'  ✓ Normal: {normal_count}/{total} ({pct_norm:.1f}%)\n')
                        self._thread_safe_log(f'  ✗ Not Normal: {not_normal_count}/{total} ({pct_not:.1f}%)\n')
                        if too_few_count > 0:
                            self._thread_safe_log(f'  ⚠ Too Few Samples: {too_few_count}/{total}\n')
                        self._thread_safe_log(f'  📈 Selected {len(selected_plots)} metabolites for QQ plots (2 normal + 2 not-normal)\n')
                # end revamped block
            
            # NOTE: Column renaming already done above before normality testing
            # norm_df already has cleaned column names at this point
            
            # Ensure deduplication key presence
            id_col = self._get_verified_id_column(norm_df)
            
            if mode == 'lipid':
                # For lipids, use verified ID column as the key
                if id_col not in norm_df.columns:
                    self._thread_safe_log(f'{label}: Warning - {id_col} column missing, cannot deduplicate.\n')
            else:
                # For metabolites, ensure Name_Key for deduplication
                if 'Name_Key' not in norm_df.columns and id_col in norm_df.columns:
                    norm_df['Name_Key'] = norm_df[id_col].astype(str).str.strip().str.lower()
            
            # Store individual polarity DataFrames (without Polarity_Source column)
            if polarity == 'positive':
                normalized_positive_df = norm_df.copy()
                # NOTE: Class storage deferred - will store after merge/split workflow for lipids
                # Preserve current ordering here; do not resort by Area so GUI
                # doesn't override the chosen representative ordering upstream.
                # (Ordering will be applied after deduplication if desired.)
            elif polarity == 'negative':
                normalized_negative_df = norm_df.copy()
                # NOTE: Class storage deferred - will store after merge/split workflow for lipids
                # Preserve current ordering here; do not resort by Area.
            
            # Add polarity indicator for tracking in merged data
            norm_df['Polarity_Source'] = label
            
            normalized_frames.append(norm_df)
        
        if not normalized_frames:
            messagebox.showwarning('Normalization', 'No data normalized.')
            return

        # Step 4b already ran per-polarity before normality; only do cross-polarity alignment here.
        if optional_4b_enabled:
            self._thread_safe_log('\n🧪 Step 4b already applied before normality; aligning polarity outputs before merge...\n')

            # Keep polarity datasets aligned if PCA outlier removal dropped samples.
            removed_samples_union = set()
            for p_report in polarity_processing_reports.values():
                try:
                    removed_samples_union.update(p_report.get('pca_outlier', {}).get('removed_samples', []) or [])
                except Exception:
                    continue

            if removed_samples_union:
                removed_samples = sorted(removed_samples_union)
                if normalized_positive_df is not None:
                    normalized_positive_df = normalized_positive_df.drop(columns=removed_samples, errors='ignore')
                if normalized_negative_df is not None:
                    normalized_negative_df = normalized_negative_df.drop(columns=removed_samples, errors='ignore')
                self._thread_safe_log(f"   • PCA alignment: dropped {len(removed_samples)} sample columns across polarity sheets.\n")

            # Rebuild merge inputs from post-4b polarity datasets.
            normalized_frames = []
            if normalized_positive_df is not None:
                pos_merge_df = normalized_positive_df.copy()
                pos_merge_df['Polarity_Source'] = 'Positive'
                normalized_frames.append(pos_merge_df)
            if normalized_negative_df is not None:
                neg_merge_df = normalized_negative_df.copy()
                neg_merge_df['Polarity_Source'] = 'Negative'
                normalized_frames.append(neg_merge_df)
        
        # Merge the normalized frames
        self._thread_safe_progress_step(step_index + 1, total_steps, "Merging normalized datasets...")
        self.stats_log.insert(tk.END, f'Merging {len(normalized_frames)} normalized datasets...\n')
        combined = pd.concat(normalized_frames, ignore_index=True)
        
        # Preserve combined ordering from upstream processing (do not resort by Area here)
        if 'Area (Max.)' in combined.columns:
            self.stats_log.insert(tk.END, 'Left combined order intact (no resort by Area).\n')
        else:
            self.stats_log.insert(tk.END, 'No Area (Max.) column found for sorting; left as concatenated.\n')
        
        # For lipids: Create Area (Max.) column as sum of all numeric columns, sort descending, then deduplicate
        if mode == 'lipid':
            # Identify numeric columns (excluding feature columns that might be numeric)
            numeric_cols = []
            for col in combined.columns:
                if pd.api.types.is_numeric_dtype(combined[col]) and not self._is_lipid_feature_col(col):
                    numeric_cols.append(col)
            
            if numeric_cols:
                combined['Area (Max.)'] = combined[numeric_cols].sum(axis=1)
                self._thread_safe_log(f'Created Area (Max.) column as sum of {len(numeric_cols)} numeric columns.\n')
                
                # Sort by Area (Max.) descending
                combined = combined.sort_values('Area (Max.)', ascending=False).reset_index(drop=True)
                self._thread_safe_log('Sorted combined data by Area (Max.) descending.\n')
            else:
                self._thread_safe_log('Warning: No numeric columns found for Area (Max.) calculation.\n')
        
        # Remove duplicates based on mode
        id_col = self._get_verified_id_column(combined)
        dedup_key = id_col if mode == 'lipid' else 'Name_Key'
        def _normalize_formula_val(v):
            try:
                return str(v).replace(' ', '').strip()
            except Exception:
                return str(v)

        def _build_formula_id_key(row):
            if mode == 'lipid':
                return None  # handled by LipidID-based dedup below
            # Priority: LipidMaps_ID > PubChem_CID > KEGG_ID > HMDB_ID
            pri = ['LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID']
            formula = _normalize_formula_val(row.get('Formula', ''))
            if not formula:
                return None
            for col in pri:
                if col in combined.columns:
                    v = row.get(col)
                    if pd.notna(v) and str(v).strip() not in ('', 'nan', 'None'):
                        return f"{formula}|{col}|{str(v).strip()}"
            return None

        def _ion_rank(ion, polarity):
            pos_order = ['[M+H]+1','[M+2H]+2','[M+H-H2O]+1','[M+H+MeOH]+1','[M+FA+H]+1']
            neg_order = ['[M-H]-1','[M-2H]-2','[M-H-H2O]-1','[M-H-MeOH]-1','[M+FA-H]-1']
            # Cross-polarity tie rule: treat '[M+H]+1' and '[M-H]-1' as equal rank
            ion_s = str(ion).strip()
            if ion_s == '[M+H]+1' or ion_s == '[M-H]-1':
                return 0
            order = pos_order if str(polarity).strip() in ['+', 'positive'] else neg_order
            return order.index(ion_s) if ion_s in order else 999

        def _ms2_rank(ms2):
            s = str(ms2).strip().lower()
            if 'dda for preferred ion' in s:
                return 0
            if 'dda for other ion' in s or 'non-preferred' in s:
                return 1
            if 'no ms2' in s:
                return 5
            return 9

        if mode != 'lipid':
            # Keep-one best per Formula+ID across combined pos/neg without summing
            try:
                combined['_fid_key'] = combined.apply(_build_formula_id_key, axis=1)
                if combined['_fid_key'].notna().any():
                    before = len(combined)
                    kept_rows = []
                    for fid, grp in combined.groupby('_fid_key'):
                        if fid is None or (isinstance(fid, float) and pd.isna(fid)):
                            for _, r in grp.iterrows():
                                kept_rows.append(r)
                            continue
                        # Do not use RT window in statistics-phase dedup; consider all rows in the group
                        cand = grp.copy()
                        # Ranking: ion → MS2 → Area
                        ion_rank = cand.apply(lambda r: _ion_rank(r.get('Reference Ion', ''), r.get('Polarity', '+')), axis=1) if 'Reference Ion' in cand.columns else 999
                        ms2_rank = cand['MS2'].apply(_ms2_rank) if 'MS2' in cand.columns else 9
                        area = pd.to_numeric(cand.get('Area (Max.)', pd.Series([0]*len(cand))), errors='coerce').fillna(0)
                        cand = cand.assign(_ion_rank=ion_rank, _ms2_rank=ms2_rank, _area=area)
                        cand = cand.sort_values(['_ion_rank','_ms2_rank','_area'], ascending=[True, True, False])
                        kept_rows.append(cand.iloc[0].drop(labels=[c for c in ['_ion_rank','_ms2_rank','_area','_rtdev'] if c in cand.columns]))
                    combined = pd.DataFrame(kept_rows).reset_index(drop=True)
                    removed = before - len(combined)
                    self.stats_log.insert(tk.END, f"De-duplicated by Formula+ID: removed {removed} rows (kept best by ion/MS2/Area).\n")
                combined = combined.drop(columns=['_fid_key'], errors='ignore')
            except Exception as e:
                self.stats_log.insert(tk.END, f"Warning: Formula+ID dedup in merge failed: {e}\n")

            # Fallback/cleanup: basic Name_Key de-dup to avoid exact duplicates
            if dedup_key in combined.columns:
                before = len(combined)
                combined = combined.drop_duplicates(subset=[dedup_key], keep='first')
                removed = before - len(combined)
                self.stats_log.insert(tk.END, f'Removed {removed} duplicate {dedup_key} rows (kept first occurrence).\n')
            else:
                self.stats_log.insert(tk.END, f'Warning: {dedup_key} column not found, skipping deduplication.\n')
        else:
            # Lipid mode: dedup by LipidID as before
            if dedup_key in combined.columns:
                before = len(combined)
                combined = combined.drop_duplicates(subset=[dedup_key], keep='first')
                removed = before - len(combined)
                self.stats_log.insert(tk.END, f'Removed {removed} duplicate {dedup_key} rows (kept first occurrence).\n')
            else:
                self.stats_log.insert(tk.END, f'Warning: {dedup_key} column not found, skipping deduplication.\n')
      
        # Remove helper columns before storing final results
        columns_to_remove = ['Polarity_Source', 'Area (Max.)']
        for col in columns_to_remove:
            if col in combined.columns:
                combined = combined.drop(columns=[col])
                self.stats_log.insert(tk.END, f'Removed helper column: {col}\n')

        # Keep merged normalized copy. Optional Step 4b is applied per-polarity before merge.
        self.raw_normalized_combined_df = combined.reset_index(drop=True)
        processed_combined_df = self.raw_normalized_combined_df

        if optional_4b_enabled:
            self.optional_processing_report = {
                'applied': True,
                'scope': 'per_polarity_pre_merge',
                'polarity_reports': polarity_processing_reports,
            }
        else:
            self.optional_processing_report = {
                'applied': False,
                'scope': 'disabled',
                'polarity_reports': {},
            }
        self.optional_processing_applied = bool(optional_4b_enabled)

        # Store final (possibly processed) results used by downstream statistical tests.
        self.normalized_combined_df = processed_combined_df.reset_index(drop=True)
        self.normalized_positive_df = normalized_positive_df.reset_index(drop=True) if normalized_positive_df is not None else None
        self.normalized_negative_df = normalized_negative_df.reset_index(drop=True) if normalized_negative_df is not None else None
        
# ...existing code...

        class_generation_report = {
            'applied': False,
            'source': 'post_normalization_polarity_classes',
            'positive': {'applied': False},
            'negative': {'applied': False},
            'merged': {'applied': False},
        }

        # Process lipid class DataFrames if in lipid mode.
        # Class creation is intentionally deferred to this stage, after polarity
        # normalization and optional post-normalization processing have completed.
        if mode == 'lipid':
            self.normalized_positive_class_df = None
            self.normalized_negative_class_df = None
            self.normalized_combined_class_df = None

            self._thread_safe_log("\n🧪 Step 7: Create final class output from deduplicated lipid ions...\n")
            self._thread_safe_log("   1) Merge raw positive and negative ions with Polarity column.\n")
            self._thread_safe_log("   2) Create Max_Area, sort descending, deduplicate by Lipid ID.\n")
            self._thread_safe_log("   3) Split deduplicated ions back into Positive/Negative.\n")
            self._thread_safe_log("   4) Group each polarity into class tables.\n")
            self._thread_safe_log("   5) Normalize each polarity's class tables.\n")
            self._thread_safe_log("   6) Merge normalized class tables using MEAN for duplicate classes.\n")

            # Capture and freeze settings used for class processing to match ion workflow.
            step7_pre_imputation_enabled = self._is_imputation_before_normalization_enabled()
            step7_post_processing_enabled = bool(optional_4b_enabled)
            step7_post_imputation_enabled = bool(self.enable_imputation_var.get()) if hasattr(self, 'enable_imputation_var') else False

            step7_pre_imp_method = self.imputation_before_method_var.get().strip().lower() if hasattr(self, 'imputation_before_method_var') else 'half_min'
            try:
                step7_pre_knn_k = int(self.imputation_before_knn_neighbors_var.get()) if hasattr(self, 'imputation_before_knn_neighbors_var') else 5
            except Exception:
                step7_pre_knn_k = 5

            step7_post_imp_method = self.imputation_method_var.get().strip().lower() if hasattr(self, 'imputation_method_var') else 'half_min'
            try:
                step7_post_knn_k = int(self.knn_neighbors_var.get()) if hasattr(self, 'knn_neighbors_var') else 5
            except Exception:
                step7_post_knn_k = 5

            self._thread_safe_log(
                "  ℹ️ Class settings inherited from ion pipeline:\n"
                f"      - Pre-normalization imputation (Step 4a): {'ENABLED' if step7_pre_imputation_enabled else 'DISABLED'}"
                + (f" [method={step7_pre_imp_method}, knn={step7_pre_knn_k}]\n" if step7_pre_imputation_enabled else "\n")
                + f"      - Post-normalization processing (Step 4b): {'ENABLED' if step7_post_processing_enabled else 'DISABLED'}\n"
                + f"      - Post-normalization imputation (Step 4b): {'ENABLED' if step7_post_imputation_enabled else 'DISABLED'}"
                + (f" [method={step7_post_imp_method}, knn={step7_post_knn_k}]\n" if step7_post_imputation_enabled else "\n")
            )

            # Step 1: Merge raw ion tables with Polarity column
            ion_frames = []
            for polarity, raw_ion_df in [('positive', raw_positive_class_df), ('negative', raw_negative_class_df)]:
                if raw_ion_df is None or raw_ion_df.empty:
                    continue
                temp_ion_df = raw_ion_df.reset_index(drop=True).copy()
                temp_ion_df['Polarity'] = polarity.capitalize()
                ion_frames.append(temp_ion_df)
                self._thread_safe_log(f"  • {polarity.capitalize()}: {len(temp_ion_df)} raw ion rows ready for merge\n")

            if not ion_frames:
                self._thread_safe_log("No ion data to process.\n")
            else:
                # Step 1: Merge ions
                combined_ions = pd.concat(ion_frames, ignore_index=True)
                
                # Step 2: Create Max_Area, sort, deduplicate at ion level
                ion_id_col = self._get_verified_id_column(combined_ions)
                if ion_id_col not in combined_ions.columns:
                    for candidate in ['LipidID', 'Lipid_ID', 'LipID', 'ID', 'Feature ID']:
                        if candidate in combined_ions.columns:
                            ion_id_col = candidate
                            break

                numeric_cols = [
                    col for col in combined_ions.columns
                    if col not in {ion_id_col, 'Polarity', 'Max_Area', 'Area (Max.)'}
                    and pd.api.types.is_numeric_dtype(combined_ions[col])
                ]

                if numeric_cols:
                    combined_ions['Max_Area'] = combined_ions[numeric_cols].sum(axis=1)
                    self._thread_safe_log(f"  • Max_Area created from {len(numeric_cols)} numeric ion columns\n")

                if 'Max_Area' in combined_ions.columns:
                    combined_ions = combined_ions.sort_values('Max_Area', ascending=False).reset_index(drop=True)
                    self._thread_safe_log(f"  • Sorted {len(combined_ions)} ions descending by Max_Area\n")

                if ion_id_col in combined_ions.columns:
                    before = len(combined_ions)
                    combined_ions = combined_ions.drop_duplicates(subset=[ion_id_col], keep='first')
                    removed = before - len(combined_ions)
                    self._thread_safe_log(f"  • Deduped ions by {ion_id_col}: removed {removed} duplicate rows\n")
                else:
                    self._thread_safe_log("  ⚠️ Could not find a Lipid ID column for ion deduplication; skipping dedup step.\n")

                # Step 3: Split back to Positive and Negative ions
                pos_ions_df = combined_ions[combined_ions['Polarity'].astype(str).str.lower() == 'positive'].copy() if 'Polarity' in combined_ions.columns else combined_ions.copy()
                neg_ions_df = combined_ions[combined_ions['Polarity'].astype(str).str.lower() == 'negative'].copy() if 'Polarity' in combined_ions.columns else pd.DataFrame(columns=combined_ions.columns)

                self._thread_safe_log(f"  • Split back into ions: Positive={len(pos_ions_df)}, Negative={len(neg_ions_df)}\n")

                # Step 4: Group each polarity into class tables
                pos_class_df = None
                neg_class_df = None

                if not pos_ions_df.empty:
                    pos_class_df, pos_class_report = self._build_lipid_class_dataframe(pos_ions_df)
                    if pos_class_df is not None and not pos_class_df.empty:
                        self._thread_safe_log(f"  • Positive ions grouped into {len(pos_class_df)} classes\n")
                    else:
                        self._thread_safe_log(f"  ⚠️ Positive class grouping failed: {pos_class_report.get('reason', 'unknown')}\n")

                if not neg_ions_df.empty:
                    neg_class_df, neg_class_report = self._build_lipid_class_dataframe(neg_ions_df)
                    if neg_class_df is not None and not neg_class_df.empty:
                        self._thread_safe_log(f"  • Negative ions grouped into {len(neg_class_df)} classes\n")
                    else:
                        self._thread_safe_log(f"  ⚠️ Negative class grouping failed: {neg_class_report.get('reason', 'unknown')}\n")

                # Step 5: Normalize each polarity's class tables (with imputation if enabled)
                def _normalize_class_frame(class_df: pd.DataFrame, polarity: str) -> pd.DataFrame | None:
                    if class_df is None or class_df.empty:
                        return None

                    # Identify sample columns (exclude n_lipids, feature columns, metadata)
                    metadata_cols = {'Lipid_Class', 'Class', 'Class_name', 'Max_Area', 'Area (Max.)', 'n_lipids', 'n_lipid', 'CalcMz', 'BaseRt', 'Rt', 'MZ', 'Mz'}
                    sample_cols = [
                        col for col in class_df.columns
                        if pd.api.types.is_numeric_dtype(class_df[col]) and col not in metadata_cols
                    ]

                    if not sample_cols:
                        return class_df.copy()

                    working_df = class_df.copy()

                    # Step 4a: APPLY pre-normalization imputation if enabled.
                    if step7_pre_imputation_enabled:
                        try:
                            from main_script.metabolite_statistics_analysis import apply_imputation
                            
                            self._thread_safe_log(f"    • {polarity.capitalize()} classes: Applying pre-normalization imputation ({step7_pre_imp_method}) to {len(class_df)} class rows...\n")
                            working_df, imp_report = apply_imputation(
                                working_df,
                                sample_cols,
                                method=step7_pre_imp_method,
                                knn_neighbors=step7_pre_knn_k,
                                debug=True,
                                log_fn=self._thread_safe_log,
                            )
                            if imp_report:
                                imputed_count = imp_report.get('imputed_count', 0)
                                self._thread_safe_log(f"    ✅ {polarity.capitalize()} classes: Imputation complete ({imputed_count} values imputed)\n")
                        except Exception as e:
                            self._thread_safe_log(f"    ⚠️ Imputation failed for {polarity} classes: {e}\n")
                    else:
                        self._thread_safe_log(f"    • {polarity.capitalize()} classes: Pre-normalization imputation is disabled (Step 4a).\n")

                    # APPLY NORMALIZATION
                    is_feature_name = None
                    qc_sample_cols = None
                    if norm_method.lower() in ('is', 'internal_standard', 'istd'):
                        if polarity == 'positive' and hasattr(self, 'verified_pos_lipid_assignments'):
                            is_feature_name = self.verified_pos_lipid_assignments.get('Internal Standard')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_lipid_assignments'):
                            is_feature_name = self.verified_neg_lipid_assignments.get('Internal Standard')
                    elif norm_method.lower() in ('loess_qc', 'loess', 'qc_correction'):
                        if polarity == 'positive' and hasattr(self, 'verified_pos_lipid_assignments'):
                            qc_sample_cols = self.verified_pos_lipid_assignments.get('_qc_columns')
                        elif polarity == 'negative' and hasattr(self, 'verified_neg_lipid_assignments'):
                            qc_sample_cols = self.verified_neg_lipid_assignments.get('_qc_columns')

                    self._thread_safe_log(f"    • {polarity.capitalize()} classes: Normalizing with method '{norm_method}' ({len(sample_cols)} sample columns)...\n")
                    normalized_df = normalize_dataframe(
                        working_df,
                        sample_cols,
                        norm_method,
                        is_feature_name=is_feature_name,
                        qc_sample_cols=qc_sample_cols,
                    )
                    self._thread_safe_log(f"    ✅ {polarity.capitalize()} classes: Normalization complete\n")

                    rename_map = class_column_mappings.get(polarity, {})
                    if rename_map:
                        normalized_df = normalized_df.rename(columns=rename_map)

                    # Step 4b: run the same optional post-normalization pipeline used for ions.
                    cleaned_sample_cols = [rename_map.get(col, col) for col in sample_cols] if rename_map else list(sample_cols)
                    cleaned_sample_cols = [col for col in cleaned_sample_cols if col in normalized_df.columns]
                    if step7_post_processing_enabled and cleaned_sample_cols:
                        self._thread_safe_log(
                            f"    • {polarity.capitalize()} classes: Applying optional post-normalization processing (Step 4b) with {len(cleaned_sample_cols)} sample columns...\n"
                        )
                        # Keep alias map aligned for group assignment resolution inside optional processing.
                        self._current_sample_column_aliases = rename_map.copy() if rename_map else {}
                        normalized_df, class_pol_report = self._apply_optional_post_normalization_processing(
                            normalized_df.reset_index(drop=True),
                            mode,
                            verified_sample_cols=cleaned_sample_cols
                        )
                        imputation_summary = ((class_pol_report or {}).get('imputation') or {})
                        if imputation_summary:
                            self._thread_safe_log(
                                f"    ✅ {polarity.capitalize()} classes: Step 4b complete; imputation filled {imputation_summary.get('imputed_cells', 0)} cells.\n"
                            )
                        else:
                            self._thread_safe_log(f"    ✅ {polarity.capitalize()} classes: Step 4b complete.\n")
                    elif step7_post_processing_enabled and not cleaned_sample_cols:
                        self._thread_safe_log(
                            f"    ⚠️ {polarity.capitalize()} classes: Step 4b enabled but no valid sample columns found after renaming.\n"
                        )

                    return normalized_df

                self.normalized_positive_class_df = _normalize_class_frame(pos_class_df, 'positive')
                self.normalized_negative_class_df = _normalize_class_frame(neg_class_df, 'negative')

                # Step 6: Merge normalized class tables with special handling for n_lipids
                final_class_frames = []
                if self.normalized_positive_class_df is not None and not self.normalized_positive_class_df.empty:
                    final_class_frames.append(self.normalized_positive_class_df)
                if self.normalized_negative_class_df is not None and not self.normalized_negative_class_df.empty:
                    final_class_frames.append(self.normalized_negative_class_df)

                if final_class_frames:
                    combined_class = pd.concat(final_class_frames, ignore_index=True)

                    # Merge duplicate classes with special aggregation rules
                    if 'Lipid_Class' in combined_class.columns:
                        # Identify all column types
                        feature_cols = {'CalcMz', 'BaseRt', 'Rt', 'MZ', 'Mz', 'Max_Area', 'Area (Max.)'}
                        count_cols = {'n_lipids', 'n_lipid'}  # Count columns that should be summed
                        
                        # Build aggregation map
                        agg_map = {}
                        for col in combined_class.columns:
                            if col == 'Lipid_Class':
                                continue  # Skip the groupby column
                            elif col in count_cols:
                                # Sum count columns
                                agg_map[col] = 'sum'
                            elif col in feature_cols or not pd.api.types.is_numeric_dtype(combined_class[col]):
                                # Use first for feature/metadata columns
                                agg_map[col] = 'first'
                            elif pd.api.types.is_numeric_dtype(combined_class[col]):
                                # Use mean for actual sample columns
                                agg_map[col] = 'mean'
                        
                        before = len(combined_class)
                        combined_class = (
                            combined_class
                            .groupby('Lipid_Class', as_index=False, sort=False)
                            .agg(agg_map)
                        )
                        merged = before - len(combined_class)
                        if merged > 0:
                            self._thread_safe_log(
                                f"  • Class merge: Merged {merged} duplicates (same Lipid_Class from pos+neg)\n"
                                f"      - MEAN aggregation for sample columns\n"
                                f"      - SUM aggregation for n_lipids count column\n"
                                f"      - FIRST aggregation for feature columns (CalcMz, BaseRt, etc.)\n"
                            )

                    # Reorder columns to put Lipid_Class first
                    if 'Lipid_Class' in combined_class.columns:
                        cols = combined_class.columns.tolist()
                        cols.insert(0, cols.pop(cols.index('Lipid_Class')))
                        combined_class = combined_class[cols]

                    self.normalized_combined_class_df = combined_class.reset_index(drop=True)
                    class_generation_report['applied'] = True
                    class_generation_report['merged']['applied'] = True
                    class_generation_report['merged']['class_count'] = len(combined_class)
                    class_generation_report['source'] = 'deduplicated_ions_grouped_to_class_then_normalized'
                    
                    # Log summary
                    pos_count = len(self.normalized_positive_class_df) if self.normalized_positive_class_df is not None and not self.normalized_positive_class_df.empty else 0
                    neg_count = len(self.normalized_negative_class_df) if self.normalized_negative_class_df is not None and not self.normalized_negative_class_df.empty else 0
                    total_before_merge = pos_count + neg_count
                    
                    # Build workflow description based on actual imputation settings
                    workflow_steps = ['Ions', 'deduplicated', 'grouped']
                    if step7_pre_imputation_enabled:
                        workflow_steps.append('imputed (4a)')
                    workflow_steps.append('normalized')
                    if step7_post_imputation_enabled:
                        workflow_steps.append('imputed (4b)')
                    workflow_steps.append('merged (MEAN)')
                    workflow_str = ' → '.join(workflow_steps)
                    
                    self._thread_safe_log(
                        f"\n✅ Step 7 complete: Class workflow finished\n"
                        f"  → Positive polarity: {pos_count} normalized classes\n"
                        f"  → Negative polarity: {neg_count} normalized classes\n"
                        f"  → Combined before merge: {total_before_merge} class rows\n"
                        f"  → Final merged output: {len(combined_class)} unique lipid classes\n"
                        f"  → Workflow: {workflow_str}\n"
                    )
                else:
                    self._thread_safe_log("No normalized class data to merge.\n")
        
        # Store the class generation report
        self.class_generation_report = class_generation_report
        if not hasattr(self, 'optional_processing_report') or not isinstance(self.optional_processing_report, dict):
            self.optional_processing_report = {}
        self.optional_processing_report['class_generation'] = self.class_generation_report
        
        # Log final summary
        sheet_summary = []
        if self.normalized_positive_df is not None:
            sheet_summary.append(f"Positive: {len(self.normalized_positive_df)} rows")
        if self.normalized_negative_df is not None:
            sheet_summary.append(f"Negative: {len(self.normalized_negative_df)} rows")
        sheet_summary.append(f"Combined: {len(self.normalized_combined_df)} rows")
        if hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
            sheet_summary.append(f"Combined Class: {len(self.normalized_combined_class_df)} rows")
        
        self._thread_safe_progress_step(total_steps - 1, total_steps, "Finalizing normalized dataframes...")
        self.stats_log.insert(tk.END, f'✅ Normalized dataframes prepared:\n')
        for summary in sheet_summary:
            self.stats_log.insert(tk.END, f'  • {summary}\n')
        self.stats_log.insert(tk.END, '✅ Normalization & merge complete with Area (Max.) priority and cleaned column names!\n')
        self.stats_log.see(tk.END)

        raw_debug_path = self._export_pre_normalization_debug_workbook(raw_pre_normalization_exports)
        if raw_debug_path:
            self.last_raw_pre_normalization_export_path = raw_debug_path
            self._thread_safe_log(f'🧪 Pre-normalization raw workbook saved to: {raw_debug_path}\n')
            try:
                backup_note = os.path.join(os.path.dirname(raw_debug_path), '..', 'logs')
                self._thread_safe_log(f'🧪 Backup copy also written to logs folder when distinct from the selected output folder.\n')
            except Exception:
                pass
        
        # Hide progress bar
        self._thread_safe_progress_step(total_steps, total_steps, "Normalization complete")
        self.hide_stats_progress()

        optional_msg = ""
        if bool(getattr(self, 'optional_processing_applied', False)):
            optional_msg = "• Optional post-normalization steps: applied\n"
        
        # Show completion popup and populate sample assignments
        messagebox.showinfo("Normalization Complete", 
                          f"Data normalization completed successfully!\n\n"
                          f"• Processed data: {len(self.normalized_combined_df)} features\n"
                          f"⚠️ IMPORTANT NEXT STEPS:\n"
                          f"Click 'Run Statistics' or 'Run Stat with Covariate Adjustment' to:\n"
                          f"'Review group assignments' & 'perform statistical tests'")
                         
        
        # Populate sample assignments with normalized data
        self.populate_sample_assignments_from_normalized_data()
    
    def _run_preprocessed_pipeline(self):
        """
        CUSTOM/PREPROCESSED: Enhanced pipeline for preprocessed/combined data.
        Skips Pos/Neg splitting and merging, but applies filtering, normalization, and normality testing.
        """
        from main_script.metabolite_statistics_analysis import normalize_dataframe, apply_min_group_size_filter
        
        total_steps = 8
        self._thread_safe_progress_step(1, total_steps, "Loading preprocessed data...")
        
        # Get preprocessed data from memory
        combined_df = self.memory_store.get('preprocessed_combined_df')
        feature_cols = self.memory_store.get('preprocessed_feature_cols', [])
        sample_cols = self.memory_store.get('preprocessed_sample_cols', [])
        mode = self.memory_store.get('preprocessed_mode', 'custom')  # Default to 'custom'
        
        # Determine if this is custom mode
        is_custom = self.memory_store.get('is_preprocessed_custom', False) or mode == 'custom'
        mode_label = "CUSTOM" if is_custom else "PREPROCESSED"
        
        if combined_df is None or not sample_cols:
            error_msg = 'Preprocessed data not properly loaded. Please re-import.'
            self._thread_safe_log(f"❌ ERROR: {error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('Invalid Preprocessed Data', error_msg))
            return
        
        self._thread_safe_log(f"📊 {mode_label} data: {len(combined_df)} rows, {len(sample_cols)} samples\n")
        self._thread_safe_log(f"📊 Feature columns: {', '.join(feature_cols[:5])}{' ...' if len(feature_cols) > 5 else ''}\n")
        
        # Store original count for comparison
        original_count = len(combined_df)
        
        # Get the user-verified Feature ID column from stored assignments
        self._thread_safe_progress_step(2, total_steps, "Preparing feature columns...")
        verified_assignments = self.memory_store.get('preprocessed_verified_assignments', {})
        
        # For custom mode, always use 'Feature ID'; for backward compatibility also check 'LipidID'
        main_feature_col = verified_assignments.get('Feature ID')
        if not main_feature_col and mode == 'lipid':
            main_feature_col = verified_assignments.get('LipidID')
        
        # If no verified assignment, try to use the first feature column from import
        if not main_feature_col and feature_cols and len(feature_cols) > 0:
            main_feature_col = feature_cols[0]
            self._thread_safe_log(f"⚠️ Using first feature column as ID: '{main_feature_col}'\n")
            # Store it for later use
            verified_assignments['Feature ID'] = main_feature_col
            self.memory_store['preprocessed_verified_assignments'] = verified_assignments
        
        if not main_feature_col:
            # If still no feature ID, log error and abort
            error_msg = 'Feature ID column not found. Please verify columns in Step 2.'
            self._thread_safe_log(f"❌ ERROR: {error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('Missing Column Assignment', error_msg))
            return
        
        self._thread_safe_log(f"✅ Using Feature ID column: '{main_feature_col}'\n")
        
        # ========== PRE-NORMALIZATION FILTERING ==========
        filter_timing = self.filter_timing_var.get() if hasattr(self, 'filter_timing_var') else 'before'
        
        if self._is_imputation_enabled() and hasattr(self, 'sample_group_vars') and self.sample_group_vars:
            self._thread_safe_progress_step(3, total_steps, "Applying imputation pre-filter before normalization...")

            # Parse group assignments
            temp_group_map = {}
            full_assignments = self._parse_group_assignments()
            for col in sample_cols:
                if col in full_assignments:
                    temp_group_map[col] = full_assignments[col]

            if temp_group_map and len(temp_group_map) == len(sample_cols):
                imp_min_pct = self._get_imputation_min_group_percent()
                imp_scope = self._get_imputation_prefilter_scope()
                scope_text = 'at least one group' if imp_scope == 'per_group' else 'all groups'
                self._thread_safe_log(
                    f"\n🔍 Applying imputation pre-filter BEFORE normalization...\n"
                    f"   Rule: keep rows with at least {imp_min_pct:.1f}% valid values in {scope_text}\n"
                )

                filtered_df, pre_imp_report = self._apply_imputation_prefilter(
                    combined_df,
                    sample_cols,
                    temp_group_map,
                    imp_min_pct,
                    scope=imp_scope,
                )

                self._thread_safe_log(
                    f"   ✅ Imputation pre-filter (pre-normalization): removed {pre_imp_report.get('removed', 0)} rows, "
                    f"kept {pre_imp_report.get('kept', 0)}.\n\n"
                )
                combined_df = filtered_df
            else:
                self._thread_safe_log(
                    "⚠️  Imputation pre-filter before normalization skipped: Not all samples assigned to groups\n"
                )
        elif filter_timing == 'before' and hasattr(self, 'sample_group_vars') and self.sample_group_vars:
            self._thread_safe_progress_step(3, total_steps, "Applying pre-normalization filter...")
            
            # Parse group assignments
            temp_group_map = {}
            full_assignments = self._parse_group_assignments()
            for col in sample_cols:
                if col in full_assignments:
                    temp_group_map[col] = full_assignments[col]
            
            if temp_group_map and len(temp_group_map) == len(sample_cols):
                # All samples have group assignments - apply filtering
                min_type = self.min_samples_type_var.get() if hasattr(self, 'min_samples_type_var') else 'absolute'
                
                try:
                    min_count_str = self.min_samples_per_group_var.get() if hasattr(self, 'min_samples_per_group_var') else '2'
                    min_count = int(min_count_str) if min_count_str else 2
                except (ValueError, Exception):
                    min_count = 2
                
                try:
                    min_percent_str = self.min_samples_percent_var.get() if hasattr(self, 'min_samples_percent_var') else '50.0'
                    min_percent = float(min_percent_str) if min_percent_str else 50.0
                except (ValueError, Exception):
                    min_percent = 50.0
                
                self._thread_safe_log(f'\n🔍 Applying pre-normalization filter...\n')
                
                from collections import Counter
                group_counts = Counter(temp_group_map.values())
                self._thread_safe_log(f'   👥 Groups: {", ".join([f"{g}={c}" for g, c in sorted(group_counts.items())])}\n')
                
                if min_type == 'percentage':
                    self._thread_safe_log(f'   Threshold: {min_percent}% of group size\n')
                else:
                    self._thread_safe_log(f'   Threshold: {min_count} samples per group\n')
                
                filtered_df, group_stats = apply_min_group_size_filter(
                    combined_df, sample_cols, temp_group_map, 
                    min_count, min_type, min_percent
                )
                
                removed = len(combined_df) - len(filtered_df)
                self._thread_safe_log(f'   ✅ Filtered: {removed} rows removed, {len(filtered_df)} rows retained\n\n')
                combined_df = filtered_df
            else:
                self._thread_safe_log(f"⚠️  Skipping pre-normalization filter: Not all samples assigned to groups\n")
        else:
            self._thread_safe_log(f"⚠️  Pre-normalization filtering: {filter_timing} (skipped in custom mode if 'after')\n")
        
        # Get normalization method
        norm_method = self.stat_norm_method.get()
        
        # Apply normalization
        self._thread_safe_progress_step(4, total_steps, f"Normalizing data ({norm_method})...")
        try:
            pre_norm_values = combined_df[sample_cols].apply(pd.to_numeric, errors='coerce')
            normalized_df = normalize_dataframe(combined_df, sample_cols, norm_method)
            post_norm_values = normalized_df[sample_cols].apply(pd.to_numeric, errors='coerce')

            # Diagnostic audit: quantify how much normalization changed the matrix.
            diff = (post_norm_values - pre_norm_values).to_numpy(dtype=float)
            finite_mask = np.isfinite(diff)
            if finite_mask.any():
                abs_diff = np.abs(diff[finite_mask])
                changed_cells = int(np.sum(abs_diff > 1e-12))
                total_cells = int(abs_diff.size)
                pct_changed = (changed_cells / total_cells * 100.0) if total_cells else 0.0
                max_abs_delta = float(abs_diff.max()) if abs_diff.size else 0.0
                mean_abs_delta = float(abs_diff.mean()) if abs_diff.size else 0.0
                self._thread_safe_log(
                    f"   Δ cells changed: {changed_cells}/{total_cells} ({pct_changed:.2f}%), "
                    f"max |Δ|={max_abs_delta:.6g}, mean |Δ|={mean_abs_delta:.6g}\n"
                )
            else:
                self._thread_safe_log("   Δ cells changed: unable to compute (all deltas non-finite).\n")

            self._thread_safe_log(f"✅ Applied {norm_method} normalization to {len(sample_cols)} sample columns.\n")
        except Exception as e:
            error_msg = f'Normalization failed: {e}'
            self._thread_safe_log(f"❌ ERROR: {error_msg}\n")
            self.root.after(0, lambda: messagebox.showerror('Normalization Error', error_msg))
            return

        # Step 4b optional processing must run before normality testing.
        self.raw_normalized_combined_df = normalized_df.reset_index(drop=True)
        processed_combined_df, processing_report = self._apply_optional_post_normalization_processing(
            self.raw_normalized_combined_df,
            self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite',
            verified_sample_cols=sample_cols
        )
        self.optional_processing_report = processing_report
        self.optional_processing_applied = bool(processing_report.get('applied', False))
        normalized_df = processed_combined_df.reset_index(drop=True)
        sample_cols = [c for c in sample_cols if c in normalized_df.columns]
        
        # ========== NORMALITY TESTING ==========
        self._thread_safe_progress_step(5, total_steps, "Testing normality...")
        if norm_method != 'none':
            try:
                from scipy.stats import shapiro
                import numpy as _np
                import pandas as _pd
                
                self._thread_safe_log(f"\n📊 Running normality tests (Shapiro-Wilk)...\n")
                
                records = []
                metabolite_id_col = main_feature_col if main_feature_col else 'metabolite_id'
                
                normal_metabolites = []
                not_normal_metabolites = []
                
                total_metabolites = len(normalized_df)
                callback_freq = max(1, total_metabolites // 100)
                
                for seq_idx, (idx, row) in enumerate(normalized_df.iterrows(), start=1):
                    if seq_idx % callback_freq == 0 or seq_idx == 1 or seq_idx == total_metabolites:
                        self._thread_safe_log(f'  Testing: {seq_idx}/{total_metabolites}\n')
                    
                    metabolite_id = row[metabolite_id_col] if metabolite_id_col in normalized_df.columns else f'Feature_{idx}'
                    
                    # Get values across all samples
                    vals_list = []
                    for col in sample_cols:
                        if col in normalized_df.columns:
                            val = _pd.to_numeric(row[col], errors='coerce')
                            if not _pd.isna(val) and val != 0:
                                vals_list.append(val)
                    
                    n = len(vals_list)
                    
                    if n < 3:
                        records.append({
                            'Metabolite': metabolite_id,
                            'N': n,
                            'Shapiro_p': _np.nan,
                            'Is_Normal': 'TooFew'
                        })
                        continue
                    
                    # Shapiro-Wilk test
                    vals_array = _np.array(vals_list)
                    shapiro_sample = vals_array if n <= 5000 else _np.random.choice(vals_array, 5000, replace=False)
                    try:
                        sh_stat, sh_p = shapiro(shapiro_sample)
                    except Exception:
                        sh_p = _np.nan
                    
                    is_norm = 'Yes' if (not _np.isnan(sh_p) and sh_p >= 0.05) else 'No' if not _np.isnan(sh_p) else 'Unknown'
                    
                    records.append({
                        'Metabolite': metabolite_id,
                        'N': n,
                        'Shapiro_p': sh_p,
                        'Is_Normal': is_norm
                    })
                    
                    if is_norm == 'Yes':
                        normal_metabolites.append((idx, metabolite_id, vals_array))
                    elif is_norm == 'No':
                        not_normal_metabolites.append((idx, metabolite_id, vals_array))
                
                normality_results = _pd.DataFrame(records) if records else _pd.DataFrame(columns=['Metabolite','N','Shapiro_p','Is_Normal'])
                
                if not hasattr(self, 'normality_test_results'):
                    self.normality_test_results = {}
                self.normality_test_results['combined'] = normality_results
                
                # Store selected metabolites for QQ plots
                if not hasattr(self, 'normality_test_targets'):
                    self.normality_test_targets = {}
                
                selected_plots = []
                if len(normal_metabolites) > 0:
                    import random
                    selected_normal = random.sample(normal_metabolites, min(2, len(normal_metabolites)))
                    selected_plots.extend([('NORMAL', m_id, vals) for _, m_id, vals in selected_normal])
                if len(not_normal_metabolites) > 0:
                    import random
                    selected_not_normal = random.sample(not_normal_metabolites, min(2, len(not_normal_metabolites)))
                    selected_plots.extend([('NOT_NORMAL', m_id, vals) for _, m_id, vals in selected_not_normal])
                
                self.normality_test_targets['combined'] = {
                    'plot_data': selected_plots,
                    'normal_count': len(normal_metabolites),
                    'not_normal_count': len(not_normal_metabolites)
                }
                
                # Log summary
                if not normality_results.empty:
                    normal_count = len(normality_results[normality_results['Is_Normal'] == 'Yes'])
                    not_normal_count = len(normality_results[normality_results['Is_Normal'] == 'No'])
                    too_few_count = len(normality_results[normality_results['Is_Normal'] == 'TooFew'])
                    total = len(normality_results)
                    pct_norm = (normal_count / total * 100.0) if total else 0.0
                    pct_not = (not_normal_count / total * 100.0) if total else 0.0
                    self._thread_safe_log(f'\n📊 Normality Test Summary:\n')
                    self._thread_safe_log(f'  ✓ Normal: {normal_count}/{total} ({pct_norm:.1f}%)\n')
                    self._thread_safe_log(f'  ✗ Not Normal: {not_normal_count}/{total} ({pct_not:.1f}%)\n')
                    if too_few_count > 0:
                        self._thread_safe_log(f'  ⚠ Too Few Samples: {too_few_count}/{total}\n')
                    self._thread_safe_log(f'  📈 Selected {len(selected_plots)} metabolites for QQ plots\n\n')
                    
            except ImportError:
                self._thread_safe_log(f"⚠️  scipy not available, skipping normality tests\n")
            except Exception as e:
                self._thread_safe_log(f"⚠️  Normality testing error: {e}\n")
        else:
            self._thread_safe_log(f"⚠️  Normality testing skipped (normalization method: none)\n")
        
        # Store normalized results (no Pos/Neg split for preprocessed data)
        self._thread_safe_progress_step(6, total_steps, "Storing normalized data...")
        self.normalized_combined_df = normalized_df.reset_index(drop=True)
        self.normalized_positive_df = None  # No polarity split in custom mode
        self.normalized_negative_df = None
        
        self._thread_safe_log(f"✅ Stored combined normalized data: {len(self.normalized_combined_df)} rows\n")
        self._thread_safe_log(f"⚠️  Note: No Positive/Negative split (custom/preprocessed mode)\n")
        
        # Populate sample assignments
        self._thread_safe_progress_step(7, total_steps, "Populating sample assignments...")
        self.root.after(0, lambda: self.populate_sample_assignments_from_normalized_data())
        
        # Completion
        self._thread_safe_progress_step(8, total_steps, "Custom mode normalization complete")
        
        filtered_count = original_count - len(self.normalized_combined_df)
        
        self._thread_safe_log(f"\n✅ {mode_label} MODE COMPLETE!\n")
        self._thread_safe_log(f"• Original features: {original_count}\n")
        if filtered_count > 0:
            self._thread_safe_log(f"• Filtered out: {filtered_count} ({filtered_count/original_count*100:.1f}%)\n")
        self._thread_safe_log(f"• Normalized features: {len(self.normalized_combined_df)}\n")
        self._thread_safe_log(f"• Method: {norm_method}\n")
        self._thread_safe_log(f"• Ready for statistical tests\n")
        self.stats_log.see(tk.END)
        
        self.hide_stats_progress()
        
        # Show completion message
        filter_msg = f"• Filtered: {filtered_count} rows removed\n" if filtered_count > 0 else ""
        self.root.after(0, lambda: messagebox.showinfo(
            "Custom Mode Normalization Complete",
            f"Preprocessed data processed successfully!\n\n"
            f"• Original data: {original_count} features\n"
            f"{filter_msg}"
            f"• Final data: {len(self.normalized_combined_df)} features\n"
            f"• Applied normalization: {norm_method}\n"
            f"• Normality tests: Completed\n"
            f"• Mode: Custom (preprocessed - no Pos/Neg split)\n\n"
            f"Ready for statistical tests!"
        ))

    def _ensure_groups_ready(self, after_config_callback=None, action_label='analysis'):
        """Ensure group assignments are confirmed before downstream actions."""
        if not hasattr(self, 'sample_group_vars') or not self.sample_group_vars:
            messagebox.showerror('Groups Not Configured', 'No sample columns detected. Complete Steps 2 and 3 before running analyses.')
            return False

        assignments = [gv.get().strip() for gv in self.sample_group_vars.values()]
        configured_once = getattr(self, '_groups_configured_once', False)

        if configured_once and all(assignments):
            return True

        try:
            self._log_stats(f'Opening Configure Groups before {action_label}...')
        except Exception:
            pass

        self.auto_assign_groups(on_done=after_config_callback)
        return False

    # --- Threaded Statistical Tests with determinate progress ---
    def run_statistical_tests(self):
        """Run statistical tests (threaded) with a determinate loading bar."""
        if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
            messagebox.showerror('No Data', 'No normalized data available. Run "Normalization & Merge" first.')
            return

        def _proceed():
            self._run_statistical_tests_core()

        if not self._ensure_groups_ready(after_config_callback=_proceed, action_label='statistical tests'):
            return

        self._run_statistical_tests_core()

    def _run_statistical_tests_core(self):
        """Run statistical tests after group configuration is confirmed."""
        if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
            messagebox.showerror('No Data', 'No normalized data available. Run "Normalization & Merge" first.')
            return

        # Clear old statistical results to prevent visualization tab from using stale data
        self._clear_statistical_results()
        
        # Sync latest user edits to group labels so defaults (Control/Disease) are not used erroneously
        try:
            if hasattr(self, 'group_id_vars'):
                for gid, var in self.group_id_vars.items():
                    val = var.get().strip()
                    if val:
                        self.group_definitions[gid] = val
        except Exception:
            pass
        # Disable button
        try:
            for widget in self.frame.winfo_children():
                if isinstance(widget, tk.Frame):
                    for child in widget.winfo_children():
                        if isinstance(child, tk.Button) and 'Statistical Tests' in child.cget('text'):
                            child.config(state='disabled')
        except Exception:
            pass
        # Reset progress context
        self._stats_total_steps = 5
        self._stats_current_step = 0
        self.show_stats_progress("Preparing statistical analysis...")
        threading.Thread(target=self._run_statistical_tests_threaded, daemon=True).start()

    def _run_two_way_anova_direct(self):
        """Execute two-way ANOVA using automated factor assignments in background thread."""
        threading.Thread(target=self._run_two_way_anova_threaded, daemon=True).start()

    def _run_two_way_anova_threaded(self):
        """Execute two-way ANOVA in background thread (prevents UI freezing)."""
        try:
            if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
                self._thread_safe_log('❌ No normalized data available.\n')
                self._thread_safe_log('🛑 Debug: normalized_combined_df is None, run normalization before Two-Way ANOVA.\n')
                return
            group_map = self._parse_group_assignments()
            if not group_map:
                self._thread_safe_log('❌ No group assignments found.\n')
                self._thread_safe_log('🛑 Debug: _parse_group_assignments returned empty; ensure samples are assigned to groups.\n')
                return
            sample_cols = list(group_map.keys())
            from main_script.metabolite_statistics_analysis import perform_two_way_anova, perform_two_way_anova_posthoc
            factor_a_map = {c: v.get() for c, v in getattr(self, 'sample_factorA_vars', {}).items() if v.get().strip()}
            factor_b_map = {c: v.get() for c, v in getattr(self, 'sample_factorB_vars', {}).items() if v.get().strip()}
            if not factor_a_map or not factor_b_map:
                error_msg = '❌ FATAL ERROR: Two-Way ANOVA factor mappings are missing or incomplete.\n\n'
                error_msg += 'Two-Way ANOVA requires ALL samples to be assigned to BOTH Factor A and Factor B.\n\n'
                error_msg += f'Factor A assignments: {len(factor_a_map)} samples\n'
                error_msg += f'Factor B assignments: {len(factor_b_map)} samples\n\n'
                error_msg += 'Please use "Configure Factor Assignments" to assign all samples to both factors.'
                self._thread_safe_log(error_msg + '\n')
                messagebox.showerror('Factor Assignments Missing', error_msg)
                return
            included = [c for c in sample_cols if c in factor_a_map and c in factor_b_map]
            if len(included) < 3:
                missing = [c for c in sample_cols if c not in factor_a_map or c not in factor_b_map]
                error_msg = f'❌ FATAL ERROR: Too few samples with both factors assigned ({len(included)}/3 minimum).\n\n'
                error_msg += f'Missing factor assignments for {len(missing)} samples:\n'
                error_msg += ', '.join(missing[:10]) + ('...' if len(missing) > 10 else '') + '\n\n'
                error_msg += 'All samples must be assigned to BOTH Factor A and Factor B for Two-Way ANOVA.'
                self._thread_safe_log(error_msg + '\n')
                messagebox.showerror('Insufficient Factor Assignments', error_msg)
                return
            self._thread_safe_log(f'▶ Running Two-Way ANOVA on {len(included)} samples...\n')
            # Determine worker count
            try:
                workers = self._get_workers_count(getattr(self, 'stats_workers', None), default=3)
            except Exception:
                workers = 3
            
            # Define progress callback for Two-Way ANOVA
            def two_way_progress(current, total, metabolite_name):
                percent = int(100 * current / total) if total > 0 else 0
                # Update every 1% or every item for small datasets
                if current % max(1, total // 100) == 0 or current == total:
                    self._thread_safe_log(f'  [{current}/{total} - {percent}%] {metabolite_name}\n')
            
            # Use NEW format function that returns complete result + individual sheets
            from main_script.two_way_anova_new_format import perform_two_way_anova_new
            
            # Get user-defined factor names
            fa_name = getattr(self, 'factorA_name', 'Factor A')
            fb_name = getattr(self, 'factorB_name', 'Factor B')
            
            # Get group order from GUI settings (respects Statistics Group Order field)
            group_order = None
            if hasattr(self, 'statistics_group_order_var'):
                order_str = self.statistics_group_order_var.get().strip()
                if order_str:
                    group_order = [g.strip() for g in order_str.split(',') if g.strip()]
            # Fallback to group_definitions order if no explicit order
            if not group_order and hasattr(self, 'group_definitions'):
                group_order = [self.group_definitions[gid] for gid in self.group_definitions.keys()]
            
            # Log group order being used
            if group_order:
                self._thread_safe_log(f'📋 Group comparison order: {" < ".join(group_order)}\n')
            
            complete_result, individual_sheets = perform_two_way_anova_new(
                self.normalized_combined_df,
                included,
                factor_a_map,
                factor_b_map,
                group_map=group_map,
                group_order=group_order,
                factor_a_name=fa_name,
                factor_b_name=fb_name,
                fdr=True,
                n_jobs=workers,
                progress_callback=two_way_progress
            )
            self._thread_safe_log(f'✅ Two-Way ANOVA computation complete (including all pairwise Tukey HSD tests)\n')
            
            # Store complete result for summary stats
            twa_df = complete_result
            self.two_way_anova_results = twa_df
            
            # NEW FORMAT: Results already in correct format from perform_two_way_anova_new()
            self._thread_safe_log(f'📊 Storing results in new format...\n')
            # Provide both keys for compatibility with Visualization tab
            combined = {
                'complete_result': complete_result,
                'enhanced_metabolites': complete_result  # visualization expects this key
            }
            combined.update(individual_sheets)  # Add all individual pairwise sheets
            self.statistical_test_results = combined
            
            # Ensure wide-format pairwise columns are in enhanced_metabolites for volcano plot compatibility
            if 'enhanced_metabolites' in self.statistical_test_results:
                enhanced_df = self.statistical_test_results['enhanced_metabolites']
                # Find id column
                id_col = None
                for candidate in ['Name', 'Protein', 'Metabolite', 'Compound', 'LipidID', 'Lipid_ID', 'metabolite']:
                    if candidate in enhanced_df.columns:
                        id_col = candidate
                        break
                if id_col is None and len(enhanced_df.columns) > 0:
                    id_col = enhanced_df.columns[0]
                
                # Check if pairwise columns are already present
                has_pairwise = any('_vs_' in col and '_log2FC' in col for col in enhanced_df.columns)
                if not has_pairwise and individual_sheets:
                    # Merge pairwise results into enhanced_df
                    enhanced_df = enhanced_df.set_index(id_col)
                    for comp_name, pairwise_df in individual_sheets.items():
                        if pairwise_df.empty:
                            continue
                        pairwise_df = pairwise_df.set_index(id_col)
                        # Add pairwise columns
                        for col in pairwise_df.columns:
                            if col != id_col and col not in enhanced_df.columns:
                                enhanced_df = enhanced_df.join(pairwise_df[[col]], how='left')
                    enhanced_df = enhanced_df.reset_index()
                    self.statistical_test_results['enhanced_metabolites'] = enhanced_df
                    self._thread_safe_log(f'✅ Merged pairwise columns into enhanced_metabolites for volcano compatibility\n')
            
            self._thread_safe_log(f'✅ Results formatted: 1 Complete Result sheet + {len(individual_sheets)} pairwise sheets\n')
            if not twa_df.empty:
                try:
                    sig_A = (twa_df.get('p_A_adj', twa_df.get('p_A')).le(0.05)).sum()
                    sig_B = (twa_df.get('p_B_adj', twa_df.get('p_B')).le(0.05)).sum()
                    sig_I = (twa_df.get('p_interaction_adj', twa_df.get('p_interaction')).le(0.05)).sum()
                    self._thread_safe_log(f'✅ Two-way ANOVA complete. Significant A:{sig_A} B:{sig_B} Interaction:{sig_I} (α=0.05)\n')
                except Exception:
                    pass
            else:
                self._thread_safe_log('⚠️ Two-way ANOVA produced no rows.\n')
            self._thread_safe_log("➡️ Use 'Export Statistical Results' to save TwoWayANOVA and Posthoc sheets.\n")
            
            # Store factor names for final completion message
            fa_name = getattr(self, 'factorA_name', 'Factor A')
            fb_name = getattr(self, 'factorB_name', 'Factor B')
            regular_pairwise_count = len(individual_sheets)
            
            try:
                # Store in shared memory for Visualization tab
                self.memory_store['statistical_test_results'] = self.statistical_test_results
                # Also push a light nudge to Visualization to refresh its data source
                viz_tab = self.get_tab_by_name("📊 Visualization")
                if viz_tab and hasattr(viz_tab, 'update_viz_data_status'):
                    try:
                        viz_tab.update_viz_data_status()
                        self._thread_safe_log("✅ Visualization updated with new statistical results\n")
                    except Exception:
                        pass
            except Exception:
                pass
            
            # Run Two-Way ANOVA on lipid class data if in lipid mode
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            if mode == 'lipid' and hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                self._thread_safe_log("\n🔬 Running Two-Way ANOVA on lipid class data...\n")
                
                # Get threshold parameters for posthoc filtering
                min_group_size_param = self._get_min_group_size()
                min_group_size_percent_param = self._get_min_group_size_percent()
                
                # Verify sample columns exist in class data
                class_missing_cols = [c for c in included if c not in self.normalized_combined_class_df.columns]
                if class_missing_cols:
                    self._thread_safe_log(f"⚠️ Warning: Some sample columns missing in class data: {class_missing_cols[:5]}...\n")
                    class_included = [c for c in included if c in self.normalized_combined_class_df.columns]
                else:
                    class_included = included
                
                if len(class_included) < 4:
                    self._thread_safe_log(f"⚠️ Insufficient samples in class data for Two-Way ANOVA (need ≥4, have {len(class_included)}). Skipping.\n")
                    self._thread_safe_log(f"🛑 Debug: class_included={class_included}\n")
                else:
                    try:
                        # Progress callback for class Two-Way ANOVA
                        def class_two_way_progress(current, total, class_name):
                            percent = int(100 * current / total) if total > 0 else 0
                            # Update every 1% or every item for small datasets
                            if current % max(1, total // 100) == 0 or current == total:
                                self._thread_safe_log(f'  [Class {current}/{total} - {percent}%] {class_name}\n')
                        
                        self._thread_safe_log(f'▶ Running Two-Way ANOVA on {len(class_included)} samples for {len(self.normalized_combined_class_df)} lipid classes...\n')
                        
                        # Use the NEW format function with group_map (same as regular lipid analysis)
                        from main_script.two_way_anova_new_format import perform_two_way_anova_new
                        
                        class_complete_result, class_individual_sheets = perform_two_way_anova_new(
                            self.normalized_combined_class_df,
                            class_included,
                            factor_a_map,
                            factor_b_map,
                            group_map=group_map,
                            group_order=group_order,
                            factor_a_name=fa_name,
                            factor_b_name=fb_name,
                            fdr=True,
                            n_jobs=workers,
                            progress_callback=class_two_way_progress
                        )
                        
                        # Extract the Two-Way ANOVA results from complete result
                        class_twa_df = class_complete_result
                        
                        self._thread_safe_log(f'✅ Lipid class Two-Way ANOVA computation complete\n')
                        
                        # Posthoc results already included in individual sheets from perform_two_way_anova_new
                        self._thread_safe_log(f'✅ Lipid class posthoc tests complete (included in pairwise sheets)\n')
                        
                        # Store results in the same format as regular analysis
                        class_posthoc = {
                            'pairwise': class_complete_result  # The complete result contains all pairwise data
                        }
                        
                        # Store class results in new format (complete result + individual pairwise sheets)
                        class_combined = {
                            'complete_result': class_complete_result,
                            'enhanced_metabolites': class_complete_result
                        }
                        class_combined.update(class_individual_sheets)  # Add all individual pairwise sheets
                        
                        self.statistical_test_results_class = class_combined
                        
                        # Rename columns for class results
                        try:
                            fa_name = getattr(self, 'factorA_name', 'Factor A')
                            fb_name = getattr(self, 'factorB_name', 'Factor B')
                            rename_map = {}
                            if 'p_A' in class_twa_df.columns:
                                rename_map['p_A'] = f'p_{fa_name}'
                            if 'p_B' in class_twa_df.columns:
                                rename_map['p_B'] = f'p_{fb_name}'
                            if 'p_A_adj' in class_twa_df.columns:
                                rename_map['p_A_adj'] = f'p_{fa_name}_adj'
                            if 'p_B_adj' in class_twa_df.columns:
                                rename_map['p_B_adj'] = f'p_{fb_name}_adj'
                            if 'F_A' in class_twa_df.columns:
                                rename_map['F_A'] = f'F_{fa_name}'
                            if 'F_B' in class_twa_df.columns:
                                rename_map['F_B'] = f'F_{fb_name}'
                            if rename_map:
                                class_twa_df = class_twa_df.rename(columns=rename_map)
                        except Exception:
                            pass
                        
                        # Store simplified class results reference
                        self.two_way_anova_results_class = class_twa_df
                        
                        if not class_twa_df.empty:
                            try:
                                sig_A = (class_twa_df.get('p_A_adj', class_twa_df.get('p_A')).le(0.05)).sum()
                                sig_B = (class_twa_df.get('p_B_adj', class_twa_df.get('p_B')).le(0.05)).sum()
                                sig_I = (class_twa_df.get('p_interaction_adj', class_twa_df.get('p_interaction')).le(0.05)).sum()
                                self._thread_safe_log(f'✅ Lipid class Two-way ANOVA complete. Significant A:{sig_A} B:{sig_B} Interaction:{sig_I} (α=0.05)\n')
                            except Exception:
                                pass
                        
                        self._thread_safe_log("✅ Lipid class statistical analysis complete!\n")
                    
                    except Exception as e:
                        self._thread_safe_log(f'⚠️ Lipid class Two-Way ANOVA error: {e}\n')
                        import traceback
                        self._thread_safe_log(f'{traceback.format_exc()}\n')
            
            # Show single completion message AFTER both regular and class analyses
            fa_name = getattr(self, 'factorA_name', 'Factor A')
            fb_name = getattr(self, 'factorB_name', 'Factor B')
            regular_pairwise_count = locals().get('regular_pairwise_count', len(individual_sheets) if 'individual_sheets' in locals() else 0)
            
            completion_msg = f"Two-Way ANOVA completed successfully!\n\nFactors: {fa_name} × {fb_name}\n"
            if mode == 'lipid' and hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
                completion_msg += f"\n✅ Regular lipid analysis: {regular_pairwise_count} pairwise comparisons\n"
                completion_msg += f"✅ Lipid class analysis: Complete\n"
            else:
                completion_msg += f"{regular_pairwise_count} pairwise comparisons\n"
            completion_msg += "\nUse 'Export Statistical Results' to save results."
            
            messagebox.showinfo("Analysis Complete", completion_msg)
        
        except Exception as e:
            self._thread_safe_log(f'❌ Two-way ANOVA error: {e}\n')
            import traceback
            self._thread_safe_log(f'{traceback.format_exc()}\n')
        finally:
            # Re-enable button and hide progress
            self._re_enable_stat_tests_button()
            self.root.after(0, self.hide_stats_progress)

    def _re_enable_stat_tests_button(self):
        def enable():
            try:
                for widget in self.frame.winfo_children():
                    if isinstance(widget, tk.Frame):
                        for child in widget.winfo_children():
                            if isinstance(child, tk.Button) and 'Statistical Tests' in child.cget('text'):
                                child.config(state='normal')
            except Exception:
                pass
        self.root.after(0, enable)

    def _run_statistical_tests_threaded(self):
        try:
            self._thread_safe_progress_step(1, 5, "Collecting group assignments...")
            group_map = self._parse_group_assignments()
            if not group_map:
                self._thread_safe_progress_step(0, 5, "No group assignments found")
                self._thread_safe_log("❌ No group assignments. Abort.\n")
                messagebox.showerror('No Groups', 'No group assignments found. Use group assignment first.')
                return
            
            # Get threshold configuration
            min_required = self._get_min_group_size()
            min_type = self._get_min_group_size_type()
            min_percent = self._get_min_group_size_percent()
            filter_timing = self.filter_timing_var.get() if hasattr(self, 'filter_timing_var') else 'before'
            imputation_enabled = self._is_imputation_enabled()
            # Current statistics data mode (metabolite or lipid)
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'

            # IMPORTANT WORKFLOW RULE:
            # If Step 4b imputation is enabled, skip original min-samples-per-group filtering.
            if imputation_enabled:
                min_required = 1
                min_type = 'absolute'
                min_percent = 0.0
            
            self._thread_safe_log(f"\n📋 Sample Filtering Configuration:\n")
            self._thread_safe_log(f"  Filter timing: {filter_timing.upper()} normalization\n")
            if imputation_enabled:
                self._thread_safe_log("  Imputation workflow: ENABLED\n")
                self._thread_safe_log("  Original min-samples-per-group filtering: SKIPPED\n")
                imp_scope = self._get_imputation_prefilter_scope()
                scope_text = 'at least one group' if imp_scope == 'per_group' else 'all groups'
                self._thread_safe_log(
                    f"  Imputation pre-filter threshold: {self._get_imputation_min_group_percent():.1f}% valid in {scope_text}\n"
                )
            elif min_type == 'percentage':
                self._thread_safe_log(f"  Minimum samples per group: {min_percent}% (percentage)\n")
            else:
                reasons = []
                if mode != 'lipid':
                    reasons.append(f"mode={mode}")
                if not hasattr(self, 'normalized_combined_class_df'):
                    reasons.append('normalized_combined_class_df attr missing')
                elif self.normalized_combined_class_df is None:
                    reasons.append('normalized_combined_class_df is None')
                msg = '; '.join(reasons) if reasons else 'unknown reason'
                self._thread_safe_log(f"\nℹ️ Skipping lipid class Two-Way ANOVA: {msg}.\n")
                self._thread_safe_log(f"  Minimum samples per group: {min_required} (absolute count)\n")
            
            if imputation_enabled:
                self._thread_safe_log("  → Statistical-test min-group filtering bypassed (handled by imputation pre-filter).\n")
            elif filter_timing == 'after':
                self._thread_safe_log(f"  → Filtering will be applied during statistical tests (per-comparison)\n")
                if min_type == 'percentage':
                    self._thread_safe_log(f"  → Metabolites will be skipped if a group has <{min_percent}% non-zero values for that comparison\n")
                    self._thread_safe_log(f"  → Example: Group with 10 samples requires ≥{int(np.ceil(10 * min_percent / 100))} non-zero values\n")
                else:
                    self._thread_safe_log(f"  → Features with <{min_required} non-zero values in a group will be skipped for that comparison\n")
            else:
                self._thread_safe_log(f"  → Filtering was already applied before normalization\n")
            
            from main_script.metabolite_statistics_analysis import get_group_metabolite_counts
            sample_cols = list(group_map.keys())
            
            # Note: Since pre-normalization filtering was already applied, 
            # the data in normalized_combined_df already has groups zeroed where insufficient
            # Just report current state
            self._thread_safe_log(f"\n📊 Data ready for statistical analysis (after pre-normalization filtering):\n")
            self._thread_safe_log(f"  Total metabolites in dataset: {len(self.normalized_combined_df)}\n")
            
            filtered_map, counts_by_group, excluded_groups = self._filter_groups_by_min_samples(group_map, min_required)
            if excluded_groups:
                excluded_summary = ', '.join(f"{grp} (n={cnt})" for grp, cnt in excluded_groups.items())
                self._thread_safe_log(f"\n⚠️ Excluding groups below minimum threshold: {excluded_summary}\n")
            included_groups = {grp: cnt for grp, cnt in counts_by_group.items() if grp not in excluded_groups}
            if len(included_groups) < 2:
                self._thread_safe_progress_step(0, 5, "Insufficient groups")
                self._thread_safe_log("❌ Not enough groups meet the minimum sample threshold.\n")
                messagebox.showerror(
                    'Insufficient Replicates',
                    f'At least two groups must have ≥ {min_required} samples. Update assignments or lower the threshold.'
                )
                return
            if hasattr(self, 'stat_base_group'):
                base_group_value = self.stat_base_group.get()
                if base_group_value and base_group_value not in included_groups:
                    self._thread_safe_log(f"⚠️ Base group '{base_group_value}' excluded due to insufficient samples.\n")
                    try:
                        self.stat_base_group.set('')
                    except Exception:
                        pass
            group_map = filtered_map
            sample_cols = list(group_map.keys())
            missing_cols = [c for c in sample_cols if c not in self.normalized_combined_df.columns]
            if missing_cols:
                self._thread_safe_log(f"❌ Missing sample columns: {missing_cols[:5]}...\n")
                messagebox.showerror('Missing Columns', f'Sample columns not found: {missing_cols[:5]}')
                return
            unique_groups = self.ordered_groups(list(dict.fromkeys(group_map.values())))
            included_summary = ', '.join(f"{grp}:{included_groups.get(grp, 0)}" for grp in unique_groups)
            self._thread_safe_log(f"Groups: {included_summary} | Samples: {len(sample_cols)} | Metabolites: {len(self.normalized_combined_df)}\n")
            # Log samples per group with user-assigned names (immediate feedback)
            from collections import Counter
            group_counts = Counter(group_map.values())
            samples_per_group = ', '.join([f"{group}: {count}" for group, count in sorted(group_counts.items())])
            self._thread_safe_log(f"📊 Samples per group: {samples_per_group}\n")
            
            self._thread_safe_progress_step(2, 5, "Preparing statistical test configuration...")
            from main_script.metabolite_statistics_analysis import perform_statistical_analysis, get_group_summary
            group_summary = get_group_summary(group_map)
            self._thread_safe_log(f"Group summary: {group_summary}\n")
            test_type = self.stat_test_type.get()
            overall_test = self.stat_overall_test.get().lower() if test_type == 'overall' else None
            pairwise_test = self.stat_pairwise_test.get().lower() if test_type == 'pairwise' else None
            self._thread_safe_progress_step(3, 5, "Running statistical tests (this may take a while)...")
            # Ordered group labels (allow explicit override)
            ordered_labels = []
            explicit_order_raw = self.statistics_group_order_var.get().strip() if hasattr(self, 'statistics_group_order_var') else ''
            if explicit_order_raw:
                provided = [p.strip() for p in explicit_order_raw.split(',') if p.strip()]
                seen = set()
                for lbl in provided:
                    if lbl in group_map.values() and lbl not in seen:
                        ordered_labels.append(lbl)
                        seen.add(lbl)
                # include any remaining groups in their default definition order
                for gid in self.group_definitions.keys():
                    lbl = self.group_definitions.get(gid, gid)
                    if lbl in group_map.values() and lbl not in seen:
                        ordered_labels.append(lbl)
                        seen.add(lbl)
                self._thread_safe_log(f"Applied explicit group order ({len(ordered_labels)}): {ordered_labels}\n")
            else:
                for gid in self.group_definitions.keys():
                    lbl = self.group_definitions.get(gid, gid)
                    if lbl in group_map.values() and lbl not in ordered_labels:
                        ordered_labels.append(lbl)
            # Store for visualization hand-off
            self.last_stats_group_order = ordered_labels.copy()
            # Parse custom comparisons if provided
            custom_comparisons = None
            if hasattr(self, 'custom_comparisons_var') and self.custom_comparisons_var.get().strip():
                custom_comp_str = self.custom_comparisons_var.get().strip()
                # Parse format: "Group1-Group2,Group3-Group4"
                custom_comparisons = []
                for pair in custom_comp_str.split(','):
                    pair = pair.strip()
                    if '-' in pair:
                        g1, g2 = pair.split('-', 1)
                        custom_comparisons.append((g1.strip(), g2.strip()))
                # Only use custom comparisons if list is non-empty; empty list = all pairs
                if custom_comparisons:
                    self._thread_safe_log(f"🎯 Custom comparisons ACTIVE: {custom_comparisons}\n")
                else:
                    custom_comparisons = None  # Treat empty list as None
                    self._thread_safe_log(f"⚠️ Custom comparisons field had text but no valid pairs parsed - using ALL pairwise\n")
            else:
                self._thread_safe_log(f"🔄 Pairwise mode: ALL combinations (custom comparisons field empty)\n")
            
            # Get use_adj_p setting
            use_adj_p = True  # Always use adjusted p-values
            
            # Prepare min_group_size parameters based on filtering type
            if filter_timing == 'after' and min_type == 'percentage':
                # Pass percentage for dynamic per-group calculation
                min_group_size_param = 1  # Use 1 as fallback for group-level filtering
                min_group_size_percent_param = min_percent
            else:
                # Use absolute count
                min_group_size_param = min_required
                min_group_size_percent_param = None
            
            # Determine whether to use adjusted p-values in enhanced output
            # If user selected 'None' for pairwise adjustment, force raw p-values in enhanced tables
            _pair_adj_method = self.pairwise_p_adjust_method.get() if hasattr(self, 'pairwise_p_adjust_method') else 'BH'
            _use_adj_for_enhanced = (use_adj_p and (_pair_adj_method != 'None'))

            # Two-Way ANOVA: Check if factor assignments are configured
            if overall_test == 'two_way_anova':
                if not hasattr(self, 'sample_factorA_vars') or not self.sample_factorA_vars:
                    self._thread_safe_log("❌ Two-Way ANOVA selected but factor assignments not configured.\n")
                    self._thread_safe_log("➡️ Open 'Configure Two-Way ANOVA' console, set up factors, click Confirm, then run again.\n")
                    self._thread_safe_progress_step(4,5,"Factor setup required")
                    return
                self._thread_safe_log("▶ Running Two-Way ANOVA (configured via console)...\n")
                self._thread_safe_progress_step(4,5,"Executing Two-Way ANOVA")
                self._run_two_way_anova_threaded()
                return
            
            # Non-Parametric Two-Way ANOVA
            elif overall_test == 'nonparametric_two_way':
                if not hasattr(self, 'sample_factorA_vars') or not self.sample_factorA_vars:
                    self._thread_safe_log("❌ Non-Parametric Two-Way ANOVA selected but factor assignments not configured.\n")
                    self._thread_safe_log("➡️ Open 'Configure Two-Way ANOVA' console, set up factors, click Confirm, then run again.\n")
                    self._thread_safe_progress_step(4,5,"Factor setup required")
                    return
                
                method = self.nonparam_method.get() if hasattr(self, 'nonparam_method') else 'art'
                self._thread_safe_log(f"▶ Running Non-Parametric Two-Way ANOVA (method: {method.upper()})...\n")
                self._thread_safe_progress_step(4,5,f"Executing {method.upper()} analysis")
                
                try:
                    from main_script.nonparametric_two_way_anova import run_nonparametric_twoway_from_gui
                    
                    # Build factor maps from GUI configuration
                    factor_a_map = {col: var.get() for col, var in self.sample_factorA_vars.items() if var.get()}
                    factor_b_map = {col: var.get() for col, var in self.sample_factorB_vars.items() if var.get()}
                    
                    # Get group map from user-defined group assignments (Control, HFD, TBI, HFD_TBI)
                    group_map = self._parse_group_assignments()
                    
                    # Get verified ID column name
                    id_col_name = self._get_verified_id_column(self.normalized_combined_df)
                    
                    # Verify column exists in dataframe
                    if id_col_name not in self.normalized_combined_df.columns:
                        # Fallback to first non-numeric column
                        for col in self.normalized_combined_df.columns:
                            if not pd.api.types.is_numeric_dtype(self.normalized_combined_df[col]):
                                id_col_name = col
                                break
                        if id_col_name not in self.normalized_combined_df.columns:
                            raise ValueError(f"Could not find verified ID column '{id_col_name}' in data")
                    
                    self._thread_safe_log(f'Using ID column: {id_col_name}\n')
                    
                    # Progress callback
                    def nonparam_progress(current, total, message):
                        percent = int(100 * current / total) if total > 0 else 0
                        if current % max(1, total // 20) == 0 or current == total:
                            self._thread_safe_log(f'  [{current}/{total} - {percent}%] {message}\n')
                    
                    results = run_nonparametric_twoway_from_gui(
                        data=self.normalized_combined_df,
                        sample_cols=sample_cols,
                        sample_factorA_map=factor_a_map,
                        sample_factorB_map=factor_b_map,
                        group_map=group_map,
                        method=method,
                        metabolite_id_col=id_col_name,
                        progress_callback=nonparam_progress
                    )
                    
                    # Store results in the same format as parametric two-way ANOVA
                    self.nonparametric_twoway_results = results
                    
                    # Format results for export (EXACT same pattern as parametric two-way ANOVA)
                    summary_df = results.get('summary')
                    posthoc_results = results.get('posthoc', {})
                    
                    if summary_df is not None and not summary_df.empty:
                        self._thread_safe_log(f"\n✅ Non-Parametric Two-Way ANOVA Complete!\n")
                        self._thread_safe_log(f"   Method: {results['method']}\n")
                        self._thread_safe_log(f"   Metabolites analyzed: {len(summary_df)}\n")
                        
                        # Count significant effects
                        sig_a = summary_df['factor_a_significant'].sum() if 'factor_a_significant' in summary_df.columns else 0
                        sig_b = summary_df['factor_b_significant'].sum() if 'factor_b_significant' in summary_df.columns else 0
                        sig_ab = summary_df['interaction_significant'].sum() if 'interaction_significant' in summary_df.columns else 0
                        
                        self._thread_safe_log(f"   Significant Factor A effects: {sig_a}\n")
                        self._thread_safe_log(f"   Significant Factor B effects: {sig_b}\n")
                        self._thread_safe_log(f"   Significant Interactions: {sig_ab}\n")
                        
                        # Get factor names and group info
                        fa_name = getattr(self, 'factorA_name', 'Factor A')
                        fb_name = getattr(self, 'factorB_name', 'Factor B')
                        
                        # Rename summary columns to match parametric format
                        # factor_a_pvalue -> p_<FactorName>, factor_a_significant -> remove
                        column_mapping = {
                            'factor_a_pvalue': f'p_{fa_name}',
                            'factor_b_pvalue': f'p_{fb_name}',
                            'interaction_pvalue': 'p_interaction',
                            'factor_a_stat': f'stat_{fa_name}',
                            'factor_b_stat': f'stat_{fb_name}',
                            'interaction_stat': 'stat_interaction'
                        }
                        
                        complete_result = summary_df.copy()
                        complete_result = complete_result.rename(columns=column_mapping)
                        
                        # Remove _significant columns and method/n_samples columns
                        cols_to_drop = [c for c in complete_result.columns if c.endswith('_significant') or c in ['method', 'n_samples']]
                        complete_result = complete_result.drop(columns=cols_to_drop, errors='ignore')
                        
                        # Add adjusted p-values using FDR correction across all metabolites
                        from statsmodels.stats.multitest import multipletests
                        for factor_col in [f'p_{fa_name}', f'p_{fb_name}', 'p_interaction']:
                            if factor_col in complete_result.columns:
                                adj_col = factor_col.replace('p_', 'adj_p_')
                                mask = complete_result[factor_col].notna()
                                if mask.any():
                                    try:
                                        _, padj_vec, _, _ = multipletests(
                                            complete_result.loc[mask, factor_col],
                                            method='fdr_bh'
                                        )
                                        complete_result.loc[mask, adj_col] = padj_vec
                                    except Exception:
                                        complete_result[adj_col] = np.nan
                                else:
                                    complete_result[adj_col] = np.nan

                        # Attach raw sample columns so Complete Result mirrors the parametric format
                        sample_cols_present = [col for col in sample_cols if col in self.normalized_combined_df.columns]
                        if sample_cols_present:
                            sample_block = self.normalized_combined_df[[id_col_name] + sample_cols_present].copy()
                            # Deduplicate in case the source table contains repeated metabolite rows
                            sample_block = sample_block.drop_duplicates(subset=[id_col_name])
                            complete_result = complete_result.merge(sample_block, on=id_col_name, how='left')
                        
                        # Use GROUP IDs from sample assignments (not concatenated factors)
                        # This ensures sheet names match user-defined groups like "Control", "HFD", "TBI", "HFD_TBI"
                        group_map = self._parse_group_assignments()
                        
                        # Get unique groups in user-defined order (respects Group1, Group2, Group3, Group4 order)
                        # This ensures comparisons are Control_vs_HFD (not HFD_vs_Control)
                        if hasattr(self, 'group_definitions'):
                            # Use group_definitions order (Group1, Group2, Group3, Group4)
                            unique_groups = [self.group_definitions[gid] for gid in sorted(self.group_definitions.keys()) if self.group_definitions[gid] in group_map.values()]
                        else:
                            unique_groups = sorted(set(group_map.values()))
                        
                        # Generate pairwise combinations respecting order (Group1 vs all others, then Group2 vs remaining, etc.)
                        pairwise_combos = []
                        for i, g1 in enumerate(unique_groups):
                            for g2 in unique_groups[i+1:]:
                                pairwise_combos.append((g1, g2))
                        
                        # Log the group order for user visibility
                        self._thread_safe_log(f"📋 Group comparison order: {' < '.join(unique_groups)}\n")
                        self._thread_safe_log(f"📊 Generating {len(pairwise_combos)} pairwise comparisons\n")
                        
                        # Build individual pairwise sheets from post-hoc results
                        # Match parametric format: Group1_vs_Group2 sheets with adj_p, FC, log2FC, neg_log10_adj_p
                        individual_sheets = {}
                        
                        # Create set of valid ordered comparison names to avoid duplicates (Group1_vs_Group2 but not Group2_vs_Group1)
                        valid_comp_names = {f"{g1}_vs_{g2}" for g1, g2 in pairwise_combos}
                        
                        if posthoc_results:
                            # Extract all pairwise comparisons from interaction post-hoc
                            for metabolite_id, posthoc_data in posthoc_results.items():
                                if 'posthoc_AB' in posthoc_data and posthoc_data['posthoc_AB'] is not None:
                                    ph_df = posthoc_data['posthoc_AB']
                                    
                                    # Get metabolite row from complete_result
                                    met_row = complete_result[complete_result[id_col_name] == metabolite_id]
                                    if met_row.empty:
                                        continue
                                    
                                    # For each pairwise comparison in post-hoc
                                    for idx, row in ph_df.iterrows():
                                        group1 = row['group1']
                                        group2 = row['group2']
                                        comp_name = f"{group1}_vs_{group2}"
                                        
                                        # Skip reverse pairs (only keep ordered pairs)
                                        if comp_name not in valid_comp_names:
                                            continue
                                        
                                        # Create columns matching parametric format
                                        if comp_name not in individual_sheets:
                                            individual_sheets[comp_name] = []
                                        
                                        # Build complete row with info + samples from both groups + stats
                                        pairwise_row = {}
                                        
                                        # Add all info columns
                                        for col in self.normalized_combined_df.columns:
                                            if col not in sample_cols:
                                                pairwise_row[col] = self.normalized_combined_df.loc[self.normalized_combined_df[id_col_name] == metabolite_id, col].iloc[0] if not self.normalized_combined_df[self.normalized_combined_df[id_col_name] == metabolite_id].empty else np.nan
                                        
                                        # Add samples from both groups
                                        samples_g1 = [s for s in sample_cols if group_map.get(s) == group1]
                                        samples_g2 = [s for s in sample_cols if group_map.get(s) == group2]
                                        
                                        # Get metabolite row from normalized data
                                        met_data = self.normalized_combined_df[self.normalized_combined_df[id_col_name] == metabolite_id]
                                        if not met_data.empty:
                                            met_data = met_data.iloc[0]
                                        
                                        # Add sample values and calculate group means
                                        g1_values = []
                                        g2_values = []
                                        for s in samples_g1:
                                            if s in self.normalized_combined_df.columns and not met_data.empty:
                                                val = met_data[s]
                                                pairwise_row[s] = val
                                                if pd.notna(val) and val != 0:
                                                    g1_values.append(val)
                                        for s in samples_g2:
                                            if s in self.normalized_combined_df.columns and not met_data.empty:
                                                val = met_data[s]
                                                pairwise_row[s] = val
                                                if pd.notna(val) and val != 0:
                                                    g2_values.append(val)
                                        
                                        # Add n per group
                                        pairwise_row[f'n_{group1}'] = len(g1_values)
                                        pairwise_row[f'n_{group2}'] = len(g2_values)
                                        
                                        # Add comparison stats (use p_adjusted from Dunn test)
                                        adj_p = row.get('p_adjusted', np.nan)
                                        pairwise_row[f'{comp_name}_adj_p'] = adj_p
                                        
                                        # Calculate fold change: FC = mean(group2) / mean(group1)
                                        if g1_values and g2_values:
                                            mean_g1 = np.mean(g1_values)
                                            mean_g2 = np.mean(g2_values)
                                            if mean_g1 > 0:
                                                fc = mean_g2 / mean_g1
                                                pairwise_row[f'{comp_name}_FC'] = fc
                                                if fc > 0:
                                                    pairwise_row[f'{comp_name}_log2FC'] = np.log2(fc)
                                                else:
                                                    pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                            else:
                                                pairwise_row[f'{comp_name}_FC'] = np.nan
                                                pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                        else:
                                            pairwise_row[f'{comp_name}_FC'] = np.nan
                                            pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                        
                                        # Calculate neg_log10 of adj_p
                                        if not np.isnan(adj_p) and adj_p > 0:
                                            pairwise_row[f'{comp_name}_neg_log10_adj_p'] = -np.log10(adj_p)
                                        else:
                                            pairwise_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
                                        
                                        individual_sheets[comp_name].append(pairwise_row)
                        
                        # Merge all pairwise stats back into complete_result (wide format)
                        # This creates columns like Control_vs_HFD_adj_p, Control_vs_HFD_log2FC, etc.
                        for comp_name, rows_list in individual_sheets.items():
                            if not rows_list:
                                continue
                            
                            # Convert list of rows to DataFrame
                            comp_df = pd.DataFrame(rows_list)
                            
                            # Extract pairwise stats for this comparison - only include columns that exist
                            pairwise_cols = [f'{comp_name}_adj_p', f'{comp_name}_FC', f'{comp_name}_log2FC', f'{comp_name}_neg_log10_adj_p']
                            existing_pairwise_cols = [c for c in pairwise_cols if c in comp_df.columns]
                            
                            if existing_pairwise_cols:
                                pairwise_df = comp_df[[id_col_name] + existing_pairwise_cols].drop_duplicates(subset=[id_col_name])
                                
                                # Merge into complete_result
                                complete_result = complete_result.merge(pairwise_df, on=id_col_name, how='left', suffixes=('', '_dup'))
                        
                        # Remove any duplicate columns created during merge
                        dup_cols = [c for c in complete_result.columns if c.endswith('_dup')]
                        if dup_cols:
                            complete_result = complete_result.drop(columns=dup_cols, errors='ignore')
                        
                        # Add n per group to complete_result (actual non-zero samples per metabolite)
                        for group in unique_groups:
                            group_samples = [s for s in sample_cols if group_map.get(s) == group]
                            valid_cols = [col for col in group_samples if col in complete_result.columns]
                            if valid_cols:
                                numeric_vals = complete_result[valid_cols].apply(pd.to_numeric, errors='coerce')
                                mask = numeric_vals.notna() & (numeric_vals != 0)
                                counts = mask.sum(axis=1).astype(int)
                                complete_result[f'n_{group}'] = counts
                            else:
                                complete_result[f'n_{group}'] = 0
                        
                        # Convert individual sheets to DataFrames and add to formatted_results
                        formatted_results = {
                            'complete_result': complete_result,
                            'enhanced_metabolites': complete_result.copy()  # For visualization compatibility
                        }
                        
                        # Get info columns (non-numeric columns)
                        info_cols = [c for c in complete_result.columns if not pd.api.types.is_numeric_dtype(complete_result[c]) or c == id_col_name]
                        
                        for comp_name, rows_list in individual_sheets.items():
                            if not rows_list:
                                continue
                            
                            # Convert list of rows to DataFrame
                            sheet_df = pd.DataFrame(rows_list)
                            
                            # Get samples for these two groups
                            g1, g2 = comp_name.split('_vs_')
                            samples_g1 = [s for s in sample_cols if group_map.get(s) == g1]
                            samples_g2 = [s for s in sample_cols if group_map.get(s) == g2]
                            
                            # Rename stat columns to simple names (without comparison prefix)
                            rename_map = {
                                f'{comp_name}_adj_p': 'adj_p',
                                f'{comp_name}_FC': 'FC',
                                f'{comp_name}_log2FC': 'log2FC',
                                f'{comp_name}_neg_log10_adj_p': 'neg_log10_adj_p'
                            }
                            sheet_df = sheet_df.rename(columns=rename_map)
                            
                            # Reorder columns: info + samples for this comparison + stats
                            cols_order = []
                            # Add info columns that exist in sheet_df
                            for col in info_cols:
                                if col in sheet_df.columns:
                                    cols_order.append(col)
                            # Add sample columns in order
                            for s in samples_g1 + samples_g2:
                                if s in sheet_df.columns:
                                    cols_order.append(s)
                            # Add stat columns (now with simple names)
                            for stat_col in ['adj_p', 'FC', 'log2FC', 'neg_log10_adj_p']:
                                if stat_col in sheet_df.columns:
                                    cols_order.append(stat_col)
                            
                            sheet_df = sheet_df[cols_order]
                            
                            # Drop rows with no p-value
                            if 'adj_p' in sheet_df.columns:
                                sheet_df = sheet_df.dropna(subset=['adj_p'])
                            
                            if not sheet_df.empty:
                                formatted_results[comp_name] = sheet_df
                        
                        # Store in statistical_test_results (same as parametric two-way ANOVA)
                        self.statistical_test_results = formatted_results
                        
                        # Store in memory for visualization tab
                        self.memory_store['statistical_test_results'] = self.statistical_test_results
                        
                        num_individual = len([k for k in formatted_results.keys() if k not in ('complete_result', 'enhanced_metabolites')])
                        self._thread_safe_log(f"✅ Results formatted: 1 Complete Result + {num_individual} pairwise sheets\n")
                        self._thread_safe_log("➡️ Use 'Export Statistical Results' to save results to Excel\n")
                        
                        # Store regular results count for completion message
                        regular_pairwise_count = num_individual
                    
                    # ================================================================================
                    # LIPID CLASS NON-PARAMETRIC TWO-WAY ANOVA (if applicable)
                    # ================================================================================
                    self._thread_safe_log(f"\n🔍 DEBUG: Checking lipid class analysis conditions...\n")
                    self._thread_safe_log(f"   mode = '{mode}'\n")
                    self._thread_safe_log(f"   has normalized_combined_class_df = {hasattr(self, 'normalized_combined_class_df')}\n")
                    self._thread_safe_log(f"   class_df is None = {not hasattr(self, 'normalized_combined_class_df') or self.normalized_combined_class_df is None}\n")
                    
                    if mode == 'lipid' and hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                        try:
                            self._thread_safe_log(f"\n{'='*80}\n")
                            self._thread_safe_log(f"🔬 Running Non-Parametric Two-Way ANOVA for LIPID CLASS data...\n")
                            self._thread_safe_log(f"{'='*80}\n")
                            
                            class_df = self.normalized_combined_class_df.copy()
                            
                            # Detect class ID column dynamically
                            class_id_col = None
                            for candidate in ['Class', 'Lipid_Class', 'Class_name']:
                                if candidate in class_df.columns:
                                    class_id_col = candidate
                                    break
                            
                            if class_id_col is None:
                                self._thread_safe_log(f"⚠️ Cannot find class ID column (Class, Lipid_Class, or Class_name). Skipping class analysis.\n")
                            else:
                                # Use verified class sample columns first; fall back to verified parent sample columns.
                                class_sample_cols = [
                                    c for c in self._get_verified_sample_columns_for_current_mode(mode, class_level=True)
                                    if c in class_df.columns
                                ]
                                if not class_sample_cols:
                                    class_sample_cols = [c for c in sample_cols if c in class_df.columns]
                                
                                # Filter to samples with group assignments
                                class_included = [s for s in class_sample_cols if s in group_map]
                                
                                if len(class_included) < 4:
                                    self._thread_safe_log(f"⚠️ Insufficient samples in class data for Non-Parametric Two-Way ANOVA (need ≥4, have {len(class_included)}). Skipping.\n")
                                else:
                                    self._thread_safe_log(f"✅ Found {len(class_included)} samples in class data with group assignments\n")
                                    
                                    # Build factor maps (same as regular lipids)
                                    factor_a_map = {col: var.get() for col, var in self.sample_factorA_vars.items() if var.get() and col in class_included}
                                    factor_b_map = {col: var.get() for col, var in self.sample_factorB_vars.items() if var.get() and col in class_included}
                                    
                                    # Progress callback for class analysis
                                    def class_nonparam_progress(current, total, message):
                                        percent = int(100 * current / total) if total > 0 else 0
                                        if current % max(1, total // 20) == 0 or current == total:
                                            self._thread_safe_log(f'  [Class {current}/{total} - {percent}%] {message}\n')
                                    
                                    # Run non-parametric two-way ANOVA on class data
                                    class_results = run_nonparametric_twoway_from_gui(
                                        data=class_df,
                                        sample_cols=class_included,
                                        sample_factorA_map=factor_a_map,
                                        sample_factorB_map=factor_b_map,
                                        group_map=group_map,
                                        method=method,
                                        metabolite_id_col=class_id_col,
                                        progress_callback=class_nonparam_progress
                                    )
                                    
                                    # Format results for export (same pattern as regular lipids)
                                    class_summary_df = class_results.get('summary')
                                    class_posthoc_results = class_results.get('posthoc', {})
                                    
                                    if class_summary_df is not None and not class_summary_df.empty:
                                        self._thread_safe_log(f"\n✅ Lipid Class Non-Parametric Two-Way ANOVA Complete!\n")
                                        self._thread_safe_log(f"   Method: {class_results['method']}\n")
                                        self._thread_safe_log(f"   Lipid classes analyzed: {len(class_summary_df)}\n")
                                        
                                        # Count significant effects
                                        sig_a = class_summary_df['factor_a_significant'].sum() if 'factor_a_significant' in class_summary_df.columns else 0
                                        sig_b = class_summary_df['factor_b_significant'].sum() if 'factor_b_significant' in class_summary_df.columns else 0
                                        sig_ab = class_summary_df['interaction_significant'].sum() if 'interaction_significant' in class_summary_df.columns else 0
                                        
                                        self._thread_safe_log(f"   Significant Factor A effects: {sig_a}\n")
                                        self._thread_safe_log(f"   Significant Factor B effects: {sig_b}\n")
                                        self._thread_safe_log(f"   Significant Interactions: {sig_ab}\n")
                                        
                                        # Get factor names
                                        fa_name = getattr(self, 'factorA_name', 'Factor A')
                                        fb_name = getattr(self, 'factorB_name', 'Factor B')
                                        
                                        # Rename summary columns to match parametric format
                                        column_mapping = {
                                            'factor_a_pvalue': f'p_{fa_name}',
                                            'factor_b_pvalue': f'p_{fb_name}',
                                            'interaction_pvalue': 'p_interaction',
                                            'factor_a_stat': f'stat_{fa_name}',
                                            'factor_b_stat': f'stat_{fb_name}',
                                            'interaction_stat': 'stat_interaction'
                                        }
                                        
                                        class_complete_result = class_summary_df.copy()
                                        class_complete_result = class_complete_result.rename(columns=column_mapping)
                                        
                                        # Remove _significant columns and method/n_samples columns
                                        cols_to_drop = [c for c in class_complete_result.columns if c.endswith('_significant') or c in ['method', 'n_samples']]
                                        class_complete_result = class_complete_result.drop(columns=cols_to_drop, errors='ignore')
                                        
                                        # Add adjusted p-values using FDR correction
                                        from statsmodels.stats.multitest import multipletests
                                        for factor_col in [f'p_{fa_name}', f'p_{fb_name}', 'p_interaction']:
                                            if factor_col in class_complete_result.columns:
                                                adj_col = factor_col.replace('p_', 'adj_p_')
                                                mask = class_complete_result[factor_col].notna()
                                                if mask.any():
                                                    try:
                                                        _, padj_vec, _, _ = multipletests(
                                                            class_complete_result.loc[mask, factor_col],
                                                            method='fdr_bh'
                                                        )
                                                        class_complete_result.loc[mask, adj_col] = padj_vec
                                                    except Exception:
                                                        class_complete_result[adj_col] = np.nan
                                                else:
                                                    class_complete_result[adj_col] = np.nan
                                        
                                        # Attach raw sample columns
                                        class_sample_cols_present = [col for col in class_included if col in class_df.columns]
                                        if class_sample_cols_present:
                                            sample_block = class_df[[class_id_col] + class_sample_cols_present].copy()
                                            sample_block = sample_block.drop_duplicates(subset=[class_id_col])
                                            class_complete_result = class_complete_result.merge(sample_block, on=class_id_col, how='left')
                                        
                                        # Get unique groups in order
                                        if hasattr(self, 'group_definitions'):
                                            unique_groups = [self.group_definitions[gid] for gid in sorted(self.group_definitions.keys()) if self.group_definitions[gid] in group_map.values()]
                                        else:
                                            unique_groups = sorted(set(group_map.values()))
                                        
                                        # Generate pairwise combinations
                                        pairwise_combos = []
                                        for i, g1 in enumerate(unique_groups):
                                            for g2 in unique_groups[i+1:]:
                                                pairwise_combos.append((g1, g2))
                                        
                                        # Build individual pairwise sheets from post-hoc results
                                        class_individual_sheets = {}
                                        
                                        # Create set of valid ordered comparison names to avoid duplicates
                                        valid_comp_names = {f"{g1}_vs_{g2}" for g1, g2 in pairwise_combos}
                                        
                                        if class_posthoc_results:
                                            for class_id, posthoc_data in class_posthoc_results.items():
                                                if 'posthoc_AB' in posthoc_data and posthoc_data['posthoc_AB'] is not None:
                                                    ph_df = posthoc_data['posthoc_AB']
                                                    
                                                    # Get class row from complete_result
                                                    class_row = class_complete_result[class_complete_result[class_id_col] == class_id]
                                                    if class_row.empty:
                                                        continue
                                                    
                                                    # For each pairwise comparison in post-hoc
                                                    for idx, row in ph_df.iterrows():
                                                        group1 = row['group1']
                                                        group2 = row['group2']
                                                        comp_name = f"{group1}_vs_{group2}"
                                                        
                                                        # Skip reverse pairs (only keep ordered pairs)
                                                        if comp_name not in valid_comp_names:
                                                            continue
                                                        
                                                        if comp_name not in class_individual_sheets:
                                                            class_individual_sheets[comp_name] = []
                                                        
                                                        # Build complete row
                                                        pairwise_row = {}
                                                        
                                                        # Add all info columns
                                                        for col in class_df.columns:
                                                            if col not in class_included:
                                                                pairwise_row[col] = class_df.loc[class_df[class_id_col] == class_id, col].iloc[0] if not class_df[class_df[class_id_col] == class_id].empty else np.nan
                                                        
                                                        # Add samples from both groups
                                                        samples_g1 = [s for s in class_included if group_map.get(s) == group1]
                                                        samples_g2 = [s for s in class_included if group_map.get(s) == group2]
                                                        
                                                        # Get class data
                                                        class_data = class_df[class_df[class_id_col] == class_id]
                                                        if not class_data.empty:
                                                            class_data = class_data.iloc[0]
                                                        
                                                        # Add sample values and calculate group means
                                                        g1_values = []
                                                        g2_values = []
                                                        for s in samples_g1:
                                                            if s in class_df.columns and not class_data.empty:
                                                                val = class_data[s]
                                                                pairwise_row[s] = val
                                                                if pd.notna(val) and val != 0:
                                                                    g1_values.append(val)
                                                        for s in samples_g2:
                                                            if s in class_df.columns and not class_data.empty:
                                                                val = class_data[s]
                                                                pairwise_row[s] = val
                                                                if pd.notna(val) and val != 0:
                                                                    g2_values.append(val)
                                                        
                                                        # Add n per group
                                                        pairwise_row[f'n_{group1}'] = len(g1_values)
                                                        pairwise_row[f'n_{group2}'] = len(g2_values)
                                                        
                                                        # Add comparison stats
                                                        adj_p = row.get('p_adjusted', np.nan)
                                                        pairwise_row[f'{comp_name}_adj_p'] = adj_p
                                                        
                                                        # Calculate fold change
                                                        if g1_values and g2_values:
                                                            mean_g1 = np.mean(g1_values)
                                                            mean_g2 = np.mean(g2_values)
                                                            if mean_g1 > 0:
                                                                fc = mean_g2 / mean_g1
                                                                pairwise_row[f'{comp_name}_FC'] = fc
                                                                if fc > 0:
                                                                    pairwise_row[f'{comp_name}_log2FC'] = np.log2(fc)
                                                                else:
                                                                    pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                                            else:
                                                                pairwise_row[f'{comp_name}_FC'] = np.nan
                                                                pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                                        else:
                                                            pairwise_row[f'{comp_name}_FC'] = np.nan
                                                            pairwise_row[f'{comp_name}_log2FC'] = np.nan
                                                        
                                                        # Calculate neg_log10 of adj_p
                                                        if not np.isnan(adj_p) and adj_p > 0:
                                                            pairwise_row[f'{comp_name}_neg_log10_adj_p'] = -np.log10(adj_p)
                                                        else:
                                                            pairwise_row[f'{comp_name}_neg_log10_adj_p'] = np.nan
                                                        
                                                        class_individual_sheets[comp_name].append(pairwise_row)
                                        
                                        # Merge all pairwise stats back into complete_result
                                        for comp_name, rows_list in class_individual_sheets.items():
                                            if not rows_list:
                                                continue
                                            
                                            comp_df = pd.DataFrame(rows_list)
                                            pairwise_cols = [
                                                f'{comp_name}_adj_p',
                                                f'{comp_name}_FC',
                                                f'{comp_name}_log2FC',
                                                f'{comp_name}_neg_log10_adj_p',
                                                f'{comp_name}_model_effect',
                                                f'{comp_name}_model_se',
                                                f'{comp_name}_ci_lower_95',
                                                f'{comp_name}_ci_upper_95',
                                            ]
                                            pairwise_df = comp_df[[class_id_col] + pairwise_cols]
                                            class_complete_result = class_complete_result.merge(pairwise_df, on=class_id_col, how='left', suffixes=('', '_dup'))
                                        
                                        # Add n per group to complete_result
                                        for group in unique_groups:
                                            group_samples = [s for s in class_included if group_map.get(s) == group]
                                            class_complete_result[f'n_{group}'] = len(group_samples)
                                        
                                        # Convert individual sheets to DataFrames
                                        class_formatted_results = {
                                            'complete_result': class_complete_result,
                                            'enhanced_metabolites': class_complete_result.copy()
                                        }
                                        
                                        # Get info columns
                                        info_cols = [c for c in class_complete_result.columns if not pd.api.types.is_numeric_dtype(class_complete_result[c]) or c == class_id_col]
                                        
                                        for comp_name, rows_list in class_individual_sheets.items():
                                            if not rows_list:
                                                continue
                                            
                                            sheet_df = pd.DataFrame(rows_list)
                                            
                                            # Get samples for these two groups
                                            g1, g2 = comp_name.split('_vs_')
                                            samples_g1 = [s for s in class_included if group_map.get(s) == g1]
                                            samples_g2 = [s for s in class_included if group_map.get(s) == g2]
                                            
                                            # Rename stat columns to simple names
                                            rename_map = {
                                                f'{comp_name}_adj_p': 'adj_p',
                                                f'{comp_name}_FC': 'FC',
                                                f'{comp_name}_log2FC': 'log2FC',
                                                f'{comp_name}_neg_log10_adj_p': 'neg_log10_adj_p'
                                            }
                                            sheet_df = sheet_df.rename(columns=rename_map)
                                            
                                            # Reorder columns
                                            cols_order = []
                                            for col in info_cols:
                                                if col in sheet_df.columns:
                                                    cols_order.append(col)
                                            for s in samples_g1 + samples_g2:
                                                if s in sheet_df.columns:
                                                    cols_order.append(s)
                                            for stat_col in ['adj_p', 'FC', 'log2FC', 'neg_log10_adj_p', 'model_effect', 'model_se', 'ci_lower_95', 'ci_upper_95']:
                                                if stat_col in sheet_df.columns:
                                                    cols_order.append(stat_col)
                                            
                                            sheet_df = sheet_df[cols_order]
                                            
                                            # Drop rows with no p-value
                                            if 'adj_p' in sheet_df.columns:
                                                sheet_df = sheet_df.dropna(subset=['adj_p'])
                                            
                                            if not sheet_df.empty:
                                                class_formatted_results[comp_name] = sheet_df
                                        
                                        # Store class results
                                        self.statistical_test_results_class = class_formatted_results
                                        
                                        num_class_individual = len([k for k in class_formatted_results.keys() if k not in ('complete_result', 'enhanced_metabolites')])
                                        self._thread_safe_log(f"✅ Lipid class results formatted: 1 Complete Result + {num_class_individual} pairwise sheets\n")
                                    else:
                                        self._thread_safe_log(f"⚠️ Non-parametric analysis returned no results for lipid class data\n")
                                        
                        except Exception as e:
                            self._thread_safe_log(f'⚠️ Lipid class Non-Parametric Two-Way ANOVA error: {e}\n')
                            import traceback
                            self._thread_safe_log(f'{traceback.format_exc()}\n')
                    
                    # Show completion message AFTER both regular and class analyses
                    self._thread_safe_progress_step(5,5,"Complete")
                    
                    completion_msg = f"Non-Parametric Two-Way ANOVA ({method.upper()}) completed successfully!\n\n"
                    if mode == 'lipid' and hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
                        completion_msg += f"✅ Regular lipid analysis: {regular_pairwise_count} pairwise comparisons\n"
                        completion_msg += f"✅ Lipid class analysis: Complete\n"
                    else:
                        completion_msg += f"{regular_pairwise_count} pairwise comparisons\n"
                    completion_msg += "\nUse 'Export Statistical Results' to save results."
                    
                    messagebox.showinfo("Complete", completion_msg)
                    
                except ImportError as ie:
                    self._thread_safe_log(f"❌ Non-parametric two-way ANOVA module not found: {ie}\n")
                    self._thread_safe_log("➡️ Ensure nonparametric_two_way_anova.py is in the project directory\n")
                    messagebox.showerror("Module Missing", f"Cannot run non-parametric analysis:\n{ie}")
                except Exception as e:
                    self._thread_safe_log(f"❌ Non-parametric two-way ANOVA failed: {e}\n")
                    import traceback
                    self._thread_safe_log(f"{traceback.format_exc()}\n")
                    messagebox.showerror("Analysis Failed", f"Non-parametric analysis error:\n{e}")
                
                return
            
            # Use new clean one-way ANOVA implementation for ANOVA and Kruskal tests
            elif overall_test in ('anova', 'kruskal'):
                # Determine worker count
                try:
                    workers = self._get_workers_count(getattr(self, 'stats_workers', None), default=3)
                except Exception:
                    workers = 3
                
                # Define progress callback
                rots_pairwise = (pairwise_test == 'rots')
                def stats_progress(current, total, metabolite_name):
                    percent = int(100 * current / total) if total > 0 else 0
                    step = 20 if rots_pairwise else 100
                    if current == 1 and rots_pairwise:
                        self._thread_safe_log("🌀 ROTS bootstrap running (serial execution). Progress updates every ~5%.\n")
                    if current % max(1, total // step) == 0 or current == total:
                        prefix = '🌀 ROTS' if rots_pairwise else '  '
                        self._thread_safe_log(f'{prefix} [{current}/{total} - {percent}%] {metabolite_name}\n')
                
                self._thread_safe_log(f"▶ Running One-Way {overall_test.upper()} using clean implementation...\n")
                self._thread_safe_progress_step(4, 5, f"Executing One-Way {overall_test.upper()}")
                
                try:
                    from main_script.one_way_anova import perform_one_way_anova
                    
                    complete_result, individual_sheets = perform_one_way_anova(
                        self.normalized_combined_df,
                        sample_cols,
                        group_map,
                        group_order=ordered_labels,
                        overall_test=overall_test,
                        pairwise_test=pairwise_test,
                        drop_zeros=True,
                        min_group_size=min_group_size_param,
                        fdr=use_adj_p,
                        n_jobs=workers,
                        progress_callback=stats_progress
                    )
                    
                    # Convert to old format for compatibility with existing GUI code
                    # The old perform_statistical_analysis returns: {'overall': df, 'pairwise': df, 'enhanced': df}
                    # Our new implementation returns: (complete_result, individual_sheets_dict)
                    
                    # Build compatible results structure
                    results = {
                        'overall': complete_result.copy(),
                        'pairwise': complete_result.copy(),  # Same data, will be filtered differently
                        'enhanced': complete_result.copy(),
                        'enhanced_metabolites': complete_result.copy(),  # Add this key for visualization compatibility
                        'individual_sheets': individual_sheets
                    }
                    
                    self._thread_safe_log(f"✅ One-Way {overall_test.upper()} complete: {len(complete_result)} metabolites\n")
                    
                    # Store regular results count for completion message (shown after class analysis)
                    regular_metabolites_count = len(complete_result)
                    regular_pairwise_count = len(individual_sheets) if individual_sheets else 0
                    
                except ImportError as ie:
                    self._thread_safe_log(f"⚠️ Clean one-way ANOVA module not found, falling back to legacy implementation: {ie}\n")
                    # Fallback to old implementation
                    from main_script.metabolite_statistics_analysis import perform_statistical_analysis
                    
                    # Define progress callback for legacy implementation
                    rots_pairwise = (pairwise_test == 'rots')
                    def stats_progress(current, total, metabolite_name):
                        percent = int(100 * current / total) if total > 0 else 0
                        step = 20 if rots_pairwise else 100
                        if current == 1 and rots_pairwise:
                            self._thread_safe_log("🌀 ROTS bootstrap running (serial execution). Progress updates every ~5%.\n")
                        if current % max(1, total // step) == 0 or current == total:
                            prefix = '🌀 ROTS' if rots_pairwise else '  '
                            self._thread_safe_log(f'{prefix} [{current}/{total} - {percent}%] {metabolite_name}\n')
                    
                    # Determine ID column name based on mode
                    mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
                    id_col_name = self._get_verified_id_column(self.normalized_combined_df)
                    
                    # Get ROTS parameters if using ROTS
                    rots_params = self._get_rots_parameters()
                    if pairwise_test == 'rots':
                        seed_display = rots_params.get('rots_seed') if rots_params.get('rots_seed') is not None else 'auto'
                        self._thread_safe_log(
                            f"🧪 ROTS parameters -> B={rots_params.get('rots_B')}, K={rots_params.get('rots_K')}, alpha={rots_params.get('rots_alpha')}, seed={seed_display}\n"
                        )
                        self._thread_safe_log("   ROTS disables multithreading; expect longer runtimes.\n")
                    
                    results = perform_statistical_analysis(
                        self.normalized_combined_df,
                        sample_cols,
                        group_map,
                        group_order=ordered_labels,
                        overall_test=overall_test,
                        pairwise_test=pairwise_test,
                        fdr=use_adj_p,
                        base_group=self.stat_base_group.get() if hasattr(self, 'stat_base_group') and self.stat_base_group.get() else None,
                        custom_comparisons=custom_comparisons,
                        fdr_scope=self.fdr_scope_var.get() if hasattr(self, 'fdr_scope_var') else 'per-comparison',
                        alpha=0.05,
                        use_adjusted_pvalues=_use_adj_for_enhanced,
                        min_group_size=min_group_size_param,
                        min_group_size_percent=min_group_size_percent_param,
                        n_jobs=workers,
                        id_column_name=id_col_name,
                        **rots_params,
                        progress_callback=stats_progress
                    )
            
            else:
                # For pairwise-only tests or other test types, use legacy implementation
                # Determine worker count
                try:
                    workers = self._get_workers_count(getattr(self, 'stats_workers', None), default=3)
                except Exception:
                    workers = 3
                
                # Define progress callback for legacy analysis
                rots_pairwise = (pairwise_test == 'rots')
                def stats_progress(current, total, metabolite_name):
                    percent = int(100 * current / total) if total > 0 else 0
                    step = 20 if rots_pairwise else 100
                    if current == 1 and rots_pairwise:
                        self._thread_safe_log("🌀 ROTS bootstrap running (serial execution). Progress updates every ~5%.\n")
                    if current % max(1, total // step) == 0 or current == total:
                        prefix = '🌀 ROTS' if rots_pairwise else '  '
                        self._thread_safe_log(f'{prefix} [{current}/{total} - {percent}%] {metabolite_name}\n')
                
                # Determine ID column name based on mode
                mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
                id_col_name = 'LipidID' if mode == 'lipid' else 'metabolite'
                
                # DEBUG: Log actual test parameters
                self._thread_safe_log(f"🐛 DEBUG: test_type='{test_type}', overall_test='{overall_test}', pairwise_test='{pairwise_test}'\n")
                
                # Get ROTS parameters if using ROTS
                rots_params = self._get_rots_parameters()
                if pairwise_test == 'rots':
                    seed_display = rots_params.get('rots_seed') if rots_params.get('rots_seed') is not None else 'auto'
                    self._thread_safe_log(
                        f"🧪 ROTS parameters -> B={rots_params.get('rots_B')}, K={rots_params.get('rots_K')}, alpha={rots_params.get('rots_alpha')}, seed={seed_display}\n"
                    )
                    self._thread_safe_log("   ROTS disables multithreading; expect longer runtimes.\n")
                
                from main_script.metabolite_statistics_analysis import perform_statistical_analysis
                results = perform_statistical_analysis(
                    self.normalized_combined_df,
                    sample_cols,
                    group_map,
                    group_order=ordered_labels,
                    overall_test=overall_test,
                    pairwise_test=pairwise_test,
                    fdr=use_adj_p,
                    base_group=self.stat_base_group.get() if hasattr(self, 'stat_base_group') and self.stat_base_group.get() else None,
                    custom_comparisons=custom_comparisons,
                    fdr_scope=self.fdr_scope_var.get() if hasattr(self, 'fdr_scope_var') else 'per-comparison',
                    alpha=0.05,
                    use_adjusted_pvalues=_use_adj_for_enhanced,
                    min_group_size=min_group_size_param,
                    min_group_size_percent=min_group_size_percent_param,
                    n_jobs=workers,
                    id_column_name=id_col_name,
                    **rots_params,
                    progress_callback=stats_progress
                )
            # Apply chosen pairwise adjustment method (override BH) AFTER base analysis
            try:
                if 'pairwise' in results and not results['pairwise'].empty and hasattr(self, 'pairwise_p_adjust_method'):
                    method = self.pairwise_p_adjust_method.get()
                    if method and method != 'BH' and method != 'None':
                        self._thread_safe_log(f"Applying pairwise p-value adjustment method: {method}\n")
                        self._apply_custom_pairwise_adjustment(results['pairwise'], method, self.fdr_scope_var.get())
                    elif method == 'None':
                        # Explicitly remove adjusted p-values for pairwise results and base Expression on raw p-values
                        pw = results['pairwise']
                        try:
                            if 'p_value_adj' in pw.columns:
                                pw.drop(columns=['p_value_adj'], inplace=True, errors='ignore')
                            if 'neg_log10_p_adj' in pw.columns:
                                pw.drop(columns=['neg_log10_p_adj'], inplace=True, errors='ignore')
                            # Recompute Expression using raw p-values at alpha=0.05
                            if 'log2_fold_change' in pw.columns and 'p_value' in pw.columns:
                                pw['Expression'] = 'Not Significant'
                                mask_sig = (pd.to_numeric(pw['p_value'], errors='coerce') <= 0.05) & pw['log2_fold_change'].notna()
                                pw.loc[mask_sig & (pw['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
                                pw.loc[mask_sig & (pw['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
                            # Ensure raw -log10 is present
                            if 'p_value' in pw.columns:
                                pw['neg_log10_p'] = -np.log10(pd.to_numeric(pw['p_value'], errors='coerce').replace(0, np.finfo(float).eps))
                            results['pairwise'] = pw
                        except Exception as _e:
                            self._thread_safe_log(f"⚠️ Failed to enforce raw pairwise p-values (None): {_e}\n")
            except Exception as e:
                self._thread_safe_log(f"⚠️ Pairwise adjustment failed ({e}); retaining BH/raw values.\n")
            # Ensure completion message metrics exist even for pairwise-only runs
            try:
                if 'regular_metabolites_count' not in locals():
                    if 'overall' in results and hasattr(results['overall'], '__len__'):
                        regular_metabolites_count = len(results['overall'])
                    elif 'pairwise' in results and hasattr(results['pairwise'], 'columns'):
                        # Use ID column if present to count unique metabolites
                        id_col = self._get_verified_id_column(self.normalized_combined_df)
                        if id_col in results['pairwise'].columns:
                            regular_metabolites_count = results['pairwise'][id_col].nunique()
                        else:
                            regular_metabolites_count = len(results['pairwise'])
                if 'regular_pairwise_count' not in locals() and 'pairwise' in results and hasattr(results['pairwise'], '__len__'):
                    regular_pairwise_count = len(results['pairwise'])
            except Exception:
                # Do not block downstream messaging if counting fails
                pass
            self.statistical_test_results = results
            self._thread_safe_progress_step(4, 5, "Summarizing results...")
            
            # DEBUG: Show first 10 rows of pairwise results (only if pairwise columns exist)
            if 'pairwise' in results and not results['pairwise'].empty:
                # Check if this is a pairwise test (has group1/group2 columns)
                if 'group1' in results['pairwise'].columns and 'group2' in results['pairwise'].columns:
                    # Try to find Control vs TBI comparison for debugging
                    try:
                        control_tbi = results['pairwise'][(results['pairwise']['group1'] == 'Control') & (results['pairwise']['group2'] == 'TBI')]
                        if not control_tbi.empty:
                            self._thread_safe_log(f"🔍 DEBUG: First 10 rows of Control_vs_TBI pairwise results:\n{control_tbi.head(10).to_string()}\n")
                    except Exception as e:
                        # Silently skip if specific comparison not found
                        pass
            
            # Report overall results
            total_metabolites = len(self.normalized_combined_df)
            self._thread_safe_log(f"\n📊 Statistical Analysis Summary:\n")
            self._thread_safe_log(f"  Total metabolites in dataset: {total_metabolites}\n")
            
            if 'overall' in results and not results['overall'].empty:
                self._thread_safe_log(f"  Overall tests: {len(results['overall'])} metabolites tested\n")
            
            # Report pairwise comparison details
            if 'pairwise' in results and not results['pairwise'].empty:
                self._thread_safe_log(f"  Pairwise comparisons: {len(results['pairwise'])} total rows\n")
            # Propagate group order to visualization preferred order if not already customized
            try:
                if hasattr(self, 'viz_preferred_group_order') and hasattr(self, 'last_stats_group_order'):
                    current_viz_val = self.viz_preferred_group_order.get()
                    # Determine if current value is default placeholder (e.g., contains 'Control,Ortho,TBI') or empty
                    if (not current_viz_val.strip()) or current_viz_val.strip() in ('Group1, Group2, Group3', 'Control, Disease, Treatment'):
                        self.viz_preferred_group_order.set(','.join(self.last_stats_group_order))
                        self._thread_safe_log(f"Visualization group order set from statistics: {self.last_stats_group_order}\n")
            except Exception as e:
                self._thread_safe_log(f"Warning: could not propagate group order to visualization: {e}\n")
            
            # If in lipid mode, also run statistics on class data
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            if mode == 'lipid' and hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                self._thread_safe_log("\n🔬 Running statistical tests on lipid class data...\n")
                
                # Verify sample columns exist in class data
                class_missing_cols = [c for c in sample_cols if c not in self.normalized_combined_class_df.columns]
                if class_missing_cols:
                    self._thread_safe_log(f"⚠️ Warning: Some sample columns missing in class data: {class_missing_cols[:5]}...\n")
                    # Use only available sample columns
                    class_sample_cols = [
                        c for c in self._get_verified_sample_columns_for_current_mode(mode, class_level=True)
                        if c in self.normalized_combined_class_df.columns
                    ]
                    if not class_sample_cols:
                        class_sample_cols = [c for c in sample_cols if c in self.normalized_combined_class_df.columns]
                    if not class_sample_cols:
                        self._thread_safe_log("❌ No matching sample columns in class data, skipping class statistics.\n")
                        self.statistical_test_results_class = {}
                    else:
                        # Create filtered group map for available samples
                        class_group_map = {k: v for k, v in group_map.items() if k in class_sample_cols}
                        
                        # Count lipid classes with sufficient data per group
                        class_initial_counts = get_group_metabolite_counts(
                            self.normalized_combined_class_df,
                            class_sample_cols,
                            class_group_map,
                            min_required
                        )
                        
                        self._thread_safe_log(f"\n📊 Lipid class counts per group (n≥{min_required}):\n")
                        for grp in sorted(class_initial_counts.keys()):
                            count = class_initial_counts[grp]
                            self._thread_safe_log(f"  • {grp}: {count} classes available\n")
                        
                        class_group_map, _class_counts, class_excluded = self._filter_groups_by_min_samples(class_group_map, min_required)
                        if class_excluded:
                            excluded_info = ', '.join(f"{grp} (n={cnt})" for grp, cnt in class_excluded.items())
                            self._thread_safe_log(f"⚠️ Skipping lipid class groups below threshold: {excluded_info}\n")
                        class_sample_cols = list(class_group_map.keys())
                        if len(set(class_group_map.values())) < 2:
                            self._thread_safe_log("⚠️ Not enough lipid class groups meet the minimum threshold. Skipping class statistics.\n")
                            self.statistical_test_results_class = {}
                        else:
                            # Determine adjusted flag for enhanced outputs (respect 'None' as raw)
                            _pair_adj_method_cls = self.pairwise_p_adjust_method.get() if hasattr(self, 'pairwise_p_adjust_method') else 'BH'
                            _use_adj_for_enhanced_cls = (use_adj_p and (_pair_adj_method_cls != 'None'))

                            # Progress callback for lipid class analysis
                            rots_pairwise = (pairwise_test == 'rots')
                            def class_progress(current, total, class_name):
                                percent = int(100 * current / total) if total > 0 else 0
                                step = 20 if rots_pairwise else 10
                                if current == 1 and rots_pairwise:
                                    self._thread_safe_log("🌀 ROTS (class data) running serially. Progress updates every ~5%.\n")
                                if current % max(1, total // step) == 0 or current == total:
                                    prefix = '🌀 ROTS Class' if rots_pairwise else '  Class'
                                    self._thread_safe_log(f'{prefix} [{current}/{total} - {percent}%] {class_name}\n')
                            
                            # Get ROTS parameters if using ROTS
                            rots_params = self._get_rots_parameters()
                            if pairwise_test == 'rots':
                                seed_display = rots_params.get('rots_seed') if rots_params.get('rots_seed') is not None else 'auto'
                                self._thread_safe_log(
                                    f"🧪 ROTS parameters (class data) -> B={rots_params.get('rots_B')}, K={rots_params.get('rots_K')}, alpha={rots_params.get('rots_alpha')}, seed={seed_display}\n"
                                )

                            # Use verified ID column detection (same as regular lipids)
                            class_id_col = self._get_verified_id_column(self.normalized_combined_class_df)
                            
                            # If verified column doesn't exist, try common class column names
                            if class_id_col not in self.normalized_combined_class_df.columns:
                                for candidate in ['Lipid_Class', 'LipidClass', 'Class', 'Class_name']:
                                    if candidate in self.normalized_combined_class_df.columns:
                                        class_id_col = candidate
                                        self._thread_safe_log(f"📝 Using '{class_id_col}' as ID column for class data\n")
                                        break
                            
                            class_results = perform_statistical_analysis(
                                self.normalized_combined_class_df,
                                class_sample_cols,
                                class_group_map,
                                group_order=ordered_labels,
                                overall_test=overall_test,
                                pairwise_test=pairwise_test,
                                fdr=use_adj_p,
                                base_group=self.stat_base_group.get() if hasattr(self, 'stat_base_group') and self.stat_base_group.get() else None,
                                custom_comparisons=custom_comparisons,
                                fdr_scope=self.fdr_scope_var.get() if hasattr(self, 'fdr_scope_var') else 'per-comparison',
                                alpha=0.05,
                                use_adjusted_pvalues=_use_adj_for_enhanced_cls,
                                min_group_size=min_required,
                                n_jobs=workers,
                                id_column_name=class_id_col,
                                **rots_params,
                                progress_callback=class_progress
                            )
                            try:
                                if 'pairwise' in class_results and not class_results['pairwise'].empty and hasattr(self, 'pairwise_p_adjust_method'):
                                    method = self.pairwise_p_adjust_method.get()
                                    pairwise_columns = list(class_results['pairwise'].columns)
                                    self._thread_safe_log(
                                        f"  🔎 Class pairwise CI debug: ci_lower_95={'ci_lower_95' in pairwise_columns}, "
                                        f"ci_upper_95={'ci_upper_95' in pairwise_columns}, model_effect={'model_effect' in pairwise_columns}, "
                                        f"model_se={'model_se' in pairwise_columns}\n"
                                    )
                                    sample_ci_cols = [c for c in ['group1', 'group2', 'ci_lower_95', 'ci_upper_95', 'model_effect', 'model_se'] if c in pairwise_columns]
                                    if sample_ci_cols:
                                        try:
                                            self._thread_safe_log(
                                                f"  🔎 Class pairwise CI sample: {class_results['pairwise'][sample_ci_cols].head(3).to_dict(orient='records')}\n"
                                            )
                                        except Exception:
                                            pass
                                    if method and method != 'BH' and method != 'None':
                                        self._thread_safe_log(f"Applying lipid class pairwise p-value adjustment: {method}\n")
                                        self._apply_custom_pairwise_adjustment(class_results['pairwise'], method, self.fdr_scope_var.get())
                                    elif method == 'None':
                                        # Enforce raw pairwise p-values for lipid class results
                                        pwc = class_results['pairwise']
                                        try:
                                            if 'p_value_adj' in pwc.columns:
                                                pwc.drop(columns=['p_value_adj'], inplace=True, errors='ignore')
                                            if 'neg_log10_p_adj' in pwc.columns:
                                                pwc.drop(columns=['neg_log10_p_adj'], inplace=True, errors='ignore')
                                            if 'log2_fold_change' in pwc.columns and 'p_value' in pwc.columns:
                                                pwc['Expression'] = 'Not Significant'
                                                mask_sig = (pd.to_numeric(pwc['p_value'], errors='coerce') <= 0.05) & pwc['log2_fold_change'].notna()
                                                pwc.loc[mask_sig & (pwc['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
                                                pwc.loc[mask_sig & (pwc['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
                                            if 'p_value' in pwc.columns:
                                                pwc['neg_log10_p'] = -np.log10(pd.to_numeric(pwc['p_value'], errors='coerce').replace(0, np.finfo(float).eps))
                                            class_results['pairwise'] = pwc
                                        except Exception as _e:
                                            self._thread_safe_log(f"⚠️ Failed to enforce raw class pairwise p-values (None): {_e}\n")
                            except Exception as e:
                                self._thread_safe_log(f"⚠️ Lipid class pairwise adjustment failed ({e}); retaining BH/raw values.\n")
                            self.statistical_test_results_class = class_results
                            
                            # Report lipid class results
                            total_classes = len(self.normalized_combined_class_df)
                            self._thread_safe_log(f"\n📊 Lipid Class Statistical Analysis Summary:\n")
                            self._thread_safe_log(f"  Total lipid classes in dataset: {total_classes}\n")
                            
                            if 'overall' in class_results and not class_results['overall'].empty:
                                self._thread_safe_log(f"  Overall tests: {len(class_results['overall'])} classes tested\n")
                            
                            if 'pairwise_summary' in class_results and not class_results['pairwise_summary'].empty:
                                self._thread_safe_log(f"\n📊 Lipid Class Pairwise Comparison Details:\n")
                                
                                for _, row in class_results['pairwise_summary'].iterrows():
                                    comp = row['comparison']
                                    tested = row['tested_metabolites']
                                    skipped = row['skipped_insufficient_n']
                                    self._thread_safe_log(f"  • {comp}: {tested} classes tested")
                                    if skipped > 0:
                                        # Get filtering settings for accurate message
                                        min_samples = getattr(self, 'min_samples_per_group_var', None)
                                        min_type = getattr(self, 'min_samples_type_var', None)
                                        if min_samples and min_type:
                                            try:
                                                threshold = min_samples.get()
                                                threshold_type = min_type.get()
                                                if threshold_type == 'percentage':
                                                    self._thread_safe_log(f" ({skipped} skipped: insufficient data <{threshold}% per group)")
                                                else:
                                                    self._thread_safe_log(f" ({skipped} skipped: insufficient data <{threshold} samples per group)")
                                            except Exception:
                                                self._thread_safe_log(f" ({skipped} skipped: insufficient data per group)")
                                        else:
                                            self._thread_safe_log(f" ({skipped} skipped: insufficient data per group)")
                                    self._thread_safe_log(f"\n")
                            elif 'pairwise' in class_results and not class_results['pairwise'].empty:
                                self._thread_safe_log(f"  Pairwise comparisons: {len(class_results['pairwise'])} total rows\n")
                            
                            self._thread_safe_log("✅ Lipid class statistical analysis complete!\n")
                else:
                    # All sample columns present
                    _pair_adj_method_cls = self.pairwise_p_adjust_method.get() if hasattr(self, 'pairwise_p_adjust_method') else 'BH'
                    _use_adj_for_enhanced_cls = (use_adj_p and (_pair_adj_method_cls != 'None'))

                    # Progress callback for lipid class analysis (all columns present path)
                    rots_pairwise = (pairwise_test == 'rots')
                    def class_progress_all(current, total, class_name):
                        percent = int(100 * current / total) if total > 0 else 0
                        step = 20 if rots_pairwise else 10
                        if current == 1 and rots_pairwise:
                            self._thread_safe_log("🌀 ROTS (class data) running serially. Progress updates every ~5%.\n")
                        if current % max(1, total // step) == 0 or current == total:
                            prefix = '🌀 ROTS Class' if rots_pairwise else '  Class'
                            self._thread_safe_log(f'{prefix} [{current}/{total} - {percent}%] {class_name}\n')
                    
                    # Get ROTS parameters if using ROTS
                    rots_params = self._get_rots_parameters()
                    if pairwise_test == 'rots':
                        seed_display = rots_params.get('rots_seed') if rots_params.get('rots_seed') is not None else 'auto'
                        self._thread_safe_log(
                            f"🧪 ROTS parameters (class data) -> B={rots_params.get('rots_B')}, K={rots_params.get('rots_K')}, alpha={rots_params.get('rots_alpha')}, seed={seed_display}\n"
                        )

                    # Use verified ID column detection (same as regular lipids)
                    class_id_col = self._get_verified_id_column(self.normalized_combined_class_df)
                    
                    # If verified column doesn't exist, try common class column names
                    if class_id_col not in self.normalized_combined_class_df.columns:
                        for candidate in ['Lipid_Class', 'LipidClass', 'Class', 'Class_name']:
                            if candidate in self.normalized_combined_class_df.columns:
                                class_id_col = candidate
                                self._thread_safe_log(f"📝 Using '{class_id_col}' as ID column for class data\n")
                                break
                    
                    class_results = perform_statistical_analysis(
                        self.normalized_combined_class_df,
                        sample_cols,
                        group_map,
                        group_order=ordered_labels,
                        overall_test=overall_test,
                        pairwise_test=pairwise_test,
                        fdr=use_adj_p,
                        base_group=self.stat_base_group.get() if hasattr(self, 'stat_base_group') and self.stat_base_group.get() else None,
                        custom_comparisons=custom_comparisons,
                        fdr_scope=self.fdr_scope_var.get() if hasattr(self, 'fdr_scope_var') else 'per-comparison',
                        alpha=0.05,
                        use_adjusted_pvalues=_use_adj_for_enhanced_cls,
                        min_group_size=min_required,
                        n_jobs=workers,
                        id_column_name=class_id_col,
                        **rots_params,
                        progress_callback=class_progress_all
                    )
                    try:
                        if 'pairwise' in class_results and not class_results['pairwise'].empty and hasattr(self, 'pairwise_p_adjust_method'):
                            method = self.pairwise_p_adjust_method.get()
                            if method and method != 'BH' and method != 'None':
                                self._thread_safe_log(f"Applying lipid class pairwise p-value adjustment: {method}\n")
                                self._apply_custom_pairwise_adjustment(class_results['pairwise'], method, self.fdr_scope_var.get())
                            elif method == 'None':
                                pwc = class_results['pairwise']
                                try:
                                    if 'p_value_adj' in pwc.columns:
                                        pwc.drop(columns=['p_value_adj'], inplace=True, errors='ignore')
                                    if 'neg_log10_p_adj' in pwc.columns:
                                        pwc.drop(columns=['neg_log10_p_adj'], inplace=True, errors='ignore')
                                    if 'log2_fold_change' in pwc.columns and 'p_value' in pwc.columns:
                                        pwc['Expression'] = 'Not Significant'
                                        mask_sig = (pd.to_numeric(pwc['p_value'], errors='coerce') <= 0.05) & pwc['log2_fold_change'].notna()
                                        pwc.loc[mask_sig & (pwc['log2_fold_change'] > 0), 'Expression'] = 'Upregulated'
                                        pwc.loc[mask_sig & (pwc['log2_fold_change'] < 0), 'Expression'] = 'Downregulated'
                                    if 'p_value' in pwc.columns:
                                        pwc['neg_log10_p'] = -np.log10(pd.to_numeric(pwc['p_value'], errors='coerce').replace(0, np.finfo(float).eps))
                                    class_results['pairwise'] = pwc
                                except Exception as _e:
                                    self._thread_safe_log(f"⚠️ Failed to enforce raw class pairwise p-values (None): {_e}\n")
                    except Exception as e:
                        self._thread_safe_log(f"⚠️ Lipid class pairwise adjustment failed ({e}); retaining BH/raw values.\n")
                    self.statistical_test_results_class = class_results
                    
                    # Report lipid class results
                    total_classes = len(self.normalized_combined_class_df)
                    self._thread_safe_log(f"\n📊 Lipid Class Statistical Analysis Summary:\n")
                    self._thread_safe_log(f"  Total lipid classes in dataset: {total_classes}\n")
                    
                    if 'overall' in class_results and not class_results['overall'].empty:
                        self._thread_safe_log(f"  Overall tests: {len(class_results['overall'])} classes tested\n")
                    
                    if 'pairwise_summary' in class_results and not class_results['pairwise_summary'].empty:
                        self._thread_safe_log(f"\n📊 Lipid Class Pairwise Comparison Details:\n")
                        
                        for _, row in class_results['pairwise_summary'].iterrows():
                            comp = row['comparison']
                            tested = row['tested_metabolites']
                            skipped = row['skipped_insufficient_n']
                            self._thread_safe_log(f"  • {comp}: {tested} classes tested")
                            if skipped > 0:
                                # Get filtering settings for accurate message
                                min_samples = getattr(self, 'min_samples_per_group_var', None)
                                min_type = getattr(self, 'min_samples_type_var', None)
                                if min_samples and min_type:
                                    try:
                                        threshold = min_samples.get()
                                        threshold_type = min_type.get()
                                        if threshold_type == 'percentage':
                                            self._thread_safe_log(f" ({skipped} skipped: insufficient data <{threshold}% per group)")
                                        else:
                                            self._thread_safe_log(f" ({skipped} skipped: insufficient data <{threshold} samples per group)")
                                    except Exception:
                                        self._thread_safe_log(f" ({skipped} skipped: insufficient data per group)")
                                else:
                                    self._thread_safe_log(f" ({skipped} skipped: insufficient data per group)")
                            self._thread_safe_log(f"\n")
                    elif 'pairwise' in class_results and not class_results['pairwise'].empty:
                        self._thread_safe_log(f"  Pairwise comparisons: {len(class_results['pairwise'])} total rows\n")
                    
                    self._thread_safe_log("✅ Lipid class statistical analysis complete!\n")
            
            # Show single completion message AFTER both regular and class analyses complete
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            # Build a safe label for completion messages
            def _test_label():
                if overall_test:
                    return overall_test.upper()
                if pairwise_test:
                    return pairwise_test.upper()
                return "ANALYSIS"

            if mode == 'lipid' and hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
                # Both regular and class analyses completed
                completion_msg = f"{_test_label()} completed successfully!\n\n"
                completion_msg += f"✅ Regular lipid analysis: {locals().get('regular_metabolites_count', 'N/A')} metabolites\n"
                completion_msg += f"✅ Lipid class analysis: Complete\n\n"
                completion_msg += "Use 'Export Statistical Results' to save results."
                messagebox.showinfo("Analysis Complete", completion_msg)
            elif locals().get('regular_metabolites_count') is not None:
                # Only regular analysis completed (no class data)
                messagebox.showinfo("Complete", 
                    f"{_test_label()} completed successfully!\n\n"
                    f"{regular_metabolites_count} metabolites analyzed\n"
                    f"{locals().get('regular_pairwise_count', 0)} pairwise comparisons\n\n"
                    f"Use 'Export Statistical Results' to save.")
            
            self._thread_safe_progress_step(5, 5, "Statistical analysis complete")
            self._thread_safe_log("✅ Statistical analysis complete!\n")
            
            # Determine test name for completion message
            test_name = "Statistical Analysis"
            if test_type == 'overall':
                if overall_test == 'anova':
                    test_name = "One-Way ANOVA"
                elif overall_test == 'kruskal':
                    test_name = "Kruskal-Wallis Test"
            elif test_type == 'pairwise':
                if pairwise_test == 'welch':
                    test_name = "Welch's t-test"
                elif pairwise_test == 'mannwhitney':
                    test_name = "Mann-Whitney U Test"
                elif pairwise_test == 'rots':
                    test_name = "ROTS Analysis"
                elif pairwise_test == 'limma':
                    test_name = "Limma (Moderated t-test)"
            
            # Count results
            # For pairwise tests, count unique metabolites from pairwise results
            # For overall tests, count from overall/enhanced results
            n_metabolites = 0
            n_comparisons = 0
            
            # First, try to get pairwise comparison count if pairwise data exists
            if 'pairwise' in results and not results['pairwise'].empty:
                pairwise_df = results['pairwise']
                
                # Check what columns actually exist in pairwise results
                if 'group1' in pairwise_df.columns and 'group2' in pairwise_df.columns:
                    n_comparisons = pairwise_df.groupby(['group1', 'group2']).ngroups
                elif 'comparison' in pairwise_df.columns:
                    # Alternative: count unique comparisons from a 'comparison' column
                    n_comparisons = pairwise_df['comparison'].nunique()
                else:
                    # No clear comparison structure - leave as 0
                    n_comparisons = 0
            
            # Now count metabolites based on test type
            if test_type == 'pairwise':
                # For pure pairwise tests, count from pairwise results
                if 'pairwise' in results and not results['pairwise'].empty:
                    pairwise_df = results['pairwise']
                    if 'metabolite' in pairwise_df.columns:
                        n_metabolites = pairwise_df['metabolite'].nunique()
                    else:
                        # Fallback: estimate from total rows divided by comparisons
                        n_metabolites = len(pairwise_df) // max(1, n_comparisons) if n_comparisons > 0 else len(pairwise_df)
            else:
                # For overall tests (ANOVA, Kruskal-Wallis), count from overall/enhanced results
                n_metabolites = len(results.get('enhanced', results.get('overall', pd.DataFrame())))
            
            # Skip showing completion message here - already shown in one-way ANOVA path above
            # messagebox.showinfo('Complete', f'{test_name} completed successfully!\n\n{n_metabolites} metabolites analyzed\n\nUse \'Export Statistical Results\' to save.')
            
            # Store results in memory_store for Visualization tab to access
            try:
                self.memory_store['statistical_test_results'] = self.statistical_test_results
                if hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
                    self.memory_store['statistical_test_results_class'] = self.statistical_test_results_class
                if hasattr(self, 'normalized_combined_df') and self.normalized_combined_df is not None:
                    self.memory_store['normalized_combined_df'] = self.normalized_combined_df
                if hasattr(self, 'sample_to_group') and self.sample_to_group:
                    self.memory_store['sample_to_group'] = self.sample_to_group
                
                # Store ID column information for visualization
                mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
                if mode == 'lipid':
                    if hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
                        # Lipid class data uses 'Class' as ID column
                        self.memory_store['id_column_class'] = 'Class'
                    # Regular lipid data uses LipidID
                    self.memory_store['id_column'] = 'LipidID'
                else:
                    # Metabolite data uses 'Name'
                    self.memory_store['id_column'] = 'Name'
                
                self._thread_safe_log("✅ Stored results in shared memory\n")
            except Exception as e:
                self._thread_safe_log(f"⚠️ Could not store results in memory: {e}\n")
            
            # Auto-notify Visualization tab without switching (switch happens on export)
            try:
                self._thread_safe_log("\n→ Auto-loading data to Visualization tab (no auto-switch)...\n")
                self.root.after(500, lambda: self._auto_load_to_visualization())
            except Exception as e:
                print(f"Warning: Could not auto-load to visualization: {e}")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self._thread_safe_log(f"❌ Statistical tests error: {e}\n{tb}\n")
            # Optionally write to a log file in output directory if available
            try:
                if hasattr(self, 'viz_output_dir') and self.viz_output_dir.get():
                    log_path = os.path.join(self.viz_output_dir.get(), 'stat_tests_error.log')
                    with open(log_path, 'a', encoding='utf-8') as lf:
                        lf.write('\n--- Statistical Test Error ---\n')
                        lf.write(tb)
                        lf.write('\n')
            except Exception:
                pass
            messagebox.showerror('Statistical Analysis Error', f"{e}\nSee log / console for traceback.")
        finally:
            self.root.after(0, self.hide_stats_progress)
            self._re_enable_stat_tests_button()

    def _generate_normality_figure(self, normality_results, polarity):
        """Generate and save a figure visualizing normality test results.
        
        Args:
            normality_results: DataFrame with normality test results
            polarity: 'positive', 'negative', or 'merged'
            
        Returns:
            Path to saved figure, or None if failed
        """
        # Guard against missing or invalid results
        try:
            if normality_results is None or not hasattr(normality_results, 'empty') or normality_results.empty:
                return None
            required_cols = {'Sample/Group', 'P_Value', 'Is_Normal'}
            if not required_cols.issubset(set(normality_results.columns)):
                return None
        except Exception:
            return None
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from datetime import datetime
        
        # Create visualizations directory if it doesn't exist
        viz_dir = os.path.join(os.getcwd(), 'visualizations')
        os.makedirs(viz_dir, exist_ok=True)
        
        # Filter out rows with insufficient data
        valid_results = normality_results[normality_results['Is_Normal'].isin(['Yes', 'No'])].copy()
        
        if valid_results.empty:
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, max(6, len(valid_results) * 0.4)))
        
        # Extract data
        samples = valid_results['Sample/Group'].values
        p_values = valid_results['P_Value'].values
        is_normal = valid_results['Is_Normal'].values
        
        # Color code bars: green if normal (p > 0.05), red if not normal
        colors = ['#2ecc71' if status == 'Yes' else '#e74c3c' for status in is_normal]
        
        # Create horizontal bar chart
        y_pos = np.arange(len(samples))
        bars = ax.barh(y_pos, p_values, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        
        # Add reference line at p = 0.05
        ax.axvline(x=0.05, color='black', linestyle='--', linewidth=2, label='Significance threshold (p = 0.05)')
        
        # Customize plot
        ax.set_yticks(y_pos)
        ax.set_yticklabels(samples, fontsize=9)
        ax.set_xlabel('P-Value', fontsize=12, fontweight='bold')
        ax.set_ylabel('Sample/Group', fontsize=12, fontweight='bold')
        ax.set_title(f'Normality Test Results - {polarity.capitalize()} Mode\n(Shapiro-Wilk or Kolmogorov-Smirnov Test)', 
                     fontsize=14, fontweight='bold', pad=20)
        
        # Add grid for easier reading
        ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
        ax.set_axisbelow(True)
        
        # Set x-axis limits
        max_p = min(1.0, max(p_values) * 1.1)
        ax.set_xlim(0, max_p)
        
        # Add legend
        normal_patch = mpatches.Patch(color='#2ecc71', label='Normal (p > 0.05)', alpha=0.7)
        not_normal_patch = mpatches.Patch(color='#e74c3c', label='Not Normal (p ≤ 0.05)', alpha=0.7)
        ax.legend(handles=[normal_patch, not_normal_patch, ax.lines[0]], 
                 loc='upper right', fontsize=10, framealpha=0.9)
        
        # Add p-value labels on bars
        for i, (bar, p_val) in enumerate(zip(bars, p_values)):
            width = bar.get_width()
            label_x = width + (max_p * 0.01)  # Slight offset to the right
            ax.text(label_x, bar.get_y() + bar.get_height()/2, 
                   f'{p_val:.4f}', 
                   ha='left', va='center', fontsize=8, fontweight='bold')
        
        # Add summary text
        normal_count = (is_normal == 'Yes').sum()
        total_count = len(is_normal)
        summary_text = f'Summary: {normal_count}/{total_count} samples/groups passed normality test'
        fig.text(0.5, 0.02, summary_text, ha='center', fontsize=11, 
                fontweight='bold', style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout(rect=[0, 0.04, 1, 1])
        
        # Save figure
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        fig_filename = f'normality_test_{polarity}_{timestamp}.png'
        fig_path = os.path.join(viz_dir, fig_filename)
        
        plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        return fig_path

    def export_normalized_results(self):
        # Export normalized data with separate sheets for merged, positive, and negative
        if not hasattr(self, 'normalized_combined_df') or self.normalized_combined_df is None:
            messagebox.showwarning('No Data', 'No normalized data to export.')
            return
        filename = filedialog.asksaveasfilename(title='Save Normalized Data', defaultextension='.xlsx', filetypes=[('Excel','*.xlsx')])
        if not filename:
            return
        
        # Show progress bar
        self.show_stats_progress("Exporting normalized data...")
        processing_report_written = False
        
        try:
            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                # Export normalized matrix after optional post-normalization steps.
                merged_df = self.normalized_combined_df.copy()
                if 'Area (Max.)' in merged_df.columns:
                    merged_df = merged_df.drop(columns=['Area (Max.)'])
                merged_df.to_excel(writer, sheet_name='Merged', index=False)

                # Add optional processing report details when post-normalization processing was enabled.
                report = getattr(self, 'optional_processing_report', {}) or {}
                class_report = getattr(self, 'class_generation_report', {}) or {}
                if bool(getattr(self, 'optional_processing_applied', False)) or bool(class_report.get('applied', False)):
                    report_rows = []

                    # New structure: optional processing is tracked per polarity before merge.
                    if report.get('scope') == 'per_polarity_pre_merge' and isinstance(report.get('polarity_reports'), dict):
                        polarity_reports = report.get('polarity_reports') or {}
                        for pol in ['positive', 'negative']:
                            pol_report = polarity_reports.get(pol) or {}
                            var_rep = pol_report.get('variability_filter') or {}
                            imp_rep = pol_report.get('imputation') or {}
                            pca_rep = pol_report.get('pca_outlier') or {}
                            report_rows.append({'Step': f'{pol.capitalize()} Variability Filter', 'Applied': bool(var_rep.get('applied', False)), 'Details': str(var_rep)})
                            report_rows.append({'Step': f'{pol.capitalize()} Imputation', 'Applied': bool(imp_rep.get('applied', False)), 'Details': str(imp_rep)})
                            report_rows.append({'Step': f'{pol.capitalize()} PCA Outlier Filter', 'Applied': bool(pca_rep.get('applied', False)), 'Details': str(pca_rep)})
                    else:
                        # Backward-compatible structure used by older/custom flow.
                        var_rep = report.get('variability_filter') or {}
                        imp_rep = report.get('imputation') or {}
                        pca_rep = report.get('pca_outlier') or {}
                        report_rows.append({'Step': 'Variability Filter', 'Applied': bool(var_rep.get('applied', False)), 'Details': str(var_rep)})
                        report_rows.append({'Step': 'Imputation', 'Applied': bool(imp_rep.get('applied', False)), 'Details': str(imp_rep)})
                        report_rows.append({'Step': 'PCA Outlier Filter', 'Applied': bool(pca_rep.get('applied', False)), 'Details': str(pca_rep)})

                    class_rep = report.get('class_generation') or class_report
                    report_rows.append({'Step': 'Class Generation', 'Applied': bool(class_rep.get('applied', False)), 'Details': str(class_rep)})
                    pd.DataFrame(report_rows).to_excel(writer, sheet_name='Processing_Report', index=False)
                    processing_report_written = True

                    if hasattr(self, 'last_pca_outlier_scores') and isinstance(self.last_pca_outlier_scores, pd.DataFrame) and not self.last_pca_outlier_scores.empty:
                        self.last_pca_outlier_scores.to_excel(writer, sheet_name='PCA_Outlier_Scores', index=False)
                
                # Save positive data if available (drop Area (Max.) if present)
                if hasattr(self, 'normalized_positive_df') and self.normalized_positive_df is not None:
                    pos_df = self.normalized_positive_df.copy()
                    if 'Area (Max.)' in pos_df.columns:
                        pos_df = pos_df.drop(columns=['Area (Max.)'])
                    pos_df.to_excel(writer, sheet_name='Positive', index=False)
                
                # Save negative data if available (drop Area (Max.) if present)
                if hasattr(self, 'normalized_negative_df') and self.normalized_negative_df is not None:
                    neg_df = self.normalized_negative_df.copy()
                    if 'Area (Max.)' in neg_df.columns:
                        neg_df = neg_df.drop(columns=['Area (Max.)'])
                    neg_df.to_excel(writer, sheet_name='Negative', index=False)
                
                # Save class data if available (drop Area (Max.) if present)
                if hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                    class_df = self.normalized_combined_class_df.copy()
                    drop_cols = [col for col in ['Area (Max.)', 'Max_Area'] if col in class_df.columns]
                    if drop_cols:
                        class_df = class_df.drop(columns=drop_cols)
                    class_df.to_excel(writer, sheet_name='Merged_Class', index=False)
                
                if hasattr(self, 'normalized_positive_class_df') and self.normalized_positive_class_df is not None:
                    pos_class_df = self.normalized_positive_class_df.copy()
                    drop_cols = [col for col in ['Area (Max.)', 'Max_Area'] if col in pos_class_df.columns]
                    if drop_cols:
                        pos_class_df = pos_class_df.drop(columns=drop_cols)
                    pos_class_df.to_excel(writer, sheet_name='Positive_Class', index=False)
                
                if hasattr(self, 'normalized_negative_class_df') and self.normalized_negative_class_df is not None:
                    neg_class_df = self.normalized_negative_class_df.copy()
                    drop_cols = [col for col in ['Area (Max.)', 'Max_Area'] if col in neg_class_df.columns]
                    if drop_cols:
                        neg_class_df = neg_class_df.drop(columns=drop_cols)
                    neg_class_df.to_excel(writer, sheet_name='Negative_Class', index=False)
                
                # Add normality test results if available (combined sheet)
                if hasattr(self, 'normality_test_results') and self.normality_test_results:
                    all_normality = []
                    for polarity, results_df in self.normality_test_results.items():
                        if results_df is not None and hasattr(results_df, 'empty') and not results_df.empty:
                            out_df = results_df.copy()
                            out_df.insert(0, 'Polarity', polarity.capitalize())
                            all_normality.append(out_df)
                    if all_normality:
                        pd.concat(all_normality, ignore_index=True).to_excel(writer, sheet_name='Normality_Test', index=False)
                
                # Generate QQ plots in same folder as export file (max 8: 2 normal + 2 not-normal per polarity)
                qq_plot_paths = []
                if hasattr(self, 'normality_test_targets') and self.normality_test_targets:
                    try:
                        from scipy.stats import norm as _norm_dist
                        import numpy as _np
                        import matplotlib.pyplot as _plt
                        from pathlib import Path as _Path
                        
                        # Create qq_plots subfolder in same directory as export file
                        export_dir = _Path(filename).parent
                        qq_dir = export_dir / 'qq_plots'
                        qq_dir.mkdir(parents=True, exist_ok=True)
                        
                        for polarity, target_data in self.normality_test_targets.items():
                            plot_data = target_data.get('plot_data', [])
                            normal_count = target_data.get('normal_count', 0)
                            not_normal_count = target_data.get('not_normal_count', 0)
                            
                            # Generate QQ plots for selected metabolites (already selected in normality test)
                            for kind, metabolite_id, values in plot_data:
                                try:
                                    all_vals = _np.array(values, dtype=float)
                                    n = len(all_vals)
                                    
                                    if n >= 3:
                                        mu, sigma = all_vals.mean(), all_vals.std(ddof=1) if n > 1 else 0.0
                                        fig, ax = _plt.subplots(figsize=(4.2, 4.2))
                                        sorted_vals = _np.sort(all_vals)
                                        quantiles = _norm_dist.ppf((_np.arange(1, n+1) - 0.5) / n) * (sigma if sigma else 1) + mu
                                        ax.scatter(quantiles, sorted_vals, s=10, alpha=0.7, edgecolor='none')
                                        min_v = min(quantiles.min(), sorted_vals.min())
                                        max_v = max(quantiles.max(), sorted_vals.max())
                                        ax.plot([min_v, max_v], [min_v, max_v], 'r--', linewidth=1)
                                        
                                        # Add status to title
                                        status = 'Normal' if kind == 'NORMAL' else 'Not Normal'
                                        ax.set_title(f'QQ Plot: {metabolite_id} ({polarity}) - {status}', fontsize=9)
                                        ax.set_xlabel('Theoretical Quantiles')
                                        ax.set_ylabel('Observed')
                                        ax.grid(alpha=0.3, linewidth=0.4)
                                        
                                        fname = qq_dir / f'qq_{polarity}_{kind}_{metabolite_id.replace(" ", "_").replace("/", "_")}.png'
                                        fig.tight_layout()
                                        fig.savefig(fname, dpi=150)
                                        _plt.close(fig)
                                        
                                        qq_plot_paths.append({
                                            'Polarity': polarity,
                                            'Status': status,
                                            'Metabolite': metabolite_id,
                                            'N': n,
                                            'Path': str(fname)
                                        })
                                except Exception:
                                    pass
                        
                        # Add summary of normality counts to QQ plots sheet
                        if qq_plot_paths:
                                qq_df = pd.DataFrame(qq_plot_paths)
                                qq_df.to_excel(writer, sheet_name='QQ_Plots', index=False)

                                # Add summary sheet with percentages
                                summary_rows = []
                                for polarity, target_data in self.normality_test_targets.items():
                                    norm_ct = target_data.get('normal_count', 0)
                                    not_ct = target_data.get('not_normal_count', 0)
                                    total_ct = norm_ct + not_ct
                                    pct_norm = (norm_ct / total_ct * 100.0) if total_ct else 0.0
                                    pct_not = (not_ct / total_ct * 100.0) if total_ct else 0.0
                                    summary_rows.append({
                                        'Polarity': polarity.capitalize(),
                                        'Normal_Count': norm_ct,
                                        'Normal_%': round(pct_norm, 2),
                                        'Not_Normal_Count': not_ct,
                                        'Not_Normal_%': round(pct_not, 2),
                                        'Total_Tested': total_ct
                                    })
                                if summary_rows:
                                    pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Normality_Summary', index=False)
                    except Exception:
                        pass
            
            # Remove read-only attribute
            self._remove_readonly_attribute(filename)
            
            # Track export path for Open Folder button
            self.last_statistics_export_path = filename
            
            # Create summary message
            sheets_saved = ['Merged']
            if processing_report_written:
                sheets_saved.append('Processing_Report')
                if hasattr(self, 'last_pca_outlier_scores') and isinstance(self.last_pca_outlier_scores, pd.DataFrame) and not self.last_pca_outlier_scores.empty:
                    sheets_saved.append('PCA_Outlier_Scores')
            if hasattr(self, 'normality_test_results') and self.normality_test_results:
                sheets_saved.append('Normality_Test')
            if hasattr(self, 'normalized_positive_df') and self.normalized_positive_df is not None:
                sheets_saved.append('Positive')
            if hasattr(self, 'normalized_negative_df') and self.normalized_negative_df is not None:
                sheets_saved.append('Negative')
            if hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                sheets_saved.append('Merged_Class')
            if hasattr(self, 'normalized_positive_class_df') and self.normalized_positive_class_df is not None:
                sheets_saved.append('Positive_Class')
            if hasattr(self, 'normalized_negative_class_df') and self.normalized_negative_class_df is not None:
                sheets_saved.append('Negative_Class')
            
            self.show_stats_progress("Normalized data saved successfully.")
            # Hide progress bar
            self.hide_stats_progress()
                
            messagebox.showinfo('Exported', f'Normalized data saved to {filename}\n\nSheets saved: {", ".join(sheets_saved)}')
        except Exception as e:
            self.hide_stats_progress()
            messagebox.showerror('Export Error', str(e))

    def open_statistics_output_folder(self):
        """Open the folder containing statistical output files"""
        import subprocess
        import platform
        
        # Try to get the last export path
        if hasattr(self, 'last_statistics_export_path') and self.last_statistics_export_path:
            folder_path = os.path.dirname(self.last_statistics_export_path)
        else:
            # No export has been made yet - inform user
            messagebox.showinfo(
                "No Export Yet", 
                "No statistical data has been exported yet.\n\n"
                "Please export your data first using:\n"
                "• 'Export Normalized Data' button, or\n"
                "• 'Export Statistical Results' button, or\n"
                "• 'Export Covariate Results' button (if using covariate adjustment)\n\n"
                "Then this button will open the folder where your data was saved."
            )
            return
        
        if not os.path.exists(folder_path):
            messagebox.showerror("Folder Not Found", f"Output folder does not exist:\n{folder_path}")
            return
        
        try:
            if platform.system() == 'Windows':
                os.startfile(folder_path)
            elif platform.system() == 'Darwin':  # macOS
                subprocess.Popen(['open', folder_path])
            else:  # Linux
                subprocess.Popen(['xdg-open', folder_path])
        except Exception as e:
            messagebox.showerror("Error", f"Could not open folder:\n{str(e)}")

    def export_statistical_results(self):
        """Threaded export of statistical test results to Excel with progress updates to avoid GUI freeze."""
        if not hasattr(self, 'statistical_test_results') or not self.statistical_test_results:
            messagebox.showwarning('No Results', 'No statistical test results available. Run "Statistical Tests" first.')
            return

        filename = filedialog.asksaveasfilename(
            title='Save Statistical Results',
            defaultextension='.xlsx',
            filetypes=[('Excel','*.xlsx')]
        )
        if not filename:
            return

        # Disable export button to prevent double clicks
        try:
            for widget in self.stats_tab.winfo_children():
                if isinstance(widget, tk.Frame):
                    for child in widget.winfo_children():
                        if isinstance(child, tk.Button) and 'Export Statistical Results' in child.cget('text'):
                            child.config(state='disabled')
        except Exception:
            pass

        # Initialize export tracking
        self._export_threads = []
        
        # Launch export in thread (non-daemon so we can track completion)
        main_thread = threading.Thread(target=self._export_statistical_results_threaded, args=(filename,), daemon=False)
        main_thread.start()
        self._export_threads.append(main_thread)
        
        # If in lipid mode and class results exist, also export class results
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        if mode == 'lipid' and hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class:
            # Generate class filename with "_class" suffix
            base_name, ext = os.path.splitext(filename)
            class_filename = f"{base_name}_class{ext}"
            # Launch class export in thread that waits for main export first
            self._thread_safe_log(f"📊 Also exporting lipid class results to: {os.path.basename(class_filename)}\n")
            class_thread = threading.Thread(target=self._export_class_statistical_results_threaded_with_wait, args=(class_filename, main_thread), daemon=False)
            class_thread.start()
            self._export_threads.append(class_thread)
            
            # Start a monitor thread to show completion popup after both exports finish
            monitor_thread = threading.Thread(target=self._monitor_export_completion, args=(self._export_threads, filename), daemon=True)
            monitor_thread.start()
        else:
            # If no class export, monitor just completes when main export finishes
            monitor_thread = threading.Thread(target=self._monitor_export_completion, args=(self._export_threads, filename), daemon=True)
            monitor_thread.start()

    def _export_statistical_results_threaded(self, filename):
        try:
            # Debug: Log what statistical results are available
            self._thread_safe_log(f"📊 Export Debug - Available results keys: {list(self.statistical_test_results.keys())}\n")
            # Determine data mode for export behavior (metabolite vs lipid)
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            
            # Determine number of steps: Two-way ANOVA (optional) + Complete Results (optional) + each pairwise comparison + original normalized + finalize
            total_steps = 3  # baseline: setup, write original normalized, finalize
            if 'two_way_anova' in self.statistical_test_results and not self.statistical_test_results['two_way_anova'].empty:
                total_steps += 1
            # Count posthoc sheets if present
            if 'two_way_posthoc_A' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_A'].empty:
                total_steps += 1
            if 'two_way_posthoc_B' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_B'].empty:
                total_steps += 1
            if 'two_way_posthoc_pairwise' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_pairwise'].empty:
                total_steps += 1  # Combined pairwise sheet
                # Count individual pairwise sheets
                phpw = self.statistical_test_results['two_way_posthoc_pairwise']
                if 'group1' in phpw.columns and 'group2' in phpw.columns:
                    num_individual_sheets = phpw.groupby(['group1','group2']).ngroups
                    total_steps += num_individual_sheets
            if 'two_way_posthoc_pairwise' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_pairwise'].empty:
                total_steps += 1
            
            # Add step for Complete Results if available
            if 'enhanced_metabolites' in self.statistical_test_results and not self.statistical_test_results['enhanced_metabolites'].empty:
                total_steps += 1  # complete results sheet
                self._thread_safe_log(f"📊 Enhanced metabolites available - will write Complete Results sheet\n")
            else:
                self._thread_safe_log(f"📊 No enhanced metabolites data available\n")
            
            # Add steps for pairwise comparisons if available
            pairwise = None
            if 'pairwise' in self.statistical_test_results and not self.statistical_test_results['pairwise'].empty:
                pairwise = self.statistical_test_results['pairwise'].copy()
                # Debug: Log pairwise columns
                self._thread_safe_log(f"🔍 DEBUG: Pairwise columns: {list(pairwise.columns)}\n")
                
                # Check if group1/group2 columns exist
                if 'group1' in pairwise.columns and 'group2' in pairwise.columns:
                    comparisons = pairwise.groupby(['group1','group2']).ngroups
                    total_steps += comparisons  # each comparison sheet
                    self._thread_safe_log(f"📊 Pairwise data available - will write {comparisons} comparison sheets\n")
                else:
                    self._thread_safe_log(f"⚠️ Pairwise data exists but lacks group1/group2 columns - skipping individual comparison sheets\n")
            else:
                self._thread_safe_log(f"📊 No pairwise comparison data available\n")
            # Show progress bar determinate
            self._stats_total_steps = total_steps
            self._stats_current_step = 0
            self.show_stats_progress("Starting export...")
            step = 1
            def prog(msg):
                self._thread_safe_progress_step(step_holder[0], total_steps, msg)
            step_holder = [step]

            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                sheets_written = []
                new_format_class_export = False
                
                # NEW FORMAT: Check if we have the new complete_result format
                if 'complete_result' in self.statistical_test_results and not self.statistical_test_results['complete_result'].empty:
                    # NEW FORMAT: Write Complete Results sheet
                    prog("Writing Complete Results sheet...")
                    complete_df = self.statistical_test_results['complete_result'].copy()
                    # Drop fully NaN pairwise stat columns to reduce clutter
                    try:
                        drop_cols = []
                        for c in complete_df.columns:
                            if ('_vs_' in c and any(suf in c for suf in ['_log2FC','_adj_p','_p_adj','_neg_log10'])):
                                if complete_df[c].isna().all():
                                    drop_cols.append(c)
                        if drop_cols:
                            complete_df = complete_df.drop(columns=drop_cols)
                            self._thread_safe_log(f"🧹 Dropped {len(drop_cols)} all-NaN pairwise columns before export\n")
                    except Exception:
                        pass
                    complete_df.to_excel(writer, sheet_name='Complete Results', index=False)
                    sheets_written.append('Complete Results')
                    step_holder[0] += 1
                    
                    # NEW FORMAT: Write individual pairwise sheets
                    for sheet_name, sheet_df in self.statistical_test_results.items():
                        # Skip duplicate enhanced_metabolites (same as complete_result)
                        if sheet_name in ('complete_result','enhanced_metabolites'):
                            continue
                        if isinstance(sheet_df, pd.DataFrame) and not sheet_df.empty:
                            # Make a copy to avoid modifying the cached results
                            sheet_df_to_write = sheet_df.copy()
                            # For each sheet, remove all-NaN pairwise stat columns
                            try:
                                drop_cols2 = []
                                for c in sheet_df_to_write.columns:
                                    if ('_vs_' in c and any(suf in c for suf in ['_log2FC','_adj_p','_p_adj','_neg_log10'])):
                                        if sheet_df_to_write[c].isna().all():
                                            drop_cols2.append(c)
                                if drop_cols2:
                                    sheet_df_to_write = sheet_df_to_write.drop(columns=drop_cols2)
                                    self._thread_safe_log(f"🧹 Dropped {len(drop_cols2)} all-NaN columns from {sheet_name} sheet\n")
                            except Exception:
                                pass
                            prog(f"Writing {sheet_name} sheet...")
                            sheet_df_to_write.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1
                    
                # OLD FORMAT: Fallback to old format if new format not available
                elif 'two_way_anova' in self.statistical_test_results and not self.statistical_test_results['two_way_anova'].empty:
                    prog("Writing Two-Way ANOVA sheet (old format)...")
                    twa_df = self.statistical_test_results['two_way_anova'].copy()

                    # Augment TwoWayANOVA sheet with Tukey pairwise adj_p columns (A and B) for compatibility
                    try:
                        import numpy as np
                        fa_name = getattr(self, 'factorA_name', 'Factor A')
                        fb_name = getattr(self, 'factorB_name', 'Factor B')
                        # Build wide pairwise columns for A
                        if 'two_way_posthoc_A' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_A'].empty:
                            pha = self.statistical_test_results['two_way_posthoc_A']
                            if 'group1' in pha.columns and 'group2' in pha.columns:
                                wide_parts = []
                                for (g1, g2), sub in pha.groupby(['group1','group2']):
                                    comp = f"{fa_name}:{g1}_vs_{g2}"
                                    part = sub[['metabolite','p_value_adj']].copy()
                                    part = part.rename(columns={'p_value_adj': f'{comp}_adj_p'})
                                    part[f'{comp}_neg_log10_adj_p'] = -np.log10(part[f'{comp}_adj_p'].replace(0, np.finfo(float).eps))
                                    wide_parts.append(part)
                                if wide_parts:
                                    from functools import reduce
                                    wideA = reduce(lambda a,b: a.merge(b, on='metabolite', how='outer'), wide_parts)
                                    if 'metabolite' in twa_df.columns:
                                        twa_df = twa_df.merge(wideA, on='metabolite', how='left')
                        # Build wide pairwise columns for B
                        if 'two_way_posthoc_B' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_B'].empty:
                            phb = self.statistical_test_results['two_way_posthoc_B']
                            if 'group1' in phb.columns and 'group2' in phb.columns:
                                wide_parts = []
                                for (g1, g2), sub in phb.groupby(['group1','group2']):
                                    comp = f"{fb_name}:{g1}_vs_{g2}"
                                    part = sub[['metabolite','p_value_adj']].copy()
                                    part = part.rename(columns={'p_value_adj': f'{comp}_adj_p'})
                                    part[f'{comp}_neg_log10_adj_p'] = -np.log10(part[f'{comp}_adj_p'].replace(0, np.finfo(float).eps))
                                    wide_parts.append(part)
                                if wide_parts:
                                    from functools import reduce
                                    wideB = reduce(lambda a,b: a.merge(b, on='metabolite', how='outer'), wide_parts)
                                    if 'metabolite' in twa_df.columns:
                                        twa_df = twa_df.merge(wideB, on='metabolite', how='left')
                    except Exception:
                        pass

                    twa_df.to_excel(writer, sheet_name='TwoWayANOVA', index=False)
                    sheets_written.append('TwoWayANOVA')
                    step_holder[0] += 1

                    # Two-Way ANOVA Posthoc sheets (if present)
                    if 'two_way_posthoc_A' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_A'].empty:
                        prog("Writing Two-Way Posthoc A sheet...")
                        pha = self.statistical_test_results['two_way_posthoc_A'].copy()
                        # Normalize names for clarity
                        import numpy as np
                        if 'p_value_adj' in pha.columns and 'adj_p' not in pha.columns:
                            pha['adj_p'] = pha['p_value_adj']
                        if 'neg_log10_adj_p' not in pha.columns and 'adj_p' in pha.columns:
                            pha['neg_log10_adj_p'] = -np.log10(pha['adj_p'].replace(0, np.finfo(float).eps))
                        pha.to_excel(writer, sheet_name='TwoWay_Posthoc_A', index=False)
                        sheets_written.append('TwoWay_Posthoc_A')
                        step_holder[0] += 1

                    if 'two_way_posthoc_B' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_B'].empty:
                        prog("Writing Two-Way Posthoc B sheet...")
                        phb = self.statistical_test_results['two_way_posthoc_B'].copy()
                        import numpy as np
                        if 'p_value_adj' in phb.columns and 'adj_p' not in phb.columns:
                            phb['adj_p'] = phb['p_value_adj']
                        if 'neg_log10_adj_p' not in phb.columns and 'adj_p' in phb.columns:
                            phb['neg_log10_adj_p'] = -np.log10(phb['adj_p'].replace(0, np.finfo(float).eps))
                        phb.to_excel(writer, sheet_name='TwoWay_Posthoc_B', index=False)
                        sheets_written.append('TwoWay_Posthoc_B')
                        step_holder[0] += 1
                    
                    # Two-Way ANOVA All Pairwise Comparisons (if present) - write combined sheet AND individual sheets
                    if 'two_way_posthoc_pairwise' in self.statistical_test_results and not self.statistical_test_results['two_way_posthoc_pairwise'].empty:
                        prog("Writing Two-Way Pairwise Comparisons sheet...")
                        phpw = self.statistical_test_results['two_way_posthoc_pairwise'].copy()
                        import numpy as np
                        if 'p_value_adj' in phpw.columns and 'adj_p' not in phpw.columns:
                            phpw['adj_p'] = phpw['p_value_adj']
                        if 'neg_log10_adj_p' not in phpw.columns and 'adj_p' in phpw.columns:
                            phpw['neg_log10_adj_p'] = -np.log10(phpw['adj_p'].replace(0, np.finfo(float).eps))
                        phpw.to_excel(writer, sheet_name='TwoWay_Pairwise', index=False)
                        sheets_written.append('TwoWay_Pairwise')
                        step_holder[0] += 1
                        
                        # Also create individual comparison sheets for Two-Way pairwise (like one-way ANOVA)
                    if 'group1' in phpw.columns and 'group2' in phpw.columns:
                        # Merge with enhanced_metabolites for feature columns
                        tw_pairwise = phpw.copy()
                        if 'enhanced_metabolites' in self.statistical_test_results and not self.statistical_test_results['enhanced_metabolites'].empty:
                            enhanced_meta = self.statistical_test_results['enhanced_metabolites'].copy()
                            exclude_patterns = [
                                'mean_', 'n_', 'p_value', 'statistic', 'fold_change', 'log2_fold_change',
                                'Expression', 'neg_log10', '_vs_', 'overall_', 'rots_', 'dunn_', 'tukey_',
                                'order_idx'
                            ]
                            feature_cols = []
                            for col in enhanced_meta.columns:
                                if col in ('metabolite_id', 'LipidID'):
                                    feature_cols.append(col)
                                    continue
                                if any(pat in col for pat in exclude_patterns):
                                    continue
                                feature_cols.append(col)
                            
                            if feature_cols:
                                feature_df = enhanced_meta[feature_cols].copy()
                                # Determine join preference
                                join_key = None
                                id_col = self._get_verified_id_column(feature_df)
                                if mode == 'lipid' and id_col in feature_df.columns and id_col in tw_pairwise.columns:
                                    join_key = id_col
                                elif 'metabolite_id' in feature_df.columns and 'metabolite_id' in tw_pairwise.columns:
                                    join_key = 'metabolite_id'
                                if join_key is None and 'metabolite' in tw_pairwise.columns:
                                    tw_pairwise['_merge_idx'] = tw_pairwise.index
                                    feature_df['_merge_idx'] = feature_df.index
                                    tw_pairwise = tw_pairwise.merge(feature_df, on='_merge_idx', how='left', suffixes=('', '_dup'))
                                    tw_pairwise.drop(columns=['_merge_idx'], inplace=True, errors='ignore')
                                    if join_key and join_key in tw_pairwise.columns and 'metabolite' in tw_pairwise.columns:
                                        tw_pairwise.drop(columns=['metabolite'], inplace=True, errors='ignore')
                                elif join_key is not None:
                                    tw_pairwise = tw_pairwise.merge(feature_df, left_on=join_key, right_on=join_key, how='left', suffixes=('', '_dup'))
                                else:
                                    tw_pairwise['_merge_idx'] = tw_pairwise.index
                                    feature_df['_merge_idx'] = feature_df.index
                                    tw_pairwise = tw_pairwise.merge(feature_df, on='_merge_idx', how='left', suffixes=('', '_dup'))
                                    tw_pairwise.drop(columns=['_merge_idx'], inplace=True, errors='ignore')
                                dup_cols = [c for c in tw_pairwise.columns if c.endswith('_dup')]
                                if dup_cols:
                                    tw_pairwise.drop(columns=dup_cols, inplace=True, errors='ignore')
                        
                        # Get sample group map
                        sample_group_map = {}
                        try:
                            if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                                for col, var in self.sample_group_vars.items():
                                    sample_group_map[col] = var.get()
                        except Exception:
                            pass
                        sample_cols_all = set(sample_group_map.keys())
                        
                        # Create individual sheets per comparison (only if group columns exist)
                        if 'group1' in tw_pairwise.columns and 'group2' in tw_pairwise.columns:
                            for (g1, g2), sub_idx in tw_pairwise.groupby(['group1','group2']).groups.items():
                                prog(f"Writing TwoWay {g1}_vs_{g2} sheet...")
                                sub = tw_pairwise.loc[sub_idx].copy()
                            if 'fold_change' in sub.columns:
                                sub.rename(columns={'fold_change':'FC', 'log2_fold_change':'log2FC'}, inplace=True)
                            sub.rename(columns={'mean_group1': f'mean_{g1}', 'mean_group2': f'mean_{g2}'}, inplace=True)
                            if 'adj_p' in sub.columns and 'neg_log10_adj_p' not in sub.columns:
                                sub['neg_log10_adj_p'] = -np.log10(sub['adj_p'].replace(0, np.finfo(float).eps))
                            
                            # Build column order
                            feature_cols = []
                            known_feature_cols = ['metabolite_id', 'Name', 'Compound', 'Metabolite', 'HMDB', 'PubChem', 
                                                'Formula', 'MW', 'KEGG', 'LipidMaps', 'InChIKey', 'SMILES',
                                                'LipidID', 'Lipid_Class', 'Class', 'Class_name']
                            for fc in known_feature_cols:
                                if fc in sub.columns and fc not in feature_cols:
                                    feature_cols.append(fc)
                            # Ensure Class_name placed after the class column
                            if 'Class_name' in feature_cols:
                                try:
                                    class_pos = None
                                    for i, col in enumerate(feature_cols):
                                        if col in ('Lipid_Class', 'Class'):
                                            class_pos = i
                                            break
                                    if class_pos is not None:
                                        feature_cols.remove('Class_name')
                                        feature_cols.insert(class_pos + 1, 'Class_name')
                                except Exception:
                                    pass
                            # In lipid mode, remove metabolite_id from output
                            if mode == 'lipid' and 'metabolite_id' in sub.columns:
                                sub.drop(columns=['metabolite_id'], inplace=True, errors='ignore')
                                if 'metabolite_id' in feature_cols:
                                    feature_cols = [c for c in feature_cols if c != 'metabolite_id']
                            
                            stat_cols = [f'mean_{g1}', f'mean_{g2}', 'FC', 'log2FC',
                                       'p_value', 'adj_p', 'neg_log10_adj_p',
                                       'model_effect', 'model_se', 'ci_lower_95', 'ci_upper_95',
                                       'statistic', 'cohen_d', 'cliffs_delta', 'Expression', 
                                       'n_group1', 'n_group2']
                            stat_present = [c for c in stat_cols if c in sub.columns]
                            
                            # Keep only sample columns for g1 and g2
                            raw_sample_cols = [c for c in sub.columns if c in sample_cols_all]
                            keep_sample_cols = [c for c in raw_sample_cols if sample_group_map.get(c) in (g1, g2)]
                            drop_samples = [c for c in raw_sample_cols if c not in keep_sample_cols]
                            if drop_samples:
                                sub.drop(columns=drop_samples, inplace=True, errors='ignore')
                            
                            remainder = [c for c in sub.columns if c not in feature_cols + stat_present + keep_sample_cols
                                       and c not in ['group1', 'group2', 'metabolite', 'factor']]
                            final_cols = feature_cols + keep_sample_cols + stat_present + remainder
                            sub = sub.loc[:, [c for c in dict.fromkeys(final_cols).keys() if c in sub.columns]]
                            
                            sheet_name = f'TwoWay_{g1}_vs_{g2}'[:31]
                            sub.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1

                # Complete Results sheet (independent of pairwise data)
                if 'enhanced_metabolites' in self.statistical_test_results and not self.statistical_test_results['enhanced_metabolites'].empty:
                    prog("Writing Complete Results sheet...")
                    complete_df = self.statistical_test_results['enhanced_metabolites'].copy()
                    complete_df = self._clean_complete_results_dataframe(complete_df)
                    # Lipid mode: remove metabolite_id from Complete Results output
                    if mode == 'lipid' and 'metabolite_id' in complete_df.columns:
                        complete_df.drop(columns=['metabolite_id'], inplace=True, errors='ignore')
                    complete_df.to_excel(writer, sheet_name='Complete Results', index=False)
                    sheets_written.append('Complete Results')
                    step_holder[0] += 1
                
                # Pairwise comparison sheets
                if pairwise is not None:
                    # Merge pairwise with enhanced_metabolites to get all feature columns
                    if 'enhanced_metabolites' in self.statistical_test_results and not self.statistical_test_results['enhanced_metabolites'].empty:
                        enhanced_meta = self.statistical_test_results['enhanced_metabolites'].copy()
                        # Build a comprehensive set of feature/ID columns by excluding known statistical/helper columns.
                        # This mirrors the class export logic so pairwise sheets carry full metadata (Formula, MW, IDs, etc.).
                        exclude_patterns = [
                            'mean_', 'n_', 'p_value', 'statistic', 'fold_change', 'log2_fold_change',
                            'Expression', 'neg_log10', '_vs_', 'overall_', 'rots_', 'dunn_', 'tukey_',
                            'order_idx'
                        ]
                        feature_cols_to_merge = []
                        for col in enhanced_meta.columns:
                            # Always keep identifier columns; other columns only if not statistical/helper
                            if col in ('metabolite_id', 'LipidID'):
                                feature_cols_to_merge.append(col)
                                continue
                            if any(pat in col for pat in exclude_patterns):
                                continue
                            # Avoid bringing raw sample columns later (we'll drop any that slip through below)
                            feature_cols_to_merge.append(col)

                        if feature_cols_to_merge:
                            feature_df = enhanced_meta[feature_cols_to_merge].copy()

                            # Choose best join key based on mode and available columns
                            join_key = None
                            id_col = self._get_verified_id_column(feature_df)
                            if mode == 'lipid' and id_col in feature_df.columns and id_col in pairwise.columns:
                                join_key = id_col
                            elif 'metabolite_id' in feature_df.columns and 'metabolite_id' in pairwise.columns:
                                join_key = 'metabolite_id'

                            if join_key is None and 'metabolite' in pairwise.columns:
                                # Fallback: merge by index when no explicit key aligns
                                pairwise['_merge_idx'] = pairwise.index
                                feature_df['_merge_idx'] = feature_df.index
                                pairwise = pairwise.merge(
                                    feature_df,
                                    on='_merge_idx',
                                    how='left',
                                    suffixes=('', '_dup')
                                )
                                pairwise.drop(columns=['_merge_idx'], inplace=True, errors='ignore')
                                # Remove placeholder metabolite column if we now have a better ID
                                if join_key and join_key in pairwise.columns and 'metabolite' in pairwise.columns:
                                    pairwise.drop(columns=['metabolite'], inplace=True, errors='ignore')
                            elif join_key is not None:
                                pairwise = pairwise.merge(
                                    feature_df,
                                    left_on=join_key,
                                    right_on=join_key,
                                    how='left',
                                    suffixes=('', '_dup')
                                )
                            else:
                                # Last resort: align on index to carry feature metadata
                                pairwise['_merge_idx'] = pairwise.index
                                feature_df['_merge_idx'] = feature_df.index
                                pairwise = pairwise.merge(
                                    feature_df,
                                    on='_merge_idx',
                                    how='left',
                                    suffixes=('', '_dup')
                                )
                                pairwise.drop(columns=['_merge_idx'], inplace=True, errors='ignore')

                            # Clean up duplicates from merge
                            dup_cols = [c for c in pairwise.columns if c.endswith('_dup')]
                            if dup_cols:
                                pairwise.drop(columns=dup_cols, inplace=True, errors='ignore')

                            # In lipid mode, remove metabolite_id from the merged pairwise table
                            if mode == 'lipid' and 'metabolite_id' in pairwise.columns:
                                pairwise.drop(columns=['metabolite_id'], inplace=True, errors='ignore')

                            self._thread_safe_log(f"✅ Merged {len(feature_cols_to_merge)} feature columns into pairwise data\n")

                            # Remove raw sample columns that may have been present in the merged DataFrame
                            # Treat any column present in the original normalized table but not explicitly chosen features
                            # and not a derived mean_/stat column as a raw sample column to drop.
                            try:
                                sample_cols_to_remove = []
                                if hasattr(self, 'normalized_combined_df') and self.normalized_combined_df is not None:
                                    norm_cols = set(self.normalized_combined_df.columns)
                                    sample_cols_to_remove = [c for c in pairwise.columns
                                                             if c in norm_cols
                                                             and c not in feature_cols_to_merge
                                                             and not c.startswith('mean_')
                                                             and c not in ('metabolite','metabolite_id','group1','group2')]
                                if sample_cols_to_remove:
                                    pairwise.drop(columns=sample_cols_to_remove, inplace=True, errors='ignore')
                                    self._thread_safe_log(f"ℹ️ Dropped {len(sample_cols_to_remove)} raw sample columns from pairwise export: {sample_cols_to_remove[:10]}\n")
                            except Exception as _:
                                pass
                    
                    # Normalize naming
                    if 'dunn_bh_pvalue' in pairwise.columns:
                        pairwise.rename(columns={'dunn_bh_pvalue':'adj_p', 'neg_log10_dunn_bh_pvalue':'neg_log10_adj_p'}, inplace=True)
                    elif 'p_value_adj' in pairwise.columns and 'adj_p' not in pairwise.columns:
                        pairwise.rename(columns={'p_value_adj':'adj_p'}, inplace=True)
                    
                    # Precompute sample → group label map for selective retention of sample intensity columns
                    sample_group_map = {}
                    try:
                        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                            for col, var in self.sample_group_vars.items():
                                sample_group_map[col] = var.get()
                    except Exception:
                        pass
                    sample_cols_all = set(sample_group_map.keys())

                    # Check if we have individual_sheets (from one-way ANOVA clean implementation)
                    if 'individual_sheets' in self.statistical_test_results and self.statistical_test_results['individual_sheets']:
                        # Use pre-split individual sheets
                        for comp_name, comp_df in self.statistical_test_results['individual_sheets'].items():
                            prog(f"Writing {comp_name} sheet...")
                            sub = comp_df.copy()
                            
                            # Normalize column names
                            if 'fold_change' in sub.columns:
                                sub.rename(columns={'fold_change':'FC', 'log2_fold_change':'log2FC'}, inplace=True)
                            
                            # Write sheet
                            sheet_name = comp_name[:31]  # Excel sheet name limit
                            sub.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1
                    
                    # Otherwise, try to split pairwise by group1/group2 columns
                    elif 'group1' in pairwise.columns and 'group2' in pairwise.columns:
                        for (g1, g2), sub_idx in pairwise.groupby(['group1','group2']).groups.items():
                            prog(f"Writing {g1}_vs_{g2} sheet...")
                            sub = pairwise.loc[sub_idx].copy()
                            if 'fold_change' in sub.columns:
                                sub.rename(columns={'fold_change':'FC', 'log2_fold_change':'log2FC'}, inplace=True)
                            sub.rename(columns={'mean_group1': f'mean_{g1}', 'mean_group2': f'mean_{g2}'}, inplace=True)
                            # Compute appropriate -log10 column based on availability and user preference
                            import numpy as np
                            if 'adj_p' in sub.columns:
                                if 'neg_log10_adj_p' not in sub.columns:
                                    sub['neg_log10_adj_p'] = -np.log10(sub['adj_p'].replace(0, np.finfo(float).eps))
                            elif 'p_value' in sub.columns:
                                # No adjusted p; ensure raw neg_log10
                                sub['neg_log10_p'] = -np.log10(sub['p_value'].replace(0, np.finfo(float).eps))
                            
                            # Core columns: Start with metabolite_id (not numeric 'metabolite'), then feature columns, then stats
                            # Identify feature columns (Name, HMDB, Formula, etc.)
                            feature_cols = []
                            known_feature_cols = ['metabolite_id', 'Name', 'Compound', 'Metabolite', 'HMDB', 'PubChem', 
                                                'Formula', 'MW', 'KEGG', 'LipidMaps', 'InChIKey', 'SMILES',
                                                'LipidID', 'Lipid_Class', 'Class', 'Class_name']
                            for fc in known_feature_cols:
                                if fc in sub.columns and fc not in feature_cols:
                                    feature_cols.append(fc)
                            # Ensure Class_name ordered after class column
                            if 'Class_name' in feature_cols:
                                try:
                                    class_pos = None
                                    for i, col in enumerate(feature_cols):
                                        if col in ('Lipid_Class', 'Class'):
                                            class_pos = i
                                            break
                                    if class_pos is not None:
                                        feature_cols.remove('Class_name')
                                        feature_cols.insert(class_pos + 1, 'Class_name')
                                except Exception:
                                    pass
                            # In lipid mode, exclude metabolite_id from outputs
                            if mode == 'lipid' and 'metabolite_id' in sub.columns:
                                sub.drop(columns=['metabolite_id'], inplace=True, errors='ignore')
                                if 'metabolite_id' in feature_cols:
                                    feature_cols = [c for c in feature_cols if c != 'metabolite_id']
                            
                            # Statistical columns: include adjusted or raw p-values accordingly
                            stat_cols = [f'mean_{g1}', f'mean_{g2}', 'FC', 'log2FC',
                                       'p_value', 'adj_p', 'p_value_adj', 
                                       'neg_log10_p', 'neg_log10_adj_p', 'neg_log10_p_adj',
                                       'statistic',
                                       'cohen_d', 'cliffs_delta', 'Expression', 
                                       'n_group1', 'n_group2']
                            stat_present = [c for c in stat_cols if c in sub.columns]
                            
                            # Remaining columns except helpers
                            remainder = [c for c in sub.columns 
                                       if c not in feature_cols + stat_present 
                                       and c not in ['group1', 'group2', 'metabolite']]
                            
                            # Separate raw sample intensity columns (retain only those belonging to g1 or g2)
                            raw_sample_cols_present = [c for c in sub.columns if c in sample_cols_all]
                            keep_sample_cols = []
                            for rc in raw_sample_cols_present:
                                grp = sample_group_map.get(rc)
                                if grp in (g1, g2):
                                    keep_sample_cols.append(rc)
                            # Remove any other sample columns
                            drop_other_samples = [c for c in raw_sample_cols_present if c not in keep_sample_cols]
                            if drop_other_samples:
                                sub.drop(columns=drop_other_samples, inplace=True, errors='ignore')
                            # Exclude kept sample columns from feature_cols if they slipped in
                            feature_cols = [c for c in feature_cols if c not in raw_sample_cols_present]

                            # Final order: features + kept sample intensities + stats + remainder
                            final_cols = feature_cols + keep_sample_cols + stat_present + remainder
                            # Drop explicitly unwanted columns if they slipped in
                            drop_unwanted = [
                                'neg_log10_p_adj',
                                'PC3_2D_2','PC3_2D_3','PC3_3D_1','PC3_3D_2','PC3_3D_3',
                                'DU145_3D_1','DU145_2D_1','DU145_2D_2','DU145_2D_3','DU145_3D_2','DU145_3D_3','PC3_2D_1',
                                'DU145_3D_vs_DU145_2D_adj_p','DU145_2D_vs_PC3_2D_adj_p','PC3_3D_vs_PC3_2D_adj_p','DU145_3D_vs_PC3_3D_adj_p'
                            ]
                            sub.drop(columns=[c for c in drop_unwanted if c in sub.columns], inplace=True, errors='ignore')
                            # Ensure uniqueness while preserving order
                            sub = sub.loc[:, [c for c in dict.fromkeys(final_cols).keys() if c in sub.columns]]
                            sheet_name = f'{g1}_vs_{g2}'[:31]
                            sub.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1
                    else:
                        # Pairwise data exists but no group1/group2 columns - write as single Combined Pairwise sheet
                        self._thread_safe_log(f"⚠️ Writing pairwise data as single 'Combined Pairwise' sheet (no group columns)\n")
                        prog("Writing Combined Pairwise sheet...")
                        pairwise.to_excel(writer, sheet_name='Combined Pairwise', index=False)
                        sheets_written.append('Combined Pairwise')
                        step_holder[0] += 1
                
                # Original normalized data
                if hasattr(self, 'normalized_combined_df') and self.normalized_combined_df is not None:
                    prog("Writing Original Normalized data...")
                    self.normalized_combined_df.to_excel(writer, sheet_name='Original Normalized', index=False)
                    sheets_written.append('Original Normalized')
                    step_holder[0] += 1
                
            # Finalize
            prog("Finalizing export...")
            # Remove read-only attribute
            self._remove_readonly_attribute(filename)
            self._thread_safe_log(f"✅ Successfully wrote {len(sheets_written)} sheets: {', '.join(sheets_written)}\n")
            self._thread_safe_log(f"✅ Exported statistical results to: {os.path.basename(filename)}\n")
            
            # Store export path for Open Folder button
            self.last_statistics_export_path = filename

            # After export completes, auto-switch to Visualization tab per user preference
            try:
                self.root.after(0, lambda: self.notify_data_ready("📊 Visualization", "statistical_results"))
                # Don't switch tab yet - wait for class export to complete if needed
                # self.root.after(300, lambda: self.switch_to_tab("📊 Visualization"))
                self._thread_safe_log(f"📊 Ready to switch to Visualization tab after class export (if any) completes.\n")
            except Exception:
                pass
        except Exception as e:
            self._thread_safe_log(f"❌ Export failed: {e}\n")
            self.root.after(0, lambda: messagebox.showerror('Export Error', f'Failed to export results: {str(e)}'))
        finally:
            # Ensure progress bar hidden and button re-enabled
            self.root.after(0, self.hide_stats_progress)
            def reenable():
                try:
                    for widget in self.stats_tab.winfo_children():
                        if isinstance(widget, tk.Frame):
                            for child in widget.winfo_children():
                                if isinstance(child, tk.Button) and 'Export Statistical Results' in child.cget('text'):
                                    child.config(state='normal')
                except Exception:
                    pass
            self.root.after(0, reenable)

    def _export_class_statistical_results_threaded(self, filename):
        """Export lipid class statistical results to Excel with progress updates."""
        from datetime import datetime
        sheets_written = []
        try:
            # Debug: Log what class statistical results are available
            self._thread_safe_log(f"📊 Class Export Debug - Available results keys: {list(self.statistical_test_results_class.keys())}\n")
            
            # Determine number of steps
            total_steps = 3  # baseline: setup, write original normalized class, finalize
            
            # Check for Two-Way ANOVA results
            if 'two_way_anova' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_anova'].empty:
                total_steps += 1
            if 'two_way_posthoc_A' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_A'].empty:
                total_steps += 1
            if 'two_way_posthoc_B' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_B'].empty:
                total_steps += 1
            if 'two_way_posthoc_pairwise' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_pairwise'].empty:
                total_steps += 1  # Combined pairwise sheet
                # Count individual pairwise sheets
                phpw = self.statistical_test_results_class['two_way_posthoc_pairwise']
                if 'group1' in phpw.columns and 'group2' in phpw.columns:
                    num_individual_sheets = phpw.groupby(['group1','group2']).ngroups
                    total_steps += num_individual_sheets
            
            # Add step for Complete Results if available
            if 'enhanced_metabolites' in self.statistical_test_results_class and not self.statistical_test_results_class['enhanced_metabolites'].empty:
                total_steps += 1
                self._thread_safe_log(f"📊 Class enhanced metabolites available - will write Complete Results sheet\n")
            else:
                self._thread_safe_log(f"📊 No class enhanced metabolites data available\n")
            
            # Add steps for pairwise comparisons if available
            pairwise = None
            if 'pairwise' in self.statistical_test_results_class and not self.statistical_test_results_class['pairwise'].empty:
                pairwise = self.statistical_test_results_class['pairwise'].copy()
                # Check if group1/group2 columns exist
                if 'group1' in pairwise.columns and 'group2' in pairwise.columns:
                    comparisons = pairwise.groupby(['group1','group2']).ngroups
                    total_steps += comparisons
                    self._thread_safe_log(f"📊 Class pairwise data available - will write {comparisons} comparison sheets\n")
                else:
                    self._thread_safe_log(f"⚠️ Class pairwise data exists but lacks group1/group2 columns\n")
            else:
                self._thread_safe_log(f"📊 No class pairwise comparison data available\n")
            
            step = 1
            step_holder = [step]
            complete_df = None  # Initialize for later use in pairwise section
            new_format_class_export = False

            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                # Write a bootstrap visible sheet first so workbook save cannot fail
                # with "At least one sheet must be visible" if any later branch errors.
                pd.DataFrame({
                    'status': ['Class export started'],
                    'timestamp': [datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
                }).to_excel(writer, sheet_name='_Export_Status', index=False)
                sheets_written.append('_Export_Status')
                
                # NEW FORMAT: Check if we have the new complete_result format (for non-parametric two-way ANOVA)
                if 'complete_result' in self.statistical_test_results_class and not self.statistical_test_results_class['complete_result'].empty:
                    # NEW FORMAT: Write Complete Results sheet
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Complete Results sheet (new format)...\n")
                    complete_df = self.statistical_test_results_class['complete_result'].copy()
                    new_format_class_export = True

                    # Ensure Complete Results carries normalized class sample columns like normal-data export.
                    try:
                        if hasattr(self, 'normalized_combined_class_df') and isinstance(self.normalized_combined_class_df, pd.DataFrame):
                            norm_class = self.normalized_combined_class_df.copy()
                            join_col = 'Lipid_Class' if 'Lipid_Class' in complete_df.columns and 'Lipid_Class' in norm_class.columns else None
                            if join_col is None and 'Class' in complete_df.columns and 'Class' in norm_class.columns:
                                join_col = 'Class'
                            if join_col is not None:
                                # Merge only missing normalized columns to avoid overwriting computed stats.
                                missing_cols = [c for c in norm_class.columns if c not in complete_df.columns]
                                if missing_cols:
                                    complete_df = complete_df.merge(
                                        norm_class[[join_col] + missing_cols],
                                        on=join_col,
                                        how='left'
                                    )
                    except Exception:
                        pass
                    
                    # Drop fully NaN pairwise stat columns to reduce clutter
                    try:
                        drop_cols = []
                        for c in complete_df.columns:
                            if ('_vs_' in c and any(suf in c for suf in ['_log2FC','_adj_p','_p_adj','_neg_log10'])):
                                if complete_df[c].isna().all():
                                    drop_cols.append(c)
                        if drop_cols:
                            complete_df = complete_df.drop(columns=drop_cols)
                            self._thread_safe_log(f"🧹 Dropped {len(drop_cols)} all-NaN pairwise columns before export\n")
                    except Exception:
                        pass
                    
                    complete_df.to_excel(writer, sheet_name='Complete Results', index=False)
                    sheets_written.append('Complete Results')
                    step_holder[0] += 1
                    
                    # NEW FORMAT: Write individual pairwise sheets
                    wrote_class_pairwise_split = False
                    for sheet_name, sheet_df in self.statistical_test_results_class.items():
                        # Skip duplicate enhanced_metabolites (same as complete_result)
                        if sheet_name in ('complete_result','enhanced_metabolites'):
                            continue
                        if isinstance(sheet_df, pd.DataFrame) and not sheet_df.empty:
                            # Make a copy to avoid modifying the cached results
                            sheet_df_to_write = sheet_df.copy()
                            # For each sheet, remove all-NaN pairwise stat columns
                            try:
                                drop_cols2 = []
                                for c in sheet_df_to_write.columns:
                                    if ('_vs_' in c and any(suf in c for suf in ['_log2FC','_adj_p','_p_adj','_neg_log10'])):
                                        if sheet_df_to_write[c].isna().all():
                                            drop_cols2.append(c)
                                if drop_cols2:
                                    sheet_df_to_write = sheet_df_to_write.drop(columns=drop_cols2)
                                    self._thread_safe_log(f"🧹 Dropped {len(drop_cols2)} all-NaN columns from class {sheet_name} sheet\n")
                            except Exception:
                                pass
                            
                            safe_sheet_name = self._safe_excel_sheet_name(sheet_name, sheets_written)
                            self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class {safe_sheet_name} sheet...\n")
                            sheet_df_to_write.to_excel(writer, sheet_name=safe_sheet_name, index=False)
                            sheets_written.append(safe_sheet_name)
                            step_holder[0] += 1

                            # If only a combined pairwise table is present, also emit per-comparison sheets.
                            if (
                                sheet_name == 'pairwise'
                                and 'group1' in sheet_df_to_write.columns
                                and 'group2' in sheet_df_to_write.columns
                            ):
                                for (g1, g2), sub_idx in sheet_df_to_write.groupby(['group1', 'group2']).groups.items():
                                    sub = sheet_df_to_write.loc[sub_idx].copy()
                                    split_sheet_name = self._safe_excel_sheet_name(f'{g1}_vs_{g2}', sheets_written)
                                    self._thread_safe_log(
                                        f"[{step_holder[0]}/{total_steps}] Writing Class {split_sheet_name} sheet...\n"
                                    )
                                    sub.to_excel(writer, sheet_name=split_sheet_name, index=False)
                                    sheets_written.append(split_sheet_name)
                                    step_holder[0] += 1
                                wrote_class_pairwise_split = True

                    # Fallback for new-format results where pairwise wasn't emitted in loop above.
                    if (
                        not wrote_class_pairwise_split
                        and 'pairwise' in self.statistical_test_results_class
                        and isinstance(self.statistical_test_results_class['pairwise'], pd.DataFrame)
                        and not self.statistical_test_results_class['pairwise'].empty
                    ):
                        pairwise_df = self.statistical_test_results_class['pairwise'].copy()
                        if 'group1' in pairwise_df.columns and 'group2' in pairwise_df.columns:
                            for (g1, g2), sub_idx in pairwise_df.groupby(['group1', 'group2']).groups.items():
                                split_sheet_name = self._safe_excel_sheet_name(f'{g1}_vs_{g2}', sheets_written)
                                self._thread_safe_log(
                                    f"[{step_holder[0]}/{total_steps}] Writing Class {split_sheet_name} sheet...\n"
                                )
                                pairwise_df.loc[sub_idx].to_excel(writer, sheet_name=split_sheet_name, index=False)
                                sheets_written.append(split_sheet_name)
                                step_holder[0] += 1
                    
                    self._thread_safe_log(f"✅ Class export complete (new format): {len(sheets_written)} sheets written\n")
                
                # OLD FORMAT: Two-Way ANOVA sheet for class data (if present)
                elif 'two_way_anova' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_anova'].empty:
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Two-Way ANOVA sheet...\n")
                    class_twa_df = self.statistical_test_results_class['two_way_anova'].copy()
                    
                    # Augment with posthoc columns if available
                    try:
                        import numpy as np
                        fa_name = getattr(self, 'factorA_name', 'Factor A')
                        fb_name = getattr(self, 'factorB_name', 'Factor B')
                        
                        # Build wide pairwise columns for Factor A
                        if 'two_way_posthoc_A' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_A'].empty:
                            pha = self.statistical_test_results_class['two_way_posthoc_A']
                            if 'group1' in pha.columns and 'group2' in pha.columns:
                                id_col = 'Class' if 'Class' in pha.columns else 'metabolite'
                                wide_parts = []
                                for (g1, g2), sub in pha.groupby(['group1','group2']):
                                    comp = f"{fa_name}:{g1}_vs_{g2}"
                                    part = sub[[id_col,'p_value_adj']].copy()
                                    part = part.rename(columns={'p_value_adj': f'{comp}_adj_p'})
                                    part[f'{comp}_neg_log10_adj_p'] = -np.log10(part[f'{comp}_adj_p'].replace(0, np.finfo(float).eps))
                                    wide_parts.append(part)
                                if wide_parts:
                                    from functools import reduce
                                    wideA = reduce(lambda a,b: a.merge(b, on=id_col, how='outer'), wide_parts)
                                    merge_col = 'Class' if 'Class' in class_twa_df.columns else 'metabolite'
                                    if merge_col in class_twa_df.columns:
                                        class_twa_df = class_twa_df.merge(wideA, on=merge_col, how='left')
                        
                        # Build wide pairwise columns for Factor B
                        if 'two_way_posthoc_B' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_B'].empty:
                            phb = self.statistical_test_results_class['two_way_posthoc_B']
                            if 'group1' in phb.columns and 'group2' in phb.columns:
                                id_col = 'Class' if 'Class' in phb.columns else 'metabolite'
                                wide_parts = []
                                for (g1, g2), sub in phb.groupby(['group1','group2']):
                                    comp = f"{fb_name}:{g1}_vs_{g2}"
                                    part = sub[[id_col,'p_value_adj']].copy()
                                    part = part.rename(columns={'p_value_adj': f'{comp}_adj_p'})
                                    part[f'{comp}_neg_log10_adj_p'] = -np.log10(part[f'{comp}_adj_p'].replace(0, np.finfo(float).eps))
                                    wide_parts.append(part)
                                if wide_parts:
                                    from functools import reduce
                                    wideB = reduce(lambda a,b: a.merge(b, on=id_col, how='outer'), wide_parts)
                                    merge_col = 'Class' if 'Class' in class_twa_df.columns else 'metabolite'
                                    if merge_col in class_twa_df.columns:
                                        class_twa_df = class_twa_df.merge(wideB, on=merge_col, how='left')
                    except Exception:
                        pass
                    
                    class_twa_df.to_excel(writer, sheet_name='TwoWayANOVA', index=False)
                    sheets_written.append('TwoWayANOVA')
                    step_holder[0] += 1
                
                # Two-Way ANOVA Posthoc sheets for class data (if present)
                if (not new_format_class_export) and 'two_way_posthoc_A' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_A'].empty:
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Two-Way Posthoc A sheet...\n")
                    pha = self.statistical_test_results_class['two_way_posthoc_A'].copy()
                    import numpy as np
                    if 'p_value_adj' in pha.columns and 'adj_p' not in pha.columns:
                        pha['adj_p'] = pha['p_value_adj']
                    if 'neg_log10_adj_p' not in pha.columns and 'adj_p' in pha.columns:
                        pha['neg_log10_adj_p'] = -np.log10(pha['adj_p'].replace(0, np.finfo(float).eps))
                    pha.to_excel(writer, sheet_name='TwoWay_Posthoc_A', index=False)
                    sheets_written.append('TwoWay_Posthoc_A')
                    step_holder[0] += 1
                
                if (not new_format_class_export) and 'two_way_posthoc_B' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_B'].empty:
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Two-Way Posthoc B sheet...\n")
                    phb = self.statistical_test_results_class['two_way_posthoc_B'].copy()
                    import numpy as np
                    if 'p_value_adj' in phb.columns and 'adj_p' not in phb.columns:
                        phb['adj_p'] = phb['p_value_adj']
                    if 'neg_log10_adj_p' not in phb.columns and 'adj_p' in phb.columns:
                        phb['neg_log10_adj_p'] = -np.log10(phb['adj_p'].replace(0, np.finfo(float).eps))
                    phb.to_excel(writer, sheet_name='TwoWay_Posthoc_B', index=False)
                    sheets_written.append('TwoWay_Posthoc_B')
                    step_holder[0] += 1
                
                # Two-Way ANOVA All Pairwise Comparisons for class data (if present) - combined sheet AND individual sheets
                if (not new_format_class_export) and 'two_way_posthoc_pairwise' in self.statistical_test_results_class and not self.statistical_test_results_class['two_way_posthoc_pairwise'].empty:
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Two-Way Pairwise Comparisons sheet...\n")
                    phpw = self.statistical_test_results_class['two_way_posthoc_pairwise'].copy()
                    import numpy as np
                    if 'p_value_adj' in phpw.columns and 'adj_p' not in phpw.columns:
                        phpw['adj_p'] = phpw['p_value_adj']
                    if 'neg_log10_adj_p' not in phpw.columns and 'adj_p' in phpw.columns:
                        phpw['neg_log10_adj_p'] = -np.log10(phpw['adj_p'].replace(0, np.finfo(float).eps))
                    phpw.to_excel(writer, sheet_name='TwoWay_Pairwise', index=False)
                    sheets_written.append('TwoWay_Pairwise')
                    step_holder[0] += 1
                    
                    # Debug: Check for group columns
                    self._thread_safe_log(f"📋 Class pairwise columns: {list(phpw.columns)}\n")
                    self._thread_safe_log(f"📋 Has group1/group2: {'group1' in phpw.columns and 'group2' in phpw.columns}\n")
                    
                    # Also create individual comparison sheets for Two-Way pairwise (like metabolite export)
                    if 'group1' in phpw.columns and 'group2' in phpw.columns:
                        self._thread_safe_log(f"📊 Creating individual pairwise comparison sheets...\n")
                        tw_pairwise_class = phpw.copy()
                        
                        # Get unique group pairs
                        unique_pairs = phpw.groupby(['group1','group2']).size().reset_index(name='count')
                        self._thread_safe_log(f"📊 Found {len(unique_pairs)} unique group pairs\n")
                        
                        for idx, (g1, g2) in enumerate(unique_pairs[['group1','group2']].values, 1):
                            self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing {g1}_vs_{g2} sheet...\n")
                            sub = phpw[(phpw['group1'] == g1) & (phpw['group2'] == g2)].copy()
                            
                            if sub.empty:
                                self._thread_safe_log(f"⚠️ No data for {g1} vs {g2}\n")
                                continue
                            
                            # Rename mean columns to group-specific names
                            if 'mean_group1' in sub.columns:
                                sub.rename(columns={'mean_group1': f'mean_{g1}'}, inplace=True)
                            if 'mean_group2' in sub.columns:
                                sub.rename(columns={'mean_group2': f'mean_{g2}'}, inplace=True)
                            
                            # Ensure -log10 p-value exists
                            if 'adj_p' in sub.columns and 'neg_log10_adj_p' not in sub.columns:
                                sub['neg_log10_adj_p'] = -np.log10(sub['adj_p'].replace(0, np.finfo(float).eps))
                            
                            # Column ordering: Class, then means, stats, etc.
                            feature_cols = [c for c in ['Lipid_Class', 'Class', 'Class_name', 'metabolite_id', 'Name'] if c in sub.columns]
                            stat_cols = [f'mean_{g1}', f'mean_{g2}', 'FC', 'log2FC',
                                       'p_value', 'adj_p', 'p_value_adj', 'neg_log10_adj_p', 'neg_log10_p_adj',
                                       'model_effect', 'model_se', 'ci_lower_95', 'ci_upper_95',
                                       'statistic', 'cohen_d', 'cliffs_delta', 'Effect_Size', 'Expression',
                                       'n_group1', 'n_group2', 'group1', 'group2']
                            stat_present = [c for c in stat_cols if c in sub.columns]
                            
                            # Ensure Lipid_Class is first, then stats
                            final_cols = feature_cols + stat_present
                            sub = sub.loc[:, [c for c in final_cols if c in sub.columns]]
                            
                            # Write sheet with just {g1}_vs_{g2} naming (not TwoWay_ prefix)
                            sheet_name = self._safe_excel_sheet_name(f'{g1}_vs_{g2}', sheets_written)
                            sub.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1
                    else:
                        self._thread_safe_log(f"⚠️ Class pairwise data has no group1/group2 columns - skipping individual sheets\n")
                
                # Complete Results sheet (independent of pairwise data)
                if (not new_format_class_export) and 'enhanced_metabolites' in self.statistical_test_results_class and not self.statistical_test_results_class['enhanced_metabolites'].empty:
                    self._thread_safe_log(f"[{step_holder[0]}/{total_steps}] Writing Class Complete Results sheet...\n")
                    complete_df = self.statistical_test_results_class['enhanced_metabolites'].copy()
                    self._thread_safe_log(f"📊 Class Complete Results - Available columns: {list(complete_df.columns)}\n")
                    # General duplicate overall_* cleanup
                    complete_df = self._clean_complete_results_dataframe(complete_df)
                    # Class-specific cleanup (deduplicate, drop helper cols, reorder)
                    complete_df = self._clean_class_complete_results_dataframe(complete_df)
                    # Lipid class: remove metabolite_id if present
                    if 'metabolite_id' in complete_df.columns:
                        complete_df.drop(columns=['metabolite_id'], inplace=True, errors='ignore')
                    complete_df.to_excel(writer, sheet_name='Complete Results', index=False)
                    sheets_written.append('Complete Results')
                    step_holder[0] += 1
                else:
                    complete_df = None
                
                # Pairwise comparison sheets for class data.
                if (not new_format_class_export) and pairwise is not None:
                    # Merge pairwise with enhanced_metabolites to get all feature columns
                    if 'enhanced_metabolites' in self.statistical_test_results_class and not self.statistical_test_results_class['enhanced_metabolites'].empty:
                        # Prefer to derive feature columns from the already-cleaned Complete Results
                        # to guarantee identical feature set/order between Complete Results and pairwise sheets.
                        if complete_df is None:
                            complete_df = self.statistical_test_results_class['enhanced_metabolites'].copy()
                            complete_df = self._clean_complete_results_dataframe(complete_df)
                            complete_df = self._clean_class_complete_results_dataframe(complete_df)

                        base_df = complete_df.copy()
                        # Identify feature columns from Complete Results: everything before mean_/n_/comparisons/overall stats
                        exclude_patterns = ['mean_', 'n_', '_vs_', 'overall_', 'p_value', 'neg_log10', 'statistic', 'FC', 'log2FC', 'Expression']
                        feature_cols_to_merge = [c for c in base_df.columns if not any(pat in c for pat in exclude_patterns)]
                        # Canonicalize: prefer Lipid_Class; drop raw 'Class' if redundant
                        if 'Class' in feature_cols_to_merge and 'Lipid_Class' in feature_cols_to_merge:
                            feature_cols_to_merge = [c for c in feature_cols_to_merge if c != 'Class']

                        # Build feature_df
                        feature_df = base_df[feature_cols_to_merge].copy()
                        # Ensure a join key
                        join_key = None
                        if 'metabolite_id' in feature_df.columns:
                            join_key = 'metabolite_id'
                        elif 'Lipid_Class' in feature_df.columns:
                            join_key = 'Lipid_Class'

                        if join_key is not None:
                            # If pairwise lacks metabolite_id but has numeric 'metabolite', align by index as fallback
                            if 'metabolite' in pairwise.columns and join_key not in pairwise.columns:
                                pairwise['_merge_idx'] = pairwise.index
                                feature_df['_merge_idx'] = feature_df.index
                                pairwise = pairwise.merge(
                                    feature_df,
                                    on='_merge_idx',
                                    how='left',
                                    suffixes=('', '_dup')
                                )
                                pairwise.drop(columns=['_merge_idx'], inplace=True, errors='ignore')
                                # Remove placeholder metabolite column if we now have a better ID
                                if join_key in pairwise.columns:
                                    pairwise.drop(columns=['metabolite'], inplace=True, errors='ignore')
                            else:
                                pairwise = pairwise.merge(
                                    feature_df,
                                    left_on=join_key,
                                    right_on=join_key,
                                    how='left',
                                    suffixes=('', '_dup')
                                )

                            # Drop merge duplicates and redundant 'Class' if Lipid_Class present
                            dup_cols = [c for c in pairwise.columns if c.endswith('_dup')]
                            if dup_cols:
                                pairwise.drop(columns=dup_cols, inplace=True, errors='ignore')
                            if 'Lipid_Class' in pairwise.columns and 'Class' in pairwise.columns:
                                # Drop 'Class' when it's identical to Lipid_Class
                                try:
                                    if pairwise['Lipid_Class'].equals(pairwise['Class']):
                                        pairwise.drop(columns=['Class'], inplace=True, errors='ignore')
                                except Exception:
                                    pass

                            self._thread_safe_log(f"✅ Merged class feature columns ({len(feature_cols_to_merge)}) into pairwise data to match Complete Results\n")
                    
                    # Normalize naming
                    if 'dunn_bh_pvalue' in pairwise.columns:
                        pairwise.rename(columns={'dunn_bh_pvalue':'adj_p', 'neg_log10_dunn_bh_pvalue':'neg_log10_adj_p'}, inplace=True)
                    elif 'p_value_adj' in pairwise.columns and 'adj_p' not in pairwise.columns:
                        pairwise.rename(columns={'p_value_adj':'adj_p'}, inplace=True)
                    # Precompute sample → group label map for class data as well
                    sample_group_map = {}
                    try:
                        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                            for col, var in self.sample_group_vars.items():
                                sample_group_map[col] = var.get()
                    except Exception:
                        pass
                    sample_cols_all = set(sample_group_map.keys())

                    # Only process individual comparison sheets if group1/group2 columns exist
                    if 'group1' in pairwise.columns and 'group2' in pairwise.columns:
                        for (g1, g2), sub_idx in pairwise.groupby(['group1','group2']).groups.items():
                            self._thread_safe_log(f"Writing Class {g1}_vs_{g2} sheet...\n")
                            sub = pairwise.loc[sub_idx].copy()
                            if 'fold_change' in sub.columns:
                                sub.rename(columns={'fold_change':'FC', 'log2_fold_change':'log2FC'}, inplace=True)
                            sub.rename(columns={'mean_group1': f'mean_{g1}', 'mean_group2': f'mean_{g2}'}, inplace=True)
                            # Compute appropriate -log10 column based on availability and user preference
                            import numpy as np
                            if 'adj_p' in sub.columns:
                                if 'neg_log10_adj_p' not in sub.columns:
                                    sub['neg_log10_adj_p'] = -np.log10(sub['adj_p'].replace(0, np.finfo(float).eps))
                            elif 'p_value' in sub.columns:
                                sub['neg_log10_p'] = -np.log10(sub['p_value'].replace(0, np.finfo(float).eps))
                            
                            # For class data: use the same feature columns as in the Complete Results
                            # Retain only raw sample intensity columns for the two groups being compared
                            raw_sample_cols_present = [c for c in sub.columns if c in sample_cols_all]
                            keep_sample_cols = []
                            for rc in raw_sample_cols_present:
                                grp = sample_group_map.get(rc)
                                if grp in (g1, g2):
                                    keep_sample_cols.append(rc)
                            drop_other_samples = [c for c in raw_sample_cols_present if c not in keep_sample_cols]
                            if drop_other_samples:
                                sub.drop(columns=drop_other_samples, inplace=True, errors='ignore')
                            feature_cols = []
                            if complete_df is not None:
                                exclude_patterns = ['mean_', 'n_', '_vs_', 'overall_', 'p_value', 'neg_log10', 'statistic', 'FC', 'log2FC', 'Expression']
                                complete_features = [c for c in complete_df.columns if not any(pat in c for pat in exclude_patterns)]
                                # Prefer Lipid_Class over Class to avoid duplication
                                if 'Lipid_Class' in complete_features and 'Class' in complete_features:
                                    complete_features = [c for c in complete_features if c != 'Class']
                                feature_cols = [c for c in complete_features if c in sub.columns]
                                
                                # Ensure Class_name is after Class/Lipid_Class in lipid mode
                                if 'Class_name' in sub.columns and 'Class_name' not in feature_cols:
                                    # Find position of Class or Lipid_Class and insert Class_name after it
                                    insert_pos = len(feature_cols)
                                    for i, col in enumerate(feature_cols):
                                        if col in ['Lipid_Class', 'Class']:
                                            insert_pos = i + 1
                                            break
                                    feature_cols.insert(insert_pos, 'Class_name')
                            else:
                                # Fallback minimal set
                                for fc in ['Lipid_Class', 'Class', 'Class_name', 'metabolite_id', 'Name', 'Compound']:
                                    if fc in sub.columns and fc not in feature_cols:
                                        feature_cols.append(fc)
                            
                            # Statistical columns
                            stat_cols = [f'mean_{g1}', f'mean_{g2}', 'FC', 'log2FC',
                                         'p_value', 'adj_p', 'p_value_adj',
                                         'neg_log10_p', 'neg_log10_adj_p', 'neg_log10_p_adj',
                                         'statistic',
                                         'cohen_d', 'cliffs_delta', 'Expression',
                                         'n_group1', 'n_group2']
                            stat_present = [c for c in stat_cols if c in sub.columns]
                            
                            # Remaining columns except helpers
                            remainder = [c for c in sub.columns 
                                       if c not in feature_cols + stat_present 
                                       and c not in ['group1', 'group2', 'metabolite']]
                            
                            # Final order and drop unwanted columns per user request
                            final_cols = feature_cols + keep_sample_cols + stat_present + remainder
                            drop_unwanted = [
                                'neg_log10_p_adj',
                                'PC3_2D_2','PC3_2D_3','PC3_3D_1','PC3_3D_2','PC3_3D_3',
                                'DU145_3D_1','DU145_2D_1','DU145_2D_2','DU145_2D_3','DU145_3D_2','DU145_3D_3','PC3_2D_1',
                                'DU145_3D_vs_DU145_2D_adj_p','DU145_2D_vs_PC3_2D_adj_p','PC3_3D_vs_PC3_2D_adj_p','DU145_3D_vs_PC3_3D_adj_p'
                            ]
                            sub.drop(columns=[c for c in drop_unwanted if c in sub.columns], inplace=True, errors='ignore')
                            sub = sub.loc[:, [c for c in dict.fromkeys(final_cols).keys() if c in sub.columns]]
                            sheet_name = self._safe_excel_sheet_name(f'{g1}_vs_{g2}', sheets_written)
                            sub.to_excel(writer, sheet_name=sheet_name, index=False)
                            sheets_written.append(sheet_name)
                            step_holder[0] += 1
                    else:
                        # Class pairwise data exists but no group1/group2 columns - write as single sheet
                        self._thread_safe_log(f"⚠️ Writing class pairwise data as single 'Combined Pairwise' sheet (no group columns)\n")
                        pairwise.to_excel(writer, sheet_name='Combined Pairwise', index=False)
                        sheets_written.append('Combined Pairwise')
                        step_holder[0] += 1

                # openpyxl requires at least one visible worksheet in a workbook.
                if not sheets_written:
                    self._thread_safe_log("⚠️ No class exportable results found; writing fallback sheet.\n")
                    pd.DataFrame({
                        'status': ['No class statistical results available for export.']
                    }).to_excel(writer, sheet_name='Summary', index=False)
                    sheets_written.append('Summary')
                
                # Original normalized class data
                if hasattr(self, 'normalized_combined_class_df') and self.normalized_combined_class_df is not None:
                    self._thread_safe_log("Writing Original Normalized Class data...\n")
                    self.normalized_combined_class_df.to_excel(writer, sheet_name='Original Normalized', index=False)
                    sheets_written.append('Original Normalized')
                    step_holder[0] += 1

                # Remove bootstrap sheet if real sheets were written.
                if '_Export_Status' in sheets_written and len(sheets_written) > 1:
                    try:
                        if '_Export_Status' in writer.book.sheetnames:
                            writer.book.remove(writer.book['_Export_Status'])
                        sheets_written = [s for s in sheets_written if s != '_Export_Status']
                    except Exception as cleanup_err:
                        self._thread_safe_log(f"⚠️ Could not remove bootstrap class sheet: {cleanup_err}\n")
                
            # Finalize
            # Remove read-only attribute
            self._remove_readonly_attribute(filename)
            self._thread_safe_log(f"✅ Successfully wrote {len(sheets_written)} class sheets: {', '.join(sheets_written)}\n")
            self._thread_safe_log(f"✅ Exported lipid class statistical results to: {os.path.basename(filename)}\n")
            # Don't show popup here - it will be shown by the monitor thread after both exports complete
        except Exception as e:
            # Capture exception message for use inside the Tkinter callback
            err_msg = str(e)
            self._thread_safe_log(f"❌ Class export failed: {err_msg}\n")
            try:
                import traceback
                self._thread_safe_log(f"❌ Class export traceback:\n{traceback.format_exc()}\n")
            except Exception:
                pass
            self.root.after(0, lambda err=err_msg: messagebox.showerror('Class Export Error', f'Failed to export class results: {err}'))

    def _export_class_statistical_results_threaded_with_wait(self, filename, main_thread):
        """Export class results after waiting for main export to complete."""
        # Wait for main export thread to finish
        main_thread.join(timeout=600)  # Wait up to 10 minutes
        # Now export class results
        self._export_class_statistical_results_threaded(filename)

    def _monitor_export_completion(self, export_threads, filename):
        """Monitor all export threads and show completion popup when all are done."""
        import time
        # Wait for all export threads to complete
        for thread in export_threads:
            thread.join(timeout=600)  # Wait up to 10 minutes for each thread
        
        # Brief delay to ensure all file writes are complete
        time.sleep(0.5)
        
        # Determine which files were exported
        base_name, ext = os.path.splitext(filename)
        class_filename = f"{base_name}_class{ext}"
        has_class_export = len(export_threads) > 1 and os.path.exists(class_filename)
        
        # Build message
        if has_class_export:
            msg = f'Statistical results exported successfully!\n\nFiles saved:\n• {os.path.basename(filename)}\n• {os.path.basename(class_filename)}'
        else:
            msg = f'Statistical results exported to:\n{filename}'
        
        # Show completion popup
        self.root.after(0, lambda m=msg: messagebox.showinfo('Export Complete', m))
        
        # Switch to Visualization tab
        try:
            self.root.after(300, lambda: self.switch_to_tab("📊 Visualization"))
            self._thread_safe_log(f"📊 Visualization tab opened after all exports complete.\n")
        except Exception:
            pass

    def _safe_excel_sheet_name(self, desired_name, existing_names=None):
        """Return an Excel-safe, unique worksheet name (max 31 chars)."""
        import re

        existing = set(existing_names or [])
        cleaned = re.sub(r'[:\\/?*\[\]]', '_', str(desired_name or 'Sheet')).strip()
        if not cleaned:
            cleaned = 'Sheet'
        base = cleaned[:31]
        candidate = base
        counter = 1
        while candidate in existing:
            suffix = f"_{counter}"
            candidate = f"{base[:31-len(suffix)]}{suffix}"
            counter += 1
        return candidate

    def _clean_complete_results_dataframe(self, df):
        """Remove duplicate overall_* columns and standardize naming for export.

        Handles cases where earlier merges produced columns like
        overall_statistic_x / overall_statistic_y or overall_p_value_x / _y.
        Preference order: _x then (no suffix) then _y. Drops unused duplicates.
        """
        import numpy as np
        # List of base names to collapse
        bases = [
            'overall_statistic',
            'overall_p_value',
            'overall_p_value_adj',
            'overall_neg_log10_p_adj'
        ]
        for base in bases:
            cols = [c for c in df.columns if c == base or c.startswith(base + '_')]
            if len(cols) <= 1:
                continue
            # choose best candidate: prefer _x if present (typical from first frame), else base, else _y
            preferred = None
            for cand in [f'{base}_x', base, f'{base}_y', f'{base}_dup']:
                if cand in cols:
                    preferred = cand
                    break
            if preferred is None:
                preferred = cols[0]
            # Move data into canonical base name
            if preferred != base:
                df[base] = df[preferred]
            # Drop others
            for c in cols:
                if c != base:
                    df.drop(columns=c, inplace=True, errors='ignore')
        # Recalculate neg_log10 if missing but adj p present
        if 'overall_p_value_adj' in df.columns and 'overall_neg_log10_p_adj' in df.columns:
            # Ensure numeric and recompute to guarantee consistency
            with np.errstate(divide='ignore'):
                df['overall_neg_log10_p_adj'] = -np.log10(
                    df['overall_p_value_adj'].replace(0, np.finfo(float).eps)
                )
        return df

    def _clean_class_complete_results_dataframe(self, df):
        """Class-specific cleanup for Complete Results export.

        Ensures:
        - Lipid_Class column present (fallback from Class/metabolite/metabolite_id)
        - Deduplicate on Lipid_Class
        - Drop helper columns: Area (Max.), duplicate metabolite_id
        - Collapse *_x/*_y duplicate comparison columns
        - Order columns: Lipid_Class, metadata, raw sample columns, mean_/n_ columns, comparisons, overall stats, leftovers
        """
        try:
            # 1. Guarantee Lipid_Class
            if 'Lipid_Class' not in df.columns:
                for cand in ['Class_x', 'Class_y', 'Class', 'metabolite_id', 'metabolite']:
                    if cand in df.columns:
                        df['Lipid_Class'] = df[cand]
                        break
            
            # If still no Lipid_Class, log and return early
            if 'Lipid_Class' not in df.columns:
                self._thread_safe_log(f"⚠️  Warning: No Lipid_Class/Class/metabolite_id/metabolite column found in class results. Columns: {list(df.columns)}\n")
                return df
            
            # 2. Deduplicate
            before = len(df)
            numeric_cols = [
                col for col in df.columns
                if col != 'Lipid_Class' and pd.api.types.is_numeric_dtype(df[col])
            ]
            non_numeric_cols = [col for col in df.columns if col not in numeric_cols + ['Lipid_Class']]
            agg_map = {col: 'sum' for col in numeric_cols}
            agg_map.update({col: 'first' for col in non_numeric_cols})
            df = df.groupby('Lipid_Class', as_index=False, sort=False).agg(agg_map)
            removed = before - len(df)
            if removed:
                self._thread_safe_log(f"Class Complete Results: merged {removed} duplicate Lipid_Class rows by summing numeric values.\n")
            # 3. Drop helpers
            if 'Area (Max.)' in df.columns:
                df.drop(columns=['Area (Max.)'], inplace=True, errors='ignore')
            if 'metabolite_id' in df.columns and 'Lipid_Class' in df.columns and df['metabolite_id'].equals(df['Lipid_Class']):
                df.drop(columns=['metabolite_id'], inplace=True, errors='ignore')
            
            # 4. Collapse duplicate suffix columns - MUST PRESERVE Lipid_Class
            for col in list(df.columns):
                if col.endswith('_x') and col != 'Lipid_Class':  # Don't process Lipid_Class itself
                    base = col[:-2]
                    if base not in df.columns and base != 'Lipid_Class':  # Don't create Lipid_Class from _x
                        df[base] = df[col]
                    ycol = f"{base}_y"
                    if ycol in df.columns:
                        mask = df[base].isna() & df[ycol].notna()
                        df.loc[mask, base] = df.loc[mask, ycol]
                        df.drop(columns=[ycol], inplace=True, errors='ignore')
                    df.drop(columns=[col], inplace=True, errors='ignore')
            
            # Ensure Lipid_Class still exists after column collapsing
            if 'Lipid_Class' not in df.columns:
                self._thread_safe_log(f"⚠️  Warning: Lipid_Class column was lost during cleanup. Available columns: {list(df.columns)}\n")
                return df
            
            # 5. Order columns
            first = ['Lipid_Class'] if 'Lipid_Class' in df.columns else []
            meta_cols = [c for c in ['Class_name', 'Class', 'LipidID', 'Name', 'Compound'] if c in df.columns and c not in first]
            overall_cols = [c for c in ['overall_statistic','overall_p_value','overall_p_value_adj','overall_neg_log10_p_adj'] if c in df.columns]
            mean_cols = [c for c in df.columns if c.startswith('mean_') or c.startswith('n_')]
            comparison_cols = [c for c in df.columns if '_vs_' in c]

            # Preserve raw class sample columns so Complete Results mirrors non-class export behavior.
            sample_cols = []
            if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                sample_cols = [c for c in self.sample_group_vars.keys() if c in df.columns]
            if not sample_cols and hasattr(self, 'normalized_combined_class_df') and isinstance(self.normalized_combined_class_df, pd.DataFrame):
                norm_cols = set(self.normalized_combined_class_df.columns)
                sample_cols = [c for c in df.columns if c in norm_cols and c not in first + meta_cols]
            
            # Statistical factor columns (from Two-Way ANOVA): F_*, p_*, _Interaction, etc.
            factor_cols = [c for c in df.columns if c.startswith(('F_', 'p_')) or 'Interaction' in c or 'Cohen' in c or 'Cliffs' in c]
            
            # Remaining columns (metadata, cleanup columns, etc.)
            ordered = first + meta_cols + sample_cols + factor_cols + mean_cols + comparison_cols + overall_cols
            remainder = [c for c in df.columns if c not in ordered]
            final_cols = ordered + remainder
            
            # Remove duplicates while preserving order
            seen = set()
            final_cols_dedup = []
            for c in final_cols:
                if c not in seen and c in df.columns:
                    final_cols_dedup.append(c)
                    seen.add(c)
            
            df = df.loc[:, final_cols_dedup]
            return df
        except Exception as e:
            self._thread_safe_log(f"❌ Class Complete Results cleaning error: {e}\n")
            import traceback
            self._thread_safe_log(f"Traceback: {traceback.format_exc()}\n")
            return df

    # ---------------- Auto-preparation for Statistics Tab ----------------
    def auto_prepare_statistics_tab(self):
        """Switch to the Statistics tab after ID annotation/Lipid cleaning and populate logs/sample columns automatically."""
        import pandas as pd  # Ensure 'pd' is defined before any use inside this function scope
        # Find statistics tab index
        stats_index = None
        try:
            for i in range(self.notebook.index('end')):
                txt = self.notebook.tab(i, 'text')
                if 'Statistics' in str(txt):
                    stats_index = i
                    break
        except Exception:
            stats_index = None
        if stats_index is not None:
            self.notebook.select(stats_index)
        # Ensure stats_log exists
        if not hasattr(self, 'stats_log'):
            return
        
        # Determine if we should auto-load lipid or metabolite data
        # Check if lipid data was just cleaned
        has_lipid_data = False
        if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
            has_lipid_data = any(key in self.memory_store and self.memory_store[key] is not None 
                                for key in ['pos_lipid_df', 'lipid_pos_df', 'neg_lipid_df', 'lipid_neg_df'])
        
        if has_lipid_data:
            # Set mode to lipid and load lipid data
            if hasattr(self, 'statistics_data_mode'):
                self.statistics_data_mode.set('lipid')
            self.stats_log.insert(tk.END, '\n===== Auto Load From Lipid Cleaning =====\n')
            self._auto_load_lipid_data()
        else:
            # Default to metabolite mode
            if hasattr(self, 'statistics_data_mode'):
                self.statistics_data_mode.set('metabolite')
            self.stats_log.insert(tk.END, '\n===== Auto Load From ID Annotation =====\n')
            self._auto_load_metabolite_data()
    
    def _auto_load_metabolite_data(self):
        """Auto-load metabolite data from ID annotation."""
        import pandas as pd
        # Summarize memory store content relevant to statistics
        if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
            try:
                self.stats_log.insert(tk.END, f"Memory store keys: {list(self.memory_store.keys())}\n")
                # Debug: show size of ALL DataFrames in memory
                for key, value in self.memory_store.items():
                    if isinstance(value, pd.DataFrame):
                        self.stats_log.insert(tk.END, f"  {key}: {len(value)} rows, {len(value.columns)} columns (empty={value.empty})\n")
            except Exception as e:
                self.stats_log.insert(tk.END, f"Debug error: {e}\n")
        else:
            self.stats_log.insert(tk.END, "⚠️ No memory_store available!\n")
        # Prefer fully annotated polarity DataFrames (Pos_id / Neg_id) if already produced
        pos_df = None
        neg_df = None
        if hasattr(self, 'memory_store'):
            # Possible keys in memory depending on pipeline stage
            cand_pos_keys = ['pos_id_df','clean_pos_id_df','pos_enhanced_df','pos_result']
            chosen_pos_key = None
            for k in cand_pos_keys:
                if k in self.memory_store and isinstance(self.memory_store[k], pd.DataFrame):
                    # Check if empty and log it
                    if self.memory_store[k].empty:
                        self.stats_log.insert(tk.END, f"⚠️ Found {k} but it's EMPTY (0 rows)\n")
                    else:
                        pos_df = self.memory_store[k]
                        chosen_pos_key = k
                        break
            cand_neg_keys = ['neg_id_df','clean_neg_id_df','neg_enhanced_df','neg_result']
            chosen_neg_key = None
            for k in cand_neg_keys:
                if k in self.memory_store and isinstance(self.memory_store[k], pd.DataFrame):
                    # Check if empty and log it
                    if self.memory_store[k].empty:
                        self.stats_log.insert(tk.END, f"⚠️ Found {k} but it's EMPTY (0 rows)\n")
                    else:
                        neg_df = self.memory_store[k]
                        chosen_neg_key = k
                        break
        # If missing, try loading annotated excel (Pos_id/Neg_id or Positive/Negative sheets) if available
        if (pos_df is None or neg_df is None):
            # Accept multiple attribute names for backward compatibility
            ann_path = None
            for attr_name in ['annotated_metabolites_excel_path','id_annotated_excel_path','annotated_ids_excel_path']:
                if hasattr(self, attr_name):
                    candidate = getattr(self, attr_name)
                    if candidate:
                        ann_path = candidate
                        break
            if ann_path and os.path.exists(ann_path):
                try:
                    xl = pd.ExcelFile(ann_path)
                    # Prefer Pos_id / Neg_id sheets when present
                    sheet_pos = None
                    sheet_neg = None
                    for cand in ['Pos_id','Positive','Pos','POS']:
                        if cand in xl.sheet_names:
                            sheet_pos = cand; break
                    for cand in ['Neg_id','Negative','Neg','NEG']:
                        if cand in xl.sheet_names:
                            sheet_neg = cand; break
                    if pos_df is None and sheet_pos:
                        pos_df = xl.parse(sheet_pos)
                        self.memory_store['pos_id_df'] = pos_df
                        self.stats_log.insert(tk.END, f'Loaded {sheet_pos} sheet: {len(pos_df)} rows.\n')
                    if neg_df is None and sheet_neg:
                        neg_df = xl.parse(sheet_neg)
                        self.memory_store['neg_id_df'] = neg_df
                        self.stats_log.insert(tk.END, f'Loaded {sheet_neg} sheet: {len(neg_df)} rows.\n')
                except Exception as e:
                    self.stats_log.insert(tk.END, f'Could not load Positive/Negative sheets from annotated file: {e}\n')
            else:
                if ann_path is None:
                    self.stats_log.insert(tk.END, 'No annotated Excel path attribute found for auto-load.\n')
                else:
                    self.stats_log.insert(tk.END, f'Annotated Excel path not found: {ann_path}\n')
        if pos_df is not None:
            if 'chosen_pos_key' in locals() and chosen_pos_key:
                self.stats_log.insert(tk.END, f'Positive dataset available ({chosen_pos_key}): {len(pos_df)} rows, {len(pos_df.columns)} columns.\n')
            else:
                self.stats_log.insert(tk.END, f'Positive dataset available: {len(pos_df)} rows, {len(pos_df.columns)} columns.\n')
        else:
            self.stats_log.insert(tk.END, 'Positive dataset: None\n')
        if neg_df is not None:
            if 'chosen_neg_key' in locals() and chosen_neg_key:
                self.stats_log.insert(tk.END, f'Negative dataset available ({chosen_neg_key}): {len(neg_df)} rows, {len(neg_df.columns)} columns.\n')
            else:
                self.stats_log.insert(tk.END, f'Negative dataset available: {len(neg_df)} rows, {len(neg_df.columns)} columns.\n')
        else:
            self.stats_log.insert(tk.END, 'Negative dataset: None\n')
        # Retry logic if nothing loaded (single retry to allow asynchronous population of memory_store)
        if pos_df is None and neg_df is None:
            retry_count = getattr(self, '_stats_auto_retry_count', 0)
            if retry_count < 1:
                self._stats_auto_retry_count = retry_count + 1
                self.stats_log.insert(tk.END, '⚠️ No polarity datasets found yet. Retrying auto-load in 1.5s...\n')
                self.stats_log.see(tk.END)
                try:
                    self.root.after(1500, self.auto_prepare_statistics_tab)
                    return
                except Exception:
                    pass
        # Auto-detect sample columns and populate assignment interface
        mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
        try:
             # Detect columns based on mode
            if mode == 'lipid':
                # For lipid data, use robust feature detection
                feature_cols = []
                sample_cols = []
                for col in df.columns:
                    # Use normalized matching for lipid features
                    if self._is_lipid_feature_col(col):
                        feature_cols.append(col)
                    elif pd.api.types.is_numeric_dtype(df[col]):
                        sample_cols.append(col)
            else:
                # For metabolite data, use existing detection function
                from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
                union_sample_cols = []
                seen = set()
                for label, df in [('Positive', pos_df), ('Negative', neg_df)]:
                    if df is None:
                        continue
                    # Treat ID columns as feature columns if present so detection does not misclassify them
                    feature_cols, sample_cols = detect_feature_and_sample_columns(df)
                    id_like = [c for c in ['LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID','ChEBI_ID','CAS','SMILES','InChI','InChIKey','IUPAC_Name','Super_Class','Class','Sub_Class','Endogenous_Source'] if c in df.columns]
                    # Ensure they are counted among features (avoid inflating sample count)
                    feature_cols = list(dict.fromkeys(list(feature_cols) + id_like))
                    self.stats_log.insert(tk.END, f'{label}: detected {len(feature_cols)} feature cols (incl. ID/meta), {len(sample_cols)} sample cols.\n')
                    for c in sample_cols:
                        if c not in seen:
                            seen.add(c)
                            union_sample_cols.append(c)
            if union_sample_cols:
                # Populate listbox & group assignment interface
                self.sample_cols_list.delete(0, tk.END)
                for c in union_sample_cols:
                    self.sample_cols_list.insert(tk.END, c)
                
                # Populate the group assignment interface
                self.populate_sample_assignments(union_sample_cols)
                
                self.stats_log.insert(tk.END, f'Union sample columns loaded: {len(union_sample_cols)} columns.\n')
                self.stats_log.insert(tk.END, 'Group assignments initialized. Modify groups in the Group Management panel as needed.\n')
            else:
                self.stats_log.insert(tk.END, 'No sample columns detected automatically. Use Detect button or load Excel.\n')
        except Exception as e:
            self.stats_log.insert(tk.END, f'Auto-detection error: {e}\n')
        self.stats_log.insert(tk.END, 'Ready for normalization. Configure groups and click "Normalization & Test Normality".\n')
        self.stats_log.see(tk.END)
    
    def _auto_load_lipid_data(self):
        """Auto-load lipid data from lipid cleaning."""
        import pandas as pd
        # Summarize memory store content
        if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
            try:
                self.stats_log.insert(tk.END, f"Memory store keys: {list(self.memory_store.keys())}\n")
            except Exception:
                pass
        
        pos_lipid_df = None
        neg_lipid_df = None
        pos_class_df = None
        neg_class_df = None
        
        if hasattr(self, 'memory_store'):
            # Try multiple naming conventions
            for key in ['pos_lipid_df', 'lipid_pos_df']:
                if key in self.memory_store:
                    pos_lipid_df = self.memory_store.get(key)
                    break
            
            for key in ['neg_lipid_df', 'lipid_neg_df']:
                if key in self.memory_store:
                    neg_lipid_df = self.memory_store.get(key)
                    break
            
            for key in ['pos_lipid_class_df', 'lipid_pos_class_df']:
                if key in self.memory_store:
                    pos_class_df = self.memory_store.get(key)
                    break
            
            for key in ['neg_lipid_class_df', 'lipid_neg_class_df']:
                if key in self.memory_store:
                    neg_class_df = self.memory_store.get(key)
                    break
        
        # Log what was found
        if pos_lipid_df is not None:
            self.stats_log.insert(tk.END, f'Positive Lipids: {len(pos_lipid_df)} rows, {len(pos_lipid_df.columns)} columns.\n')
        else:
            self.stats_log.insert(tk.END, 'Positive Lipids: None\n')
        
        if neg_lipid_df is not None:
            self.stats_log.insert(tk.END, f'Negative Lipids: {len(neg_lipid_df)} rows, {len(neg_lipid_df.columns)} columns.\n')
        else:
            self.stats_log.insert(tk.END, 'Negative Lipids: None\n')
        
        if pos_class_df is not None:
            self.stats_log.insert(tk.END, f'Positive Lipid Class: {len(pos_class_df)} rows\n')
        if neg_class_df is not None:
            self.stats_log.insert(tk.END, f'Negative Lipid Class: {len(neg_class_df)} rows\n')
        
        # Auto-detect sample columns from lipid data
        union_sample_cols = []
        seen = set()
        
        for label, df in [('Positive Lipids', pos_lipid_df), ('Negative Lipids', neg_lipid_df)]:
            if df is None:
                continue
            sample_cols = []
            for col in df.columns:
                # Use robust lipid feature detection
                if self._is_lipid_feature_col(col):
                    continue
                # Numeric columns that are not features are sample intensity columns
                if pd.api.types.is_numeric_dtype(df[col]):
                    sample_cols.append(col)
            self.stats_log.insert(tk.END, f'{label}: detected {len(sample_cols)} sample cols.\n')
            for c in sample_cols:
                if c not in seen:
                    seen.add(c)
                    union_sample_cols.append(c)
        
        # Also check class DataFrames
        for label, df in [('Positive Class', pos_class_df), ('Negative Class', neg_class_df)]:
            if df is None:
                continue
            sample_cols = []
            for col in df.columns:
                if col is None:
                    continue
                # Use normalized matching for 'Class' column
                if self._normalize_col(col) == 'class':
                    continue
                try:
                    if self._is_lipid_feature_col(col):
                        continue
                except Exception:
                    pass
                if pd.api.types.is_numeric_dtype(df[col]):
                    sample_cols.append(col)
            self.stats_log.insert(tk.END, f'{label}: detected {len(sample_cols)} sample cols.\n')
            for c in sample_cols:
                if c not in seen:
                    seen.add(c)
                    union_sample_cols.append(c)
        
        if union_sample_cols:
            # Populate listbox & group assignment interface
            self.sample_cols_list.delete(0, tk.END)
            for c in union_sample_cols:
                self.sample_cols_list.insert(tk.END, c)
            
            # Populate the group assignment interface
            self.populate_sample_assignments(union_sample_cols)
            
            self.stats_log.insert(tk.END, f'Union sample columns loaded: {len(union_sample_cols)} columns.\n')
            self.stats_log.insert(tk.END, 'Group assignments initialized. Modify groups in the Group Management panel as needed.\n')
        else:
            self.stats_log.insert(tk.END, 'No sample columns detected automatically.\n')
        
        self.stats_log.insert(tk.END, 'Ready for normalization. Configure groups and click "Normalization & Test Normality".\n')
        self.stats_log.see(tk.END)
    
    def load_lipid_data_from_memory(self):
        """Public method to load lipid data from memory - called from Data Cleaning tab"""
        # Switch to lipid mode first
        if hasattr(self, 'data_mode') and self.data_mode.get() != 'lipid':
            self.data_mode.set('lipid')
            self.on_mode_change()  # Trigger mode change to update UI
        
        # Clear existing log
        if hasattr(self, 'stats_log'):
            self.stats_log.delete(1.0, tk.END)
            self.stats_log.insert(tk.END, "📊 Loading lipid data from memory...\n")
        
        # Call the internal auto-load method
        self._auto_load_lipid_data()
        
        # Log success
        if hasattr(self, 'stats_log'):
            self.stats_log.insert(tk.END, "\n✅ Lipid data loaded successfully!\n")
            self.stats_log.see(tk.END)
    
    def _save_nonparametric_twoway_results(self, results, method):
        """Save non-parametric two-way ANOVA results to Excel file."""
        try:
            import pandas as pd
            from datetime import datetime
            import os
            
            # Get output directory
            output_dir = self.stats_results_folder.get() if hasattr(self, 'stats_results_folder') else os.getcwd()
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            
            # Create timestamped filename
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f'NonParametric_TwoWay_{method.upper()}_{timestamp}.xlsx'
            filepath = os.path.join(output_dir, filename)
            
            # Write results to Excel
            with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
                # Write summary sheet
                summary_df = results.get('summary')
                if summary_df is not None and not summary_df.empty:
                    summary_df.to_excel(writer, sheet_name='Summary', index=False)
                    self._thread_safe_log(f"   ✓ Wrote Summary sheet ({len(summary_df)} metabolites)\n")
                
                # Write post-hoc results
                posthoc_results = results.get('posthoc', {})
                
                # Compile Factor A post-hoc results
                if posthoc_results:
                    factor_a_posthoc = []
                    factor_b_posthoc = []
                    interaction_posthoc = []
                    
                    for metabolite_id, posthoc_data in posthoc_results.items():
                        # Factor A
                        if 'posthoc_A' in posthoc_data and posthoc_data['posthoc_A'] is not None:
                            df_a = posthoc_data['posthoc_A'].copy()
                            df_a.insert(0, 'metabolite', metabolite_id)
                            factor_a_posthoc.append(df_a)
                        
                        # Factor B
                        if 'posthoc_B' in posthoc_data and posthoc_data['posthoc_B'] is not None:
                            df_b = posthoc_data['posthoc_B'].copy()
                            df_b.insert(0, 'metabolite', metabolite_id)
                            factor_b_posthoc.append(df_b)
                        
                        # Interaction
                        if 'posthoc_AB' in posthoc_data and posthoc_data['posthoc_AB'] is not None:
                            df_ab = posthoc_data['posthoc_AB'].copy()
                            df_ab.insert(0, 'metabolite', metabolite_id)
                            interaction_posthoc.append(df_ab)
                    
                    # Write Factor A post-hoc
                    if factor_a_posthoc:
                        combined_a = pd.concat(factor_a_posthoc, ignore_index=True)
                        combined_a.to_excel(writer, sheet_name='Posthoc_FactorA', index=False)
                        self._thread_safe_log(f"   ✓ Wrote Posthoc_FactorA sheet ({len(combined_a)} comparisons)\n")
                    
                    # Write Factor B post-hoc
                    if factor_b_posthoc:
                        combined_b = pd.concat(factor_b_posthoc, ignore_index=True)
                        combined_b.to_excel(writer, sheet_name='Posthoc_FactorB', index=False)
                        self._thread_safe_log(f"   ✓ Wrote Posthoc_FactorB sheet ({len(combined_b)} comparisons)\n")
                    
                    # Write Interaction post-hoc
                    if interaction_posthoc:
                        combined_ab = pd.concat(interaction_posthoc, ignore_index=True)
                        combined_ab.to_excel(writer, sheet_name='Posthoc_Interaction', index=False)
                        self._thread_safe_log(f"   ✓ Wrote Posthoc_Interaction sheet ({len(combined_ab)} comparisons)\n")
            
            self._thread_safe_log(f"\n💾 Results saved to:\n   {filepath}\n")
            
        except Exception as e:
            self._thread_safe_log(f"⚠️ Failed to save results: {e}\n")
            import traceback
            self._thread_safe_log(f"{traceback.format_exc()}\n")


    def open_covariate_dialog(self):
        """Open the covariate adjustment console in a separate dialog window"""
        # First ensure groups are configured before opening the covariate dialog
        def _proceed_to_open_dialog():
            self._open_covariate_dialog_core()
        
        if not self._ensure_groups_ready(after_config_callback=_proceed_to_open_dialog, action_label='covariate adjustment'):
            return
        
        # Groups already configured, open dialog immediately
        self._open_covariate_dialog_core()
    
    def _open_covariate_dialog_core(self):
        """Open the covariate adjustment console (after groups are confirmed)"""
        try:
            from gui.tabs.covariate_adjustment_section import CovariateAdjustmentDialog
            
            # Data provider function
            def _provide_data_for_covariates():
                """Provide current data to covariate adjustment section"""
                try:
                    # Determine current mode
                    mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
                    
                    # 🔓 BACKDOOR MODE CHECK: Use preprocessed sample columns if in backdoor mode
                    is_backdoor = self.memory_store.get('is_preprocessed_backdoor', False) if hasattr(self, 'memory_store') else False
                    
                    # Get current dataframes - check normalized first, then raw
                    # Use the EXACT same logic as statistics pipeline
                    df_pos = getattr(self, 'normalized_positive_df', None)
                    if df_pos is None:
                        df_pos = getattr(self, 'pos_df', None)
                    
                    df_neg = getattr(self, 'normalized_negative_df', None)
                    if df_neg is None:
                        df_neg = getattr(self, 'neg_df', None)
                    
                    # Use combined normalized if available (preferred), otherwise use positive if available, otherwise negative
                    df_intensities = getattr(self, 'normalized_combined_df', None)
                    if df_intensities is None:
                        df_intensities = df_pos if df_pos is not None else df_neg
                    
                    if df_intensities is None:
                        self._log_stats('🔍 COVARIATE DEBUG: No df_intensities available (import data in Step 2 first)')
                        return {}
                    
                    # Get sample columns and group map from sample_group_vars
                    # Configure Groups is done before covariate analysis and shows the correct 
                    # (cleaned) column names from the normalized data, so use those directly
                    sample_cols = []
                    group_map = {}
                    
                    if is_backdoor:
                        # 🔓 BACKDOOR MODE: Use preprocessed sample columns
                        self._log_stats('🔓 COVARIATE: Using backdoor mode sample columns\n')
                        preprocessed_cols = self.memory_store.get('preprocessed_sample_cols', [])
                        sample_cols = [col for col in preprocessed_cols if col in df_intensities.columns]
                        # Build group map from sample_group_vars
                        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                            for col, var in self.sample_group_vars.items():
                                group_label = var.get()
                                if group_label and col in sample_cols:
                                    group_map[col] = group_label
                    else:
                        # Normal mode - Use sample_group_vars directly since Configure Groups
                        # already shows the correct column names from the normalized dataframe
                        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
                            for col, var in self.sample_group_vars.items():
                                group_label = var.get()
                                if group_label:
                                    # Verify the column exists in the dataframe
                                    if col in df_intensities.columns:
                                        sample_cols.append(col)
                                        group_map[col] = group_label
                            
                            self._log_stats(f'🔍 COVARIATE DEBUG: Got {len(sample_cols)} sample columns from Configure Groups')
                    
                    if not sample_cols:
                        self._log_stats('🔍 COVARIATE DEBUG: No sample_cols available (verify columns in Step 2)')
                        return {}
                    
                    self._log_stats(f'🔍 COVARIATE DEBUG: sample_cols={len(sample_cols)}, group_map={len(group_map)} entries')
                    if len(group_map) < len(sample_cols):
                        missing = len(sample_cols) - len(group_map)
                        self._log_stats(f'⚠️ COVARIATE: {missing} samples missing group assignments. Use Step 3 "Auto-Assign Groups".')
                    
                    # Get metabolite ID column from verified assignments
                    metabolite_id_col = None
                    
                    if is_backdoor:
                        # 🔓 BACKDOOR MODE: Get from preprocessed verified assignments
                        verified_assignments = self.memory_store.get('preprocessed_verified_assignments', {})
                        if verified_assignments:
                            feature_id_key = 'LipidID' if mode == 'lipid' else 'Feature ID'
                            metabolite_id_col = verified_assignments.get(feature_id_key)
                            
                            # Debug: Log dataframe columns
                            self._log_stats(f'🔍 COVARIATE DEBUG (BACKDOOR): df_intensities columns: {list(df_intensities.columns)[:10]}{"..." if len(df_intensities.columns) > 10 else ""}\n')
                            self._log_stats(f'🔍 COVARIATE DEBUG (BACKDOOR): Looking for verified column: "{metabolite_id_col}"\n')
                            
                            # Verify the column actually exists in the dataframe
                            if metabolite_id_col and metabolite_id_col not in df_intensities.columns:
                                self._log_stats(f'⚠️ WARNING: Verified column "{metabolite_id_col}" not found in dataframe\n')
                                self._log_stats(f'⚠️ Available non-sample columns: {[col for col in df_intensities.columns if col not in sample_cols]}\n')
                                metabolite_id_col = None
                        
                        if not metabolite_id_col:
                            self._log_stats(f'❌ ERROR: No Feature ID column verified in backdoor mode. Cannot proceed.\n')
                            return {}
                    else:
                        # NORMAL MODE: Get from regular verified assignments
                        if mode == 'lipid':
                            # Try lipid assignments first - use 'LipidID' key for lipid mode
                            metabolite_id_col = getattr(self, 'verified_pos_lipid_assignments', {}).get('LipidID')
                            if not metabolite_id_col:
                                metabolite_id_col = getattr(self, 'verified_neg_lipid_assignments', {}).get('LipidID')
                        else:
                            # Try metabolite assignments - use 'Feature ID' or 'Name' key for metabolite mode
                            metabolite_id_col = getattr(self, 'verified_pos_assignments', {}).get('Feature ID')
                            if not metabolite_id_col:
                                metabolite_id_col = getattr(self, 'verified_pos_assignments', {}).get('Name')
                            if not metabolite_id_col:
                                metabolite_id_col = getattr(self, 'verified_neg_assignments', {}).get('Feature ID')
                            if not metabolite_id_col:
                                metabolite_id_col = getattr(self, 'verified_neg_assignments', {}).get('Name')
                        
                        # Verify the column exists in dataframe
                        if metabolite_id_col and metabolite_id_col not in df_intensities.columns:
                            self._log_stats(f'⚠️ WARNING: Verified column "{metabolite_id_col}" not found in combined dataframe\n')
                            # Try to find it in the original dataframes
                            if not metabolite_id_col:
                                id_col_name = 'LipidID' if mode == 'lipid' else 'Feature ID/Name'
                                self._log_stats(f'❌ ERROR: No {id_col_name} column found. Please verify columns in Step 2.\n')
                                return {}
                        
                        if not metabolite_id_col:
                            id_col_name = 'LipidID' if mode == 'lipid' else 'Feature ID/Name'
                            self._log_stats(f'❌ ERROR: No {id_col_name} column verified. Please complete Step 2: Verify Columns.\n')
                            return {}
                    
                    self._log_stats(f'✅ COVARIATE: Using Metabolite ID column: {metabolite_id_col}\n')
                    
                    # Get lipid class data if in lipid mode
                    df_class = None
                    class_id_col = None
                    if mode == 'lipid':
                        df_class = getattr(self, 'normalized_combined_class_df', None)
                        if df_class is not None:
                            # Find class ID column - typically 'Class' or first non-sample column
                            non_sample_cols = [c for c in df_class.columns if c not in sample_cols]
                            for candidate in ['Class', 'LipidClass', 'Lipid_Class', 'class']:
                                if candidate in non_sample_cols:
                                    class_id_col = candidate
                                    break
                            if not class_id_col and non_sample_cols:
                                class_id_col = non_sample_cols[0]
                            self._log_stats(f'✅ COVARIATE: Found lipid class data ({len(df_class)} classes), ID column: {class_id_col}\n')
                    
                    return {
                        'df_intensities': df_intensities,
                        'sample_cols': sample_cols,
                        'group_map': group_map,
                        'metabolite_id_col': metabolite_id_col,
                        'mode': mode,
                        'df_class': df_class,
                        'class_id_col': class_id_col
                    }
                except Exception as e:
                    self._log_stats(f'❌ Error providing data for covariates: {e}')
                    import traceback
                    self._log_stats(f'   Traceback: {traceback.format_exc()}')
                    return {}
            
            # Open dialog
            dialog = CovariateAdjustmentDialog(
                parent=self.root,
                log_callback=self._log_stats,
                data_provider=_provide_data_for_covariates,
                parent_tab=self  # Pass the Statistics tab instance
            )
            
        except ImportError as e:
            messagebox.showerror(
                'Missing Dependency',
                f'Covariate adjustment requires statsmodels.\n\n'
                f'Install with: pip install statsmodels\n\n'
                f'Error: {e}'
            )
        except Exception as e:
            messagebox.showerror('Error', f'Could not open covariate console:\n{str(e)}')
            logger.error(f'Error opening covariate dialog: {e}', exc_info=True)






