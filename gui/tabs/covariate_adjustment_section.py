"""
Covariate Adjustment UI Section for Statistics Tab

This module provides a GUI section for covariate adjustment in metabolite analysis.
Can be embedded into the statistics tab.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import os
from typing import Optional, List, Dict, Callable
import logging
import threading
import contextvars

logger = logging.getLogger(__name__)


class CovariateAdjustmentSection:
    """
    UI section for covariate adjustment configuration and execution.
    
    Features:
    - Load covariates from separate file or detect from main data
    - Preview covariate data
    - Select covariates to include in model
    - Run covariate-adjusted analysis
    - Export results
    """
    
    def __init__(
        self,
        parent_frame: tk.Frame,
        log_callback: Optional[Callable[[str], None]] = None,
        data_provider: Optional[Callable[[], Dict]] = None,
        parent_tab = None
    ):
        """
        Initialize covariate adjustment section.
        
        Parameters
        ----------
        parent_frame : tk.Frame
            Parent frame to place widgets in
        log_callback : Optional[Callable]
            Function to call for logging messages (signature: log_callback(message))
        data_provider : Optional[Callable]
            Function that returns dict with keys:
            - 'df_intensities': metabolite intensity dataframe
            - 'sample_cols': list of sample column names
            - 'group_map': dict mapping sample -> group
            - 'metabolite_id_col': name of metabolite ID column
        parent_tab : Optional
            Parent statistics tab for accessing attributes like last_statistics_export_path
        """
        self.parent_frame = parent_frame
        self.log_callback = log_callback or print
        self.data_provider = data_provider
        self.parent_tab = parent_tab
        
        # State variables
        self.covariate_file_path: Optional[str] = None
        self.df_covariates: Optional[pd.DataFrame] = None
        self.detected_continuous_covars: List[str] = []
        self.detected_categorical_covars: List[str] = []
        self.selected_covariates: List[str] = []
        self.covariate_checkboxes: Dict[str, tk.BooleanVar] = {}
        self.mapped_sample_id_col: Optional[str] = None
        self.mapped_covariate_cols: List[str] = []
        self._analysis_run_id: int = 0
        
        # Create UI
        self.create_ui()
    
    def _log(self, message: str):
        """Log a message using the callback"""
        self.log_callback(message)

    def _widget_exists(self, widget) -> bool:
        """Safely check whether a Tk widget still exists."""
        try:
            return widget is not None and bool(widget.winfo_exists())
        except Exception:
            return False

    def _safe_widget_config(self, widget, **kwargs):
        """Configure a widget only if it still exists."""
        if self._widget_exists(widget):
            try:
                widget.config(**kwargs)
            except Exception:
                pass

    def _ui_after(self, delay_ms: int, callback):
        """Schedule UI callbacks only when the owning frame is alive."""
        if self._widget_exists(self.parent_frame):
            try:
                self.parent_frame.after(delay_ms, callback)
            except Exception:
                pass
    
    def _auto_load_to_visualization(self):
        """Auto-load covariate results to Visualization tab without switching."""
        try:
            self._log('🔍 VIZ DEBUG: Starting auto-load to visualization...')
            
            if not self.parent_tab:
                self._log('❌ VIZ DEBUG: No parent tab found')
                return
            
            self._log(f'🔍 VIZ DEBUG: Parent tab exists: {type(self.parent_tab).__name__}')
            
            # Get the sample-to-group mapping from parent tab
            group_map = None
            if hasattr(self.parent_tab, 'sample_group_vars') and self.parent_tab.sample_group_vars:
                group_map = {}
                for sample_col, group_var in self.parent_tab.sample_group_vars.items():
                    group_val = group_var.get()
                    if group_val:
                        group_map[sample_col] = group_val
                self._log(f'🔍 VIZ DEBUG: Built group_map with {len(group_map)} samples')
            else:
                self._log('⚠️ VIZ DEBUG: No sample_group_vars available')
            
            # Store sample-to-group mapping in memory_store for visualization
            if group_map and hasattr(self.parent_tab, 'memory_store'):
                self.parent_tab.memory_store['sample_to_group'] = group_map
                self._log(f'✅ VIZ DEBUG: Stored group mappings: {len(group_map)} samples across {len(set(group_map.values()))} groups')
            else:
                self._log('⚠️ VIZ DEBUG: Could not store group mappings in memory_store')
            
            # Check what's in memory_store
            if hasattr(self.parent_tab, 'memory_store'):
                self._log(f'🔍 VIZ DEBUG: memory_store keys: {list(self.parent_tab.memory_store.keys())}')
            
            # Check what's in statistical_test_results
            if hasattr(self.parent_tab, 'statistical_test_results'):
                self._log(f'🔍 VIZ DEBUG: statistical_test_results keys: {list(self.parent_tab.statistical_test_results.keys())}')
            else:
                self._log('⚠️ VIZ DEBUG: No statistical_test_results attribute on parent tab')
            
            # Notify Visualization tab that results are ready
            if hasattr(self.parent_tab, 'notify_data_ready'):
                self._log('🔍 VIZ DEBUG: Calling notify_data_ready...')
                self.parent_tab.notify_data_ready("📊 Visualization", "covariate_results")
            else:
                self._log('⚠️ VIZ DEBUG: Parent tab missing notify_data_ready method')
            
            # Get Visualization tab and trigger its data loading
            viz_tab = None
            if hasattr(self.parent_tab, 'get_tab_by_name'):
                self._log('🔍 VIZ DEBUG: Getting visualization tab by name...')
                viz_tab = self.parent_tab.get_tab_by_name("📊 Visualization")
                if viz_tab:
                    self._log(f'✅ VIZ DEBUG: Found viz_tab: {type(viz_tab).__name__}')
                else:
                    self._log('❌ VIZ DEBUG: get_tab_by_name returned None')
            else:
                self._log('⚠️ VIZ DEBUG: Parent tab missing get_tab_by_name method')
            
            if viz_tab:
                # Call visualization tab's data update method if available
                if hasattr(viz_tab, 'update_viz_data_status'):
                    self._log('🔍 VIZ DEBUG: Calling update_viz_data_status on viz_tab...')
                    viz_tab.update_viz_data_status()
                    self._log('✅ Covariate results auto-loaded to Visualization tab')
                else:
                    self._log('⚠️ VIZ DEBUG: viz_tab missing update_viz_data_status method')
            else:
                self._log('❌ VIZ DEBUG: Visualization tab not found for auto-load')
        except Exception as e:
            self._log(f'❌ VIZ DEBUG: Exception in auto-load: {e}')
            import traceback
            self._log(f'❌ VIZ DEBUG: Traceback:\n{traceback.format_exc()}')
    
    def create_ui(self):
        """Create the covariate adjustment UI section"""
        # Main container - use parent_frame directly for dialog mode
        self.main_frame = self.parent_frame
        
        # Step 1: Load covariates
        step1_frame = tk.LabelFrame(
            self.main_frame, text='Step 1: Load Covariate Data',
            bg='#f0f0f0', font=('Arial', 9, 'bold')
        )
        step1_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        btn_frame = tk.Frame(step1_frame, bg='#f0f0f0')
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Button(
            btn_frame, text='📂 Load Covariate File',
            command=self.load_covariate_file,
            bg='#2980b9', fg='white', font=('Arial', 9, 'bold'),
            relief='raised', bd=2, pady=3
        ).pack(side='left', fill='x', expand=True, padx=2)
        
        tk.Button(
            btn_frame, text='🔍 Detect from Main Data',
            command=self.detect_covariates_from_main_data,
            bg='#16a085', fg='white', font=('Arial', 9, 'bold'),
            relief='raised', bd=2, pady=3
        ).pack(side='left', fill='x', expand=True, padx=2)

        tk.Button(
            btn_frame, text='🧾 Preview Columns',
            command=self.preview_covariate_columns,
            bg='#7f8c8d', fg='white', font=('Arial', 9, 'bold'),
            relief='raised', bd=2, pady=3
        ).pack(side='left', fill='x', expand=True, padx=2)
        
        # File status label
        self.covariate_file_label = tk.Label(
            step1_frame, text='No covariate data loaded',
            bg='#f0f0f0', font=('Arial', 8), fg='#666'
        )
        self.covariate_file_label.pack(padx=5, pady=(0, 5))
        
        # Step 2: Select covariates
        step2_frame = tk.LabelFrame(
            self.main_frame, text='Step 2: Select Covariates',
            bg='#f0f0f0', font=('Arial', 9, 'bold')
        )
        step2_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        # Scrollable covariate list
        list_container = tk.Frame(step2_frame, bg='#f0f0f0')
        list_container.pack(fill='both', expand=True, padx=5, pady=5)
        
        canvas = tk.Canvas(list_container, bg='white', height=120, highlightthickness=1)
        scrollbar = ttk.Scrollbar(list_container, orient='vertical', command=canvas.yview)
        self.covariate_list_frame = tk.Frame(canvas, bg='white')
        
        self.covariate_list_frame.bind(
            '<Configure>',
            lambda e: canvas.configure(scrollregion=canvas.bbox('all'))
        )
        
        canvas.create_window((0, 0), window=self.covariate_list_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        scrollbar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        
        # Initially empty
        self.empty_label = tk.Label(
            self.covariate_list_frame,
            text='Load covariates to select variables',
            bg='white', fg='#999', font=('Arial', 9, 'italic')
        )
        self.empty_label.pack(pady=20)
        
        # Select/deselect buttons
        select_btn_frame = tk.Frame(step2_frame, bg='#f0f0f0')
        select_btn_frame.pack(fill='x', padx=5, pady=(0, 5))
        
        tk.Button(
            select_btn_frame, text='Select All',
            command=self.select_all_covariates,
            bg='#27ae60', fg='white', font=('Arial', 8),
            relief='raised', bd=1, pady=2
        ).pack(side='left', fill='x', expand=True, padx=2)
        
        tk.Button(
            select_btn_frame, text='Deselect All',
            command=self.deselect_all_covariates,
            bg='#e74c3c', fg='white', font=('Arial', 8),
            relief='raised', bd=1, pady=2
        ).pack(side='left', fill='x', expand=True, padx=2)
        
        # Step 3: Run analysis
        step3_frame = tk.LabelFrame(
            self.main_frame, text='Step 3: Run Analysis',
            bg='#f0f0f0', font=('Arial', 9, 'bold')
        )
        step3_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        # Settings
        settings_frame = tk.Frame(step3_frame, bg='#f0f0f0')
        settings_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Label(
            settings_frame, text='Reference Group:',
            bg='#f0f0f0', font=('Arial', 9)
        ).grid(row=0, column=0, sticky='w', padx=5, pady=2)
        
        self.ref_group_var = tk.StringVar(value='')
        self.ref_group_combo = ttk.Combobox(
            settings_frame, textvariable=self.ref_group_var,
            state='readonly', width=15
        )
        self.ref_group_combo.grid(row=0, column=1, sticky='w', padx=5, pady=2)
        
        tk.Label(
            settings_frame, text='(leave blank for alphabetical first)',
            bg='#f0f0f0', font=('Arial', 7), fg='#666'
        ).grid(row=0, column=2, sticky='w', padx=5, pady=2)
        
        # Method selection (Linear Model or Limma)
        tk.Label(
            settings_frame, text='Analysis Method:',
            bg='#f0f0f0', font=('Arial', 9)
        ).grid(row=1, column=0, sticky='w', padx=5, pady=2)
        
        self.analysis_method_var = tk.StringVar(value='Linear Model (OLS)')
        ttk.Combobox(
            settings_frame,
            values=['Linear Model (OLS)', 'Limma'],
            textvariable=self.analysis_method_var,
            state='readonly',
            width=15
        ).grid(row=1, column=1, sticky='w', padx=5, pady=2)
        
        tk.Label(
            settings_frame, text='(Empirical Bayes moderated statistics)',
            bg='#f0f0f0', font=('Arial', 7), fg='#666'
        ).grid(row=1, column=2, sticky='w', padx=5, pady=2)
        
        # Pairwise p-value adjustment method (same as statistics tab)
        tk.Label(
            settings_frame, text='Pairwise p-value adjustment:',
            bg='#f0f0f0', font=('Arial', 9)
        ).grid(row=2, column=0, sticky='w', padx=5, pady=2)
        
        # Get correction method from parent tab if available, otherwise default to BH
        if self.parent_tab and hasattr(self.parent_tab, 'pairwise_p_adjust_method'):
            self.pairwise_p_adjust_method = self.parent_tab.pairwise_p_adjust_method
        else:
            self.pairwise_p_adjust_method = tk.StringVar(value='BH')
        
        ttk.Combobox(
            settings_frame, 
            values=['BH','Bonferroni','Holm','Hochberg','BY','None'], 
            textvariable=self.pairwise_p_adjust_method, 
            state='readonly',
            width=15
        ).grid(row=2, column=1, sticky='w', padx=5, pady=2)
        
        # Alpha threshold
        tk.Label(
            settings_frame, text='Significance α:',
            bg='#f0f0f0', font=('Arial', 9)
        ).grid(row=3, column=0, sticky='w', padx=5, pady=2)
        
        self.alpha_var = tk.StringVar(value='0.05')
        tk.Entry(
            settings_frame, textvariable=self.alpha_var, width=8
        ).grid(row=3, column=1, sticky='w', padx=5, pady=2)
        
        # Run button
        self.run_covariate_btn = tk.Button(
            step3_frame, text='▶️ Run Analysis',
            command=self.run_covariate_analysis,
            bg='#8e44ad', fg='white', font=('Arial', 9, 'bold'),
            relief='raised', bd=2, pady=5
        )
        self.run_covariate_btn.pack(fill='x', padx=5, pady=(5, 5))
        
        # Export button
        self.export_covariate_btn = tk.Button(
            step3_frame, text='💾 Export Covariate Results',
            command=self.export_covariate_results,
            bg='#2c3e50', fg='white', font=('Arial', 9, 'bold'),
            relief='raised', bd=2, pady=3, state='disabled'
        )
        self.export_covariate_btn.pack(fill='x', padx=5, pady=(0, 5))
        
        # Results storage
        self.covariate_results = None
        self.covariate_class_results = None  # For lipid class results
    
    def load_covariate_file(self):
        """Load covariate data from a file"""
        file_path = filedialog.askopenfilename(
            title='Select Covariate File',
            filetypes=[
                ('Excel files', '*.xlsx *.xls'),
                ('CSV files', '*.csv'),
                ('Text files', '*.txt *.tsv'),
                ('All files', '*.*')
            ]
        )
        
        if not file_path:
            return
        
        try:
            raw_df = self._read_covariate_file_raw(file_path)

            sample_cols = []
            if self.data_provider:
                try:
                    sample_cols = self.data_provider().get('sample_cols', []) or []
                except Exception:
                    sample_cols = []

            mapping = self._show_covariate_mapping_dialog(raw_df, sample_cols)
            if mapping is None:
                self._log('Covariate load cancelled by user during column mapping.')
                return

            sample_id_col = mapping['sample_id_col']
            selected_covariates = mapping['covariate_cols']

            if sample_id_col not in raw_df.columns:
                raise ValueError(f"Selected sample ID column '{sample_id_col}' not found in file")

            df_cov = raw_df.copy()
            df_cov[sample_id_col] = df_cov[sample_id_col].astype(str).str.strip()

            # Drop blank sample IDs and collapse duplicates (first occurrence kept)
            df_cov = df_cov[df_cov[sample_id_col].str.len() > 0].copy()
            if df_cov[sample_id_col].duplicated().any():
                dup_count = int(df_cov[sample_id_col].duplicated().sum())
                self._log(f'⚠️ Found {dup_count} duplicate sample IDs in covariate file; keeping first occurrence per sample.')
                df_cov = df_cov.drop_duplicates(subset=[sample_id_col], keep='first')

            keep_covars = [c for c in selected_covariates if c in df_cov.columns]
            self.df_covariates = df_cov.set_index(sample_id_col)[keep_covars]
            self.covariate_file_path = file_path
            self.mapped_sample_id_col = sample_id_col
            self.mapped_covariate_cols = keep_covars
            
            # Update UI
            filename = os.path.basename(file_path)
            n_samples = len(self.df_covariates)
            n_covars = len(self.df_covariates.columns)
            
            self.covariate_file_label.config(
                text=f'✓ Loaded: {filename} ({n_samples} samples, {n_covars} variables)',
                fg='#27ae60'
            )
            
            # Detect covariate types
            self._detect_and_display_covariates()
            
            self._log(f'Loaded covariates from {filename}')
            self._log(f'🔗 Mapping used -> sample ID column: {sample_id_col}; selected covariates: {", ".join(keep_covars) if keep_covars else "None"}')
            cols_preview = ', '.join([str(c) for c in self.df_covariates.columns[:20]])
            if len(self.df_covariates.columns) > 20:
                cols_preview += ', ...'
            self._log(f'📋 Covariate columns ({n_covars}): {cols_preview}')
            idx_preview = ', '.join([str(i) for i in self.df_covariates.index[:5]])
            self._log(f'🧬 Sample ID preview from covariate file index: {idx_preview}')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load covariate file:\n{str(e)}')
            logger.error(f'Error loading covariate file: {e}', exc_info=True)

    def _read_covariate_file_raw(self, file_path: str) -> pd.DataFrame:
        """Read covariate file without forcing any index column."""
        if file_path.endswith(('.xlsx', '.xls')):
            return pd.read_excel(file_path)
        if file_path.endswith('.csv'):
            return pd.read_csv(file_path)
        if file_path.endswith(('.txt', '.tsv')):
            return pd.read_csv(file_path, sep='\t')
        raise ValueError(f'Unsupported file format: {file_path}')

    def _show_covariate_mapping_dialog(self, raw_df: pd.DataFrame, sample_cols: List[str]):
        """Popup for mapping sample ID and covariate columns from uploaded file."""
        if raw_df is None or raw_df.empty:
            messagebox.showerror('Error', 'Uploaded covariate file is empty.')
            return None

        columns = [str(c) for c in raw_df.columns]
        if not columns:
            messagebox.showerror('Error', 'No columns found in covariate file.')
            return None

        normalize = lambda x: str(x).strip().casefold()
        expected_norm = {normalize(s) for s in sample_cols} if sample_cols else set()

        # Heuristic default for sample ID column
        sample_candidates = [c for c in columns if any(k in normalize(c) for k in ['sample', 'sampleid', 'sample_id', 'id', 'subject'])]
        default_sample_col = sample_candidates[0] if sample_candidates else columns[0]

        result = {'confirmed': False, 'sample_id_col': None, 'covariate_cols': []}

        dialog = tk.Toplevel(self.parent_frame.winfo_toplevel())
        dialog.title('Covariate Column Mapping')
        dialog.geometry('760x560')
        dialog.minsize(700, 520)
        dialog.transient(self.parent_frame.winfo_toplevel())

        tk.Label(
            dialog,
            text='Map covariate file columns',
            font=('Arial', 12, 'bold')
        ).pack(anchor='w', padx=12, pady=(12, 4))

        tk.Label(
            dialog,
            text='1) Choose Sample ID column. 2) Select covariate columns (Age, Sex, BMI, PMI, etc.).',
            font=('Arial', 9),
            fg='#333'
        ).pack(anchor='w', padx=12, pady=(0, 8))

        top_frame = tk.Frame(dialog)
        top_frame.pack(fill='x', padx=12, pady=6)

        tk.Label(top_frame, text='Sample ID column:', font=('Arial', 9, 'bold')).pack(side='left')
        sample_col_var = tk.StringVar(value=default_sample_col)
        sample_col_combo = ttk.Combobox(top_frame, textvariable=sample_col_var, values=columns, state='readonly', width=40)
        sample_col_combo.pack(side='left', padx=8)

        status_var = tk.StringVar(value='')
        tk.Label(dialog, textvariable=status_var, font=('Arial', 9), fg='#1565c0').pack(anchor='w', padx=12, pady=(0, 6))

        middle_frame = tk.Frame(dialog)
        middle_frame.pack(fill='both', expand=True, padx=12, pady=6)

        tk.Label(middle_frame, text='Covariate columns to include in model:', font=('Arial', 9, 'bold')).pack(anchor='w')

        list_frame = tk.Frame(middle_frame)
        list_frame.pack(fill='both', expand=True, pady=6)

        covar_listbox = tk.Listbox(list_frame, selectmode='extended', exportselection=False)
        covar_scroll = ttk.Scrollbar(list_frame, orient='vertical', command=covar_listbox.yview)
        covar_listbox.configure(yscrollcommand=covar_scroll.set)
        covar_listbox.pack(side='left', fill='both', expand=True)
        covar_scroll.pack(side='right', fill='y')

        for col in columns:
            covar_listbox.insert('end', col)

        btns = tk.Frame(dialog)
        btns.pack(fill='x', padx=12, pady=(4, 2))

        def _select_all_covars():
            covar_listbox.selection_clear(0, 'end')
            for i, col in enumerate(columns):
                if col != sample_col_var.get():
                    covar_listbox.selection_set(i)

        def _clear_covars():
            covar_listbox.selection_clear(0, 'end')

        tk.Button(btns, text='Select All (except Sample ID)', command=_select_all_covars).pack(side='left')
        tk.Button(btns, text='Clear', command=_clear_covars).pack(side='left', padx=6)

        preview_text = tk.Text(dialog, height=10, wrap='word')
        preview_text.pack(fill='x', padx=12, pady=(8, 4))

        def _update_status(*_):
            sample_col = sample_col_var.get()
            if sample_col not in raw_df.columns:
                status_var.set('Select a valid sample ID column.')
                return

            sample_values = raw_df[sample_col].dropna().astype(str).str.strip()
            sample_norm = {normalize(v) for v in sample_values if str(v).strip()}
            matched = len(sample_norm & expected_norm) if expected_norm else 0

            if expected_norm:
                status_var.set(f'Matched samples with current analysis set: {matched}/{len(expected_norm)}')
            else:
                status_var.set(f'Loaded sample IDs in selected column: {len(sample_norm)} (no active sample set to compare)')

            preview_lines = []
            preview_lines.append(f'Sample ID column: {sample_col}')
            preview_lines.append(f'Sample values preview: {", ".join(sample_values.head(10).tolist())}')
            if expected_norm:
                missing = [s for s in sample_cols if normalize(s) not in sample_norm]
                if missing:
                    preview_lines.append(f'Missing from covariate file (first 10): {", ".join(missing[:10])}')
                else:
                    preview_lines.append('All analysis samples are present in covariate file.')

            preview_text.config(state='normal')
            preview_text.delete('1.0', 'end')
            preview_text.insert('1.0', '\n'.join(preview_lines))
            preview_text.config(state='disabled')

            # Keep sample column deselected from covariate list.
            for i, col in enumerate(columns):
                if col == sample_col:
                    covar_listbox.selection_clear(i)

        sample_col_combo.bind('<<ComboboxSelected>>', _update_status)
        _select_all_covars()
        _update_status()

        footer = tk.Frame(dialog)
        footer.pack(fill='x', padx=12, pady=10)

        def _on_ok():
            sample_col = sample_col_var.get()
            if sample_col not in columns:
                messagebox.showerror('Mapping Error', 'Please choose a valid Sample ID column.', parent=dialog)
                return

            selected_idxs = covar_listbox.curselection()
            selected_covars = [columns[i] for i in selected_idxs if columns[i] != sample_col]

            if not selected_covars:
                proceed = messagebox.askyesno(
                    'No Covariates Selected',
                    'No covariate columns were selected. Continue anyway?',
                    parent=dialog
                )
                if not proceed:
                    return

            if expected_norm:
                sample_values = raw_df[sample_col].dropna().astype(str).str.strip()
                sample_norm = {normalize(v) for v in sample_values if str(v).strip()}
                matched = len(sample_norm & expected_norm)
                if matched == 0:
                    messagebox.showerror(
                        'Mapping Error',
                        'Selected Sample ID column does not match any currently assigned samples.\n'
                        'Please choose the correct Sample ID column.',
                        parent=dialog
                    )
                    return

            result['confirmed'] = True
            result['sample_id_col'] = sample_col
            result['covariate_cols'] = selected_covars
            dialog.destroy()

        def _on_cancel():
            dialog.destroy()

        tk.Button(footer, text='Cancel', command=_on_cancel, width=10).pack(side='right')
        tk.Button(footer, text='Apply Mapping', command=_on_ok, width=14, bg='#27ae60', fg='white').pack(side='right', padx=8)

        dialog.grab_set()
        dialog.wait_window()

        if result['confirmed']:
            return {'sample_id_col': result['sample_id_col'], 'covariate_cols': result['covariate_cols']}
        return None
    
    def detect_covariates_from_main_data(self):
        """Detect covariates from the main metabolite data"""
        if not self.data_provider:
            messagebox.showwarning(
                'Warning',
                'No data provider configured. Load covariate file instead.'
            )
            return
        
        try:
            # Get data from provider
            data = self.data_provider()
            df_intensities = data.get('df_intensities')
            sample_cols = data.get('sample_cols', [])
            
            if df_intensities is None or not sample_cols:
                messagebox.showwarning(
                    'Warning',
                    'No metabolite data loaded. Import data first in Step 2.'
                )
                return
            
            # Import detection function
            from main_script.covariate_adjustment import detect_covariate_columns
            
            # Detect covariates
            continuous, categorical = detect_covariate_columns(
                df_intensities, sample_cols
            )
            
            if not continuous and not categorical:
                messagebox.showinfo(
                    'No Covariates Found',
                    'No covariate columns detected in the main data.\n'
                    'Load covariates from a separate file instead.'
                )
                return
            
            # Create covariate dataframe (transposed if needed)
            covar_cols = continuous + categorical
            self.df_covariates = df_intensities[covar_cols].T
            self.df_covariates.columns = sample_cols
            
            self.covariate_file_path = None  # Embedded, not from file
            
            # Update UI
            self.covariate_file_label.config(
                text=f'✓ Detected from main data ({len(covar_cols)} variables)',
                fg='#27ae60'
            )
            
            # Display covariates
            self._detect_and_display_covariates()
            
            self._log(f'Detected {len(continuous)} continuous and {len(categorical)} categorical covariates')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to detect covariates:\n{str(e)}')
            logger.error(f'Error detecting covariates: {e}', exc_info=True)
    
    def _detect_and_display_covariates(self):
        """Detect covariate types and display in UI"""
        if self.df_covariates is None:
            return
        
        from main_script.covariate_adjustment import detect_covariate_columns
        
        # The df_covariates structure after load_covariate_file:
        # - Index: Sample IDs
        # - Columns: Covariate names (Sex, Age, PMI, etc.)
        # So we detect on the columns directly, no need to transpose
        
        # Get mapped covariates when user explicitly selected them during upload.
        if self.mapped_covariate_cols:
            all_covar_cols = [c for c in self.mapped_covariate_cols if c in self.df_covariates.columns]
        else:
            all_covar_cols = list(self.df_covariates.columns)
        
        # Detect which are continuous vs categorical
        continuous = []
        categorical = []
        
        for col in all_covar_cols:
            series = self.df_covariates[col].dropna()
            if len(series) == 0:
                continue
            
            # Try to convert to numeric
            try:
                numeric_series = pd.to_numeric(series, errors='coerce')
                numeric_ratio = numeric_series.notna().sum() / len(series)
                
                if numeric_ratio > 0.8:  # Mostly numeric
                    unique_count = numeric_series.nunique()
                    if unique_count <= 10:
                        # Few unique values -> categorical
                        categorical.append(col)
                    else:
                        # Many unique values -> continuous
                        continuous.append(col)
                else:
                    # Not numeric -> categorical
                    unique_count = series.nunique()
                    if unique_count <= 20:
                        categorical.append(col)
            except:
                # If conversion fails, treat as categorical if unique count is low
                unique_count = series.nunique()
                if unique_count <= 20:
                    categorical.append(col)
        
        self.detected_continuous_covars = continuous
        self.detected_categorical_covars = categorical
        
        # Display in UI
        self._populate_covariate_list()
        
        # Enable run button if covariates available
        if continuous or categorical:
            self._safe_widget_config(self.run_covariate_btn, state='normal')
            
        # Try to populate reference group dropdown
        self._try_populate_reference_groups()
    
    def _populate_covariate_list(self):
        """Populate the covariate selection list"""
        # Clear existing
        for widget in self.covariate_list_frame.winfo_children():
            widget.destroy()
        
        self.covariate_checkboxes = {}
        
        if not self.detected_continuous_covars and not self.detected_categorical_covars:
            self.empty_label = tk.Label(
                self.covariate_list_frame,
                text='No covariates available',
                bg='white', fg='#999', font=('Arial', 9, 'italic')
            )
            self.empty_label.pack(pady=20)
            return
        
        # Header
        tk.Label(
            self.covariate_list_frame, text='Select covariates to include in model:',
            bg='white', font=('Arial', 9, 'bold')
        ).pack(anchor='w', padx=5, pady=(5, 2))
        
        # Continuous covariates
        if self.detected_continuous_covars:
            tk.Label(
                self.covariate_list_frame, text='Continuous:',
                bg='white', font=('Arial', 9, 'bold'), fg='#2980b9'
            ).pack(anchor='w', padx=10, pady=(5, 2))
            
            for covar in self.detected_continuous_covars:
                var = tk.BooleanVar(value=True)  # Selected by default
                self.covariate_checkboxes[covar] = var
                
                tk.Checkbutton(
                    self.covariate_list_frame, text=covar,
                    variable=var, bg='white', font=('Arial', 9)
                ).pack(anchor='w', padx=20)
        
        # Categorical covariates
        if self.detected_categorical_covars:
            tk.Label(
                self.covariate_list_frame, text='Categorical:',
                bg='white', font=('Arial', 9, 'bold'), fg='#16a085'
            ).pack(anchor='w', padx=10, pady=(5, 2))
            
            for covar in self.detected_categorical_covars:
                var = tk.BooleanVar(value=True)  # Selected by default
                self.covariate_checkboxes[covar] = var
                
                tk.Checkbutton(
                    self.covariate_list_frame, text=covar,
                    variable=var, bg='white', font=('Arial', 9)
                ).pack(anchor='w', padx=20)
    
    def select_all_covariates(self):
        """Select all covariates"""
        for var in self.covariate_checkboxes.values():
            var.set(True)
        self._log('Selected all covariates')
    
    def deselect_all_covariates(self):
        """Deselect all covariates"""
        for var in self.covariate_checkboxes.values():
            var.set(False)
        self._log('Deselected all covariates')
    
    def _get_selected_covariates(self) -> List[str]:
        """Get list of selected covariate names"""
        return [
            covar for covar, var in self.covariate_checkboxes.items()
            if var.get()
        ]
    
    def _update_reference_group_options(self, groups: List[str]):
        """Update reference group dropdown with available groups"""
        self.ref_group_combo['values'] = [''] + groups
        if not self.ref_group_var.get() and groups:
            self.ref_group_var.set('')  # Default to auto-select
    
    def _try_populate_reference_groups(self):
        """Try to populate reference group dropdown from current data"""
        if not self.data_provider:
            return
        
        try:
            data = self.data_provider()
            group_map = data.get('group_map', {})
            
            if group_map:
                groups = sorted(list(set(group_map.values())))
                self._update_reference_group_options(groups)
                self._log(f'📋 Available groups for reference: {", ".join(groups)}')
            else:
                self._log('⚠️ No groups assigned yet. Complete Step 3 "Auto-Assign Groups" first.')
        except Exception as e:
            pass  # Silently fail, will retry when running analysis

    def _prepare_covariate_dataframe(self, sample_cols: List[str]):
        """Align covariate dataframe to sample columns and standardize labels."""
        if self.df_covariates is None:
            return None, []

        df = self.df_covariates.copy()

        # Normalize labels to improve matching (sample IDs only; keep covariate names intact)
        normalize = lambda x: str(x).strip().casefold()
        idx_map = {normalize(idx): idx for idx in df.index}
        col_map = {normalize(col): col for col in df.columns}

        # Build mappings for sample names
        normalized_samples = {s: normalize(s) for s in sample_cols}
        rows_hit = [s for s, norm in normalized_samples.items() if norm in idx_map]
        cols_hit = [s for s, norm in normalized_samples.items() if norm in col_map]

        # Decide orientation: prefer samples-as-rows when available
        use_rows = len(rows_hit) >= len(cols_hit)
        missing = []

        if use_rows and rows_hit:
            aligned_rows = [idx_map[normalized_samples[s]] for s in sample_cols if normalized_samples[s] in idx_map]
            df_aligned = df.loc[aligned_rows]
            # Restore original sample names as index
            df_aligned.index = [s for s in sample_cols if normalized_samples[s] in idx_map]
            missing = [s for s, norm in normalized_samples.items() if norm not in idx_map]
            orientation = 'samples_as_rows'
        elif cols_hit:
            aligned_cols = [col_map[normalized_samples[s]] for s in sample_cols if normalized_samples[s] in col_map]
            df_aligned = df[aligned_cols].T
            df_aligned.index = [s for s in sample_cols if normalized_samples[s] in col_map]
            missing = [s for s, norm in normalized_samples.items() if norm not in col_map]
            orientation = 'samples_as_columns_transposed'
        else:
            return None, sample_cols

        self._log(
            f'🔗 Covariate mapping: orientation={orientation}, matched_samples={len(df_aligned.index)}/{len(sample_cols)}, covariate_columns={len(df_aligned.columns)}'
        )
        mapped_cols_preview = ', '.join([str(c) for c in df_aligned.columns[:20]])
        if len(df_aligned.columns) > 20:
            mapped_cols_preview += ', ...'
        self._log(f'📋 Mapped covariate columns: {mapped_cols_preview}')

        return df_aligned, missing

    def preview_covariate_columns(self):
        """Show all detected covariate columns and sample ID preview for mapping checks."""
        if self.df_covariates is None:
            messagebox.showinfo('Covariate Preview', 'No covariate file loaded yet.')
            return

        n_rows = len(self.df_covariates.index)
        n_cols = len(self.df_covariates.columns)
        col_lines = [f'{i + 1}. {col}' for i, col in enumerate(self.df_covariates.columns)]
        sample_preview = ', '.join([str(i) for i in self.df_covariates.index[:10]])
        if len(self.df_covariates.index) > 10:
            sample_preview += ', ...'

        messagebox.showinfo(
            'Covariate Columns',
            f'Rows (sample IDs): {n_rows}\n'
            f'Columns (covariates): {n_cols}\n\n'
            f'Sample ID preview:\n{sample_preview}\n\n'
            f'Covariate columns:\n' + '\n'.join(col_lines)
        )
    
    def run_covariate_analysis(self):
        """Run covariate-adjusted analysis."""
        self._log('🔍 COVARIATE DEBUG: run_covariate_analysis called')
        # Groups are already confirmed when dialog was opened, proceed directly
        self._run_covariate_analysis_core()

    def _run_covariate_analysis_core(self):
        """Execute covariate-adjusted analysis once group configuration is confirmed."""
        if not self.data_provider:
            messagebox.showerror('Error', 'No data provider configured')
            return
        
        # Get selected covariates (can be empty now)
        selected_covars = self._get_selected_covariates()
        
        try:
            # Get data
            data = self.data_provider()
            df_intensities = data.get('df_intensities')
            sample_cols = data.get('sample_cols')
            group_map = data.get('group_map')
            metabolite_id_col = data.get('metabolite_id_col')
            # Get lipid class data if available
            mode = data.get('mode', 'metabolite')
            df_class = data.get('df_class')
            class_id_col = data.get('class_id_col')
            
            # Enhanced error checking with specific messages
            if df_intensities is None:
                messagebox.showerror(
                    'Error',
                    'No metabolite data loaded.\n\nPlease complete Step 2: Import Data first.'
                )
                return
                
            if not sample_cols:
                messagebox.showerror(
                    'Error',
                    'No sample columns detected.\n\nPlease complete Step 2: Import Data first.'
                )
                return
                
            if not group_map:
                messagebox.showerror(
                    'Error',
                    'No group assignments found.\n\n'
                    'Please complete Step 3:\n'
                    '1. Set your group labels (Control, PD, etc.)\n'
                    '2. Click "Auto-Assign Groups" to assign samples to groups'
                )
                return

            # Ensure every sample has a group assignment
            missing_groups = [s for s in sample_cols if s not in group_map or not group_map[s]]
            if missing_groups:
                messagebox.showerror(
                    'Error',
                    'Some samples are missing group assignments: '\
                    f'{", ".join(missing_groups)}. Please set groups in Step 3.'
                )
                return
            
            # Get settings
            ref_group = self.ref_group_var.get() or None
            # Get p-value adjustment method
            p_adjust_method = self.pairwise_p_adjust_method.get()
            apply_fdr = (p_adjust_method != 'None')
            
            try:
                alpha = float(self.alpha_var.get())
            except ValueError:
                alpha = 0.05
            
            # Update reference group options
            groups = list(set(group_map.values()))
            self._update_reference_group_options(groups)
            
            # Determine analysis description based on covariates
            if len(selected_covars) > 0:
                self._log(f'Running analysis with {len(selected_covars)} covariate(s)...')
            else:
                self._log('Running analysis without covariates...')
            
            self._safe_widget_config(self.run_covariate_btn, state='disabled', text='⏳ Running...')
            self.parent_frame.update()

            # Freeze a deterministic snapshot for this run so method switches or UI edits
            # after clicking Run do not leak into the executing analysis.
            self._analysis_run_id += 1
            current_run_id = self._analysis_run_id
            analysis_method = str(self.analysis_method_var.get())
            p_adjust_method = str(self.pairwise_p_adjust_method.get())
            selected_covars_snapshot = list(selected_covars)
            sample_cols_snapshot = list(sample_cols)
            group_map_snapshot = dict(group_map)
            df_intensities_snapshot = df_intensities.copy(deep=True)
            df_class_snapshot = df_class.copy(deep=True) if isinstance(df_class, pd.DataFrame) else None

            # Align covariate dataframe to sample order (only if covariates are selected)
            aligned_covars = None
            if len(selected_covars) > 0:
                aligned_covars, missing_samples = self._prepare_covariate_dataframe(sample_cols)
                if aligned_covars is None:
                    messagebox.showerror(
                        'Error',
                        'Could not align covariate file to your sample IDs. Please verify Sample_ID names.'
                    )
                    self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis')
                    return

                if missing_samples:
                    self._log(
                        f'⚠️ Covariate file missing these samples (not used): {", ".join(missing_samples)}'
                    )

                aligned_covars_norm = {str(c).strip().casefold() for c in aligned_covars.columns}
                missing_covars = [c for c in selected_covars if str(c).strip().casefold() not in aligned_covars_norm]
                if missing_covars:
                    self._log(
                        f'⚠️ Selected covariates not found after mapping (will be ignored): {", ".join(missing_covars)}'
                    )
            aligned_covars_snapshot = aligned_covars.copy(deep=True) if isinstance(aligned_covars, pd.DataFrame) else None
            
            # Log data dimensions for debugging
            self._log(f'📊 Data dimensions: df_intensities shape={df_intensities.shape}, metabolite_col={metabolite_id_col}')
            self._log(f'📊 Sample columns: {len(sample_cols)} samples')
            if aligned_covars is not None:
                self._log(f'📊 Covariate data: {aligned_covars.shape[0]} samples × {aligned_covars.shape[1]} available variables')
            else:
                self._log('📊 Covariate data: none (running without covariates)')
            self._log(f'📊 Group map: {len(group_map)} assignments')

            # Get group order from parent tab if available
            group_order = None
            if self.parent_tab and hasattr(self.parent_tab, 'group_definitions'):
                group_order = [self.parent_tab.group_definitions[gid] 
                              for gid in sorted(self.parent_tab.group_definitions.keys())]
            
            # Get filtering parameters from parent tab
            min_samples_per_group = 2
            min_samples_type = 'absolute'
            if self.parent_tab:
                # First get the type to know which variable to read
                if hasattr(self.parent_tab, 'min_samples_type_var'):
                    min_samples_type = self.parent_tab.min_samples_type_var.get()
                    self._log(f'🔍 COVARIATE DEBUG: Read min_samples_type={min_samples_type} from parent tab')
                else:
                    self._log(f'⚠️ COVARIATE DEBUG: Parent tab missing min_samples_type_var')
                
                # Now read the appropriate value based on type
                if min_samples_type == 'percentage':
                    if hasattr(self.parent_tab, 'min_samples_percent_var'):
                        try:
                            min_samples_per_group = float(self.parent_tab.min_samples_percent_var.get())
                            self._log(f'🔍 COVARIATE DEBUG: Read min_samples_percent={min_samples_per_group}% from parent tab')
                        except Exception as e:
                            self._log(f'⚠️ COVARIATE DEBUG: Failed to read min_samples_percent_var: {e}')
                            min_samples_per_group = 50.0
                    else:
                        self._log(f'⚠️ COVARIATE DEBUG: Parent tab missing min_samples_percent_var')
                        min_samples_per_group = 50.0
                else:  # absolute
                    if hasattr(self.parent_tab, 'min_samples_per_group_var'):
                        try:
                            min_samples_per_group = int(self.parent_tab.min_samples_per_group_var.get())
                            self._log(f'🔍 COVARIATE DEBUG: Read min_samples_per_group={min_samples_per_group} samples from parent tab')
                        except Exception as e:
                            self._log(f'⚠️ COVARIATE DEBUG: Failed to read min_samples_per_group: {e}')
                            min_samples_per_group = 2
                    else:
                        self._log(f'⚠️ COVARIATE DEBUG: Parent tab missing min_samples_per_group_var')
                        min_samples_per_group = 2
            else:
                self._log(f'⚠️ COVARIATE DEBUG: No parent tab available')
            
            # Log filtering settings
            if min_samples_type == 'percentage':
                self._log(f'📊 Filtering: ≥{min_samples_per_group}% valid samples per group')
            else:
                self._log(f'📊 Filtering: ≥{min_samples_per_group} valid samples per group')
            
            # Log lipid class data availability
            if mode == 'lipid' and df_class is not None:
                self._log(f'📊 Lipid class data: {len(df_class)} classes will also be analyzed')
            
            # Run analysis in background thread to prevent UI hanging
            def run_analysis_thread():
                try:
                    # Select the appropriate analysis function
                    if analysis_method == 'Limma':
                        # Pure Python Limma implementation using empirical Bayes moderated statistics
                        from main_script.covariate_adjustment import run_limma_covariate_analysis, SCIPY_AVAILABLE
                        if not SCIPY_AVAILABLE:
                            error_msg = (
                                "❌ Limma requires scipy\n\n"
                                "Please install scipy:\n"
                                "   pip install scipy\n\n"
                                "Alternatively, use 'Linear Model (OLS)' method instead."
                            )
                            self._ui_after(0, lambda: self._log(f'\n❌ ERROR: {error_msg}\n'))
                            self._ui_after(0, lambda: messagebox.showerror("Limma Not Available", error_msg))
                            self._ui_after(0, lambda: self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis') if current_run_id == self._analysis_run_id else None)
                            return
                        analysis_func = run_limma_covariate_analysis
                        self._ui_after(0, lambda: self._log(f'🔬 Using Limma (empirical Bayes) method for analysis'))
                    else:
                        from main_script.covariate_adjustment import run_covariate_adjusted_analysis
                        analysis_func = run_covariate_adjusted_analysis
                        self._ui_after(0, lambda: self._log(f'🔬 Using Linear Model (OLS) for analysis'))
                    
                    # Prepare covariate parameters (can be None/empty for no covariates)
                    covariate_params = {
                        'covariate_data': aligned_covars_snapshot if len(selected_covars_snapshot) > 0 else None,
                        'covariate_cols': selected_covars_snapshot if len(selected_covars_snapshot) > 0 else None
                    }
                    
                    # Run main analysis on individual lipids/metabolites
                    local_covariate_results = analysis_func(
                        df_intensities=df_intensities_snapshot,
                        sample_cols=sample_cols_snapshot,
                        group_map=group_map_snapshot,
                        **covariate_params,
                        reference_group=ref_group,
                        apply_fdr=apply_fdr,
                        fdr_method=p_adjust_method,
                        alpha=alpha,
                        metabolite_id_col=metabolite_id_col,
                        return_adjusted_intensities=False,
                        group_order=group_order,
                        min_samples_per_group=min_samples_per_group,
                        min_samples_type=min_samples_type
                    )
                    
                    # Run analysis on lipid class data if available
                    local_covariate_class_results = None
                    if mode == 'lipid' and df_class_snapshot is not None and class_id_col is not None:
                        try:
                            self._ui_after(0, lambda: self._log('🔬 Running analysis on lipid class data...'))
                            local_covariate_class_results = analysis_func(
                                df_intensities=df_class_snapshot,
                                sample_cols=sample_cols_snapshot,
                                group_map=group_map_snapshot,
                                **covariate_params,
                                reference_group=ref_group,
                                apply_fdr=apply_fdr,
                                fdr_method=p_adjust_method,
                                alpha=alpha,
                                metabolite_id_col=class_id_col,
                                return_adjusted_intensities=False,
                                group_order=group_order,
                                min_samples_per_group=min_samples_per_group,
                                min_samples_type=min_samples_type
                            )
                            self._ui_after(0, lambda: self._log(f'✅ Lipid class analysis complete: {local_covariate_class_results.summary_stats.get("n_metabolites_tested", 0)} classes tested'))
                        except Exception as e:
                            self._ui_after(0, lambda: self._log(f'⚠️ Lipid class analysis failed: {e}'))
                            logger.warning(f'Lipid class covariate analysis failed: {e}')
                    
                    # Display results summary (on main thread)
                    def _finish_success():
                        if current_run_id != self._analysis_run_id:
                            self._log(f'ℹ️ Ignoring stale covariate run #{current_run_id}; newer run #{self._analysis_run_id} is active.')
                            return
                        self.covariate_results = local_covariate_results
                        self.covariate_class_results = local_covariate_class_results
                        self._show_analysis_results(alpha, selected_covars_snapshot)
                    self._ui_after(0, _finish_success)
                    
                except ImportError as e:
                    error_msg = str(e)
                    if 'scipy' in error_msg:
                        self._ui_after(0, lambda: messagebox.showerror(
                            'Missing Dependency',
                            'Limma analysis requires scipy.\n\n'
                            'Install with: pip install scipy\n\n'
                            'Or use "Linear Model (OLS)" method instead.'
                        ))
                    else:
                        self._ui_after(0, lambda: messagebox.showerror(
                            'Missing Dependency',
                            'Analysis requires statsmodels and scipy.\n\n'
                            'Install with: pip install statsmodels scipy'
                        ))
                    logger.error(f'Import error: {e}')
                    self._ui_after(0, lambda: self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis') if current_run_id == self._analysis_run_id else None)
                
                except Exception as e:
                    self._ui_after(0, lambda: messagebox.showerror('Error', f'Analysis failed:\n{str(e)}'))
                    logger.error(f'Covariate analysis error: {e}', exc_info=True)
                    self._ui_after(0, lambda: self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis') if current_run_id == self._analysis_run_id else None)
            
            # Start analysis thread
            # rpy2 stores conversion rules in a contextvar; ensure the new thread inherits it
            ctx = contextvars.copy_context()
            analysis_thread = threading.Thread(target=ctx.run, args=(run_analysis_thread,), daemon=True)
            analysis_thread.start()
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to start analysis:\n{str(e)}')
            logger.error(f'Analysis startup error: {e}', exc_info=True)
            self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis')
    
    def _show_analysis_results(self, alpha, selected_covars):
        """Show analysis results and re-enable UI (called on main thread)"""
        try:
            if not self._widget_exists(self.parent_frame):
                return

            summary = self.covariate_results.summary_stats
            n_sig_per_comp = summary.get('n_significant_per_comparison', {})
            n_tested = summary.get('n_metabolites_tested', 0)
            n_excluded = summary.get('n_metabolites_excluded', 0)
            missing_covars = summary.get('covariates_missing', 'None')
            used_covars = summary.get('covariates', 'None')
            requested_covars = summary.get('covariates_requested', 'None')
            ref_group_used = summary.get('reference_group', 'None')
            
            # Format significance results per comparison
            sig_msg = []
            for comp, count in n_sig_per_comp.items():
                sig_msg.append(f'{comp}: {count}')
            
            self._log(
                f'✓ Covariate analysis complete: {n_tested} metabolites tested '
                f'({n_excluded} excluded due to insufficient data, α={alpha}, p-adjust={self.pairwise_p_adjust_method.get()})'
            )
            for msg in sig_msg:
                self._log(f'  Significant in {msg}')
            
            # Log lipid class results if available
            class_msg = ''
            if self.covariate_class_results is not None:
                class_summary = self.covariate_class_results.summary_stats
                class_tested = class_summary.get('n_metabolites_tested', 0)
                class_excluded = class_summary.get('n_metabolites_excluded', 0)
                class_sig = class_summary.get('n_significant_per_comparison', {})
                
                self._log(f'✓ Lipid class analysis: {class_tested} classes tested ({class_excluded} excluded)')
                for comp, count in class_sig.items():
                    self._log(f'  Class significant in {comp}: {count}')
                
                class_msg = f'\n\nLipid Classes tested: {class_tested}\nClass significant per comparison:\n'
                class_msg += '\n'.join(f'  {comp}: {count}' for comp, count in class_sig.items())
            
            # Prepare message based on whether covariates were used
            if selected_covars:
                covariate_info = (
                    f'Covariates requested: {requested_covars}\n'
                    f'Covariates used after mapping: {used_covars}'
                )
            else:
                covariate_info = 'No covariates used'
            analysis_type = summary.get('method', 'Linear Model')
            if missing_covars and missing_covars != 'None':
                self._log(f'⚠️ Covariates requested but not used after mapping: {missing_covars}')
            
            messagebox.showinfo(
                'Analysis Complete',
                f'Analysis finished using {analysis_type}!\n\n'
                f'Metabolites tested: {n_tested}\n'
                f'Metabolites excluded: {n_excluded} (insufficient valid data)\n'
                f'Reference group used: {ref_group_used}\n'
                f'Significant per comparison:\n' + '\n'.join(f'  {msg}' for msg in sig_msg) + 
                class_msg + '\n\n'
                f'{covariate_info}\n'
                f'Covariates not used after mapping: {missing_covars}\n\n'
                f'Use "Export Results" to save.'
            )
            
            # Enable export button and reset run button
            self._safe_widget_config(self.export_covariate_btn, state='normal')
            self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Run Analysis')
            
            # Auto-load results to Visualization tab
            try:
                self._log('🔍 VIZ DEBUG: Starting result storage for visualization...')
                if self.parent_tab:
                    complete_df = self.covariate_results.metabolite_results
                    self._log(f'🔍 VIZ DEBUG: complete_df shape: {complete_df.shape if complete_df is not None else "None"}')
                    
                    if complete_df is not None:
                        self._log(f'🔍 VIZ DEBUG: complete_df columns: {list(complete_df.columns)[:10]}...')
                        
                        # Store covariate results in parent tab's statistical_test_results for visualization access
                        # Use 'enhanced_metabolites' key which is what visualization tab expects
                        if not hasattr(self.parent_tab, 'statistical_test_results'):
                            self.parent_tab.statistical_test_results = {}
                            self._log('🔍 VIZ DEBUG: Created new statistical_test_results dict')
                        
                        # Store as 'enhanced_metabolites' - the key that visualization tab looks for
                        self.parent_tab.statistical_test_results['enhanced_metabolites'] = complete_df
                        self._log('✅ VIZ DEBUG: Stored as enhanced_metabolites')
                        
                        # Also store under covariate-specific key for backward compatibility
                        self.parent_tab.statistical_test_results['covariate_complete'] = complete_df
                        self._log('✅ VIZ DEBUG: Stored as covariate_complete')
                        
                        # Add pairwise sheets if available
                        if hasattr(self.covariate_results, 'pairwise_sheets'):
                            self._log(f'🔍 VIZ DEBUG: Storing {len(self.covariate_results.pairwise_sheets)} pairwise sheets')
                            for comp_name, comp_df in self.covariate_results.pairwise_sheets.items():
                                self.parent_tab.statistical_test_results[f'covariate_{comp_name}'] = comp_df
                        
                        # Store lipid class results if available
                        if self.covariate_class_results is not None:
                            class_complete_df = self.covariate_class_results.metabolite_results
                            if class_complete_df is not None:
                                self.parent_tab.statistical_test_results['covariate_class_complete'] = class_complete_df
                                self._log('✅ VIZ DEBUG: Stored lipid class results as covariate_class_complete')
                                
                                if hasattr(self.covariate_class_results, 'pairwise_sheets'):
                                    for comp_name, comp_df in self.covariate_class_results.pairwise_sheets.items():
                                        self.parent_tab.statistical_test_results[f'covariate_class_{comp_name}'] = comp_df
                        
                        # Store in memory_store for Visualization tab to access
                        if hasattr(self.parent_tab, 'memory_store'):
                            # Store in the format that visualization tab expects
                            self.parent_tab.memory_store['statistical_test_results'] = {
                                'enhanced_metabolites': complete_df
                            }
                            self.parent_tab.memory_store['complete_df'] = complete_df
                            self.parent_tab.memory_store['statistical_results'] = complete_df
                            self.parent_tab.memory_store['id_column'] = 'Name'
                            self._log('✅ Stored covariate results in shared memory for visualization')
                            self._log(f'🔍 VIZ DEBUG: memory_store now has keys: {list(self.parent_tab.memory_store.keys())}')
                        else:
                            self._log('⚠️ VIZ DEBUG: Parent tab has no memory_store attribute')
                        
                        # Schedule auto-load to visualization tab (delayed to ensure data is ready)
                        if hasattr(self.parent_tab, 'root'):
                            self._log('🔍 VIZ DEBUG: Scheduling _auto_load_to_visualization with 500ms delay')
                            self.parent_tab.root.after(500, self._auto_load_to_visualization)
                        else:
                            # Fallback: try immediate load
                            self._log('🔍 VIZ DEBUG: No root, calling _auto_load_to_visualization immediately')
                            self._auto_load_to_visualization()
                    else:
                        self._log('❌ VIZ DEBUG: complete_df is None, cannot store results')
                else:
                    self._log('❌ VIZ DEBUG: No parent_tab available for storing results')
            except Exception as e:
                import traceback
                self._log(f'❌ VIZ DEBUG: Exception during result storage: {e}')
                self._log(f'❌ VIZ DEBUG: Traceback:\n{traceback.format_exc()}')
            
        except Exception as e:
            logger.error(f'Error displaying results: {e}', exc_info=True)
            self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Statistics with Covariate')
        
        finally:
            self._safe_widget_config(self.run_covariate_btn, state='normal', text='▶️ Statistics with Covariate')
    
    def export_covariate_results(self):
        """Export covariate analysis results"""
        if self.covariate_results is None:
            messagebox.showwarning('Warning', 'No results to export. Run analysis first.')
            return
        
        # Ask for output file
        file_path = filedialog.asksaveasfilename(
            title='Save Covariate Results',
            defaultextension='.xlsx',
            filetypes=[
                ('Excel files', '*.xlsx'),
                ('All files', '*.*')
            ]
        )
        
        if not file_path:
            return
        
        # Disable export button during export
        self._safe_widget_config(self.export_covariate_btn, state='disabled', text='⏳ Exporting...')
        self.parent_frame.update()
        
        # Get lipid class results if available
        class_results = getattr(self, 'covariate_class_results', None)
        
        # Run export in background thread to prevent UI hanging
        def export_thread():
            try:
                from main_script.covariate_adjustment import export_covariate_results
                
                export_covariate_results(
                    self.covariate_results,
                    file_path,
                    include_diagnostics=True,
                    include_coefficients=True,
                    class_results=class_results  # Pass lipid class results
                )
                
                # Track export path in parent tab for Open Folder button (on main thread)
                def update_export_path():
                    if self.parent_tab and hasattr(self.parent_tab, 'last_statistics_export_path'):
                        self.parent_tab.last_statistics_export_path = file_path
                    
                    self._log(f'✓ Exported covariate results to {os.path.basename(file_path)}')
                    
                    messagebox.showinfo(
                        'Export Complete',
                        f'Covariate results exported successfully!\n\n{file_path}'
                    )
                    
                    # Re-enable export button
                    self._safe_widget_config(self.export_covariate_btn, state='normal', text='💾 Export Covariate Results')
                    
                    # Auto-switch to Visualization tab after export (like other stats)
                    try:
                        if self.parent_tab and hasattr(self.parent_tab, 'switch_to_tab'):
                            self._ui_after(300, lambda: self.parent_tab.switch_to_tab("📊 Visualization"))
                            self._log('✓ Switched to Visualization tab')
                    except Exception as e:
                        logger.warning(f'Could not auto-switch to visualization: {e}')
                
                self._ui_after(0, update_export_path)
                
            except Exception as e:
                def show_error():
                    messagebox.showerror('Export Error', f'Export failed:\n{str(e)}')
                    logger.error(f'Export error: {e}', exc_info=True)
                    self._safe_widget_config(self.export_covariate_btn, state='normal', text='💾 Export Covariate Results')
                
                self._ui_after(0, show_error)
        
        # Start export thread
        threading.Thread(target=export_thread, daemon=True).start()


class CovariateAdjustmentDialog:
    """
    Standalone dialog window for covariate adjustment.
    Opens in a separate window with all covariate adjustment features.
    """
    
    def __init__(
        self,
        parent: tk.Tk,
        log_callback: Optional[Callable[[str], None]] = None,
        data_provider: Optional[Callable[[], Dict]] = None,
        parent_tab = None
    ):
        """
        Initialize covariate adjustment dialog.
        
        Parameters
        ----------
        parent : tk.Tk
            Parent window
        log_callback : Optional[Callable]
            Function to call for logging messages
        data_provider : Optional[Callable]
            Function that provides metabolite data and group assignments
        parent_tab : Optional
            The parent Statistics tab instance for accessing group configuration
        """
        self.parent = parent
        self.parent_tab_instance = parent_tab  # Store the actual parent tab
        self.log_callback = log_callback or print
        self.data_provider = data_provider
        
        # Create dialog window
        self.dialog = tk.Toplevel(parent)
        self.dialog.title('🎯 Covariate Adjustment Console')
        
        # Set initial size and make resizable
        self.dialog.geometry('900x700')
        self.dialog.minsize(750, 600)  # Set minimum size
        self.dialog.configure(bg='#f0f0f0')
        
        # Make it NOT modal so it has minimize/maximize buttons
        self.dialog.transient(parent)
        # Don't grab_set() - this removes minimize/maximize buttons
        
        # Center on screen
        self.dialog.update_idletasks()
        width = self.dialog.winfo_width()
        height = self.dialog.winfo_height()
        x = (self.dialog.winfo_screenwidth() // 2) - (width // 2)
        y = (self.dialog.winfo_screenheight() // 2) - (height // 2)
        self.dialog.geometry(f'{width}x{height}+{x}+{y}')
        
        # Create main content
        self._create_dialog_ui()
    
    def _create_dialog_ui(self):
        """Create the dialog UI"""
        # Header
        header_frame = tk.Frame(self.dialog, bg='#2c3e50', height=60)
        header_frame.pack(fill='x', side='top')
        header_frame.pack_propagate(False)
        
        tk.Label(
            header_frame,
            text='🎯 Covariate Adjustment Console',
            bg='#2c3e50', fg='white',
            font=('Arial', 16, 'bold')
        ).pack(pady=15)
        
        # Info banner
        info_frame = tk.Frame(self.dialog, bg='#e3f2fd', relief='solid', borderwidth=1)
        info_frame.pack(fill='x', padx=10, pady=10)
        
        info_text = (
            "Adjust for covariates (Age, Sex, BMI, PMI, etc.) using linear regression (ANCOVA)."
            "This removes the effect of confounding variables while testing group differences.\n"
            "📋 Requirements: Complete Steps 1-3 in the main Statistics tab first."
        )
        tk.Label(
            info_frame,
            text=info_text,
            bg='#e3f2fd', fg='#1565c0',
            font=('Arial', 9),
            justify='left',
            wraplength=700
        ).pack(padx=10, pady=10)
        
        # Main content with scrollbar
        canvas_frame = tk.Frame(self.dialog, bg='#f0f0f0')
        canvas_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        canvas = tk.Canvas(canvas_frame, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(canvas_frame, orient='vertical', command=canvas.yview)
        content_frame = tk.Frame(canvas, bg='#f0f0f0')
        
        def on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox('all'))
        
        def on_canvas_configure(event):
            # Update the window width to match canvas width
            canvas.itemconfig(canvas_window, width=event.width)
        
        content_frame.bind('<Configure>', on_frame_configure)
        canvas.bind('<Configure>', on_canvas_configure)
        
        canvas_window = canvas.create_window((0, 0), window=content_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mousewheel scrolling
        def on_mousewheel(event):
            try:
                if canvas.winfo_exists():
                    canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except:
                pass
        
        canvas.bind("<MouseWheel>", on_mousewheel)
        content_frame.bind("<MouseWheel>", on_mousewheel)
        
        # Cleanup binding when dialog closes
        def on_close():
            try:
                canvas.unbind("<MouseWheel>")
                content_frame.unbind("<MouseWheel>")
            except:
                pass
        
        self.dialog.protocol("WM_DELETE_WINDOW", lambda: [on_close(), self.dialog.destroy()])
        
        scrollbar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        
        # Create covariate adjustment section inside the content frame
        self.covariate_section = CovariateAdjustmentSection(
            parent_frame=content_frame,
            log_callback=self.log_callback,
            data_provider=self.data_provider,
            parent_tab=self.parent_tab_instance  # Pass the actual parent Statistics tab
        )
        
        # Footer with close button
        footer_frame = tk.Frame(self.dialog, bg='#f0f0f0', height=50)
        footer_frame.pack(fill='x', side='bottom', padx=10, pady=10)
        footer_frame.pack_propagate(False)
        
        tk.Button(
            footer_frame,
            text='❌ Close',
            command=self.dialog.destroy,
            bg='#e74c3c', fg='white',
            font=('Arial', 10, 'bold'),
            relief='raised', bd=2,
            padx=20, pady=5
        ).pack(side='right', padx=5)
        
        tk.Label(
            footer_frame,
            text='Tip: Results are logged to the main Statistics tab',
            bg='#f0f0f0', fg='#666',
            font=('Arial', 8, 'italic')
        ).pack(side='left', padx=5)

