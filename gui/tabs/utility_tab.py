"""
Utility Tab - Contains various utility tools for data visualization
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import os
import json
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
import threading

# Import shared components
from gui.shared.base_tab import BaseTab
from gui.shared.column_assignment import show_column_assignment_dialog
from main_script.glycan_classification import process_glycan_dataframe

logger = logging.getLogger(__name__)


class UtilityTab(BaseTab):
    """Utility Tab - Contains Chord Diagram, Venn Diagram, Pie Chart, and Heatmap tools"""
    
    def __init__(self, parent, data_manager):
        """Initialize Utility tab"""
        super().__init__(parent, data_manager)
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Setup memory_store as reference to data_manager's memory store
        self.memory_store = self.data_manager.memory_store
        
        # Config file path for persisting group configurations
        self.config_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.utility_config_file = os.path.join(self.config_dir, 'utility_group_configs.json')
        
        # Create UI
        self.setup_ui()
        logger.info("[OK] Utility Tab initialized")
    
    def _load_group_configs(self):
        """Load saved group configurations from file."""
        try:
            if os.path.exists(self.utility_config_file):
                with open(self.utility_config_file, 'r') as f:
                    return json.load(f)
        except Exception as e:
            logger.warning(f"Failed to load group configs: {e}")
        return {}
    
    def _save_group_configs(self, configs):
        """Save group configurations to file."""
        try:
            with open(self.utility_config_file, 'w') as f:
                json.dump(configs, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save group configs: {e}")
    
    def setup_ui(self):
        """Setup the UI for the Utility tab with internal sub-tabs"""
        # Create internal notebook for sub-tabs
        self.notebook = ttk.Notebook(self.frame)
        self.notebook.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create sub-tabs
        self.setup_chord_diagram_tab()
        self.setup_venn_diagram_tab()
        self.setup_pie_chart_tab()
        self.setup_heatmap_tab()
        self.setup_effect_size_plot_tab()
        self.setup_covariate_category_stats_tab()
        self.setup_linear_regression_plot_tab()
        self.setup_glycan_classification_tab()

    # ==================== COVARIATE CATEGORY STATS TAB ====================
    def setup_covariate_category_stats_tab(self):
        """Setup flexible category/covariate statistics utility tab."""
        stats_frame = ttk.Frame(self.notebook)
        self.notebook.add(stats_frame, text='Regression Utility')

        # Scrollable host so the full utility remains accessible on smaller windows.
        scroll_host = tk.Frame(stats_frame, bg='#f0f0f0')
        scroll_host.pack(fill='both', expand=True)

        h_scroll = ttk.Scrollbar(scroll_host, orient='horizontal')
        h_scroll.pack(side='bottom', fill='x')

        v_scroll = ttk.Scrollbar(scroll_host, orient='vertical')
        v_scroll.pack(side='right', fill='y')

        canvas = tk.Canvas(
            scroll_host,
            bg='#f0f0f0',
            highlightthickness=0,
            yscrollcommand=v_scroll.set,
            xscrollcommand=h_scroll.set
        )
        canvas.pack(side='left', fill='both', expand=True)
        v_scroll.configure(command=canvas.yview)
        h_scroll.configure(command=canvas.xview)

        container = tk.Frame(canvas, bg='#f0f0f0')
        container_window = canvas.create_window((0, 0), window=container, anchor='nw')

        def _update_scroll_region(_event=None):
            canvas.configure(scrollregion=canvas.bbox('all'))

        def _fit_container_width(event):
            required = container.winfo_reqwidth()
            canvas.itemconfigure(container_window, width=max(event.width, required))

        container.bind('<Configure>', _update_scroll_region)
        canvas.bind('<Configure>', _fit_container_width)

        tk.Label(
            container,
            text='Flexible Outcome Statistics Utility',
            font=('Arial', 14, 'bold'),
            bg='#f0f0f0'
        ).pack(anchor='w', pady=(0, 8))

        tk.Label(
            container,
            text=(
                'Load a table with sample IDs, covariates/outcomes, and metabolite columns. '
                'Select mappings manually and run per-outcome analysis with optional covariate adjustment.'
            ),
            bg='#f0f0f0',
            font=('Arial', 9),
            wraplength=1000,
            justify='left'
        ).pack(anchor='w', pady=(0, 8))

        # Step 1: data
        load_frame = tk.LabelFrame(container, text='Step 1: Load Data', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        load_frame.pack(fill='x', pady=4)

        row = tk.Frame(load_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=6, pady=6)

        tk.Button(
            row,
            text='Upload Data File',
            command=self.upload_outcome_stats_data,
            bg='#2980b9',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=(0, 8))

        self.outcome_stats_file_label = tk.Label(row, text='No file loaded', bg='#f0f0f0', font=('Arial', 9), fg='#555')
        self.outcome_stats_file_label.pack(side='left')

        # Step 2: mapping
        map_frame = tk.LabelFrame(container, text='Step 2: Map Columns', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        map_frame.pack(fill='both', expand=False, pady=4)

        sample_row = tk.Frame(map_frame, bg='#f0f0f0')
        sample_row.pack(fill='x', padx=6, pady=5)
        tk.Label(sample_row, text='Sample ID Column:', bg='#f0f0f0', width=20, anchor='w').pack(side='left')

        self.outcome_stats_sample_col_var = tk.StringVar(value='')
        self.outcome_stats_sample_col_combo = ttk.Combobox(
            sample_row,
            textvariable=self.outcome_stats_sample_col_var,
            state='readonly',
            width=40
        )
        self.outcome_stats_sample_col_combo.pack(side='left', padx=5)

        tk.Button(
            sample_row,
            text='Auto-detect',
            command=self._auto_detect_outcome_stats_columns,
            bg='#16a085',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=4)

        list_holder = tk.Frame(map_frame, bg='#f0f0f0')
        list_holder.pack(fill='x', padx=6, pady=(0, 6))
        list_holder.grid_columnconfigure(0, weight=1)
        list_holder.grid_columnconfigure(1, weight=1)
        list_holder.grid_columnconfigure(2, weight=1)

        self.outcome_stats_categories_list = self._create_multiselect_listbox(
            list_holder,
            0,
            'Outcome Category Columns (analyzed independently)'
        )
        self.outcome_stats_covariates_list = self._create_multiselect_listbox(
            list_holder,
            1,
            'Covariate Columns (optional)'
        )
        self.outcome_stats_metabolites_list = self._create_multiselect_listbox(
            list_holder,
            2,
            'Metabolite Value Columns'
        )

        map_actions = tk.Frame(map_frame, bg='#f0f0f0')
        map_actions.pack(fill='x', padx=6, pady=(0, 6))
        tk.Button(
            map_actions,
            text='Auto-select Remaining as Metabolites',
            command=self._auto_select_remaining_metabolites,
            bg='#6c5ce7',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left')
        tk.Label(
            map_actions,
            text='Metabolite options exclude Sample ID + selected Category/Covariate columns.',
            bg='#f0f0f0',
            fg='#555',
            font=('Arial', 8)
        ).pack(side='left', padx=8)

        self.outcome_stats_sample_col_combo.bind('<<ComboboxSelected>>', self._on_outcome_mapping_changed)
        self.outcome_stats_categories_list.bind('<<ListboxSelect>>', self._on_outcome_mapping_changed)
        self.outcome_stats_covariates_list.bind('<<ListboxSelect>>', self._on_outcome_mapping_changed)

        # Step 3: analysis settings
        settings_frame = tk.LabelFrame(container, text='Step 3: Analysis Settings', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=4)

        srow = tk.Frame(settings_frame, bg='#f0f0f0')
        srow.pack(fill='x', padx=6, pady=6)

        tk.Label(srow, text='Method:', bg='#f0f0f0').pack(side='left')
        self.outcome_stats_method_var = tk.StringVar(value='Linear Model (OLS)')
        ttk.Combobox(
            srow,
            textvariable=self.outcome_stats_method_var,
            values=['Linear Model (OLS)', 'Limma', 'Ordinal Regression'],
            state='readonly',
            width=20
        ).pack(side='left', padx=5)

        tk.Label(srow, text='P-adjust:', bg='#f0f0f0').pack(side='left', padx=(10, 0))
        self.outcome_stats_padj_var = tk.StringVar(value='BH')
        ttk.Combobox(
            srow,
            textvariable=self.outcome_stats_padj_var,
            values=['BH', 'Bonferroni', 'Holm', 'Hochberg', 'BY', 'None'],
            state='readonly',
            width=12
        ).pack(side='left', padx=5)

        tk.Label(srow, text='Alpha:', bg='#f0f0f0').pack(side='left', padx=(10, 0))
        self.outcome_stats_alpha_var = tk.StringVar(value='0.05')
        tk.Entry(srow, textvariable=self.outcome_stats_alpha_var, width=8).pack(side='left', padx=5)

        cluster_row = tk.Frame(settings_frame, bg='#f0f0f0')
        cluster_row.pack(fill='x', padx=6, pady=(0, 6))
        tk.Label(cluster_row, text='Numeric outcome handling:', bg='#f0f0f0', width=24, anchor='w').pack(side='left')
        self.outcome_stats_numeric_mode_var = tk.StringVar(value='Continuous (OLS)')
        ttk.Combobox(
            cluster_row,
            textvariable=self.outcome_stats_numeric_mode_var,
            values=[
                'Continuous (OLS)',
                'Binary (0 vs 1)',
                'Cluster to 2 groups',
                'Cluster to 3 groups',
                'Auto-cluster (2 or 3)'
            ],
            state='readonly',
            width=24
        ).pack(side='left', padx=5)
        tk.Label(
            cluster_row,
            text='Data-driven bins (k-means) for numeric outcomes',
            bg='#f0f0f0',
            fg='#555',
            font=('Arial', 8)
        ).pack(side='left', padx=8)

        strict_row = tk.Frame(settings_frame, bg='#f0f0f0')
        strict_row.pack(fill='x', padx=6, pady=(0, 6))

        self.outcome_stats_strict_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            strict_row,
            text='Strict outcome handling (drop unknown + force binary 0/1 for two-level text outcomes)',
            variable=self.outcome_stats_strict_var,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(side='left')

        unknown_row = tk.Frame(settings_frame, bg='#f0f0f0')
        unknown_row.pack(fill='x', padx=6, pady=(0, 6))
        tk.Label(unknown_row, text='Unknown tokens to drop:', bg='#f0f0f0', width=20, anchor='w').pack(side='left')
        self.outcome_stats_unknown_tokens_var = tk.StringVar(value='unknown,unk,n/a,na,missing,null')
        tk.Entry(unknown_row, textvariable=self.outcome_stats_unknown_tokens_var, width=50).pack(side='left', padx=5)

        out_row = tk.Frame(settings_frame, bg='#f0f0f0')
        out_row.pack(fill='x', padx=6, pady=(0, 6))
        tk.Label(out_row, text='Output Folder:', bg='#f0f0f0', width=20, anchor='w').pack(side='left')

        self.outcome_stats_output_dir_var = tk.StringVar(value='')
        tk.Entry(out_row, textvariable=self.outcome_stats_output_dir_var, width=60).pack(side='left', padx=5)
        tk.Button(
            out_row,
            text='Browse',
            command=self._browse_outcome_stats_output_dir,
            bg='#7f8c8d',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=4)

        # Run
        run_row = tk.Frame(container, bg='#f0f0f0')
        run_row.pack(fill='x', pady=6)
        self.outcome_stats_run_btn = tk.Button(
            run_row,
            text='Run Per-Category Analysis and Export',
            command=self.run_outcome_stats_analysis,
            bg='#8e44ad',
            fg='white',
            font=('Arial', 10, 'bold'),
            pady=6
        )
        self.outcome_stats_run_btn.pack(fill='x')

        # Log
        log_frame = tk.LabelFrame(container, text='Analysis Log', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        log_frame.pack(fill='both', expand=True, pady=4)

        self.outcome_stats_log = tk.Text(log_frame, height=10, wrap='word', font=('Consolas', 9), bg='white')
        self.outcome_stats_log.pack(side='left', fill='both', expand=True, padx=(6, 0), pady=6)
        log_scroll = ttk.Scrollbar(log_frame, orient='vertical', command=self.outcome_stats_log.yview)
        log_scroll.pack(side='right', fill='y', padx=(0, 6), pady=6)
        self.outcome_stats_log.configure(yscrollcommand=log_scroll.set)

        self.outcome_stats_df = None
        self.outcome_stats_all_columns: List[str] = []
        self._outcome_stats_log('Ready. Upload data and map columns.')

    # ==================== LINEAR REGRESSION PLOT TAB ====================
    def setup_linear_regression_plot_tab(self):
        """Setup sub-tab for continuous scatter + regression line plots."""
        plot_frame = ttk.Frame(self.notebook)
        self.notebook.add(plot_frame, text='Linear Regression Plot')

        scroll_host = tk.Frame(plot_frame, bg='#f0f0f0')
        scroll_host.pack(fill='both', expand=True)

        v_scroll = ttk.Scrollbar(scroll_host, orient='vertical')
        v_scroll.pack(side='right', fill='y')

        canvas = tk.Canvas(
            scroll_host,
            bg='#f0f0f0',
            highlightthickness=0,
            yscrollcommand=v_scroll.set
        )
        canvas.pack(side='left', fill='both', expand=True)
        v_scroll.configure(command=canvas.yview)

        container = tk.Frame(canvas, bg='#f0f0f0')
        container_window = canvas.create_window((0, 0), window=container, anchor='nw')

        def _update_scroll_region(_event=None):
            canvas.configure(scrollregion=canvas.bbox('all'))

        def _fit_container_width(event):
            required = container.winfo_reqwidth()
            canvas.itemconfigure(container_window, width=max(event.width, required))

        container.bind('<Configure>', _update_scroll_region)
        canvas.bind('<Configure>', _fit_container_width)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), 'units')

        canvas.bind('<Enter>', lambda _e: canvas.bind('<MouseWheel>', _on_mousewheel))
        canvas.bind('<Leave>', lambda _e: canvas.unbind('<MouseWheel>'))

        inner = tk.Frame(container, bg='#f0f0f0')
        inner.pack(fill='both', expand=True, padx=10, pady=10)

        tk.Label(
            inner,
            text='Continuous Linear Regression Plot Utility',
            font=('Arial', 14, 'bold'),
            bg='#f0f0f0'
        ).pack(anchor='w', pady=(0, 8))

        tk.Label(
            inner,
            text=(
                'Create scatter plots for continuous outcomes with OLS fit, 95% confidence band, '
                'and optional Y-jitter for overlapping values.'
            ),
            bg='#f0f0f0',
            font=('Arial', 9),
            fg='#333',
            justify='left',
            wraplength=1100
        ).pack(anchor='w', pady=(0, 8))

        load_frame = tk.LabelFrame(inner, text='Step 1: Data Source', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        load_frame.pack(fill='x', pady=4)

        load_row = tk.Frame(load_frame, bg='#f0f0f0')
        load_row.pack(fill='x', padx=6, pady=6)

        tk.Button(
            load_row,
            text='Upload Data File',
            command=self.upload_linear_plot_data,
            bg='#2980b9',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=(0, 6))

        tk.Button(
            load_row,
            text='Use Regression Utility Data',
            command=self.use_outcome_stats_data_for_linear_plot,
            bg='#16a085',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=(0, 8))

        self.linear_plot_file_label = tk.Label(
            load_row,
            text='No file loaded',
            bg='#f0f0f0',
            font=('Arial', 9),
            fg='#555'
        )
        self.linear_plot_file_label.pack(side='left')

        map_frame = tk.LabelFrame(inner, text='Step 2: Map Variables', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        map_frame.pack(fill='x', pady=4)

        mode_row = tk.Frame(map_frame, bg='#f0f0f0')
        mode_row.pack(fill='x', padx=6, pady=5)
        tk.Label(mode_row, text='Plot mode:', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_mode_var = tk.StringVar(value='Single metabolite (X column)')
        ttk.Combobox(
            mode_row,
            textvariable=self.linear_plot_mode_var,
            state='readonly',
            values=[
                'Single metabolite (X column)',
                'Selected metabolites - individual plots',
                'Selected metabolites - combined panel'
            ],
            width=45
        ).pack(side='left', padx=5)

        x_row = tk.Frame(map_frame, bg='#f0f0f0')
        x_row.pack(fill='x', padx=6, pady=5)
        tk.Label(x_row, text='Single X column:', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_x_var = tk.StringVar(value='')
        self.linear_plot_x_combo = ttk.Combobox(x_row, textvariable=self.linear_plot_x_var, state='readonly', width=50)
        self.linear_plot_x_combo.pack(side='left', padx=5)

        y_row = tk.Frame(map_frame, bg='#f0f0f0')
        y_row.pack(fill='x', padx=6, pady=5)
        tk.Label(y_row, text='Y-axis (Continuous outcome):', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_y_var = tk.StringVar(value='')
        self.linear_plot_y_combo = ttk.Combobox(y_row, textvariable=self.linear_plot_y_var, state='readonly', width=50)
        self.linear_plot_y_combo.pack(side='left', padx=5)

        id_row = tk.Frame(map_frame, bg='#f0f0f0')
        id_row.pack(fill='x', padx=6, pady=5)
        tk.Label(id_row, text='Sample ID (optional):', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_sample_var = tk.StringVar(value='')
        self.linear_plot_sample_combo = ttk.Combobox(id_row, textvariable=self.linear_plot_sample_var, state='readonly', width=50)
        self.linear_plot_sample_combo.pack(side='left', padx=5)

        metabolite_frame = tk.LabelFrame(
            map_frame,
            text='Metabolite Selection (for individual/combined modes)',
            bg='#f0f0f0',
            font=('Arial', 9, 'bold')
        )
        metabolite_frame.pack(fill='x', padx=6, pady=(3, 6))

        met_actions = tk.Frame(metabolite_frame, bg='#f0f0f0')
        met_actions.pack(fill='x', padx=4, pady=4)
        tk.Button(
            met_actions,
            text='Refresh Numeric Candidates',
            command=self._refresh_linear_plot_metabolite_candidates,
            bg='#3498db',
            fg='white',
            font=('Arial', 8, 'bold')
        ).pack(side='left', padx=(0, 5))
        tk.Button(
            met_actions,
            text='Clear Selection',
            command=lambda: self.linear_plot_metabolites_list.selection_clear(0, 'end'),
            bg='#95a5a6',
            fg='white',
            font=('Arial', 8, 'bold')
        ).pack(side='left', padx=(0, 8))
        tk.Label(
            met_actions,
            text='Use Ctrl/Shift to select specific metabolites only.',
            bg='#f0f0f0',
            fg='#555',
            font=('Arial', 8)
        ).pack(side='left')

        list_holder = tk.Frame(metabolite_frame, bg='#f0f0f0')
        list_holder.pack(fill='x', padx=4, pady=(0, 5))
        self.linear_plot_metabolites_list = tk.Listbox(
            list_holder,
            selectmode='extended',
            exportselection=False,
            height=8
        )
        self.linear_plot_metabolites_list.pack(side='left', fill='x', expand=True, padx=(0, 4))
        met_scroll = ttk.Scrollbar(list_holder, orient='vertical', command=self.linear_plot_metabolites_list.yview)
        met_scroll.pack(side='right', fill='y')
        self.linear_plot_metabolites_list.configure(yscrollcommand=met_scroll.set)

        settings_frame = tk.LabelFrame(inner, text='Step 3: Plot Settings', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=4)

        row1 = tk.Frame(settings_frame, bg='#f0f0f0')
        row1.pack(fill='x', padx=6, pady=5)

        self.linear_plot_logx_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row1,
            text='Log10-transform X (drop non-positive values)',
            variable=self.linear_plot_logx_var,
            bg='#f0f0f0'
        ).pack(side='left', padx=(0, 12))

        self.linear_plot_ci_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            row1,
            text='Show 95% confidence band',
            variable=self.linear_plot_ci_var,
            bg='#f0f0f0'
        ).pack(side='left')

        row2 = tk.Frame(settings_frame, bg='#f0f0f0')
        row2.pack(fill='x', padx=6, pady=5)

        tk.Label(row2, text='Y jitter (optional):', bg='#f0f0f0').pack(side='left')
        self.linear_plot_jitter_var = tk.StringVar(value='0.0')
        tk.Entry(row2, textvariable=self.linear_plot_jitter_var, width=10).pack(side='left', padx=5)

        tk.Label(row2, text='Point size:', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.linear_plot_point_size_var = tk.StringVar(value='18')
        tk.Entry(row2, textvariable=self.linear_plot_point_size_var, width=8).pack(side='left', padx=5)

        tk.Label(row2, text='Point alpha (0-1):', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.linear_plot_alpha_var = tk.StringVar(value='0.75')
        tk.Entry(row2, textvariable=self.linear_plot_alpha_var, width=8).pack(side='left', padx=5)

        row3 = tk.Frame(settings_frame, bg='#f0f0f0')
        row3.pack(fill='x', padx=6, pady=5)

        tk.Label(row3, text='Figure Width:', bg='#f0f0f0').pack(side='left')
        self.linear_plot_fig_w_var = tk.StringVar(value='7.0')
        tk.Entry(row3, textvariable=self.linear_plot_fig_w_var, width=8).pack(side='left', padx=5)

        tk.Label(row3, text='Figure Height:', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.linear_plot_fig_h_var = tk.StringVar(value='5.0')
        tk.Entry(row3, textvariable=self.linear_plot_fig_h_var, width=8).pack(side='left', padx=5)

        out_frame = tk.LabelFrame(inner, text='Step 4: Export', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        out_frame.pack(fill='x', pady=4)

        out_row = tk.Frame(out_frame, bg='#f0f0f0')
        out_row.pack(fill='x', padx=6, pady=6)
        tk.Label(out_row, text='Output Folder:', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_output_dir_var = tk.StringVar(value='')
        tk.Entry(out_row, textvariable=self.linear_plot_output_dir_var, width=60).pack(side='left', padx=5)
        tk.Button(
            out_row,
            text='Browse',
            command=self._browse_linear_plot_output_dir,
            bg='#7f8c8d',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=4)

        file_row = tk.Frame(out_frame, bg='#f0f0f0')
        file_row.pack(fill='x', padx=6, pady=(0, 6))
        tk.Label(file_row, text='File Prefix:', width=22, anchor='w', bg='#f0f0f0').pack(side='left')
        self.linear_plot_prefix_var = tk.StringVar(value='regression_plot')
        tk.Entry(file_row, textvariable=self.linear_plot_prefix_var, width=30).pack(side='left', padx=5)

        action_row = tk.Frame(inner, bg='#f0f0f0')
        action_row.pack(fill='x', pady=6)

        self.linear_plot_generate_btn = tk.Button(
            action_row,
            text='Generate Plot and Save (PNG)',
            command=self.generate_linear_regression_plot,
            bg='#8e44ad',
            fg='white',
            font=('Arial', 10, 'bold'),
            pady=6
        )
        self.linear_plot_generate_btn.pack(fill='x')

        log_frame = tk.LabelFrame(inner, text='Plot Log', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        log_frame.pack(fill='both', expand=True, pady=4)

        self.linear_plot_log = tk.Text(log_frame, height=8, wrap='word', font=('Consolas', 9), bg='white')
        self.linear_plot_log.pack(side='left', fill='both', expand=True, padx=(6, 0), pady=6)
        plot_log_scroll = ttk.Scrollbar(log_frame, orient='vertical', command=self.linear_plot_log.yview)
        plot_log_scroll.pack(side='right', fill='y', padx=(0, 6), pady=6)
        self.linear_plot_log.configure(yscrollcommand=plot_log_scroll.set)

        self.linear_plot_df = None
        self.linear_plot_all_columns: List[str] = []
        self.linear_plot_numeric_candidates: List[str] = []

        self.linear_plot_y_combo.bind('<<ComboboxSelected>>', lambda _e: self._refresh_linear_plot_metabolite_candidates())
        self.linear_plot_sample_combo.bind('<<ComboboxSelected>>', lambda _e: self._refresh_linear_plot_metabolite_candidates())
        self._linear_plot_log('Ready. Load data and choose X/Y columns.')

    def _linear_plot_log(self, message: str):
        """Append message to linear plot log panel."""
        timestamp = datetime.now().strftime('%H:%M:%S')
        try:
            self.linear_plot_log.insert('end', f'[{timestamp}] {message}\n')
            self.linear_plot_log.see('end')
        except Exception:
            pass

    def _populate_linear_plot_columns(self):
        """Populate combobox options after data load."""
        if self.linear_plot_df is None or self.linear_plot_df.empty:
            return

        cols = list(self.linear_plot_df.columns)
        self.linear_plot_all_columns = cols
        self.linear_plot_x_combo['values'] = cols
        self.linear_plot_y_combo['values'] = cols
        self.linear_plot_sample_combo['values'] = [''] + cols

        # Prefer previously selected columns when valid.
        if self.linear_plot_x_var.get() not in cols:
            self.linear_plot_x_var.set(cols[0] if cols else '')
        if self.linear_plot_y_var.get() not in cols:
            self.linear_plot_y_var.set(cols[1] if len(cols) > 1 else (cols[0] if cols else ''))
        if self.linear_plot_sample_var.get() not in cols:
            guess_sample = ''
            for c in cols:
                key = str(c).strip().casefold()
                if any(tok in key for tok in ['sample', 'sampleid', 'sample_id', 'subject', 'id']):
                    guess_sample = c
                    break
            self.linear_plot_sample_var.set(guess_sample)

        self._refresh_linear_plot_metabolite_candidates()

    def _refresh_linear_plot_metabolite_candidates(self):
        """Refresh multi-select metabolite candidates using numeric columns."""
        if self.linear_plot_df is None or self.linear_plot_df.empty:
            return

        y_col = self.linear_plot_y_var.get().strip()
        sample_col = self.linear_plot_sample_var.get().strip()

        keep_selected = set(self._get_selected_listbox_values(self.linear_plot_metabolites_list))
        candidates = []
        for col in self.linear_plot_df.columns:
            if col == y_col or col == sample_col:
                continue
            numeric = pd.to_numeric(self.linear_plot_df[col], errors='coerce')
            if float(numeric.notna().mean()) >= 0.5 and int(numeric.nunique(dropna=True)) >= 4:
                candidates.append(col)

        self.linear_plot_numeric_candidates = candidates
        self._populate_listbox(self.linear_plot_metabolites_list, candidates)
        self._set_listbox_selection(
            self.linear_plot_metabolites_list,
            candidates,
            [c for c in candidates if c in keep_selected]
        )

    def upload_linear_plot_data(self):
        """Load table for linear regression plotting."""
        file_path = filedialog.askopenfilename(
            title='Select Data File for Linear Plot',
            filetypes=[('Excel/CSV', '*.xlsx *.xls *.csv *.tsv *.txt'), ('All Files', '*.*')]
        )
        if not file_path:
            return

        try:
            if file_path.lower().endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            elif file_path.lower().endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                df = pd.read_csv(file_path, sep='\t')

            self.linear_plot_df = df
            self.linear_plot_file_label.config(
                text=f"Loaded: {os.path.basename(file_path)} ({df.shape[0]} rows x {df.shape[1]} cols)",
                fg='green'
            )
            self._populate_linear_plot_columns()
            self._linear_plot_log(f'Data loaded from {os.path.basename(file_path)}')
        except Exception as e:
            messagebox.showerror('Load Error', f'Failed to load file:\n{e}')
            logger.error(f'Linear plot data load failed: {e}', exc_info=True)

    def use_outcome_stats_data_for_linear_plot(self):
        """Reuse Regression Utility loaded data for quick plotting."""
        if getattr(self, 'outcome_stats_df', None) is None or self.outcome_stats_df.empty:
            messagebox.showwarning('No Data', 'No data is currently loaded in Regression Utility.')
            return

        self.linear_plot_df = self.outcome_stats_df.copy()
        self.linear_plot_file_label.config(
            text=(
                f'Using Regression Utility data '
                f'({self.linear_plot_df.shape[0]} rows x {self.linear_plot_df.shape[1]} cols)'
            ),
            fg='green'
        )
        self._populate_linear_plot_columns()

        # Carry over sample column when available.
        sample_col = getattr(self, 'outcome_stats_sample_col_var', tk.StringVar(value='')).get().strip()
        if sample_col in self.linear_plot_df.columns:
            self.linear_plot_sample_var.set(sample_col)

        self._linear_plot_log('Imported data from Regression Utility tab.')

    def _browse_linear_plot_output_dir(self):
        folder = filedialog.askdirectory(title='Select Output Folder for Linear Plot')
        if folder:
            self.linear_plot_output_dir_var.set(folder)

    def generate_linear_regression_plot(self):
        """Generate regression plots in single, individual, or combined mode."""
        if self.linear_plot_df is None or self.linear_plot_df.empty:
            messagebox.showwarning('No Data', 'Please load plotting data first.')
            return

        mode = self.linear_plot_mode_var.get().strip()
        x_col = self.linear_plot_x_var.get().strip()
        y_col = self.linear_plot_y_var.get().strip()
        sample_col = self.linear_plot_sample_var.get().strip()

        if not y_col or y_col not in self.linear_plot_df.columns:
            messagebox.showwarning('Missing Mapping', 'Please select a valid Y-axis column.')
            return

        if mode == 'Single metabolite (X column)':
            if not x_col or x_col not in self.linear_plot_df.columns:
                messagebox.showwarning('Missing Mapping', 'Please select a valid X-axis column.')
                return
            selected_metabolites = [x_col]
        else:
            selected_metabolites = [
                c for c in self._get_selected_listbox_values(self.linear_plot_metabolites_list)
                if c in self.linear_plot_df.columns
            ]
            if not selected_metabolites:
                messagebox.showwarning(
                    'Missing Metabolites',
                    'Select one or more metabolites in the list for individual/combined plotting modes.'
                )
                return

        output_dir = self.linear_plot_output_dir_var.get().strip()
        if not output_dir:
            messagebox.showwarning('Missing Output Folder', 'Please choose an output folder.')
            return

        try:
            jitter = float(self.linear_plot_jitter_var.get())
        except Exception:
            jitter = 0.0
        jitter = max(0.0, jitter)

        try:
            point_size = float(self.linear_plot_point_size_var.get())
        except Exception:
            point_size = 18.0
        point_size = max(1.0, point_size)

        try:
            point_alpha = float(self.linear_plot_alpha_var.get())
        except Exception:
            point_alpha = 0.75
        point_alpha = min(1.0, max(0.05, point_alpha))

        try:
            fig_w = float(self.linear_plot_fig_w_var.get())
            fig_h = float(self.linear_plot_fig_h_var.get())
        except Exception:
            fig_w, fig_h = 7.0, 5.0
        fig_w = max(3.0, fig_w)
        fig_h = max(3.0, fig_h)

        os.makedirs(output_dir, exist_ok=True)

        data = self.linear_plot_df.copy()
        if sample_col and sample_col in data.columns:
            data[sample_col] = data[sample_col].astype(str).str.strip()
            data = data[data[sample_col] != ''].drop_duplicates(subset=[sample_col], keep='first')

        try:
            import statsmodels.api as sm
        except Exception as e:
            messagebox.showerror('Missing Dependency', f'statsmodels is required:\n{e}')
            return

        y_raw = pd.to_numeric(data[y_col], errors='coerce')

        def _prepare_metabolite_result(metabolite_col: str):
            x_raw = pd.to_numeric(data[metabolite_col], errors='coerce')
            kept = x_raw.notna() & y_raw.notna()
            if bool(self.linear_plot_logx_var.get()):
                kept = kept & (x_raw > 0)

            x_vals = x_raw[kept].astype(float)
            y_vals = y_raw[kept].astype(float)
            if len(x_vals) < 4:
                return None

            x_arr = np.log10(x_vals.to_numpy()) if bool(self.linear_plot_logx_var.get()) else x_vals.to_numpy()
            y_arr = y_vals.to_numpy()
            design = sm.add_constant(x_arr, has_constant='add')
            model = sm.OLS(y_arr, design).fit()

            x_grid = np.linspace(float(np.min(x_arr)), float(np.max(x_arr)), 200)
            pred = model.get_prediction(sm.add_constant(x_grid, has_constant='add')).summary_frame(alpha=0.05)

            y_scatter = y_arr.copy()
            if jitter > 0:
                y_scatter = y_scatter + np.random.normal(loc=0.0, scale=jitter, size=len(y_scatter))

            slope_v = float(model.params[1]) if len(model.params) > 1 else np.nan
            pval_v = float(model.pvalues[1]) if len(model.pvalues) > 1 else np.nan
            r2_v = float(model.rsquared)
            n_obs_v = int(model.nobs)
            dropped_v = int((~kept).sum())

            return {
                'metabolite': metabolite_col,
                'x': x_arr,
                'y': y_scatter,
                'y_raw': y_arr,
                'x_grid': x_grid,
                'pred': pred,
                'slope': slope_v,
                'pval': pval_v,
                'r2': r2_v,
                'n': n_obs_v,
                'dropped': dropped_v,
            }

        fit_results = []
        for metabolite in selected_metabolites:
            res = _prepare_metabolite_result(metabolite)
            if res is None:
                self._linear_plot_log(f'Skipped {metabolite}: fewer than 4 valid paired values.')
                continue
            fit_results.append(res)

        if not fit_results:
            messagebox.showwarning('Insufficient Data', 'No selected metabolite had enough valid values to fit regression.')
            return

        def _spearman_rho(x_values, y_values):
            try:
                return float(pd.Series(x_values).corr(pd.Series(y_values), method='spearman'))
            except Exception:
                return float('nan')

        def _style_axis(ax):
            ax.tick_params(axis='both', labelsize=11, width=1.2)
            for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
                tick_label.set_fontweight('bold')
            ax.xaxis.label.set_size(13)
            ax.yaxis.label.set_size(13)
            ax.xaxis.label.set_fontweight('bold')
            ax.yaxis.label.set_fontweight('bold')
            ax.title.set_size(14)
            ax.title.set_fontweight('bold')

        def _annotate_axis(ax, res):
            rho = _spearman_rho(res['x'], res['y_raw'])
            ax.text(
                0.02,
                0.98,
                (
                    f'slope={res["slope"]:.4g}  '
                    f'p={res["pval"]:.3g}  '
                    f'$R^2$={res["r2"]:.3f}\n'
                    f'Spearman ρ={rho:.3f}'
                ),
                transform=ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.75, edgecolor='none')
            )

        prefix = self._sanitize_filename(self.linear_plot_prefix_var.get().strip() or 'regression_plot')
        y_name = self._sanitize_filename(y_col)
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        x_label_prefix = 'log10' if bool(self.linear_plot_logx_var.get()) else 'raw'
        saved_paths = []

        if mode == 'Selected metabolites - combined panel':
            n = len(fit_results)
            n_cols = min(3, n)
            n_rows = int(np.ceil(n / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w * n_cols, fig_h * n_rows))
            axes_arr = np.array(axes).reshape(-1)

            for idx, res in enumerate(fit_results):
                ax = axes_arr[idx]
                ax.scatter(res['x'], res['y'], s=point_size, alpha=point_alpha, color='#1f4e79', edgecolor='none')
                ax.plot(res['x_grid'], res['pred']['mean'].to_numpy(dtype=float), color='#c0392b', linewidth=2.0)
                if bool(self.linear_plot_ci_var.get()):
                    ax.fill_between(
                        res['x_grid'],
                        res['pred']['mean_ci_lower'].to_numpy(dtype=float),
                        res['pred']['mean_ci_upper'].to_numpy(dtype=float),
                        color='#c0392b',
                        alpha=0.2
                    )
                ax.set_xlabel(f'{x_label_prefix}({res["metabolite"]})' if bool(self.linear_plot_logx_var.get()) else res['metabolite'])
                ax.set_ylabel(y_col)
                ax.set_title(res['metabolite'])
                _style_axis(ax)
                ax.grid(alpha=0.2)
                _annotate_axis(ax, res)

            for j in range(len(fit_results), len(axes_arr)):
                axes_arr[j].axis('off')

            fig.suptitle(f'{y_col} vs Selected Metabolites (combined panel)', fontsize=13, fontweight='bold')
            fig.tight_layout(rect=[0, 0, 1, 0.98])

            out_base = os.path.join(output_dir, f'{prefix}_{y_name}_combined_{ts}')
            out_png = f'{out_base}.png'
            fig.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close(fig)
            saved_paths.append(out_png)

            self._linear_plot_log(f'Combined panel saved for {len(fit_results)} metabolite(s).')
        else:
            for res in fit_results:
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                ax.scatter(res['x'], res['y'], s=point_size, alpha=point_alpha, color='#1f4e79', edgecolor='none')
                ax.plot(res['x_grid'], res['pred']['mean'].to_numpy(dtype=float), color='#c0392b', linewidth=2.0, label='OLS fit')

                if bool(self.linear_plot_ci_var.get()):
                    ax.fill_between(
                        res['x_grid'],
                        res['pred']['mean_ci_lower'].to_numpy(dtype=float),
                        res['pred']['mean_ci_upper'].to_numpy(dtype=float),
                        color='#c0392b',
                        alpha=0.2,
                        label='95% CI'
                    )

                x_label = f'log10({res["metabolite"]})' if bool(self.linear_plot_logx_var.get()) else res['metabolite']
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_col)
                ax.set_title(f'{y_col} vs {res["metabolite"]}')
                _style_axis(ax)
                ax.grid(alpha=0.2)
                ax.legend(loc='best')
                _annotate_axis(ax, res)

                x_name = self._sanitize_filename(res['metabolite'])
                out_base = os.path.join(output_dir, f'{prefix}_{y_name}_vs_{x_name}_{ts}')
                out_png = f'{out_base}.png'
                fig.tight_layout()
                fig.savefig(out_png, dpi=300, bbox_inches='tight')
                plt.close(fig)
                saved_paths.append(out_png)

                self._linear_plot_log(
                    f'Generated {res["metabolite"]}: n={res["n"]}, dropped={res["dropped"]}, '
                    f'slope={res["slope"]:.4g}, p={res["pval"]:.3g}, R2={res["r2"]:.3f}'
                )

        shown_files = '\n'.join(saved_paths[:6])
        extra = '' if len(saved_paths) <= 6 else f'\n... plus {len(saved_paths) - 6} more file(s)'
        self._linear_plot_log(f'Saved {len(saved_paths)} file(s) to {output_dir}')
        messagebox.showinfo(
            'Plot Generated',
            f'Saved {len(saved_paths)} file(s):\n{shown_files}{extra}'
        )

    def _create_multiselect_listbox(self, parent, col_idx: int, title: str):
        """Create a labeled multi-select listbox in a grid column."""
        frame = tk.LabelFrame(parent, text=title, bg='#f0f0f0', font=('Arial', 9, 'bold'))
        frame.grid(row=0, column=col_idx, sticky='nsew', padx=4, pady=2)

        listbox = tk.Listbox(frame, selectmode='extended', exportselection=False, height=10)
        listbox.pack(side='left', fill='both', expand=True, padx=(5, 0), pady=5)
        scroll = ttk.Scrollbar(frame, orient='vertical', command=listbox.yview)
        scroll.pack(side='right', fill='y', padx=(0, 5), pady=5)
        listbox.configure(yscrollcommand=scroll.set)
        return listbox

    def _outcome_stats_log(self, message: str):
        """Append message to utility analysis log."""
        timestamp = datetime.now().strftime('%H:%M:%S')
        def _append():
            try:
                self.outcome_stats_log.insert('end', f'[{timestamp}] {message}\n')
                self.outcome_stats_log.see('end')
            except Exception:
                pass

        try:
            if threading.current_thread() is threading.main_thread():
                _append()
            else:
                self.root.after(0, _append)
        except Exception:
            pass

    def upload_outcome_stats_data(self):
        """Load tabular data for flexible outcome statistics."""
        file_path = filedialog.askopenfilename(
            title='Select Data File',
            filetypes=[('Excel/CSV', '*.xlsx *.xls *.csv *.tsv *.txt'), ('All Files', '*.*')]
        )
        if not file_path:
            return

        try:
            if file_path.lower().endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            elif file_path.lower().endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                df = pd.read_csv(file_path, sep='\t')

            self.outcome_stats_df = df
            self.outcome_stats_file_label.config(
                text=f"Loaded: {os.path.basename(file_path)} ({df.shape[0]} rows x {df.shape[1]} cols)",
                fg='green'
            )

            cols = list(df.columns)
            self.outcome_stats_all_columns = cols
            self.outcome_stats_sample_col_combo['values'] = cols
            self._populate_listbox(self.outcome_stats_categories_list, cols)
            self._populate_listbox(self.outcome_stats_covariates_list, cols)
            self._populate_listbox(self.outcome_stats_metabolites_list, cols)
            self._auto_detect_outcome_stats_columns()

            self._outcome_stats_log(f'Data loaded successfully from {os.path.basename(file_path)}')
        except Exception as e:
            messagebox.showerror('Load Error', f'Failed to load file:\n{e}')
            logger.error(f'Outcome stats data load failed: {e}', exc_info=True)

    def _populate_listbox(self, listbox: tk.Listbox, values: List[str]):
        listbox.delete(0, 'end')
        for val in values:
            listbox.insert('end', val)

    def _set_listbox_selection(self, listbox: tk.Listbox, all_values: List[str], selected_values: List[str]):
        listbox.selection_clear(0, 'end')
        selected_set = set(selected_values)
        for idx, name in enumerate(all_values):
            if name in selected_set:
                listbox.selection_set(idx)

    def _get_selected_listbox_values(self, listbox: tk.Listbox) -> List[str]:
        return [listbox.get(i) for i in listbox.curselection()]

    def _get_outcome_stats_assigned_non_metabolites(self) -> set:
        """Return columns currently assigned as sample/category/covariate."""
        assigned = set()

        sample_col = self.outcome_stats_sample_col_var.get().strip()
        if sample_col:
            assigned.add(sample_col)

        assigned.update(self._get_selected_listbox_values(self.outcome_stats_categories_list))
        assigned.update(self._get_selected_listbox_values(self.outcome_stats_covariates_list))
        return assigned

    def _refresh_outcome_metabolite_options(self, auto_select_remaining: bool = False):
        """Refresh metabolite options by removing columns assigned elsewhere."""
        if self.outcome_stats_df is None or self.outcome_stats_df.empty:
            return

        all_cols = self.outcome_stats_all_columns or list(self.outcome_stats_df.columns)
        assigned = self._get_outcome_stats_assigned_non_metabolites()

        current_selected = set(self._get_selected_listbox_values(self.outcome_stats_metabolites_list))
        available = [c for c in all_cols if c not in assigned]

        if auto_select_remaining:
            selected_after = available
        else:
            selected_after = [c for c in available if c in current_selected]

        self._populate_listbox(self.outcome_stats_metabolites_list, available)
        self._set_listbox_selection(self.outcome_stats_metabolites_list, available, selected_after)

    def _on_outcome_mapping_changed(self, _event=None):
        """Keep metabolite options synced when sample/category/covariate selections change."""
        self._refresh_outcome_metabolite_options(auto_select_remaining=False)

    def _auto_select_remaining_metabolites(self):
        """Select all currently available metabolite columns after exclusions."""
        self._refresh_outcome_metabolite_options(auto_select_remaining=True)
        count = len(self._get_selected_listbox_values(self.outcome_stats_metabolites_list))
        self._outcome_stats_log(f'Auto-selected {count} remaining metabolite columns.')

    def _auto_detect_outcome_stats_columns(self):
        """Auto-suggest sample/category/covariate/metabolite columns; user can adjust manually."""
        if self.outcome_stats_df is None or self.outcome_stats_df.empty:
            return

        df = self.outcome_stats_df
        cols = list(df.columns)

        # Sample ID guess
        normalized = {c: str(c).strip().casefold() for c in cols}
        sample_candidates = [
            c for c in cols
            if any(k in normalized[c] for k in ['sample', 'sampleid', 'sample_id', 'subject', 'id'])
        ]
        sample_col = sample_candidates[0] if sample_candidates else cols[0]
        self.outcome_stats_sample_col_var.set(sample_col)

        cat_cols = []
        covar_cols = []
        metabolite_cols = []

        for c in cols:
            if c == sample_col:
                continue
            series = df[c]
            non_na = series.dropna()
            if non_na.empty:
                continue

            numeric = pd.to_numeric(non_na, errors='coerce')
            numeric_ratio = float(numeric.notna().mean()) if len(non_na) else 0.0
            n_unique = int(non_na.nunique())

            # Outcome/category candidates: low-cardinality or explicit clinical terms
            lc = str(c).lower()
            is_clinical_name = any(k in lc for k in ['symptom', 'loc', 'amnesia', 'status', 'outcome', 'diagnosis'])

            if is_clinical_name or n_unique <= 12:
                cat_cols.append(c)
                # Numeric low-cardinality can also serve as covariate if user wants.
                if numeric_ratio > 0.95:
                    covar_cols.append(c)
                continue

            if numeric_ratio > 0.90:
                # Heuristic: high-cardinality numeric columns are metabolites.
                if n_unique > 15:
                    metabolite_cols.append(c)
                else:
                    covar_cols.append(c)
            else:
                # Non-numeric with higher-cardinality often covariate metadata.
                covar_cols.append(c)

        # Ensure disjoint defaults (used for log/suggestions only).
        covar_cols = [c for c in covar_cols if c not in cat_cols]
        metabolite_cols = [c for c in metabolite_cols if c not in cat_cols and c not in covar_cols]

        # IMPORTANT: Do not auto-select categories/covariates.
        # User selections must explicitly control what is treated as outcome/covariate.
        self._set_listbox_selection(self.outcome_stats_categories_list, cols, [])
        self._set_listbox_selection(self.outcome_stats_covariates_list, cols, [])

        # Start with all remaining columns (excluding sample/category/covariate assignments).
        self._refresh_outcome_metabolite_options(auto_select_remaining=True)

        self._outcome_stats_log(
            f'Auto-detected sample={sample_col}. Suggestions (not auto-selected): '
            f'categories={len(cat_cols)}, covariates={len(covar_cols)}. '
            f'Auto-selected remaining metabolites={len(self._get_selected_listbox_values(self.outcome_stats_metabolites_list))}.'
        )

    def _browse_outcome_stats_output_dir(self):
        folder = filedialog.askdirectory(title='Select Output Folder')
        if folder:
            self.outcome_stats_output_dir_var.set(folder)

    def _sanitize_filename(self, text: str) -> str:
        clean = ''.join(ch if ch.isalnum() or ch in ('-', '_') else '_' for ch in str(text))
        return clean.strip('_') or 'category'

    def _bh_adjust(self, pvals: np.ndarray) -> np.ndarray:
        """Benjamini-Hochberg adjustment."""
        pvals = np.asarray(pvals, dtype=float)
        n = len(pvals)
        if n == 0:
            return pvals
        order = np.argsort(pvals)
        ranked = pvals[order]
        adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
        adj = np.clip(adj, 0, 1)
        out = np.empty_like(adj)
        out[order] = adj
        return out

    def _normalize_binary_values(self, series: pd.Series):
        """Try mapping common binary text labels to 0/1. Returns mapped series or None."""
        cleaned = series.copy()
        cleaned = cleaned.where(cleaned.notna(), np.nan)
        cleaned = cleaned.astype(str).str.strip().str.casefold()
        cleaned = cleaned.where(cleaned != '', np.nan)
        unique_vals = sorted(v for v in cleaned.dropna().unique())
        if len(unique_vals) != 2:
            return None, None

        map_yes = {'yes', 'y', 'true', 'present', 'positive', '1', 'any'}
        map_no = {'no', 'n', 'false', 'absent', 'negative', '0', 'none'}

        mapping = {}
        for val in unique_vals:
            if val in map_yes:
                mapping[val] = 1.0
            elif val in map_no:
                mapping[val] = 0.0

        if len(mapping) < 2:
            # Fallback deterministic mapping by lexical order.
            mapping = {unique_vals[0]: 0.0, unique_vals[1]: 1.0}

        mapped = cleaned.map(mapping)
        return mapped, mapping

    def _encode_numeric_as_binary(self, numeric_series: pd.Series) -> Tuple[pd.Series, str]:
        """Encode numeric outcome as binary: min value -> 0, all others -> 1."""
        clean = numeric_series.dropna().copy()
        
        if len(clean) < 2:
            raise ValueError('Insufficient valid values for binary encoding')
        
        min_val = float(clean.min())
        max_val = float(clean.max())
        
        if min_val == max_val:
            raise ValueError('All values are identical; cannot create binary groups')
        
        # Create binary encoding: min -> 0, others -> 1
        encoded = pd.Series(index=numeric_series.index, dtype=int)
        encoded[numeric_series == min_val] = 0
        encoded[numeric_series > min_val] = 1
        
        encoding_note = f'Binary: {min_val} vs >{min_val}'
        return encoded, encoding_note

    def _cluster_numeric_outcome(self, numeric_series: pd.Series, mode: str) -> Tuple[pd.Series, Dict[str, Any]]:
        """Cluster a numeric outcome into 2/3 ordered groups using k-means."""
        try:
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_score
        except Exception as e:
            raise RuntimeError(f'scikit-learn is required for numeric clustering mode: {e}')

        clean = pd.to_numeric(numeric_series, errors='coerce').dropna()
        if clean.shape[0] < 6:
            raise ValueError('Need at least 6 numeric samples for clustering.')
        if int(clean.nunique()) < 2:
            raise ValueError('Numeric outcome has fewer than 2 unique values; clustering is not possible.')

        if mode == 'cluster2':
            k_candidates = [2]
        elif mode == 'cluster3':
            k_candidates = [3]
        else:
            k_candidates = [2, 3]

        values = clean.to_numpy(dtype=float).reshape(-1, 1)
        best_model = None
        best_k = None
        best_score = -np.inf

        for k in k_candidates:
            if values.shape[0] <= k or int(clean.nunique()) < k:
                continue

            model = KMeans(n_clusters=k, random_state=42, n_init=20)
            labels = model.fit_predict(values)
            if len(np.unique(labels)) < 2:
                continue

            score = silhouette_score(values, labels)
            if score > best_score:
                best_score = float(score)
                best_model = model
                best_k = k

        if best_model is None or best_k is None:
            raise ValueError('Could not derive stable clusters from numeric outcome values.')

        labels_raw = best_model.predict(values)
        centers = best_model.cluster_centers_.flatten()
        order = np.argsort(centers)
        rank_map = {int(cluster_id): rank for rank, cluster_id in enumerate(order)}

        if best_k == 2:
            rank_names = ['Low', 'High']
        elif best_k == 3:
            rank_names = ['Low', 'Mid', 'High']
        else:
            rank_names = [f'G{i+1}' for i in range(best_k)]

        label_series = pd.Series(index=clean.index, dtype=object)
        intervals: Dict[str, str] = {}
        center_map: Dict[str, float] = {}

        for cluster_id in range(best_k):
            rank = rank_map[cluster_id]
            grp_name = rank_names[rank]
            mask = labels_raw == cluster_id
            cluster_vals = clean.iloc[np.where(mask)[0]]
            label_series.loc[cluster_vals.index] = grp_name

            center_map[grp_name] = float(centers[cluster_id])
            intervals[grp_name] = f"{float(cluster_vals.min()):.3g}-{float(cluster_vals.max()):.3g}"

        details = {
            'k': int(best_k),
            'silhouette': best_score,
            'intervals': intervals,
            'centers': center_map,
            'ordered_groups': rank_names,
        }
        return label_series, details

    def _run_continuous_outcome_ols(
        self,
        data_df: pd.DataFrame,
        sample_col: str,
        outcome_col: str,
        metabolite_cols: List[str],
        covariate_cols: List[str],
        alpha: float,
        p_adjust_method: str,
        preencoded_outcome: Optional[pd.Series] = None,
        encoding_note: Optional[str] = None
    ) -> pd.DataFrame:
        """Run outcome ~ metabolite (+ covariates) OLS for each metabolite."""
        try:
            import statsmodels.api as sm
        except Exception as e:
            raise RuntimeError(f'statsmodels is required for OLS outcome analysis: {e}')

        working = data_df.copy()
        working[sample_col] = working[sample_col].astype(str).str.strip()
        working = working[working[sample_col] != '']
        working = working.drop_duplicates(subset=[sample_col], keep='first')
        working = working.set_index(sample_col)

        # Outcome preparation: use provided encoding when available; otherwise auto-detect.
        outcome_raw = working[outcome_col]
        if preencoded_outcome is not None:
            y = preencoded_outcome.reindex(working.index)
            y = pd.to_numeric(y, errors='coerce')
            mapping_note = encoding_note or 'preencoded'
        else:
            y = pd.to_numeric(outcome_raw, errors='coerce')
            mapping_note = 'numeric'

            if y.notna().sum() < 3:
                mapped, mapping = self._normalize_binary_values(outcome_raw)
                if mapped is not None:
                    y = mapped
                    mapping_note = f'binary mapped {mapping}'

        # Build numeric covariate matrix (categorical covariates are one-hot encoded).
        cov_df = pd.DataFrame(index=working.index)
        for cov in covariate_cols:
            if cov not in working.columns:
                continue
            cov_series = working[cov]
            numeric_cov = pd.to_numeric(cov_series, errors='coerce')
            if numeric_cov.notna().mean() > 0.9:
                cov_df[cov] = numeric_cov
            else:
                dummies = pd.get_dummies(cov_series.astype(str), prefix=cov, drop_first=True, dtype=float)
                for dcol in dummies.columns:
                    cov_df[dcol] = dummies[dcol]

        rows = []
        for metabolite in metabolite_cols:
            if metabolite not in working.columns:
                continue

            x_met = pd.to_numeric(working[metabolite], errors='coerce')
            model_df = pd.DataFrame({'y': y, 'metabolite': x_met}, index=working.index)
            if not cov_df.empty:
                model_df = pd.concat([model_df, cov_df], axis=1)
            model_df = model_df.dropna(axis=0)

            if len(model_df) < max(4, model_df.shape[1] + 1):
                continue

            try:
                X = sm.add_constant(model_df.drop(columns=['y']), has_constant='add')
                model = sm.OLS(model_df['y'], X).fit()
                if 'metabolite' not in model.params.index:
                    continue

                ci = model.conf_int(alpha=0.05)
                rows.append({
                    'Metabolite': metabolite,
                    'Outcome': outcome_col,
                    'Outcome_Encoding': mapping_note,
                    'n_samples': int(model.nobs),
                    'beta': float(model.params['metabolite']),
                    'std_err': float(model.bse['metabolite']),
                    't_statistic': float(model.tvalues['metabolite']),
                    'pvalue': float(model.pvalues['metabolite']),
                    'ci_lower_95': float(ci.loc['metabolite', 0]),
                    'ci_upper_95': float(ci.loc['metabolite', 1]),
                    'r_squared': float(model.rsquared),
                    'adj_r_squared': float(model.rsquared_adj)
                })
            except Exception:
                continue

        result_df = pd.DataFrame(rows)
        if result_df.empty:
            return result_df

        pvals = result_df['pvalue'].to_numpy(dtype=float)
        if p_adjust_method == 'None':
            result_df['adj_p'] = pvals
        elif p_adjust_method == 'Bonferroni':
            result_df['adj_p'] = np.minimum(pvals * len(pvals), 1.0)
        else:
            result_df['adj_p'] = self._bh_adjust(pvals)

        result_df['neg_log10_adj_p'] = -np.log10(np.maximum(result_df['adj_p'], np.finfo(float).eps))
        result_df['significant'] = result_df['adj_p'] < alpha
        return result_df.sort_values('adj_p', ascending=True)

    def _run_ordinal_outcome_regression(
        self,
        data_df: pd.DataFrame,
        sample_col: str,
        outcome_col: str,
        metabolite_cols: List[str],
        covariate_cols: List[str],
        alpha: float,
        p_adjust_method: str,
        preencoded_outcome: Optional[pd.Series] = None,
        encoding_note: Optional[str] = None,
        progress_callback: Optional[Any] = None,
        progress_every: int = 25,
        maxiter: int = 120
    ) -> pd.DataFrame:
        """Run ordinal logistic regression: outcome ~ metabolite (+ covariates) for each metabolite."""
        try:
            import statsmodels.api as sm
            from statsmodels.miscmodels.ordinal_model import OrderedModel
        except Exception as e:
            raise RuntimeError(f'statsmodels is required for ordinal regression: {e}')

        working = data_df.copy()
        working[sample_col] = working[sample_col].astype(str).str.strip()
        working = working[working[sample_col] != '']
        working = working.drop_duplicates(subset=[sample_col], keep='first')
        working = working.set_index(sample_col)

        # Outcome preparation
        outcome_raw = working[outcome_col]
        if preencoded_outcome is not None:
            y = preencoded_outcome.reindex(working.index)
            y = pd.to_numeric(y, errors='coerce')
        else:
            y = pd.to_numeric(outcome_raw, errors='coerce')

            if y.notna().sum() < 3:
                mapped, mapping = self._normalize_binary_values(outcome_raw)
                if mapped is not None:
                    y = mapped

        if y.notna().sum() < 3:
            raise ValueError('Outcome has fewer than 3 valid values for ordinal regression')

        # Ensure y is categorical with ordered levels
        y_valid = y.dropna()
        if len(y_valid) < 3:
            raise ValueError('Insufficient valid outcome values')

        unique_vals = sorted(y_valid.unique())
        if len(unique_vals) < 2:
            raise ValueError('Outcome must have at least 2 unique values')

        if progress_callback is not None:
            try:
                progress_callback(
                    f'{outcome_col}: ordinal setup complete - n={len(y_valid)}, '
                    f'levels={len(unique_vals)}, metabolites={len(metabolite_cols)}'
                )
            except Exception:
                pass

        y_cat = pd.Categorical(y, categories=unique_vals, ordered=True)

        # Build metabolite data
        met_data = working[metabolite_cols].apply(pd.to_numeric, errors='coerce')

        # Build covariate matrix
        cov_df = pd.DataFrame(index=working.index)
        for cov in covariate_cols:
            if cov not in working.columns:
                continue
            cov_series = working[cov]
            numeric_cov = pd.to_numeric(cov_series, errors='coerce')
            if numeric_cov.notna().sum() > 0:
                numeric_cov = numeric_cov.fillna(numeric_cov.mean())
                cov_df[cov] = numeric_cov
            else:
                cat_cov = pd.Categorical(cov_series).codes
                cov_df[cov] = cat_cov

        results = []

        total_mets = len(metabolite_cols)
        for idx, metabolite in enumerate(metabolite_cols, start=1):
            if progress_callback is not None and (idx == 1 or idx % max(1, progress_every) == 0 or idx == total_mets):
                try:
                    progress_callback(f'{outcome_col}: ordinal progress {idx}/{total_mets}')
                except Exception:
                    pass
            try:
                x_raw = met_data[metabolite].copy()
                valid = y_cat.notna() & x_raw.notna()

                if valid.sum() < 4:
                    results.append({
                        'metabolite': metabolite,
                        'estimate': np.nan,
                        'std_err': np.nan,
                        'z_value': np.nan,
                        'p': 1.0,
                        'adj_p': 1.0,
                        'n_samples': int(valid.sum()),
                        'note': 'Too few valid samples'
                    })
                    continue

                x = x_raw[valid].copy()
                y_subset = y_cat[valid]

                # Standardize metabolite
                x_mean = x.mean()
                x_std = x.std()
                if x_std > 0:
                    x = (x - x_mean) / x_std
                else:
                    results.append({
                        'metabolite': metabolite,
                        'estimate': np.nan,
                        'std_err': np.nan,
                        'z_value': np.nan,
                        'p': 1.0,
                        'adj_p': 1.0,
                        'n_samples': int(valid.sum()),
                        'note': 'Zero variance'
                    })
                    continue

                # Prepare design matrix
                design = pd.DataFrame({'metabolite': x}, index=x.index)

                if not cov_df.empty and len(covariate_cols) > 0:
                    cov_subset = cov_df.loc[valid, covariate_cols]
                    for col in cov_subset.columns:
                        if cov_subset[col].notna().sum() > 0:
                            col_data = cov_subset[col].fillna(cov_subset[col].mean())
                            col_std = col_data.std()
                            if col_std > 0:
                                design[col] = (col_data - col_data.mean()) / col_std
                            else:
                                design[col] = col_data

                # Fit ordinal regression model
                try:
                    model = OrderedModel(y_subset, design, distr='logit')
                    result = model.fit(method='lbfgs', maxiter=maxiter, disp=False)

                    # Extract coefficient for metabolite
                    est = float(result.params.get('metabolite', np.nan))
                    se = float(result.bse.get('metabolite', np.nan))
                    if not np.isnan(se) and se > 0:
                        z_val = est / se
                        from scipy.stats import norm
                        pval = 2 * (1 - norm.cdf(abs(z_val)))
                    else:
                        z_val = np.nan
                        pval = 1.0

                    results.append({
                        'metabolite': metabolite,
                        'estimate': est,
                        'std_err': se,
                        'z_value': z_val,
                        'p': pval,
                        'adj_p': pval,
                        'n_samples': int(valid.sum()),
                        'note': 'Ordinal logistic'
                    })
                except Exception as fit_err:
                    results.append({
                        'metabolite': metabolite,
                        'estimate': np.nan,
                        'std_err': np.nan,
                        'z_value': np.nan,
                        'p': 1.0,
                        'adj_p': 1.0,
                        'n_samples': int(valid.sum()),
                        'note': f'Model fit failed: {str(fit_err)[:50]}'
                    })

            except Exception as e:
                results.append({
                    'metabolite': metabolite,
                    'estimate': np.nan,
                    'std_err': np.nan,
                    'z_value': np.nan,
                    'p': 1.0,
                    'adj_p': 1.0,
                    'n_samples': 0,
                    'note': f'Error: {str(e)[:50]}'
                })

        result_df = pd.DataFrame(results)

        if len(result_df) == 0:
            raise ValueError('No valid models produced')

        # Adjust p-values
        pvals = result_df['p'].to_numpy()
        if p_adjust_method == 'BH':
            result_df['adj_p'] = self._bh_adjust(pvals)
        elif p_adjust_method == 'Bonferroni':
            result_df['adj_p'] = np.minimum(pvals * len(pvals), 1.0)
        elif p_adjust_method == 'Holm':
            from scipy import stats
            result_df['adj_p'] = stats.holm(pvals, alpha=1.0)
        elif p_adjust_method == 'Hochberg':
            from scipy import stats
            result_df['adj_p'] = stats.hochberg(pvals, alpha=1.0)
        elif p_adjust_method == 'BY':
            result_df['adj_p'] = self._by_adjust(pvals)
        else:
            result_df['adj_p'] = pvals

        result_df['neg_log10_adj_p'] = -np.log10(np.maximum(result_df['adj_p'], np.finfo(float).eps))
        result_df['significant'] = result_df['adj_p'] < alpha

        return result_df.sort_values('adj_p', ascending=True)

    def _plot_ordinal_regression_results(
        self,
        data_df: pd.DataFrame,
        outcome_col: str,
        metabolite_cols: List[str],
        output_dir: str,
        significant_metabolites: Optional[List[str]] = None,
        n_top: int = 6
    ) -> List[str]:
        """Create boxplots for ordinal regression results showing metabolite distributions across outcome categories."""
        saved_paths = []

        try:
            data_df = data_df.copy()
            outcome_values = pd.to_numeric(data_df[outcome_col], errors='coerce')
            y_unique = sorted([v for v in outcome_values.unique() if pd.notna(v)])

            if len(y_unique) < 2:
                return saved_paths

            # Select metabolites to plot
            if significant_metabolites and len(significant_metabolites) > 0:
                plot_mets = [m for m in significant_metabolites if m in metabolite_cols][:n_top]
            else:
                plot_mets = metabolite_cols[:n_top]

            if not plot_mets:
                return saved_paths

            # Create grid of subplots
            n = len(plot_mets)
            n_cols = min(3, n)
            n_rows = int(np.ceil(n / n_cols))

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.5 * n_rows))
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1 or n_cols == 1:
                axes = np.array(axes).reshape(n_rows, n_cols)

            ts = datetime.now().strftime('%Y%m%d_%H%M%S')
            prefix = self._sanitize_filename(outcome_col)

            for idx, met in enumerate(plot_mets):
                row = idx // n_cols
                col = idx % n_cols
                ax = axes[row, col]

                met_vals = pd.to_numeric(data_df[met], errors='coerce')
                valid = outcome_values.notna() & met_vals.notna()

                if valid.sum() < 2:
                    ax.text(0.5, 0.5, f'{met}\n(no data)', ha='center', va='center', fontsize=9, color='red')
                    ax.axis('off')
                    continue

                plot_data = []
                labels = []
                for y_val in y_unique:
                    mask = (outcome_values == y_val) & valid
                    if mask.sum() > 0:
                        plot_data.append(met_vals[mask].to_numpy())
                        labels.append(str(int(y_val)) if float(y_val).is_integer() else f'{y_val:.2f}')

                if len(plot_data) < 2:
                    ax.text(0.5, 0.5, f'{met}\n(insufficient groups)', ha='center', va='center', fontsize=9, color='orange')
                    ax.axis('off')
                    continue

                bp = ax.boxplot(plot_data, labels=labels, patch_artist=True, showmeans=True)

                for patch in bp['boxes']:
                    patch.set_facecolor('#3498db')
                    patch.set_alpha(0.7)

                for whisker in bp['whiskers']:
                    whisker.set(color='#2c3e50', linewidth=1)

                for median in bp['medians']:
                    median.set(color='#c0392b', linewidth=2)

                for mean in bp['means']:
                    mean.set(marker='D', markerfacecolor='#2ecc71', markeredgecolor='#27ae60', markersize=6)

                ax.set_xlabel(outcome_col, fontsize=9)
                ax.set_ylabel(met, fontsize=9)
                ax.set_title(met, fontsize=10, fontweight='bold')
                ax.grid(axis='y', alpha=0.3)

            # Hide unused subplots
            for idx in range(len(plot_mets), n_rows * n_cols):
                row = idx // n_cols
                col = idx % n_cols
                axes[row, col].axis('off')

            fig.suptitle(f'Ordinal Regression: {outcome_col}', fontsize=12, fontweight='bold')
            fig.tight_layout(rect=[0, 0, 1, 0.97])

            out_base = os.path.join(output_dir, f'ordinal_regression_{prefix}_{ts}')
            out_png = f'{out_base}.png'
            fig.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close(fig)
            saved_paths.append(out_png)

        except Exception as e:
            logger.warning(f'Could not generate ordinal regression plots: {e}')

        return saved_paths

    def run_outcome_stats_analysis(self):
        """Run per-category analyses and export one workbook per category in background."""
        if self.outcome_stats_df is None or self.outcome_stats_df.empty:
            messagebox.showwarning('No Data', 'Please upload data first.')
            return

        sample_col = self.outcome_stats_sample_col_var.get().strip()
        if not sample_col or sample_col not in self.outcome_stats_df.columns:
            messagebox.showwarning('Missing Mapping', 'Please select a valid Sample ID column.')
            return

        category_cols = self._get_selected_listbox_values(self.outcome_stats_categories_list)
        covariate_cols = self._get_selected_listbox_values(self.outcome_stats_covariates_list)
        metabolite_cols = self._get_selected_listbox_values(self.outcome_stats_metabolites_list)

        # Ensure disjoint sets and valid columns.
        category_cols = [c for c in category_cols if c in self.outcome_stats_df.columns and c != sample_col]
        covariate_cols = [c for c in covariate_cols if c in self.outcome_stats_df.columns and c not in category_cols and c != sample_col]
        metabolite_cols = [
            c for c in metabolite_cols
            if c in self.outcome_stats_df.columns and c not in category_cols and c not in covariate_cols and c != sample_col
        ]

        if not category_cols:
            messagebox.showwarning('Missing Mapping', 'Select at least one outcome category column.')
            return
        if not metabolite_cols:
            messagebox.showwarning('Missing Mapping', 'Select at least one metabolite value column.')
            return

        output_dir = self.outcome_stats_output_dir_var.get().strip()
        if not output_dir:
            messagebox.showwarning('Missing Output Folder', 'Please choose an output folder.')
            return
        os.makedirs(output_dir, exist_ok=True)

        method = self.outcome_stats_method_var.get().strip()
        p_adjust_method = self.outcome_stats_padj_var.get().strip()
        strict_var = getattr(self, 'outcome_stats_strict_var', None)
        strict_mode = bool(strict_var.get()) if strict_var is not None else True
        unknown_var = getattr(self, 'outcome_stats_unknown_tokens_var', None)
        unknown_text = unknown_var.get() if unknown_var is not None else 'unknown'
        unknown_tokens = {
            token.strip().casefold()
            for token in str(unknown_text).split(',')
            if token.strip()
        }
        if strict_mode:
            unknown_tokens.update({'unknown'})
        try:
            alpha = float(self.outcome_stats_alpha_var.get())
        except Exception:
            alpha = 0.05

        if hasattr(self, 'outcome_stats_run_btn'):
            try:
                self.outcome_stats_run_btn.config(state='disabled', text='Running Analysis...')
            except Exception:
                pass

        # Copy data/parameters now (main thread) so worker does not read Tk variables.
        params = {
            'df': self.outcome_stats_df.copy(),
            'sample_col': sample_col,
            'category_cols': category_cols,
            'covariate_cols': covariate_cols,
            'metabolite_cols': metabolite_cols,
            'output_dir': output_dir,
            'method': method,
            'p_adjust_method': p_adjust_method,
            'strict_mode': strict_mode,
            'unknown_tokens': unknown_tokens,
            'alpha': alpha,
            'numeric_mode': (getattr(self, 'outcome_stats_numeric_mode_var', tk.StringVar(value='Continuous (OLS)')).get().strip()),
        }

        self._outcome_stats_log('Launching analysis in background thread...')
        threading.Thread(
            target=self._run_outcome_stats_analysis_worker,
            args=(params,),
            daemon=True
        ).start()

    def _run_outcome_stats_analysis_worker(self, params: Dict[str, Any]):
        """Background worker for regression utility analysis/export."""
        try:
            df = params['df']
            sample_col = params['sample_col']
            category_cols = params['category_cols']
            covariate_cols = params['covariate_cols']
            metabolite_cols = params['metabolite_cols']
            output_dir = params['output_dir']
            method = params['method']
            p_adjust_method = params['p_adjust_method']
            strict_mode = params['strict_mode']
            unknown_tokens = params['unknown_tokens']
            alpha = params['alpha']
            numeric_mode_raw = str(params.get('numeric_mode', 'Continuous (OLS)'))

            if numeric_mode_raw == 'Binary (0 vs 1)':
                numeric_mode = 'binary'
            elif numeric_mode_raw == 'Cluster to 2 groups':
                numeric_mode = 'cluster2'
            elif numeric_mode_raw == 'Cluster to 3 groups':
                numeric_mode = 'cluster3'
            elif numeric_mode_raw == 'Auto-cluster (2 or 3)':
                numeric_mode = 'cluster_auto'
            else:
                numeric_mode = 'continuous'

            self._outcome_stats_log('Starting per-category analysis...')
            self._outcome_stats_log(f'Sample column: {sample_col}')
            self._outcome_stats_log(f'Categories: {category_cols}')
            self._outcome_stats_log(f'Covariates: {covariate_cols if covariate_cols else "None"}')
            self._outcome_stats_log(f'Metabolite columns: {len(metabolite_cols)} selected')
            self._outcome_stats_log(f'Numeric outcome mode: {numeric_mode_raw}')
            self._outcome_stats_log(
                f'Strict mode: {"ON" if strict_mode else "OFF"}; '
                f'unknown tokens: {sorted(unknown_tokens) if unknown_tokens else "None"}'
            )

            from main_script.covariate_adjustment import (
                run_covariate_adjusted_analysis,
                run_limma_covariate_analysis,
                export_covariate_results,
                CovariateAnalysisResult
            )

            exported_files = []
            ts = datetime.now().strftime('%Y%m%d_%H%M%S')

            working = df.copy()
            working[sample_col] = working[sample_col].astype(str).str.strip()
            working = working[working[sample_col] != ''].copy()
            if working[sample_col].duplicated().any():
                dup_count = int(working[sample_col].duplicated().sum())
                self._outcome_stats_log(f'Warning: {dup_count} duplicate sample IDs found; using first occurrence.')
                working = working.drop_duplicates(subset=[sample_col], keep='first')
            working = working.set_index(sample_col)

            for category in category_cols:
                try:
                    self._outcome_stats_log(f'Analyzing category: {category}')

                    cat_series = working[category].copy()
                    cat_clean = cat_series.astype(str).str.strip()
                    cat_clean = cat_clean.where(cat_series.notna(), np.nan)

                    dropped_unknown = 0
                    if strict_mode and unknown_tokens:
                        unknown_mask = cat_clean.notna() & cat_clean.str.casefold().isin(unknown_tokens)
                        dropped_unknown = int(unknown_mask.sum())
                        cat_clean = cat_clean.where(~unknown_mask, np.nan)
                        if dropped_unknown > 0:
                            self._outcome_stats_log(
                                f'{category}: dropped {dropped_unknown} samples matching unknown tokens.'
                            )

                    valid_samples = cat_clean[cat_clean.notna()].index.tolist()
                    if len(valid_samples) < 4:
                        self._outcome_stats_log(f'Skipped {category}: too few samples with values.')
                        continue

                    cat_non_na = cat_clean.loc[valid_samples]
                    numeric_cat = pd.to_numeric(cat_non_na, errors='coerce')
                    numeric_ratio = float(numeric_cat.notna().mean()) if len(cat_non_na) else 0.0
                    n_unique = int(cat_non_na.nunique())
                    max_ordinal_levels = 15

                    out_name = self._sanitize_filename(category)
                    out_file = os.path.join(output_dir, f'{out_name}_{ts}.xlsx')
                    cluster_assignment_df = None

                    encoded_binary = None
                    encoding_note = None
                    if strict_mode and numeric_ratio <= 0.90 and n_unique == 2:
                        encoded_binary, mapping = self._normalize_binary_values(cat_non_na)
                        if encoded_binary is not None:
                            encoded_binary = encoded_binary.reindex(working.index)
                            encoding_note = f'strict binary mapped {mapping}'
                            self._outcome_stats_log(
                                f'{category}: strict binary encoding applied ({mapping}).'
                            )

                    is_numeric_outcome = (numeric_ratio > 0.90 and n_unique >= 2) or (encoded_binary is not None)
                    is_continuous_numeric_for_ordinal = (
                        numeric_ratio > 0.90 and n_unique >= 3 and encoded_binary is None and n_unique <= max_ordinal_levels
                    )

                    if method == 'Ordinal Regression' and numeric_mode == 'continuous' and is_continuous_numeric_for_ordinal:
                        try:
                            self._outcome_stats_log(f'Running ordinal regression for {category}...')

                            reg_df = self._run_ordinal_outcome_regression(
                                data_df=working.reset_index(),
                                sample_col=sample_col,
                                outcome_col=category,
                                metabolite_cols=metabolite_cols,
                                covariate_cols=[c for c in covariate_cols if c != category],
                                alpha=alpha,
                                p_adjust_method=p_adjust_method,
                                preencoded_outcome=encoded_binary,
                                encoding_note=encoding_note,
                                progress_callback=self._outcome_stats_log,
                                progress_every=25,
                                maxiter=120
                            )

                            if reg_df.empty:
                                self._outcome_stats_log(f'Skipped {category}: no valid models after filtering.')
                                continue

                            summary_df = pd.DataFrame([
                                ('category', category),
                                ('analysis_type', 'ordinal_logistic_regression'),
                                ('method', 'Ordinal Regression'),
                                ('n_samples', len(valid_samples)),
                                ('n_unknown_dropped', dropped_unknown),
                                ('n_metabolites_tested', len(reg_df)),
                                ('n_significant', int((reg_df['adj_p'] < alpha).sum())),
                                ('covariates', ', '.join([c for c in covariate_cols if c != category]) if covariate_cols else 'None'),
                                ('strict_mode', strict_mode),
                                ('alpha', alpha),
                                ('p_adjust_method', p_adjust_method)
                            ], columns=['Parameter', 'Value'])

                            with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                                reg_df.to_excel(writer, sheet_name='Results', index=False)
                                summary_df.to_excel(writer, sheet_name='Summary', index=False)

                            # Generate plots for significant metabolites
                            sig_mets = reg_df[reg_df['adj_p'] < alpha]['metabolite'].tolist()
                            plot_paths = self._plot_ordinal_regression_results(
                                data_df=working.reset_index(),
                                outcome_col=category,
                                metabolite_cols=metabolite_cols,
                                output_dir=output_dir,
                                significant_metabolites=sig_mets,
                                n_top=6
                            )
                            exported_files.extend(plot_paths)

                            exported_files.append(out_file)
                            self._outcome_stats_log(
                                f'Exported {category} (Ordinal Regression, {int((reg_df["adj_p"] < alpha).sum())} significant) -> {os.path.basename(out_file)}'
                            )
                            continue

                        except Exception as ord_e:
                            self._outcome_stats_log(
                                f'Warning: Ordinal regression failed for {category} ({str(ord_e)[:80]}); '
                                f'falling back to OLS outcome regression'
                            )
                            logger.warning(f'Ordinal regression failed: {ord_e}', exc_info=True)

                    if is_numeric_outcome and numeric_mode == 'continuous':
                        if method == 'Limma':
                            self._outcome_stats_log(
                                f'Notice: {category} detected as numeric outcome. Using OLS outcome regression '
                                f'(Limma applies to group-comparison workflow).'
                            )
                        elif method == 'Ordinal Regression' and not is_continuous_numeric_for_ordinal:
                            if n_unique > max_ordinal_levels:
                                self._outcome_stats_log(
                                    f'Notice: Ordinal regression skipped for {category}: '
                                    f'{n_unique} outcome levels (> {max_ordinal_levels}) are near-continuous. '
                                    f'Using OLS outcome regression as optimal fallback.'
                                )
                            else:
                                self._outcome_stats_log(
                                    f'Notice: Ordinal regression skipped for {category} because outcome is not '
                                    f'continuous numeric (likely binary-like). Using OLS outcome regression.'
                                )

                        reg_df = self._run_continuous_outcome_ols(
                            data_df=working.reset_index(),
                            sample_col=sample_col,
                            outcome_col=category,
                            metabolite_cols=metabolite_cols,
                            covariate_cols=[c for c in covariate_cols if c != category],
                            alpha=alpha,
                            p_adjust_method=p_adjust_method,
                            preencoded_outcome=encoded_binary,
                            encoding_note=encoding_note
                        )

                        if reg_df.empty:
                            self._outcome_stats_log(f'Skipped {category}: no valid models after filtering.')
                            continue

                        summary_df = pd.DataFrame([
                            ('category', category),
                            ('analysis_type', 'continuous_outcome_ols'),
                            ('method', 'Linear Model (OLS)'),
                            ('n_samples', len(valid_samples)),
                            ('n_unknown_dropped', dropped_unknown),
                            ('n_metabolites_tested', len(reg_df)),
                            ('n_significant', int((reg_df['adj_p'] < alpha).sum())),
                            ('covariates', ', '.join([c for c in covariate_cols if c != category]) if covariate_cols else 'None'),
                            ('strict_mode', strict_mode),
                            ('alpha', alpha),
                            ('p_adjust_method', p_adjust_method)
                        ], columns=['Parameter', 'Value'])

                        with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                            reg_df.to_excel(writer, sheet_name='Results', index=False)
                            summary_df.to_excel(writer, sheet_name='Summary', index=False)

                        exported_files.append(out_file)
                        self._outcome_stats_log(f'Exported {category} -> {os.path.basename(out_file)}')
                        continue

                    if is_numeric_outcome and numeric_mode in {'cluster2', 'cluster3', 'cluster_auto'} and encoded_binary is None:
                        try:
                            numeric_for_cluster = pd.to_numeric(cat_non_na, errors='coerce')
                            clustered_labels, cluster_info = self._cluster_numeric_outcome(
                                numeric_for_cluster,
                                numeric_mode
                            )
                            cat_non_na = clustered_labels.dropna()
                            valid_samples = cat_non_na.index.tolist()

                            cluster_assignment_df = pd.DataFrame({
                                'Sample': valid_samples,
                                'Original_Numeric_Outcome': [float(numeric_for_cluster.loc[s]) for s in valid_samples],
                                'Cluster_Label': [str(cat_non_na.loc[s]) for s in valid_samples],
                                'Cluster_Range': [cluster_info['intervals'].get(str(cat_non_na.loc[s]), '') for s in valid_samples],
                                'Cluster_Center': [cluster_info['centers'].get(str(cat_non_na.loc[s]), np.nan) for s in valid_samples],
                                'Cluster_K': int(cluster_info['k']),
                                'Silhouette_Score': float(cluster_info['silhouette']),
                            })

                            grp_txt = ', '.join(
                                [f"{g}~[{cluster_info['intervals'].get(g, '?')}]" for g in cluster_info['ordered_groups'] if g in set(cat_non_na.values)]
                            )
                            self._outcome_stats_log(
                                f"{category}: clustered numeric outcome into {cluster_info['k']} groups "
                                f"(silhouette={cluster_info['silhouette']:.3f}): {grp_txt}"
                            )
                        except Exception as ce:
                            fallback_reason = str(ce)
                            self._outcome_stats_log(
                                f"{category}: clustering failed ({fallback_reason}); falling back to continuous OLS mode."
                            )

                            reg_df = self._run_continuous_outcome_ols(
                                data_df=working.reset_index(),
                                sample_col=sample_col,
                                outcome_col=category,
                                metabolite_cols=metabolite_cols,
                                covariate_cols=[c for c in covariate_cols if c != category],
                                alpha=alpha,
                                p_adjust_method=p_adjust_method,
                                preencoded_outcome=encoded_binary,
                                encoding_note=encoding_note
                            )

                            if reg_df.empty:
                                self._outcome_stats_log(f'Skipped {category}: no valid models after filtering.')
                                continue

                            summary_df = pd.DataFrame([
                                ('category', category),
                                ('analysis_type', 'continuous_outcome_ols_fallback'),
                                ('method', 'Linear Model (OLS)'),
                                ('requested_numeric_mode', numeric_mode_raw),
                                ('n_samples', len(valid_samples)),
                                ('n_unknown_dropped', dropped_unknown),
                                ('n_metabolites_tested', len(reg_df)),
                                ('n_significant', int((reg_df['adj_p'] < alpha).sum())),
                                ('covariates', ', '.join([c for c in covariate_cols if c != category]) if covariate_cols else 'None'),
                                ('strict_mode', strict_mode),
                                ('alpha', alpha),
                                ('p_adjust_method', p_adjust_method),
                                ('fallback_reason', fallback_reason)
                            ], columns=['Parameter', 'Value'])

                            fallback_assignment_df = pd.DataFrame({
                                'Sample': valid_samples,
                                'Original_Numeric_Outcome': [pd.to_numeric(cat_non_na.loc[s], errors='coerce') for s in valid_samples],
                                'Cluster_Label': ['N/A'] * len(valid_samples),
                                'Cluster_Range': ['N/A'] * len(valid_samples),
                                'Cluster_Center': [np.nan] * len(valid_samples),
                                'Cluster_K': [np.nan] * len(valid_samples),
                                'Silhouette_Score': [np.nan] * len(valid_samples),
                                'Status': ['Clustering_Failed_Fallback_to_OLS'] * len(valid_samples),
                                'Fallback_Reason': [fallback_reason] * len(valid_samples),
                            })

                            with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                                reg_df.to_excel(writer, sheet_name='Results', index=False)
                                summary_df.to_excel(writer, sheet_name='Summary', index=False)
                                fallback_assignment_df.to_excel(writer, sheet_name='Outcome_Cluster_Assignments', index=False)

                            exported_files.append(out_file)
                            self._outcome_stats_log(f'Exported {category} -> {os.path.basename(out_file)}')
                            continue

                    if is_numeric_outcome and numeric_mode == 'binary' and encoded_binary is None:
                        try:
                            numeric_for_binary = pd.to_numeric(cat_non_na, errors='coerce')
                            binary_encoded, binary_note = self._encode_numeric_as_binary(numeric_for_binary)
                            cat_non_na = binary_encoded.dropna()
                            valid_samples = cat_non_na.index.tolist()

                            self._outcome_stats_log(
                                f"{category}: numeric outcome encoded as binary ({binary_note})"
                            )
                        except Exception as be:
                            self._outcome_stats_log(
                                f"{category}: binary encoding failed ({str(be)[:80]}); "
                                f"falling back to continuous OLS mode."
                            )

                            reg_df = self._run_continuous_outcome_ols(
                                data_df=working.reset_index(),
                                sample_col=sample_col,
                                outcome_col=category,
                                metabolite_cols=metabolite_cols,
                                covariate_cols=[c for c in covariate_cols if c != category],
                                alpha=alpha,
                                p_adjust_method=p_adjust_method,
                                preencoded_outcome=encoded_binary,
                                encoding_note=encoding_note
                            )

                            if reg_df.empty:
                                self._outcome_stats_log(f'Skipped {category}: no valid models after filtering.')
                                continue

                            summary_df = pd.DataFrame([
                                ('category', category),
                                ('analysis_type', 'continuous_outcome_ols_fallback'),
                                ('method', 'Linear Model (OLS)'),
                                ('requested_numeric_mode', numeric_mode_raw),
                                ('n_samples', len(valid_samples)),
                                ('n_unknown_dropped', dropped_unknown),
                                ('n_metabolites_tested', len(reg_df)),
                                ('n_significant', int((reg_df['adj_p'] < alpha).sum())),
                                ('covariates', ', '.join([c for c in covariate_cols if c != category]) if covariate_cols else 'None'),
                                ('strict_mode', strict_mode),
                                ('alpha', alpha),
                                ('p_adjust_method', p_adjust_method),
                                ('fallback_reason', str(be))
                            ], columns=['Parameter', 'Value'])

                            with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                                reg_df.to_excel(writer, sheet_name='Results', index=False)
                                summary_df.to_excel(writer, sheet_name='Summary', index=False)

                            exported_files.append(out_file)
                            self._outcome_stats_log(f'Exported {category} -> {os.path.basename(out_file)}')
                            continue

                    if method == 'Ordinal Regression' and (numeric_mode != 'continuous' or not is_continuous_numeric_for_ordinal):
                        if n_unique > max_ordinal_levels and numeric_mode == 'continuous' and is_numeric_outcome:
                            self._outcome_stats_log(
                                f'Notice: {category} has {n_unique} unique numeric levels. '
                                f'Using OLS fallback instead of ordinal for efficiency/stability.'
                            )
                        else:
                            self._outcome_stats_log(
                                f'Notice: Ordinal regression requires numeric continuous outcome. '
                                f'(mode={numeric_mode_raw}, outcome_type={"numeric" if is_numeric_outcome else "categorical"}). '
                                f'Using standard analysis instead.'
                            )

                    sample_cols = []
                    group_map = {}
                    for s in valid_samples:
                        val = str(cat_non_na.loc[s]).strip()
                        if val == '' or val.lower() == 'nan':
                            continue
                        sample_cols.append(s)
                        group_map[s] = val

                    unique_groups = sorted(set(group_map.values()))
                    if len(unique_groups) < 2:
                        self._outcome_stats_log(f'Skipped {category}: fewer than 2 groups after cleaning.')
                        continue

                    if strict_mode and any(str(g).strip().casefold() in unknown_tokens for g in unique_groups):
                        self._outcome_stats_log(
                            f'Skipped {category}: unknown labels remain after strict filtering.'
                        )
                        continue

                    met_df = working.loc[sample_cols, metabolite_cols].apply(pd.to_numeric, errors='coerce')
                    df_intensities = met_df.T.reset_index().rename(columns={'index': 'Name'})
                    ordered_cols = ['Name'] + sample_cols
                    df_intensities = df_intensities[ordered_cols]

                    cov_df = pd.DataFrame(index=sample_cols)
                    use_covariates = [c for c in covariate_cols if c != category and c in working.columns]
                    if use_covariates:
                        cov_df = working.loc[sample_cols, use_covariates].copy()

                    analysis_func = run_limma_covariate_analysis if method == 'Limma' else run_covariate_adjusted_analysis
                    result = analysis_func(
                        df_intensities=df_intensities,
                        sample_cols=sample_cols,
                        group_map=group_map,
                        covariate_data=cov_df if use_covariates else None,
                        covariate_cols=use_covariates if use_covariates else None,
                        reference_group=None,
                        apply_fdr=(p_adjust_method != 'None'),
                        fdr_method=p_adjust_method,
                        alpha=alpha,
                        metabolite_id_col='Name',
                        return_adjusted_intensities=False,
                        group_order=unique_groups,
                        min_samples_per_group=2,
                        min_samples_type='absolute'
                    )

                    export_covariate_results(result, out_file, include_diagnostics=True, include_coefficients=True, class_results=None)

                    if cluster_assignment_df is not None and not cluster_assignment_df.empty:
                        with pd.ExcelWriter(out_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                            cluster_assignment_df.to_excel(writer, sheet_name='Outcome_Cluster_Assignments', index=False)

                    exported_files.append(out_file)
                    self._outcome_stats_log(
                        f'Exported {category} ({len(unique_groups)} groups: {", ".join(unique_groups)}) -> {os.path.basename(out_file)}'
                    )

                except Exception as e:
                    self._outcome_stats_log(f'Error in category {category}: {e}')
                    logger.error(f'Outcome stats category failed ({category}): {e}', exc_info=True)

            def _finish_ui():
                try:
                    if hasattr(self, 'outcome_stats_run_btn'):
                        self.outcome_stats_run_btn.config(state='normal', text='Run Per-Category Analysis and Export')
                except Exception:
                    pass

                if exported_files:
                    self._outcome_stats_log(f'Completed. Exported {len(exported_files)} file(s).')
                    messagebox.showinfo(
                        'Analysis Complete',
                        f'Exported {len(exported_files)} category-specific file(s) to:\n{output_dir}'
                    )
                else:
                    self._outcome_stats_log('No outputs generated (check mappings/data validity).')
                    messagebox.showwarning('No Output', 'No output files were generated. Check log for details.')

            try:
                self.root.after(0, _finish_ui)
            except Exception:
                pass
        except Exception as e:
            logger.error(f'Outcome stats worker crashed: {e}', exc_info=True)

            def _fatal_ui():
                try:
                    if hasattr(self, 'outcome_stats_run_btn'):
                        self.outcome_stats_run_btn.config(state='normal', text='Run Per-Category Analysis and Export')
                except Exception:
                    pass
                self._outcome_stats_log(f'Fatal error: {e}')
                messagebox.showerror('Analysis Error', f'Analysis failed:\n{e}')

            try:
                self.root.after(0, _fatal_ui)
            except Exception:
                pass
    
    # ==================== CHORD DIAGRAM TAB ====================
    def setup_chord_diagram_tab(self):
        """Setup Chord Diagram sub-tab"""
        chord_frame = ttk.Frame(self.notebook)
        self.notebook.add(chord_frame, text='🔄 Chord Diagram')

        # Scrollable container so controls remain reachable on smaller screens.
        scroll_host = tk.Frame(chord_frame, bg='#f0f0f0')
        scroll_host.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(scroll_host, bg='#f0f0f0', highlightthickness=0)
        v_scroll = ttk.Scrollbar(scroll_host, orient='vertical', command=canvas.yview)
        container = tk.Frame(canvas, bg='#f0f0f0')
        container.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))

        canvas_window = canvas.create_window((0, 0), window=container, anchor='nw')
        canvas.configure(yscrollcommand=v_scroll.set)

        def _fit_chord_container_width(event):
            canvas.itemconfigure(canvas_window, width=event.width)

        canvas.bind('<Configure>', _fit_chord_container_width)
        canvas.pack(side='left', fill='both', expand=True)
        v_scroll.pack(side='right', fill='y')
        
        # Title
        title_label = tk.Label(container, text='Chord Diagram Generator', 
                               font=('Arial', 14, 'bold'), bg='#f0f0f0')
        title_label.pack(pady=10)
        
        # Upload section
        upload_frame = tk.LabelFrame(container, text='Step 1: Upload Data', 
                                     bg='#f0f0f0', font=('Arial', 10, 'bold'))
        upload_frame.pack(fill='x', pady=5)
        
        btn_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Button(btn_frame, text='📂 Upload Excel File', 
                 command=self.upload_chord_data, bg='#4CAF50', fg='white',
                 font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        
        self.chord_file_label = tk.Label(btn_frame, text='No file uploaded', 
                                         bg='#f0f0f0', font=('Arial', 9))
        self.chord_file_label.pack(side='left', padx=5)
        
        # Column mapping section
        mapping_frame = tk.LabelFrame(container, text='Step 2: Map Columns', 
                                      bg='#f0f0f0', font=('Arial', 10, 'bold'))
        mapping_frame.pack(fill='x', pady=5)
        
        # Compound column
        compound_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        compound_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(compound_frame, text='Compound Column:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_compound_var = tk.StringVar()
        self.chord_compound_combo = ttk.Combobox(compound_frame, 
                                                  textvariable=self.chord_compound_var, 
                                                  state='readonly', width=30)
        self.chord_compound_combo.pack(side='left', padx=5)
        
        # Classification column
        class_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        class_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(class_frame, text='Classification Column:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_class_var = tk.StringVar()
        self.chord_class_combo = ttk.Combobox(class_frame, 
                                               textvariable=self.chord_class_var, 
                                               state='readonly', width=30)
        self.chord_class_combo.pack(side='left', padx=5)
        
        # Optional: Log2FC column for coloring
        fc_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        fc_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(fc_frame, text='Log2FC Column (optional):', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_fc_var = tk.StringVar()
        self.chord_fc_combo = ttk.Combobox(fc_frame, 
                                            textvariable=self.chord_fc_var, 
                                            state='readonly', width=30)
        self.chord_fc_combo.pack(side='left', padx=5)
        
        # Optional: Value/Weight column
        value_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        value_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(value_frame, text='Value/Weight (optional):', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_value_var = tk.StringVar()
        self.chord_value_combo = ttk.Combobox(value_frame, 
                                               textvariable=self.chord_value_var, 
                                               state='readonly', width=30)
        self.chord_value_combo.pack(side='left', padx=5)
        
        # Settings section
        settings_frame = tk.LabelFrame(container, text='Step 3: Configure Settings', 
                                       bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=5)
        
        # Upper case option
        upper_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        upper_frame.pack(fill='x', padx=5, pady=5)
        self.chord_uppercase_var = tk.BooleanVar(value=True)
        tk.Checkbutton(upper_frame, text='Uppercase Labels', 
                      variable=self.chord_uppercase_var, 
                      bg='#f0f0f0').pack(side='left')
        
        # Label orientation
        orient_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        orient_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(orient_frame, text='Label Orientation:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_orientation_var = tk.StringVar(value='outward')
        orient_combo = ttk.Combobox(orient_frame, 
                                     textvariable=self.chord_orientation_var,
                                     values=['outward', 'inward', 'tangent'],
                                     state='readonly', width=15)
        orient_combo.pack(side='left', padx=5)
        
        # Font size controls
        fontsize_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        fontsize_frame.pack(fill='x', padx=5, pady=5)
        
        # Compound font size
        tk.Label(fontsize_frame, text='Compound Font Size:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_compound_fontsize_var = tk.IntVar(value=8)
        compound_fs_spin = tk.Spinbox(fontsize_frame, from_=6, to=20, 
                                       textvariable=self.chord_compound_fontsize_var,
                                       width=5)
        compound_fs_spin.pack(side='left', padx=5)
        
        # Classification font size
        tk.Label(fontsize_frame, text='Classification Font Size:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.chord_class_fontsize_var = tk.IntVar(value=9)
        class_fs_spin = tk.Spinbox(fontsize_frame, from_=6, to=20,
                                    textvariable=self.chord_class_fontsize_var,
                                    width=5)
        class_fs_spin.pack(side='left', padx=5)
        
        # Max character control
        maxchar_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        maxchar_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(maxchar_frame, text='Max Characters Before Line Break:', bg='#f0f0f0', 
                width=30, anchor='w').pack(side='left')
        self.chord_max_chars_var = tk.IntVar(value=16)
        maxchar_spin = tk.Spinbox(maxchar_frame, from_=10, to=50,
                                   textvariable=self.chord_max_chars_var,
                                   width=5)
        maxchar_spin.pack(side='left', padx=5)
        
        # Legend options
        legend_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        legend_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(legend_frame, text='Legends:', bg='#f0f0f0', 
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        
        self.chord_show_compound_legend_var = tk.BooleanVar(value=True)
        tk.Checkbutton(legend_frame, text='Show Compound Legend (bottom)', 
                      variable=self.chord_show_compound_legend_var, 
                      bg='#f0f0f0').pack(side='left', padx=5)
        
        self.chord_show_class_legend_var = tk.BooleanVar(value=True)
        tk.Checkbutton(legend_frame, text='Show Classification Legend (right)', 
                      variable=self.chord_show_class_legend_var, 
                      bg='#f0f0f0').pack(side='left', padx=5)
        
        # Figure size controls
        size_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        size_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(size_frame, text='Figure Size:', bg='#f0f0f0', 
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        
        tk.Label(size_frame, text='Width:', bg='#f0f0f0', 
                width=8, anchor='w').pack(side='left', padx=2)
        self.chord_fig_width_var = tk.IntVar(value=10)
        tk.Spinbox(size_frame, from_=6, to=20, textvariable=self.chord_fig_width_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(size_frame, text='Height:', bg='#f0f0f0', 
                width=8, anchor='w').pack(side='left', padx=2)
        self.chord_fig_height_var = tk.IntVar(value=10)
        tk.Spinbox(size_frame, from_=6, to=20, textvariable=self.chord_fig_height_var,
                  width=5).pack(side='left', padx=2)
        
        # Generate button
        gen_frame = tk.Frame(container, bg='#f0f0f0')
        gen_frame.pack(fill='x', pady=10)
        
        tk.Button(gen_frame, text='🎨 Generate Chord Diagram', 
                 command=self.generate_chord_diagram, 
                 bg='#2196F3', fg='white',
                 font=('Arial', 12, 'bold'), height=2).pack(pady=5)
        
        # Status label
        self.chord_status_label = tk.Label(container, text='', 
                                           bg='#f0f0f0', font=('Arial', 9), 
                                           fg='blue')
        self.chord_status_label.pack(pady=5)
        
        # Initialize data storage
        self.chord_data = None
    
    def upload_chord_data(self):
        """Upload Excel file for chord diagram"""
        file_path = filedialog.askopenfilename(
            title='Select Excel File',
            filetypes=[('Excel Files', '*.xlsx *.xls'), ('All Files', '*.*')]
        )
        
        if not file_path:
            return
        
        try:
            # Read Excel file
            df = pd.read_excel(file_path)
            self.chord_data = df
            
            # Update file label
            filename = os.path.basename(file_path)
            self.chord_file_label.config(text=f'✓ {filename} ({len(df)} rows)', fg='green')
            
            # Update column dropdowns
            columns = [''] + list(df.columns)
            self.chord_compound_combo['values'] = columns
            self.chord_class_combo['values'] = columns
            self.chord_fc_combo['values'] = columns
            self.chord_value_combo['values'] = columns
            
            # Try to auto-detect columns
            cols_lower = {c.lower(): c for c in df.columns}
            
            # Compound column
            for candidate in ['metabolite', 'compound', 'name', 'molecule', 'glycan', 'glycopeptide']:
                if candidate in cols_lower:
                    self.chord_compound_var.set(cols_lower[candidate])
                    break
            
            # Classification column
            for candidate in ['pathway', 'pathways', 'classification', 'class', 'enzyme', 'gene', 'site']:
                if candidate in cols_lower:
                    self.chord_class_var.set(cols_lower[candidate])
                    break
            
            # Log2FC column
            for candidate in ['log2fc', 'log2_fc', 'logfc', 'fc', 'foldchange']:
                if candidate in cols_lower:
                    self.chord_fc_var.set(cols_lower[candidate])
                    break
            
            # Value column
            for candidate in ['value', 'weight', 'count', 'n', 'hits']:
                if candidate in cols_lower:
                    self.chord_value_var.set(cols_lower[candidate])
                    break
            
            self.chord_status_label.config(text='✓ File uploaded successfully!', fg='green')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load file:\n{str(e)}')
            logger.error(f'Error loading chord data: {e}')
    
    def generate_chord_diagram(self):
        """Generate chord diagram from uploaded data"""
        if self.chord_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return
        
        compound_col = self.chord_compound_var.get()
        class_col = self.chord_class_var.get()
        
        if not compound_col or not class_col:
            messagebox.showwarning('Missing Columns', 
                                   'Please select both Compound and Classification columns.')
            return
        
        try:
            # Import pycirclize
            try:
                import pycirclize
            except ImportError:
                messagebox.showerror('Missing Library', 
                                     'pycirclize library is required. Install with:\npip install pycirclize')
                return
            
            df = self.chord_data.copy()
            
            # Get optional columns
            fc_col = self.chord_fc_var.get()
            value_col = self.chord_value_var.get()
            
            # Prepare data
            df_clean = pd.DataFrame({
                'Compound': df[compound_col].astype(str).str.strip(),
                'Classification': df[class_col].astype(str).str.strip()
            })
            
            if fc_col:
                df_clean['log2FC'] = pd.to_numeric(df[fc_col], errors='coerce').fillna(0)
            else:
                df_clean['log2FC'] = 0.0
            
            if value_col:
                df_clean['value'] = pd.to_numeric(df[value_col], errors='coerce').fillna(1)
            else:
                df_clean['value'] = 1.0
            
            df_clean = df_clean.dropna(subset=['Compound', 'Classification'])
            
            if len(df_clean) == 0:
                messagebox.showwarning('No Data', 'No valid data after cleaning.')
                return
            
            # Generate chord diagram using the existing logic
            self.chord_status_label.config(text='⏳ Generating chord diagram...', fg='blue')
            self.frame.update()
            
            fig = self._create_chord_diagram(
                df_clean,
                uppercase=self.chord_uppercase_var.get(),
                orientation=self.chord_orientation_var.get(),
                compound_fontsize=self.chord_compound_fontsize_var.get(),
                class_fontsize=self.chord_class_fontsize_var.get(),
                max_chars=self.chord_max_chars_var.get(),
                show_compound_legend=self.chord_show_compound_legend_var.get(),
                show_class_legend=self.chord_show_class_legend_var.get(),
                fig_width=self.chord_fig_width_var.get(),
                fig_height=self.chord_fig_height_var.get()
            )
            
            # Save option
            save_path = filedialog.asksaveasfilename(
                defaultextension='.png',
                filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All Files', '*.*')]
            )
            
            if save_path:
                fig.savefig(save_path, dpi=300, bbox_inches='tight')
                self.chord_status_label.config(text=f'✓ Saved to {os.path.basename(save_path)}', fg='green')
            else:
                self.chord_status_label.config(text='✓ Diagram generated successfully!', fg='green')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to generate chord diagram:\n{str(e)}')
            logger.error(f'Error generating chord diagram: {e}', exc_info=True)
            self.chord_status_label.config(text='✗ Generation failed', fg='red')
    
    def _create_chord_diagram(self, df, uppercase=True, orientation='outward',
                             compound_fontsize=9, class_fontsize=9, max_chars=25,
                             show_compound_legend=True, show_class_legend=True,
                             fig_width=10, fig_height=10):
        """Create chord diagram using pycirclize with optional legends"""
        import random
        import pycirclize
        from matplotlib.gridspec import GridSpec
        from matplotlib.patches import Rectangle
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        import matplotlib.cm as cm
        
        # Helper functions from chord_diagram.py
        def hex_to_rgb01(h):
            h = h.lstrip("#")
            return (int(h[0:2], 16)/255, int(h[2:4], 16)/255, int(h[4:6], 16)/255)
        
        def is_too_red_or_green(hex_color):
            r, g, b = hex_to_rgb01(hex_color)
            if r > 0.70 and g < 0.45:
                return True
            if g > 0.70 and r < 0.45:
                return True
            if (r - g) > 0.35 and r > 0.6:
                return True
            if (g - r) > 0.35 and g > 0.6:
                return True
            return False
        
        def random_distinct_colors(items, seed=1):
            random.seed(seed)
            colors = {}
            while len(colors) < len(items):
                item = items[len(colors)]
                r = random.randint(40, 235)
                g = random.randint(40, 235)
                b = random.randint(40, 235)
                hexc = f"#{r:02x}{g:02x}{b:02x}"
                if is_too_red_or_green(hexc):
                    continue
                colors[item] = hexc
            return colors
        
        def log2fc_to_hex(log2fc, vmin, vmax):
            if vmax == vmin:
                t = 0.5
            else:
                t = (log2fc - vmin) / (vmax - vmin)
                t = max(0.0, min(1.0, t))
            
            green = np.array([0x18, 0xA5, 0x58]) / 255.0
            gray = np.array([0xD9, 0xD9, 0xD9]) / 255.0
            red = np.array([0xE5, 0x39, 0x35]) / 255.0
            
            if t < 0.5:
                a = t * 2.0
                rgb = (1 - a) * green + a * gray
            else:
                a = (t - 0.5) * 2.0
                rgb = (1 - a) * gray + a * red
            
            return "#{:02x}{:02x}{:02x}".format(
                int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
            )
        
        def wrap_label(label, max_chars):
            """Wrap label text at max_chars boundary"""
            if len(label) <= max_chars:
                return label
            words = label.split()
            lines = []
            current_line = []
            for word in words:
                if len(" ".join(current_line + [word])) <= max_chars:
                    current_line.append(word)
                else:
                    if current_line:
                        lines.append(" ".join(current_line))
                    current_line = [word]
            if current_line:
                lines.append(" ".join(current_line))
            return "\n".join(lines) if lines else label
        
        # Calculate mean log2FC for each compound
        compound_fc = df.groupby('Compound', as_index=False)['log2FC'].mean().fillna(0.0)
        compound_fc_map = dict(zip(compound_fc['Compound'], compound_fc['log2FC']))
        
        # Get log2FC range
        vals = compound_fc['log2FC'].values
        max_abs = np.nanmax(np.abs(vals)) if len(vals) else 1.0
        max_abs = max(0.5, float(max_abs))
        vmin, vmax = -max_abs, max_abs
        
        # Get unique classifications and compounds
        classifications = sorted(df['Classification'].unique().tolist())
        compounds = sorted(df['Compound'].unique().tolist())
        
        # Assign colors
        classification_colors = random_distinct_colors(classifications, seed=1)
        compound_colors = {c: log2fc_to_hex(compound_fc_map.get(c, 0.0), vmin, vmax) 
                          for c in compounds}
        
        # Build edges
        edges = []
        for _, r in df.iterrows():
            edges.append({
                'classification': r['Classification'],
                'compound': r['Compound'],
                'log2fc': float(r['log2FC']) if np.isfinite(r['log2FC']) else 0.0,
                'value': float(r['value']) if np.isfinite(r['value']) else 1.0
            })
        
        # Calculate sector sizes
        sectors = {}
        for c in classifications:
            sectors[c] = int(sum(1 for e in edges if e['classification'] == c)) or 1
        for c in compounds:
            sectors[c] = int(sum(1 for e in edges if e['compound'] == c)) or 1
        
        # Create circos plot
        circos = pycirclize.Circos(sectors, space=5)
        
        # Add tracks and labels
        for sector in circos.sectors:
            name = sector.name
            
            # Sector color
            if name in compound_colors:
                sec_color = compound_colors[name]
            elif name in classification_colors:
                sec_color = classification_colors[name]
            else:
                sec_color = "#CCCCCC"
            
            # Outer ring
            tr = sector.add_track((95, 100))
            tr.rect(sector.start, sector.end, fc=sec_color, ec="black", lw=0.6)
            
            # Label
            label = name.upper() if uppercase else name
            label = wrap_label(label, max_chars)
            
            # Determine font size based on label type
            is_compound = name in compound_colors
            font_size = compound_fontsize if is_compound else class_fontsize
            
            # Compute mid angle
            deg_start, deg_end = sector.deg_lim
            mid_deg = (deg_start + deg_end) / 2.0
            
            # Rotation and alignment
            if 90 < mid_deg < 270:
                rotation = 90 - mid_deg + 180
                ha = 'right'
            else:
                rotation = 90 - mid_deg
                ha = 'left'
            
            # Normalize rotation
            rotation = ((rotation + 180) % 360) - 180
            if rotation < -90:
                rotation += 180
                ha = 'right' if ha == 'left' else 'left'
            elif rotation > 90:
                rotation -= 180
                ha = 'right' if ha == 'left' else 'left'
            
            sector.text(
                label,
                r=112,
                size=font_size,
                weight="bold",
                rotation=rotation,
                ha=ha,
                va="center",
                rotation_mode="anchor",
                adjust_rotation=False
            )
        
        # Add links
        for e in edges:
            c = circos.get_sector(e['classification'])
            m = circos.get_sector(e['compound'])
            
            lw = max(0.6, min(3.0, 0.6 + 0.35 * abs(e['value'])))
            circos.link(
                (c.name, c.start, c.end),
                (m.name, m.start, m.end),
                color=classification_colors.get(e['classification'], "#999999"),
                alpha=0.55,
                linewidth=lw
            )
        
        # Create figure with circos plot and legends using simple approach
        num_sectors = len(sectors)
        
        # Let pycirclize create its own figure
        fig = circos.plotfig()
        if fig is None:
            fig = plt.gcf()
        
        # Set size and colors
        adjusted_width = fig_width + (2 if show_class_legend else 0)
        adjusted_height = fig_height + (1.5 if show_compound_legend else 0)
        fig.set_size_inches(adjusted_width, adjusted_height)
        fig.patch.set_facecolor("white")
        
        # Get the existing axes from circos plot
        ax_circos = circos.ax
        
        # Create compound legend (bottom)
        if show_compound_legend and len(compounds) > 0:
            # Use text annotations instead of separate axes
            sorted_compounds = sorted(compounds)
            legend_y_start = -0.15
            legend_x_start = 0.05
            items_per_row = 4
            
            for idx, comp in enumerate(sorted_compounds):
                row = idx // items_per_row
                col = idx % items_per_row
                x = legend_x_start + (col * 0.22)
                y = legend_y_start - (row * 0.05)
                
                color = compound_colors.get(comp, "#CCCCCC")
                display_label = comp.upper() if uppercase else comp
                
                # Add colored square marker
                fig.text(x, y, '■', fontsize=12, color=color, transform=fig.transFigure)
                # Add label
                fig.text(x + 0.015, y, display_label, fontsize=compound_fontsize, 
                        va='center', transform=fig.transFigure)
        
        # Create classification legend (right)
        if show_class_legend and len(classifications) > 0:
            sorted_classifications = sorted(classifications)
            legend_x = 0.88
            legend_y_start = 0.95
            
            for idx, clf in enumerate(sorted_classifications):
                y = legend_y_start - (idx * 0.05)
                
                color = classification_colors.get(clf, "#CCCCCC")
                display_label = clf.upper() if uppercase else clf
                
                # Add colored square marker
                fig.text(legend_x, y, '■', fontsize=10, color=color, transform=fig.transFigure)
                # Add label
                fig.text(legend_x + 0.015, y, display_label, fontsize=class_fontsize, 
                        va='center', transform=fig.transFigure)
        
        plt.tight_layout()
        
        fig.patch.set_facecolor("white")
        
        return fig
    
    # ==================== VENN DIAGRAM TAB ====================
    def setup_venn_diagram_tab(self):
        """Setup Venn Diagram sub-tab"""
        venn_frame = ttk.Frame(self.notebook)
        self.notebook.add(venn_frame, text='⭕ Venn Diagram')
        
        # Main container
        container = tk.Frame(venn_frame, bg='#f0f0f0')
        container.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Title
        title_label = tk.Label(container, text='Venn Diagram Generator', 
                               font=('Arial', 14, 'bold'), bg='#f0f0f0')
        title_label.pack(pady=10)
        
        # Number of sets
        num_frame = tk.Frame(container, bg='#f0f0f0')
        num_frame.pack(fill='x', pady=5)
        tk.Label(num_frame, text='Number of Sets:', bg='#f0f0f0', 
                font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        self.venn_num_sets_var = tk.IntVar(value=2)
        for i in range(2, 7):  # Support up to 6 sets
            tk.Radiobutton(num_frame, text=str(i), variable=self.venn_num_sets_var, 
                          value=i, bg='#f0f0f0', 
                          command=self.update_venn_sets).pack(side='left', padx=5)
        
        # Scrollable frame for sets
        canvas = tk.Canvas(container, bg='#f0f0f0', height=400)
        scrollbar = ttk.Scrollbar(container, orient='vertical', command=canvas.yview)
        self.venn_sets_frame = tk.Frame(canvas, bg='#f0f0f0')
        
        self.venn_sets_frame.bind(
            '<Configure>',
            lambda e: canvas.configure(scrollregion=canvas.bbox('all'))
        )
        
        canvas.create_window((0, 0), window=self.venn_sets_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side='left', fill='both', expand=True, pady=10)
        scrollbar.pack(side='right', fill='y')
        
        # Settings frame
        settings_frame = tk.LabelFrame(container, text='⚙️ Settings', 
                                       bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=5)
        
        # Figure size controls
        size_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        size_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(size_frame, text='Width:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.venn_fig_width_var = tk.IntVar(value=10)
        tk.Spinbox(size_frame, from_=6, to=18, textvariable=self.venn_fig_width_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(size_frame, text='Height:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.venn_fig_height_var = tk.IntVar(value=8)
        tk.Spinbox(size_frame, from_=6, to=16, textvariable=self.venn_fig_height_var,
                  width=5).pack(side='left', padx=2)
        
        # Generate button
        tk.Button(container, text='🎨 Generate Venn Diagram', 
                 command=self.generate_venn_diagram, 
                 bg='#2196F3', fg='white',
                 font=('Arial', 12, 'bold'), height=2).pack(pady=10)
        
        # Status label
        self.venn_status_label = tk.Label(container, text='', 
                                          bg='#f0f0f0', font=('Arial', 9), 
                                          fg='blue')
        self.venn_status_label.pack(pady=5)
        
        # Initialize set UI
        self.venn_set_widgets = []
        self.update_venn_sets()
    
    def update_venn_sets(self):
        """Update the Venn sets UI based on number of sets"""
        # Clear existing widgets
        for widget in self.venn_sets_frame.winfo_children():
            widget.destroy()
        self.venn_set_widgets = []
        
        num_sets = self.venn_num_sets_var.get()
        
        for i in range(num_sets):
            set_frame = tk.LabelFrame(self.venn_sets_frame, 
                                      text=f'Set {i+1}', 
                                      bg='#f0f0f0', 
                                      font=('Arial', 10, 'bold'))
            set_frame.pack(fill='x', padx=5, pady=5)
            
            # Name
            name_frame = tk.Frame(set_frame, bg='#f0f0f0')
            name_frame.pack(fill='x', padx=5, pady=5)
            tk.Label(name_frame, text='Name:', bg='#f0f0f0', width=10, 
                    anchor='w').pack(side='left')
            name_var = tk.StringVar(value=f'Set {i+1}')
            name_entry = tk.Entry(name_frame, textvariable=name_var, width=30)
            name_entry.pack(side='left', padx=5)
            
            # Color
            color_frame = tk.Frame(set_frame, bg='#f0f0f0')
            color_frame.pack(fill='x', padx=5, pady=5)
            tk.Label(color_frame, text='Color:', bg='#f0f0f0', width=10, 
                    anchor='w').pack(side='left')
            
            default_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#95E1D3', '#F38181']
            color_var = tk.StringVar(value=default_colors[i % len(default_colors)])
            color_entry = tk.Entry(color_frame, textvariable=color_var, width=15)
            color_entry.pack(side='left', padx=5)
            
            def choose_color(cv=color_var):
                from tkinter import colorchooser
                color = colorchooser.askcolor(initialcolor=cv.get())
                if color[1]:
                    cv.set(color[1])
            
            tk.Button(color_frame, text='🎨', command=choose_color).pack(side='left', padx=2)
            
            # Data input method
            method_frame = tk.Frame(set_frame, bg='#f0f0f0')
            method_frame.pack(fill='x', padx=5, pady=5)
            method_var = tk.StringVar(value='paste')
            tk.Radiobutton(method_frame, text='Paste Data', variable=method_var, 
                          value='paste', bg='#f0f0f0').pack(side='left', padx=5)
            tk.Radiobutton(method_frame, text='Upload Excel', variable=method_var, 
                          value='upload', bg='#f0f0f0').pack(side='left', padx=5)
            
            # Text area for pasting
            text_frame = tk.Frame(set_frame, bg='#f0f0f0')
            text_frame.pack(fill='both', expand=True, padx=5, pady=5)
            tk.Label(text_frame, text='Data (one item per line):', 
                    bg='#f0f0f0', anchor='w').pack(anchor='w')
            text_widget = tk.Text(text_frame, height=5, width=40)
            text_widget.pack(fill='both', expand=True)
            
            # Upload button
            upload_btn_frame = tk.Frame(set_frame, bg='#f0f0f0')
            upload_btn_frame.pack(fill='x', padx=5, pady=5)
            
            file_label_var = tk.StringVar(value='No file')
            file_label = tk.Label(upload_btn_frame, textvariable=file_label_var, 
                                  bg='#f0f0f0', fg='gray')
            file_label.pack(side='right', padx=5)
            
            def upload_set_data(tw=text_widget, flv=file_label_var):
                file_path = filedialog.askopenfilename(
                    title='Select Excel File',
                    filetypes=[('Excel Files', '*.xlsx *.xls'), ('All Files', '*.*')]
                )
                if file_path:
                    try:
                        df = pd.read_excel(file_path)
                        if len(df.columns) > 0:
                            # Use first column
                            data = df.iloc[:, 0].dropna().astype(str).tolist()
                            tw.delete('1.0', tk.END)
                            tw.insert('1.0', '\n'.join(data))
                            flv.set(f'✓ {os.path.basename(file_path)}')
                    except Exception as e:
                        messagebox.showerror('Error', f'Failed to load file:\n{str(e)}')
            
            tk.Button(upload_btn_frame, text='📂 Upload Excel', 
                     command=upload_set_data).pack(side='left', padx=5)
            
            self.venn_set_widgets.append({
                'name': name_var,
                'color': color_var,
                'method': method_var,
                'text': text_widget,
                'file_label': file_label_var
            })
    
    def generate_venn_diagram(self):
        """Generate Venn diagram using venn package (supports 2-6 sets)"""
        try:
            # Import venn package
            try:
                import venn
            except ImportError:
                messagebox.showerror('Missing Library', 
                                     'venn library is required. Install with:\npip install venn')
                return
            
            num_sets = self.venn_num_sets_var.get()
            
            # Collect data from each set
            sets_data = {}  # venn package uses dict
            colors_dict = {}
            
            for widget_dict in self.venn_set_widgets[:num_sets]:
                # Get data
                text_data = widget_dict['text'].get('1.0', tk.END).strip()
                if not text_data:
                    messagebox.showwarning('Empty Set', 
                                           f"Set '{widget_dict['name'].get()}' is empty.")
                    return
                
                # Split by newline and clean
                items = set([item.strip() for item in text_data.split('\n') 
                            if item.strip()])
                
                set_name = widget_dict['name'].get()
                sets_data[set_name] = items
                colors_dict[set_name] = widget_dict['color'].get()
            
            # Generate Venn diagram
            self.venn_status_label.config(text='⏳ Generating Venn diagram...', fg='blue')
            self.frame.update()
            
            # Create figure with user-specified size
            fig_width = self.venn_fig_width_var.get()
            fig_height = self.venn_fig_height_var.get()
            
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            
            # Convert colors dict to list in same order as sets_data
            set_names = list(sets_data.keys())
            colors_list = [colors_dict.get(name, '#CCCCCC') for name in set_names]
            
            # Some venn package builds do not expose a typed `colors` parameter.
            # Use a broadly compatible call signature.
            venn.venn(sets_data, ax=ax, fontsize=10)
            
            plt.title('Venn Diagram', fontsize=14, fontweight='bold', pad=20)
            plt.tight_layout()
            
            # Save option
            save_path = filedialog.asksaveasfilename(
                defaultextension='.png',
                filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All Files', '*.*')]
            )
            
            if save_path:
                fig.savefig(save_path, dpi=300, bbox_inches='tight')
                self.venn_status_label.config(text=f'✓ Saved to {os.path.basename(save_path)}', fg='green')
            else:
                self.venn_status_label.config(text='✓ Diagram generated successfully!', fg='green')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to generate Venn diagram:\n{str(e)}')
            logger.error(f'Error generating Venn diagram: {e}', exc_info=True)
            self.venn_status_label.config(text='✗ Generation failed', fg='red')
    
    # ==================== PIE CHART TAB ====================
    def setup_pie_chart_tab(self):
        """Setup Pie Chart sub-tab"""
        pie_frame = ttk.Frame(self.notebook)
        self.notebook.add(pie_frame, text='🥧 Pie Chart')

        # Two-column layout: controls on left, live log on right.
        body = tk.Frame(pie_frame, bg='#f0f0f0')
        body.pack(fill='both', expand=True, padx=10, pady=10)

        # Scrollable container so all controls remain reachable on smaller windows.
        scroll_host = tk.Frame(body, bg='#f0f0f0')
        scroll_host.pack(side='left', fill='both', expand=True)

        canvas = tk.Canvas(scroll_host, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(scroll_host, orient='vertical', command=canvas.yview)
        container = tk.Frame(canvas, bg='#f0f0f0')
        container.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))

        canvas_window = canvas.create_window((0, 0), window=container, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        def _fit_pie_container_width(event):
            canvas.itemconfigure(canvas_window, width=event.width)

        canvas.bind('<Configure>', _fit_pie_container_width)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        pie_log_frame = tk.LabelFrame(body, text='Pie Chart Log', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        pie_log_frame.pack(side='right', fill='y', padx=(10, 0))
        self.pie_log_text = tk.Text(pie_log_frame, width=45, wrap='word', font=('Consolas', 9), bg='white')
        self.pie_log_text.pack(side='left', fill='both', expand=True, padx=(6, 0), pady=6)
        pie_log_scroll = ttk.Scrollbar(pie_log_frame, orient='vertical', command=self.pie_log_text.yview)
        pie_log_scroll.pack(side='right', fill='y', padx=(0, 6), pady=6)
        self.pie_log_text.configure(yscrollcommand=pie_log_scroll.set)
        
        # Title
        title_label = tk.Label(container, text='Pie Chart Generator', 
                               font=('Arial', 14, 'bold'), bg='#f0f0f0')
        title_label.pack(pady=10)
        
        # Upload section
        upload_frame = tk.LabelFrame(container, text='Step 1: Upload Data', 
                                     bg='#f0f0f0', font=('Arial', 10, 'bold'))
        upload_frame.pack(fill='x', pady=5)
        
        btn_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Button(btn_frame, text='📂 Upload Excel File', 
                 command=self.upload_pie_data, bg='#4CAF50', fg='white',
                 font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        
        self.pie_file_label = tk.Label(btn_frame, text='No file uploaded', 
                                       bg='#f0f0f0', font=('Arial', 9))
        self.pie_file_label.pack(side='left', padx=5)
        
        # Column mapping
        mapping_frame = tk.LabelFrame(container, text='Step 2: Map Columns', 
                                      bg='#f0f0f0', font=('Arial', 10, 'bold'))
        mapping_frame.pack(fill='x', pady=5)
        
        # Class column
        class_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        class_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(class_frame, text='Class Column:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.pie_class_var = tk.StringVar()
        self.pie_class_combo = ttk.Combobox(class_frame, 
                                             textvariable=self.pie_class_var, 
                                             state='readonly', width=30)
        self.pie_class_combo.pack(side='left', padx=5)
        
        # Sample columns selection
        sample_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        sample_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(sample_frame, text='Sample Columns:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        tk.Button(sample_frame, text='🔍 Verify Columns', 
                 command=self.select_pie_sample_columns, bg='#4CAF50', fg='white').pack(side='left', padx=5)
        self.pie_sample_label = tk.Label(sample_frame, text='No columns selected', 
                                         bg='#f0f0f0', font=('Arial', 9))
        self.pie_sample_label.pack(side='left', padx=5)
        
        # Group configuration
        group_frame = tk.LabelFrame(container, text='Step 3: Configure Groups', 
                                    bg='#f0f0f0', font=('Arial', 10, 'bold'))
        group_frame.pack(fill='x', pady=5)
        
        tk.Button(group_frame, text='⚙️ Configure Groups & Patterns', 
                 command=self.configure_pie_groups,
                 bg='#FF9800', fg='white',
                 font=('Arial', 10, 'bold')).pack(padx=5, pady=5)
        
        self.pie_group_label = tk.Label(group_frame, text='No groups configured', 
                                        bg='#f0f0f0', font=('Arial', 9))
        self.pie_group_label.pack(padx=5, pady=5)
        
        # Settings frame (collapsible)
        settings_frame = tk.LabelFrame(container, text='⚙️ Settings', 
                                       bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=5)
        
        # Label display options
        label_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        label_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(label_frame, text='Pie Labels:', bg='#f0f0f0', 
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.pie_label_style_var = tk.StringVar(value='inside')
        for style in ['inside', 'legend', 'remove']:
            tk.Radiobutton(label_frame, text=style.capitalize(), 
                          variable=self.pie_label_style_var, value=style, 
                          bg='#f0f0f0').pack(side='left', padx=5)
        
        # Font size controls
        font_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        font_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(font_frame, text='Title Font Size:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.pie_title_fontsize_var = tk.IntVar(value=12)
        tk.Spinbox(font_frame, from_=8, to=24, textvariable=self.pie_title_fontsize_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(font_frame, text='Label Font Size:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.pie_label_fontsize_var = tk.IntVar(value=10)
        tk.Spinbox(font_frame, from_=6, to=20, textvariable=self.pie_label_fontsize_var,
                  width=5).pack(side='left', padx=2)
        
        # Image size
        size_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        size_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(size_frame, text='Width:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.pie_fig_width_var = tk.IntVar(value=6)
        tk.Spinbox(size_frame, from_=4, to=16, textvariable=self.pie_fig_width_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(size_frame, text='Height:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.pie_fig_height_var = tk.IntVar(value=6)
        tk.Spinbox(size_frame, from_=4, to=16, textvariable=self.pie_fig_height_var,
                  width=5).pack(side='left', padx=2)
        
        # Save options
        save_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        save_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(save_frame, text='Save:', bg='#f0f0f0', 
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.pie_save_by_group_var = tk.BooleanVar(value=False)
        tk.Checkbutton(save_frame, text='Save Separately by Group', 
                      variable=self.pie_save_by_group_var, bg='#f0f0f0').pack(side='left', padx=5)
        
        # Generate button
        tk.Button(container, text='🎨 Generate Pie Charts', 
                 command=self.generate_pie_chart, 
                 bg='#2196F3', fg='white',
                 font=('Arial', 12, 'bold'), height=2).pack(pady=10)
        
        # Status label
        self.pie_status_label = tk.Label(container, text='', 
                                         bg='#f0f0f0', font=('Arial', 9), 
                                         fg='blue')
        self.pie_status_label.pack(pady=5)
        
        # Initialize data storage
        self.pie_data = None
        self.pie_sample_columns = []
        self.pie_group_definitions = {}
        self.pie_sample_group_mapping = {}
        self.pie_group_patterns = {}
        self._pie_log('Ready. Upload data, verify columns, configure groups, then generate charts.')

    def _pie_log(self, message: str):
        """Append message to pie chart log panel."""
        try:
            timestamp = datetime.now().strftime('%H:%M:%S')
            self.pie_log_text.insert('end', f'[{timestamp}] {message}\n')
            self.pie_log_text.see('end')
            self.pie_log_text.update_idletasks()
        except Exception:
            pass
    
    def upload_pie_data(self):
        """Upload Excel file for pie chart"""
        file_path = filedialog.askopenfilename(
            title='Select Excel File',
            filetypes=[('Excel Files', '*.xlsx *.xls'), ('All Files', '*.*')]
        )
        
        if not file_path:
            return
        
        try:
            df = pd.read_excel(file_path)
            self.pie_data = df
            
            filename = os.path.basename(file_path)
            self.pie_file_label.config(text=f'✓ {filename} ({len(df)} rows)', fg='green')
            
            # Update column dropdowns
            columns = [''] + list(df.columns)
            self.pie_class_combo['values'] = columns
            
            # Try to auto-detect class column
            cols_lower = {c.lower(): c for c in df.columns}
            for candidate in ['class', 'classification', 'category', 'type', 'pathway']:
                if candidate in cols_lower:
                    self.pie_class_var.set(cols_lower[candidate])
                    break
            
            # Reset column-specific configuration when new file is uploaded
            self.pie_sample_columns = []
            self.pie_sample_label.config(text='No columns selected', fg='red')
            self.pie_sample_group_mapping.clear()  # Clear old mappings
            self.pie_group_label.config(text='No groups configured', fg='gray')
            
            self.pie_status_label.config(text='✓ File uploaded successfully!', fg='green')
            self._pie_log(f'Loaded file: {filename} ({len(df)} rows, {len(df.columns)} columns).')
            if self.pie_class_var.get():
                self._pie_log(f'Auto-detected class column: {self.pie_class_var.get()}')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load file:\n{str(e)}')
            logger.error(f'Error loading pie data: {e}')
            self._pie_log(f'Failed to load file: {e}')
    
    def _verify_sample_columns_with_statistics_dialog(self, df: pd.DataFrame, dialog_title: str) -> Optional[List[str]]:
        """Use the same shared Verify Columns dialog used in Statistics tab."""
        result = show_column_assignment_dialog(
            parent=self.root,
            df=df,
            tab_type='statistics_metabolite',
            auto_calculate=False,
            dialog_title=dialog_title,
            detected_sample_cols=None,
            allow_skip=False,
        )

        if not result:
            return None

        sample_cols = result.get('sample_cols', [])
        if sample_cols:
            return sample_cols

        assignments = result.get('assignments', {})
        if isinstance(assignments, dict):
            return [col for col, col_type in assignments.items() if col_type == 'Sample Column']
        return []

    def _open_statistics_style_group_dialog(
        self,
        sample_columns: List[str],
        group_definitions: Dict[str, str],
        group_patterns: Dict[str, List[str]],
        sample_group_mapping: Dict[str, str],
        status_label: tk.Label,
        on_save: Optional[callable] = None,
    ):
        """Open Statistics-style 'Auto-Assign Groups by Pattern' dialog."""
        if not sample_columns:
            messagebox.showwarning('No Columns', 'Please verify columns first.')
            return

        if not group_definitions:
            group_definitions.update({'Group1': 'Group 1', 'Group2': 'Group 2'})

        pattern_window = tk.Toplevel(self.root)
        pattern_window.title('Auto-Assign Groups by Pattern')
        pattern_window.geometry('900x700')
        pattern_window.configure(bg='#f0f0f0')

        main_frame = tk.Frame(pattern_window, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=10, pady=5)

        canvas = tk.Canvas(main_frame, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(main_frame, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#f0f0f0')
        scrollable_frame.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))

        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        def configure_scroll(event):
            canvas.configure(scrollregion=canvas.bbox('all'))
            canvas.itemconfig(canvas_window, width=event.width)

        canvas.bind('<Configure>', configure_scroll)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        pattern_frame = tk.LabelFrame(scrollable_frame, text='Define Patterns', bg='#f0f0f0')
        pattern_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(pattern_frame, text='Enter keywords/patterns for each group (one per line):',
                 bg='#f0f0f0').pack(anchor='w', padx=5, pady=5)
        pattern_group_holder = tk.Frame(pattern_frame, bg='#f0f0f0')
        pattern_group_holder.pack(fill='x', padx=0, pady=0)

        label_frame = tk.LabelFrame(scrollable_frame, text='Edit Group Labels', bg='#f0f0f0')
        label_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(label_frame, text='Change the visible group names here, then refresh the assignment list.',
                 bg='#f0f0f0').pack(anchor='w', padx=5, pady=5)

        group_label_vars = {}
        label_row_holder = tk.Frame(label_frame, bg='#f0f0f0')
        label_row_holder.pack(fill='x', padx=5, pady=(0, 5))

        pattern_vars = {}

        group_manage_row = tk.Frame(label_frame, bg='#f0f0f0')
        group_manage_row.pack(fill='x', padx=5, pady=(0, 5))

        remove_group_var = tk.StringVar()

        def _next_group_id():
            n = 1
            while f'Group{n}' in group_definitions:
                n += 1
            return f'Group{n}', f'Group {n}'

        def _persist_current_group_widgets():
            for gid, txt in pattern_vars.items():
                group_patterns[gid] = [p.strip() for p in txt.get('1.0', tk.END).split('\n') if p.strip()]
            for gid, label_var in group_label_vars.items():
                group_definitions[gid] = (label_var.get().strip() or gid)

        def _refresh_remove_group_choices():
            options = [f"{gid}: {group_definitions.get(gid, gid)}" for gid in group_definitions.keys()]
            remove_group_combo['values'] = options
            if options:
                if remove_group_var.get() not in options:
                    remove_group_var.set(options[0])
            else:
                remove_group_var.set('')

        def _render_group_editors():
            for w in label_row_holder.winfo_children():
                w.destroy()
            for w in pattern_group_holder.winfo_children():
                w.destroy()
            group_label_vars.clear()
            pattern_vars.clear()

            for group_id in group_definitions.keys():
                current_label = group_definitions[group_id]
                label_row = tk.Frame(label_row_holder, bg='#f0f0f0')
                label_row.pack(fill='x', pady=2)
                tk.Label(label_row, text=f'{group_id}:', bg='#f0f0f0', width=10, anchor='w').pack(side='left')
                label_var = tk.StringVar(value=current_label)
                group_label_vars[group_id] = label_var
                tk.Entry(label_row, textvariable=label_var, width=32).pack(side='left', padx=5)

                group_frame = tk.LabelFrame(pattern_group_holder, text=f'{group_id}: {current_label}', bg='#f0f0f0')
                group_frame.pack(fill='x', padx=5, pady=5)
                pattern_text = tk.Text(group_frame, height=3, font=('Arial', 9))
                pattern_text.pack(fill='x', padx=5, pady=5)
                saved_patterns = group_patterns.get(group_id, [])
                if isinstance(saved_patterns, list) and saved_patterns:
                    pattern_text.insert('1.0', '\n'.join(saved_patterns))
                pattern_vars[group_id] = pattern_text

            _refresh_remove_group_choices()

        def _add_group():
            _persist_current_group_widgets()
            gid, label = _next_group_id()
            group_definitions[gid] = label
            group_patterns.setdefault(gid, [])
            _render_group_editors()

        def _remove_group():
            if len(group_definitions) <= 1:
                messagebox.showwarning('Cannot Remove', 'At least one group is required.')
                return

            _persist_current_group_widgets()
            raw = remove_group_var.get().strip()
            if not raw or ':' not in raw:
                messagebox.showwarning('No Group', 'Please choose a group to remove.')
                return

            gid = raw.split(':', 1)[0].strip()
            if gid not in group_definitions:
                return

            removed_label = group_definitions.get(gid, gid)
            group_definitions.pop(gid, None)
            group_patterns.pop(gid, None)
            for sample_col, group_var in sample_group_vars.items():
                if group_var.get().strip() == removed_label:
                    group_var.set('')

            _render_group_editors()
            refresh_group_labels()

        tk.Button(
            group_manage_row,
            text='Add Group',
            command=_add_group,
            bg='#27ae60',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=(0, 5))
        remove_group_combo = ttk.Combobox(group_manage_row, textvariable=remove_group_var, state='readonly', width=30)
        remove_group_combo.pack(side='left', padx=5)
        tk.Button(
            group_manage_row,
            text='Remove Group',
            command=_remove_group,
            bg='#c0392b',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=5)

        assignment_frame = tk.LabelFrame(scrollable_frame, text='Current Group Assignments', bg='#f0f0f0')
        assignment_frame.pack(fill='both', expand=True, padx=5, pady=5)

        assignment_header = tk.Frame(assignment_frame, bg='#dfeef5')
        assignment_header.pack(fill='x', padx=5, pady=(5, 2))
        tk.Label(assignment_header, text='Sample Column', bg='#dfeef5', font=('Arial', 9, 'bold'), width=30, anchor='w').pack(side='left', padx=5)
        tk.Label(assignment_header, text='Assigned Group', bg='#dfeef5', font=('Arial', 9, 'bold'), width=22, anchor='w').pack(side='left', padx=5)

        assignment_canvas = tk.Canvas(assignment_frame, bg='white', height=250)
        assignment_scrollbar = ttk.Scrollbar(assignment_frame, orient='vertical', command=assignment_canvas.yview)
        assignment_scrollable = tk.Frame(assignment_canvas, bg='white')
        assignment_scrollable.bind('<Configure>', lambda e: assignment_canvas.configure(scrollregion=assignment_canvas.bbox('all')))
        assignment_canvas.create_window((0, 0), window=assignment_scrollable, anchor='nw')
        assignment_canvas.configure(yscrollcommand=assignment_scrollbar.set)
        assignment_canvas.pack(side='left', fill='both', expand=True, padx=(5, 0), pady=(0, 5))
        assignment_scrollbar.pack(side='right', fill='y', padx=(0, 5), pady=(0, 5))

        group_labels = [group_definitions[group_id] for group_id in group_definitions.keys()]
        label_to_gid = {label: gid for gid, label in group_definitions.items()}
        group_combo_widgets = []

        sample_group_vars = {}
        for sample_col in sorted(sample_columns):
            row = tk.Frame(assignment_scrollable, bg='white')
            row.pack(fill='x', pady=1)
            tk.Label(row, text=sample_col, bg='white', font=('Arial', 8), width=30, anchor='w').pack(side='left', padx=5)

            existing_gid = sample_group_mapping.get(sample_col, '')
            existing_label = group_definitions.get(existing_gid, '')
            group_var = tk.StringVar(value=existing_label)
            sample_group_vars[sample_col] = group_var

            ttk.Combobox(
                row,
                values=[''] + group_labels,
                textvariable=group_var,
                state='readonly',
                width=22,
                font=('Arial', 8),
            ).pack(side='left', padx=5)

            combo_widget = row.winfo_children()[-1]
            group_combo_widgets.append(combo_widget)

        _render_group_editors()

        def refresh_group_labels():
            old_label_to_gid = dict(label_to_gid)
            for gid, label_var in group_label_vars.items():
                new_label = label_var.get().strip() or gid
                group_definitions[gid] = new_label

            new_group_labels = [group_definitions[group_id] for group_id in group_definitions.keys()]
            label_to_gid.clear()
            label_to_gid.update({label: gid for gid, label in group_definitions.items()})

            # Refresh the pattern section titles and assignment combo values.
            for child in pattern_frame.winfo_children():
                if isinstance(child, tk.LabelFrame) and ': ' in child.cget('text'):
                    group_id = child.cget('text').split(':', 1)[0].strip()
                    if group_id in group_definitions:
                        child.config(text=f'{group_id}: {group_definitions[group_id]}')

            for sample_col, group_var in sample_group_vars.items():
                selected_label = group_var.get().strip()
                selected_gid = old_label_to_gid.get(selected_label, '')
                new_value = group_definitions.get(selected_gid, '') if selected_gid else selected_label
                group_var.set(new_value if new_value in new_group_labels else '')

            for combo_widget in group_combo_widgets:
                combo_widget['values'] = [''] + new_group_labels

            _refresh_remove_group_choices()

            status_label.config(text=f'✓ {len(group_definitions)} groups updated', fg='green')

        bottom_btn_frame = tk.Frame(pattern_window, bg='#f0f0f0')
        bottom_btn_frame.pack(fill='x', padx=10, pady=(0, 10))

        def apply_patterns():
            for gid, pattern_widget in pattern_vars.items():
                patterns = [
                    p.strip() for p in pattern_widget.get('1.0', tk.END).strip().split('\n') if p.strip()
                ]
                target_label = group_definitions.get(gid, gid)
                for pattern in patterns:
                    for sample_col, group_var in sample_group_vars.items():
                        if pattern.lower() in sample_col.lower():
                            group_var.set(target_label)

        def save_and_close():
            # Ensure latest edited labels are synchronized before persisting.
            refresh_group_labels()

            for gid, pattern_widget in pattern_vars.items():
                group_patterns[gid] = [
                    p.strip() for p in pattern_widget.get('1.0', tk.END).split('\n') if p.strip()
                ]

            sample_group_mapping.clear()
            for sample_col, group_var in sample_group_vars.items():
                selected_label = group_var.get().strip()
                if selected_label and selected_label in label_to_gid:
                    sample_group_mapping[sample_col] = label_to_gid[selected_label]

            assigned = len(sample_group_mapping)
            status_label.config(
                text=f'✓ {len(group_definitions)} groups, {assigned}/{len(sample_columns)} samples assigned',
                fg='green'
            )
            if callable(on_save):
                try:
                    on_save({
                        'groups': len(group_definitions),
                        'assigned': assigned,
                        'total_samples': len(sample_columns),
                    })
                except Exception:
                    pass
            pattern_window.destroy()

        tk.Button(bottom_btn_frame, text='Apply Patterns', command=apply_patterns,
                 bg='#4CAF50', fg='white', font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        tk.Button(bottom_btn_frame, text='Refresh Labels', command=refresh_group_labels,
             bg='#8e44ad', fg='white', font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        tk.Button(bottom_btn_frame, text='Done', command=save_and_close,
                 bg='#3498db', fg='white', font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        tk.Button(bottom_btn_frame, text='Cancel', command=pattern_window.destroy,
                 bg='#e74c3c', fg='white', font=('Arial', 10, 'bold')).pack(side='right', padx=5)

    def select_pie_sample_columns(self):
        """Use shared Statistics Verify Columns dialog to select pie sample columns."""
        if self.pie_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return

        selected_columns = self._verify_sample_columns_with_statistics_dialog(
            self.pie_data,
            'Pie Chart - Verify Columns'
        )
        if selected_columns is None:
            return

        # Check if columns have changed
        columns_changed = set(selected_columns) != set(self.pie_sample_columns)
        
        self.pie_sample_columns = selected_columns
        if self.pie_sample_columns:
            self.pie_sample_label.config(text=f'✓ {len(self.pie_sample_columns)} columns verified', fg='green')
            preview = ', '.join(self.pie_sample_columns[:8])
            suffix = ' ...' if len(self.pie_sample_columns) > 8 else ''
            self._pie_log(f'Verified sample columns: {len(self.pie_sample_columns)} selected ({preview}{suffix}).')
            
            # If columns changed, update group mappings to only include valid columns
            if columns_changed and self.pie_sample_group_mapping:
                valid_samples = set(self.pie_sample_columns)
                # Remove mappings for columns that no longer exist
                removed = [s for s in self.pie_sample_group_mapping.keys() if s not in valid_samples]
                for sample in removed:
                    del self.pie_sample_group_mapping[sample]
                
                if removed:
                    self._pie_log(f'⚠ Removed {len(removed)} old group mappings (columns no longer exist).')
                    # Update status to show configuration needs refresh
                    assigned = len(self.pie_sample_group_mapping)
                    self.pie_group_label.config(
                        text=f'⚠ {len(self.pie_group_definitions)} groups, {assigned}/{len(self.pie_sample_columns)} samples assigned',
                        fg='orange'
                    )
        else:
            self.pie_sample_label.config(text='No columns verified', fg='red')
            self.pie_sample_group_mapping.clear()
            self.pie_group_label.config(text='No groups configured', fg='gray')
            self._pie_log('Verify Columns returned no sample columns.')

    def configure_pie_groups(self):
        """Configure pie groups using Statistics-style Auto-Assign Groups dialog."""
        self._pie_log('Opening group configuration dialog...')
        self._open_statistics_style_group_dialog(
            sample_columns=self.pie_sample_columns,
            group_definitions=self.pie_group_definitions,
            group_patterns=self.pie_group_patterns,
            sample_group_mapping=self.pie_sample_group_mapping,
            status_label=self.pie_group_label,
            on_save=lambda info: self._pie_log(
                f"Group configuration saved: {info.get('groups', 0)} groups, "
                f"{info.get('assigned', 0)}/{info.get('total_samples', 0)} samples assigned."
            ),
        )
    
    def generate_pie_chart(self):
        """Generate pie charts"""
        if self.pie_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return
        
        class_col = self.pie_class_var.get()
        if not class_col:
            messagebox.showwarning('Missing Column', 'Please select a class column.')
            return
        
        if not self.pie_sample_columns:
            messagebox.showwarning('No Samples', 'Please select sample columns.')
            return
        
        if not self.pie_sample_group_mapping:
            messagebox.showwarning('No Groups', 'Please configure groups first.')
            return
        
        try:
            self.pie_status_label.config(text='⏳ Generating pie charts...', fg='blue')
            self.frame.update()
            self._pie_log('Starting pie chart generation...')
            
            df = self.pie_data.copy()
            self._pie_log(f'Using class column: {class_col}')
            self._pie_log(f'Using {len(self.pie_sample_columns)} verified sample columns.')
            
            # Group samples by group
            groups = {}
            for sample, gid in self.pie_sample_group_mapping.items():
                if gid not in groups:
                    groups[gid] = []
                groups[gid].append(sample)
            self._pie_log(f'Detected {len(groups)} configured groups in mapping.')
            
            # Calculate percentages for each group
            group_data = {}
            for gid, samples in groups.items():
                # Get numeric columns
                numeric_samples = []
                for sample in samples:
                    if sample in df.columns:
                        numeric_samples.append(sample)
                
                if not numeric_samples:
                    continue
                
                # Calculate mean across samples
                mean_values = df[numeric_samples].mean(axis=1)
                
                # Group by class
                class_totals = df.groupby(class_col).apply(
                    lambda x: x[numeric_samples].mean(axis=1).mean()
                )
                
                # Convert to percentage
                total = class_totals.sum()
                if total > 0:
                    percentages = (class_totals / total * 100).to_dict()
                    group_data[gid] = percentages
                    self._pie_log(f'Group {self.pie_group_definitions.get(gid, gid)}: {len(percentages)} class slices.')
            
            # Generate pie charts
            num_groups = len(group_data)
            fig_width = self.pie_fig_width_var.get()
            fig_height = self.pie_fig_height_var.get()
            fig, axes = plt.subplots(1, num_groups, figsize=(fig_width*num_groups, fig_height))
            
            if num_groups == 1:
                axes = [axes]
            
            # Get label style and font sizes
            label_style = self.pie_label_style_var.get()
            title_fontsize = self.pie_title_fontsize_var.get()
            label_fontsize = self.pie_label_fontsize_var.get()
            
            # DEBUG: Log the label style being used
            self._pie_log(f'DEBUG: Generating pie charts with label_style={repr(label_style)}, autopct will be {"None" if label_style in ["legend", "remove"] else "%1.1f%%"}')
            
            for idx, (gid, data) in enumerate(group_data.items()):
                ax = axes[idx]
                
                classes = list(data.keys())
                values = list(data.values())
                
                # Determine pie parameters based on label style
                if label_style == 'remove':
                    labels = None
                    autopct = None
                elif label_style == 'legend':
                    # No labels on pie, percentages only in legend
                    labels = None
                    autopct = None
                else:  # inside (default)
                    labels = classes
                    autopct = '%1.1f%%'
                
                # Create pie chart
                pie_result = ax.pie(
                    values,
                    labels=labels,
                    autopct=autopct,
                    startangle=90,
                    textprops={'fontsize': label_fontsize, 'weight': 'bold'}
                )
                wedges = pie_result[0]
                
                # Style wedges without hatch patterns (improves readability)
                for wedge in wedges:
                    wedge.set_edgecolor('black')
                    wedge.set_linewidth(1.5)
                
                # Add legend if requested
                if label_style == 'legend':
                    # Legend format: "ClassName X.X%"
                    legend_labels = [f'{cls} {val:.1f}%' for cls, val in zip(classes, values)]
                    ax.legend(legend_labels, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=label_fontsize, frameon=True, ncol=1)
                
                # Title
                group_label = self.pie_group_definitions.get(gid, gid)
                ax.set_title(group_label, fontsize=title_fontsize, fontweight='bold')
            
            # Adjust layout to reserve space for legends on the right (25% of figure width)
            if label_style == 'legend':
                fig.subplots_adjust(right=0.75)
                plt.tight_layout(rect=(0, 0, 0.75, 1))
            else:
                plt.tight_layout()
            
            # Save option
            save_by_group = self.pie_save_by_group_var.get()
            
            if save_by_group:
                # Save each group separately
                base_dir = filedialog.askdirectory(title='Select Directory to Save Group Charts')
                if base_dir:
                    try:
                        for idx, (gid, data) in enumerate(group_data.items()):
                            group_label = self.pie_group_definitions.get(gid, gid)
                            filename = f'pie_chart_{group_label.replace(" ", "_")}.png'
                            filepath = os.path.join(base_dir, filename)
                            
                            # Create individual figure for this group
                            fig_single, ax_single = plt.subplots(figsize=(self.pie_fig_width_var.get(), 
                                                                          self.pie_fig_height_var.get()))
                            classes = list(data.keys())
                            values = list(data.values())
                            
                            if label_style == 'remove':
                                labels = None
                                autopct = None
                            elif label_style == 'legend':
                                # No labels or percentages on pie - legend shows "ClassName X.X%"
                                labels = None
                                autopct = None
                            else:  # 'inside' (default)
                                labels = classes
                                autopct = '%1.1f%%'
                            
                            pie_result = ax_single.pie(
                                values,
                                labels=labels,
                                autopct=autopct,
                                startangle=90,
                                textprops={'fontsize': label_fontsize, 'weight': 'bold'}
                            )
                            wedges = pie_result[0]
                            
                            for wedge in wedges:
                                wedge.set_edgecolor('black')
                                wedge.set_linewidth(1.5)
                            
                            if label_style == 'legend':
                                # Legend format: "ClassName X.X%" - must match combined chart format
                                legend_labels = [f'{cls} {val:.1f}%' for cls, val in zip(classes, values)]
                                ax_single.legend(legend_labels, loc='upper left', bbox_to_anchor=(1.05, 1), 
                                              fontsize=label_fontsize, frameon=True, ncol=1)
                                # Adjust layout to reserve space for legend (right 25% of figure width)
                                fig_single.subplots_adjust(right=0.75)
                            
                            ax_single.set_title(group_label, fontsize=title_fontsize, fontweight='bold')
                            fig_single.tight_layout(rect=(0, 0, 0.75, 1) if label_style == 'legend' else None)
                            fig_single.savefig(filepath, dpi=300, bbox_inches='tight')
                            plt.close(fig_single)
                            self._pie_log(f'Saved group pie chart: {os.path.basename(filepath)}')
                        
                        self.pie_status_label.config(text=f'✓ Saved {len(group_data)} charts by group', fg='green')
                        self._pie_log(f'Completed pie chart export: {len(group_data)} files.')
                    except Exception as e:
                        messagebox.showerror('Error', f'Failed to save charts:\n{str(e)}')
                        self._pie_log(f'Failed to save group pie charts: {e}')
            else:
                # Save all as one figure
                save_path = filedialog.asksaveasfilename(
                    defaultextension='.png',
                    filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All Files', '*.*')]
                )
                
                if save_path:
                    fig.savefig(save_path, dpi=300, bbox_inches='tight')
                    self.pie_status_label.config(text=f'✓ Saved to {os.path.basename(save_path)}', fg='green')
                    self._pie_log(f'Saved combined pie chart figure: {os.path.basename(save_path)}')
                else:
                    self.pie_status_label.config(text='✓ Charts generated successfully!', fg='green')
                    self._pie_log('Pie charts generated (not saved).')
            
        except Exception as e:
            messagebox.showerror('Error', f'Failed to generate pie charts:\n{str(e)}')
            logger.error(f'Error generating pie charts: {e}', exc_info=True)
            self.pie_status_label.config(text='✗ Generation failed', fg='red')
            self._pie_log(f'Pie chart generation failed: {e}')
    
    # ==================== HEATMAP TAB ====================
    def setup_heatmap_tab(self):
        """Setup Heatmap sub-tab (similar to pie chart but generates heatmap)"""
        heatmap_frame = ttk.Frame(self.notebook)
        self.notebook.add(heatmap_frame, text='🔥 Heatmap')

        body = tk.Frame(heatmap_frame, bg='#f0f0f0')
        body.pack(fill='both', expand=True, padx=10, pady=10)

        # Scrollable container so the heatmap controls can be reached via the right scrollbar.
        scroll_host = tk.Frame(body, bg='#f0f0f0')
        scroll_host.pack(side='left', fill='both', expand=True)

        canvas = tk.Canvas(scroll_host, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(scroll_host, orient='vertical', command=canvas.yview)
        container = tk.Frame(canvas, bg='#f0f0f0')
        container.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))

        canvas_window = canvas.create_window((0, 0), window=container, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)

        def _fit_heatmap_container_width(event):
            canvas.itemconfigure(canvas_window, width=event.width)

        canvas.bind('<Configure>', _fit_heatmap_container_width)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        heatmap_log_frame = tk.LabelFrame(body, text='Heatmap Log', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        heatmap_log_frame.pack(side='right', fill='y', padx=(10, 0))
        self.heatmap_log_text = tk.Text(heatmap_log_frame, width=45, wrap='word', font=('Consolas', 9), bg='white')
        self.heatmap_log_text.pack(side='left', fill='both', expand=True, padx=(6, 0), pady=6)
        heatmap_log_scroll = ttk.Scrollbar(heatmap_log_frame, orient='vertical', command=self.heatmap_log_text.yview)
        heatmap_log_scroll.pack(side='right', fill='y', padx=(0, 6), pady=6)
        self.heatmap_log_text.configure(yscrollcommand=heatmap_log_scroll.set)
        
        # Title
        title_label = tk.Label(container, text='Heatmap Generator', 
                               font=('Arial', 14, 'bold'), bg='#f0f0f0')
        title_label.pack(pady=10)
        
        # Upload section
        upload_frame = tk.LabelFrame(container, text='Step 1: Upload Data', 
                                     bg='#f0f0f0', font=('Arial', 10, 'bold'))
        upload_frame.pack(fill='x', pady=5)
        
        btn_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Button(btn_frame, text='📂 Upload Excel File', 
                 command=self.upload_heatmap_data, bg='#4CAF50', fg='white',
                 font=('Arial', 10, 'bold')).pack(side='left', padx=5)
        
        self.heatmap_file_label = tk.Label(btn_frame, text='No file uploaded', 
                                           bg='#f0f0f0', font=('Arial', 9))
        self.heatmap_file_label.pack(side='left', padx=5)
        
        # Column mapping
        mapping_frame = tk.LabelFrame(container, text='Step 2: Map Columns', 
                                      bg='#f0f0f0', font=('Arial', 10, 'bold'))
        mapping_frame.pack(fill='x', pady=5)
        
        # Class column
        class_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        class_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(class_frame, text='Class Column:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        self.heatmap_class_var = tk.StringVar()
        self.heatmap_class_combo = ttk.Combobox(class_frame, 
                                                 textvariable=self.heatmap_class_var, 
                                                 state='readonly', width=30)
        self.heatmap_class_combo.pack(side='left', padx=5)
        
        # Sample columns selection
        sample_frame = tk.Frame(mapping_frame, bg='#f0f0f0')
        sample_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(sample_frame, text='Sample Columns:', bg='#f0f0f0', 
                width=20, anchor='w').pack(side='left')
        tk.Button(sample_frame, text='🔍 Verify Columns', 
                 command=self.select_heatmap_sample_columns, bg='#4CAF50', fg='white').pack(side='left', padx=5)
        self.heatmap_sample_label = tk.Label(sample_frame, text='No columns selected', 
                                             bg='#f0f0f0', font=('Arial', 9))
        self.heatmap_sample_label.pack(side='left', padx=5)
        
        # Group configuration
        group_frame = tk.LabelFrame(container, text='Step 3: Configure Groups', 
                                    bg='#f0f0f0', font=('Arial', 10, 'bold'))
        group_frame.pack(fill='x', pady=5)
        
        tk.Button(group_frame, text='⚙️ Configure Groups & Patterns', 
                 command=self.configure_heatmap_groups,
                 bg='#FF9800', fg='white',
                 font=('Arial', 10, 'bold')).pack(padx=5, pady=5)
        
        self.heatmap_group_label = tk.Label(group_frame, text='No groups configured', 
                                            bg='#f0f0f0', font=('Arial', 9))
        self.heatmap_group_label.pack(padx=5, pady=5)
        
        # Settings frame (collapsible)
        settings_frame = tk.LabelFrame(container, text='⚙️ Settings', 
                                       bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=5)
        
        # Font size controls
        font_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        font_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(font_frame, text='Title Font Size:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.heatmap_title_fontsize_var = tk.IntVar(value=14)
        tk.Spinbox(font_frame, from_=10, to=24, textvariable=self.heatmap_title_fontsize_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(font_frame, text='Label Font Size:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.heatmap_label_fontsize_var = tk.IntVar(value=10)
        tk.Spinbox(font_frame, from_=6, to=20, textvariable=self.heatmap_label_fontsize_var,
                  width=5).pack(side='left', padx=2)
        
        # Image size
        size_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        size_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(size_frame, text='Width:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.heatmap_fig_width_var = tk.IntVar(value=10)
        tk.Spinbox(size_frame, from_=6, to=20, textvariable=self.heatmap_fig_width_var,
                  width=5).pack(side='left', padx=2)
        
        tk.Label(size_frame, text='Height:', bg='#f0f0f0', 
                font=('Arial', 9)).pack(side='left', padx=5)
        self.heatmap_fig_height_var = tk.IntVar(value=8)
        tk.Spinbox(size_frame, from_=4, to=20, textvariable=self.heatmap_fig_height_var,
                  width=5).pack(side='left', padx=2)
        
        # Save options
        save_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        save_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(save_frame, text='Save:', bg='#f0f0f0', 
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.heatmap_save_by_group_var = tk.BooleanVar(value=False)
        tk.Checkbutton(save_frame, text='Save Separately by Group', 
                      variable=self.heatmap_save_by_group_var, bg='#f0f0f0').pack(side='left', padx=5)

        # Normalization options
        norm_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        norm_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(norm_frame, text='Normalization:', bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.heatmap_normalization_var = tk.StringVar(value='zscore')
        for norm in ['Raw', 'Z-Score']:
            tk.Radiobutton(norm_frame, text=norm, variable=self.heatmap_normalization_var,
                          value=norm.lower().replace('-', ''), bg='#f0f0f0').pack(side='left', padx=5)

        # Color scheme options
        color_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        color_frame.pack(fill='x', padx=5, pady=5)
        tk.Label(color_frame, text='Color Scheme:', bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(side='left', padx=5)
        self.heatmap_colormap_var = tk.StringVar(value='RdYlGn_r')
        colormap_options = ['RdYlGn_r', 'RdBu_r', 'coolwarm', 'seismic', 'viridis', 'plasma', 'inferno']
        colormap_combo = ttk.Combobox(color_frame, textvariable=self.heatmap_colormap_var,
                                      values=colormap_options, state='readonly', width=15)
        colormap_combo.pack(side='left', padx=5)

        # Generate button
        tk.Button(container, text='🎨 Generate Heatmap',
                 command=self.generate_heatmap,
                 bg='#2196F3', fg='white',
                 font=('Arial', 12, 'bold'), height=2).pack(pady=10)

        # Status label
        self.heatmap_status_label = tk.Label(container, text='',
                                             bg='#f0f0f0', font=('Arial', 9),
                                             fg='blue')
        self.heatmap_status_label.pack(pady=5)

        # Initialize data storage
        self.heatmap_data = None
        self.heatmap_sample_columns = []
        self.heatmap_group_definitions = {}
        self.heatmap_sample_group_mapping = {}
        self.heatmap_group_patterns = {}
        self._heatmap_log('Ready. Upload data, verify columns, configure groups, then generate heatmap.')

    def _heatmap_log(self, message: str):
        """Append message to heatmap log panel."""
        try:
            timestamp = datetime.now().strftime('%H:%M:%S')
            self.heatmap_log_text.insert('end', f'[{timestamp}] {message}\n')
            self.heatmap_log_text.see('end')
            self.heatmap_log_text.update_idletasks()
        except Exception:
            pass
    
    def upload_heatmap_data(self):
        """Upload Excel file for heatmap"""
        file_path = filedialog.askopenfilename(
            title='Select Excel File',
            filetypes=[('Excel Files', '*.xlsx *.xls'), ('All Files', '*.*')]
        )

        if not file_path:
            return

        try:
            df = pd.read_excel(file_path)
            self.heatmap_data = df

            filename = os.path.basename(file_path)
            self.heatmap_file_label.config(text=f'✓ {filename} ({len(df)} rows)', fg='green')

            columns = [''] + list(df.columns)
            self.heatmap_class_combo['values'] = columns

            cols_lower = {c.lower(): c for c in df.columns}
            for candidate in ['class', 'classification', 'category', 'type', 'pathway', 'name', 'metabolite']:
                if candidate in cols_lower:
                    self.heatmap_class_var.set(cols_lower[candidate])
                    break

            # Reset column-specific configuration when new file is uploaded
            self.heatmap_sample_columns = []
            self.heatmap_sample_label.config(text='No columns selected', fg='red')
            self.heatmap_sample_group_mapping.clear()  # Clear old mappings
            self.heatmap_group_label.config(text='No groups configured', fg='gray')

            self.heatmap_status_label.config(text='✓ File uploaded successfully!', fg='green')
            self._heatmap_log(f'Loaded file: {filename} ({len(df)} rows, {len(df.columns)} columns).')
            if self.heatmap_class_var.get():
                self._heatmap_log(f'Auto-detected class column: {self.heatmap_class_var.get()}')

        except Exception as e:
            messagebox.showerror('Error', f'Failed to load file:\n{str(e)}')
            logger.error(f'Error loading heatmap data: {e}')
            self._heatmap_log(f'Failed to load heatmap data: {e}')

    def select_heatmap_sample_columns(self):
        """Use shared Statistics Verify Columns dialog to select heatmap sample columns."""
        if self.heatmap_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return

        selected_columns = self._verify_sample_columns_with_statistics_dialog(
            self.heatmap_data,
            'Heatmap - Verify Columns'
        )
        if selected_columns is None:
            return

        # Check if columns have changed
        columns_changed = set(selected_columns) != set(self.heatmap_sample_columns)
        
        self.heatmap_sample_columns = selected_columns
        if self.heatmap_sample_columns:
            self.heatmap_sample_label.config(text=f'✓ {len(self.heatmap_sample_columns)} columns verified', fg='green')
            preview = ', '.join(self.heatmap_sample_columns[:8])
            suffix = ' ...' if len(self.heatmap_sample_columns) > 8 else ''
            self._heatmap_log(f'Verified sample columns: {len(self.heatmap_sample_columns)} selected ({preview}{suffix}).')
            
            # If columns changed, update group mappings to only include valid columns
            if columns_changed and self.heatmap_sample_group_mapping:
                valid_samples = set(self.heatmap_sample_columns)
                # Remove mappings for columns that no longer exist
                removed = [s for s in list(self.heatmap_sample_group_mapping.keys()) if s not in valid_samples]
                for sample in removed:
                    del self.heatmap_sample_group_mapping[sample]
                
                if removed:
                    self._heatmap_log(f'⚠ Removed {len(removed)} old group mappings (columns no longer exist).')
                    # Update status to show configuration needs refresh
                    assigned = len(self.heatmap_sample_group_mapping)
                    if self.heatmap_group_definitions:
                        self.heatmap_group_label.config(
                            text=f'⚠ {len(self.heatmap_group_definitions)} groups, {assigned}/{len(self.heatmap_sample_columns)} samples assigned',
                            fg='orange'
                        )
        else:
            self.heatmap_sample_label.config(text='No columns verified', fg='red')
            self.heatmap_sample_group_mapping.clear()
            self.heatmap_group_label.config(text='No groups configured', fg='gray')
            self._heatmap_log('Verify Columns returned no sample columns for heatmap.')

    def configure_heatmap_groups(self):
        """Configure heatmap groups using Statistics-style Auto-Assign Groups dialog."""
        self._heatmap_log('Opening group configuration dialog...')
        self._open_statistics_style_group_dialog(
            sample_columns=self.heatmap_sample_columns,
            group_definitions=self.heatmap_group_definitions,
            group_patterns=self.heatmap_group_patterns,
            sample_group_mapping=self.heatmap_sample_group_mapping,
            status_label=self.heatmap_group_label,
            on_save=lambda info: self._heatmap_log(
                f"Group configuration saved: {info.get('groups', 0)} groups, "
                f"{info.get('assigned', 0)}/{info.get('total_samples', 0)} samples assigned."
            ),
        )

    def generate_heatmap(self):
        """Generate heatmap"""
        if self.heatmap_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return

        class_col = self.heatmap_class_var.get()
        if not class_col:
            messagebox.showwarning('Missing Column', 'Please select a class column.')
            return

        if not self.heatmap_sample_columns:
            messagebox.showwarning('No Samples', 'Please select sample columns.')
            return

        if not self.heatmap_sample_group_mapping:
            messagebox.showwarning('No Groups', 'Please configure groups first.')
            return

        try:
            import seaborn as sns

            self.heatmap_status_label.config(text='⏳ Generating heatmap...', fg='blue')
            self.frame.update()
            self._heatmap_log('Starting heatmap generation...')

            df = self.heatmap_data.copy()
            self._heatmap_log(f'Class column: {class_col}')

            # Enforce sample ordering by configured groups (prevents scattered columns).
            valid_samples = [c for c in self.heatmap_sample_columns if c in df.columns]
            group_order = list(self.heatmap_group_definitions.keys())
            ordered_samples = []
            for gid in group_order:
                ordered_samples.extend([c for c in valid_samples if self.heatmap_sample_group_mapping.get(c) == gid])
            # Keep any mapped samples for unknown groups, then any unassigned leftovers.
            ordered_samples.extend([c for c in valid_samples if c not in ordered_samples and c in self.heatmap_sample_group_mapping])
            ordered_samples.extend([c for c in valid_samples if c not in ordered_samples])
            if not ordered_samples:
                raise ValueError('No valid sample columns available for heatmap after verification.')

            self._heatmap_log(f'Using {len(ordered_samples)} ordered samples based on group configuration.')

            group_block_sizes = []
            for gid in group_order:
                block = [c for c in ordered_samples if self.heatmap_sample_group_mapping.get(c) == gid]
                if block:
                    group_block_sizes.append((gid, len(block)))
            if group_block_sizes:
                block_text = ', '.join(
                    f"{self.heatmap_group_definitions.get(gid, gid)}={size}" for gid, size in group_block_sizes
                )
                self._heatmap_log(f'Group sample blocks: {block_text}')

            df_heat = df[[class_col] + ordered_samples].copy()

            for col in ordered_samples:
                df_heat[col] = pd.to_numeric(df_heat[col], errors='coerce')

            # Aggregate duplicate class rows while preserving first-seen class order.
            df_grouped = df_heat.set_index(class_col).groupby(level=0, sort=False).agg(['mean', 'sem'])
            df_heat_agg = df_grouped.xs('mean', level=1, axis=1)
            if isinstance(df_heat_agg, pd.Series):
                df_heat_agg = df_heat_agg.to_frame()
            df_heat_agg.columns = [col[0] if isinstance(col, tuple) else col for col in df_heat_agg.columns]
            df_heat_agg = df_heat_agg.dropna(how='all')
            self._heatmap_log(f'Class rows after aggregation: {len(df_heat_agg)}')

            # Apply selected normalization while keeping configured sample order.
            norm_mode = (self.heatmap_normalization_var.get() or 'zscore').strip().lower()
            if norm_mode == 'raw':
                df_heat_plot = df_heat_agg.copy()
                cbar_label = 'Group Mean (Raw)'
                center_value = None
            else:
                df_heat_plot = df_heat_agg.apply(
                    lambda row: (row - row.mean()) / row.std() if row.std() > 0 else row,
                    axis=1
                )
                if isinstance(df_heat_plot, pd.Series):
                    df_heat_plot = df_heat_plot.to_frame()
                cbar_label = 'Z-Score (Group Mean)'
                center_value = 0
            self._heatmap_log(f'Normalization mode: {norm_mode}')

            fig_width = self.heatmap_fig_width_var.get()
            fig_height = self.heatmap_fig_height_var.get()
            colormap = self.heatmap_colormap_var.get()

            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(
                df_heat_plot,
                cmap=colormap,
                center=center_value,
                cbar_kws={'label': cbar_label},
                linewidths=0.5,
                linecolor='gray',
                ax=ax
            )

            title_fontsize = self.heatmap_title_fontsize_var.get()
            plt.title('Heatmap (Samples grouped by configured groups)', fontsize=title_fontsize, fontweight='bold')
            plt.ylabel('Glycan Class')
            plt.xlabel('Sample Groups')

            # Add separators between group blocks for easier visual comparison.
            boundaries = []
            running = 0
            for gid, size in group_block_sizes:
                running += size
                boundaries.append(running)
            y0, y1 = ax.get_ylim()
            for boundary in boundaries[:-1]:
                ax.vlines(boundary, y0, y1, colors='black', linewidth=2.0, alpha=0.9)

            # Add subtle separators for each class row.
            x0, x1 = ax.get_xlim()
            for row_boundary in range(1, len(df_heat_plot.index)):
                ax.hlines(row_boundary, x0, x1, colors='white', linewidth=0.8, alpha=0.85)

            self._heatmap_log('Applied group separators and class row separators to heatmap.')
            plt.tight_layout()

            save_path = filedialog.asksaveasfilename(
                defaultextension='.png',
                filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All Files', '*.*')]
            )

            if save_path:
                fig.savefig(save_path, dpi=300, bbox_inches='tight')
                self.heatmap_status_label.config(text=f'✓ Saved to {os.path.basename(save_path)}', fg='green')
                self._heatmap_log(f'Saved heatmap: {os.path.basename(save_path)}')
            else:
                self.heatmap_status_label.config(text='✓ Heatmap generated successfully!', fg='green')
                self._heatmap_log('Heatmap generated (not saved).')

        except Exception as e:
            messagebox.showerror('Error', f'Failed to generate heatmap:\n{str(e)}')
            logger.error(f'Error generating heatmap: {e}', exc_info=True)
            self.heatmap_status_label.config(text='✗ Generation failed', fg='red')
            self._heatmap_log(f'Heatmap generation failed: {e}')

    # ==================== EFFECT SIZE TAB ====================
    def setup_effect_size_plot_tab(self):
        """Setup single-tab effect-size plot utility with optional row selection."""
        effect_frame = ttk.Frame(self.notebook)
        self.notebook.add(effect_frame, text='📉 Effect Plot')

        scroll_host = tk.Frame(effect_frame, bg='#f0f0f0')
        scroll_host.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(scroll_host, bg='#f0f0f0', highlightthickness=0)
        v_scroll = ttk.Scrollbar(scroll_host, orient='vertical', command=canvas.yview)
        container = tk.Frame(canvas, bg='#f0f0f0')
        container.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas_window = canvas.create_window((0, 0), window=container, anchor='nw')
        canvas.configure(yscrollcommand=v_scroll.set)

        def _fit_effect_container_width(event):
            canvas.itemconfigure(canvas_window, width=event.width)

        canvas.bind('<Configure>', _fit_effect_container_width)
        canvas.pack(side='left', fill='both', expand=True)
        v_scroll.pack(side='right', fill='y')

        tk.Label(container, text='Effect Size Plot Generator', font=('Arial', 14, 'bold'), bg='#f0f0f0').pack(pady=8)

        upload_frame = tk.LabelFrame(container, text='Step 1: Upload Data', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        upload_frame.pack(fill='x', pady=5)
        row = tk.Frame(upload_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=6, pady=6)
        tk.Button(
            row,
            text='Upload Excel/CSV',
            command=self.upload_effect_size_data,
            bg='#4CAF50',
            fg='white',
            font=('Arial', 10, 'bold'),
        ).pack(side='left', padx=5)
        self.effect_file_label = tk.Label(row, text='No file uploaded', bg='#f0f0f0', font=('Arial', 9))
        self.effect_file_label.pack(side='left', padx=6)

        map_frame = tk.LabelFrame(container, text='Step 2: Map Required Columns', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        map_frame.pack(fill='x', pady=5)

        self.effect_group_var = tk.StringVar()
        self.effect_estimate_var = tk.StringVar()
        self.effect_sig_var = tk.StringVar()
        self.effect_ci_low_var = tk.StringVar()
        self.effect_ci_high_var = tk.StringVar()

        def _add_combo(parent, label, var):
            r = tk.Frame(parent, bg='#f0f0f0')
            r.pack(fill='x', padx=6, pady=3)
            tk.Label(r, text=label, bg='#f0f0f0', width=28, anchor='w').pack(side='left')
            combo = ttk.Combobox(r, textvariable=var, state='readonly', width=46)
            combo.pack(side='left', padx=5)
            return combo

        self.effect_group_combo = _add_combo(map_frame, 'Grouping variable (X):', self.effect_group_var)
        self.effect_estimate_combo = _add_combo(map_frame, 'Effect size / Estimate (Y):', self.effect_estimate_var)
        self.effect_sig_combo = _add_combo(map_frame, 'Significance (adj p):', self.effect_sig_var)
        self.effect_ci_low_combo = _add_combo(map_frame, 'CI lower (95%):', self.effect_ci_low_var)
        self.effect_ci_high_combo = _add_combo(map_frame, 'CI upper (95%):', self.effect_ci_high_var)
        self.effect_group_combo.bind('<<ComboboxSelected>>', lambda _e: self._refresh_effect_row_selector())

        settings_frame = tk.LabelFrame(container, text='Step 3: Plot Options', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        settings_frame.pack(fill='x', pady=5)

        row1 = tk.Frame(settings_frame, bg='#f0f0f0')
        row1.pack(fill='x', padx=6, pady=4)
        tk.Label(row1, text='X-axis label:', bg='#f0f0f0', width=16, anchor='w').pack(side='left')
        self.effect_xlabel_var = tk.StringVar(value='Lipid Class')
        tk.Entry(row1, textvariable=self.effect_xlabel_var, width=34).pack(side='left', padx=5)
        tk.Label(row1, text='Y-axis label:', bg='#f0f0f0', width=16, anchor='w').pack(side='left', padx=(12, 0))
        self.effect_ylabel_var = tk.StringVar(value='Effect Size')
        tk.Entry(row1, textvariable=self.effect_ylabel_var, width=34).pack(side='left', padx=5)

        row2 = tk.Frame(settings_frame, bg='#f0f0f0')
        row2.pack(fill='x', padx=6, pady=4)
        tk.Label(row2, text='Significance alpha:', bg='#f0f0f0', width=16, anchor='w').pack(side='left')
        self.effect_alpha_var = tk.StringVar(value='0.05')
        tk.Entry(row2, textvariable=self.effect_alpha_var, width=8).pack(side='left', padx=5)
        self.effect_show_pvalues_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row2,
            text='Show p-values',
            variable=self.effect_show_pvalues_var,
            bg='#f0f0f0'
        ).pack(side='left', padx=(10, 0))
        tk.Label(row2, text='Figure width:', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.effect_fig_w_var = tk.IntVar(value=11)
        tk.Spinbox(row2, from_=6, to=30, textvariable=self.effect_fig_w_var, width=5).pack(side='left', padx=4)
        tk.Label(row2, text='Figure height:', bg='#f0f0f0').pack(side='left', padx=(8, 0))
        self.effect_fig_h_var = tk.IntVar(value=8)
        tk.Spinbox(row2, from_=4, to=40, textvariable=self.effect_fig_h_var, width=5).pack(side='left', padx=4)

        row3 = tk.Frame(settings_frame, bg='#f0f0f0')
        row3.pack(fill='x', padx=6, pady=4)
        tk.Label(row3, text='Sort by:', bg='#f0f0f0', width=16, anchor='w').pack(side='left')
        self.effect_sort_var = tk.StringVar(value='Effect Size')
        sort_combo = ttk.Combobox(row3, textvariable=self.effect_sort_var, state='readonly', width=18)
        sort_combo['values'] = ('Effect Size', 'Name (Alphabetical)', 'P-value')
        sort_combo.pack(side='left', padx=5)
        self.effect_limit_top_p_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row3,
            text='Show only Top X by p-value',
            variable=self.effect_limit_top_p_var,
            command=self._toggle_effect_top_n_state,
            bg='#f0f0f0'
        ).pack(side='left', padx=(10, 0))
        self.effect_top_n_p_var = tk.IntVar(value=10)
        self.effect_top_n_spin = tk.Spinbox(row3, from_=1, to=500, textvariable=self.effect_top_n_p_var, width=5)
        self.effect_top_n_spin.pack(side='left', padx=4)
        self._toggle_effect_top_n_state()
        self.effect_show_vlines_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            row3,
            text='Background vertical lines',
            variable=self.effect_show_vlines_var,
            bg='#f0f0f0'
        ).pack(side='left', padx=(10, 0))
        tk.Label(row3, text='X-label font:', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.effect_xlabel_fontsize_var = tk.IntVar(value=12)
        tk.Spinbox(row3, from_=8, to=24, textvariable=self.effect_xlabel_fontsize_var, width=4).pack(side='left', padx=3)
        tk.Label(row3, text='Y-label font:', bg='#f0f0f0').pack(side='left', padx=(8, 0))
        self.effect_ylabel_fontsize_var = tk.IntVar(value=12)
        tk.Spinbox(row3, from_=8, to=24, textvariable=self.effect_ylabel_fontsize_var, width=4).pack(side='left', padx=3)

        row4 = tk.Frame(settings_frame, bg='#f0f0f0')
        row4.pack(fill='x', padx=6, pady=4)
        tk.Label(row4, text='X-tick font:', bg='#f0f0f0', width=16, anchor='w').pack(side='left')
        self.effect_xtick_fontsize_var = tk.IntVar(value=10)
        tk.Spinbox(row4, from_=7, to=20, textvariable=self.effect_xtick_fontsize_var, width=4).pack(side='left', padx=5)
        tk.Label(row4, text='Y-tick font:', bg='#f0f0f0').pack(side='left', padx=(8, 0))
        self.effect_ytick_fontsize_var = tk.IntVar(value=10)
        tk.Spinbox(row4, from_=7, to=20, textvariable=self.effect_ytick_fontsize_var, width=4).pack(side='left', padx=3)
        tk.Label(row4, text='X-tick thickness:', bg='#f0f0f0').pack(side='left', padx=(12, 0))
        self.effect_xtick_width_var = tk.DoubleVar(value=1.5)
        tk.Spinbox(row4, from_=0.5, to=4.0, increment=0.25, textvariable=self.effect_xtick_width_var, width=4).pack(side='left', padx=3)
        tk.Label(row4, text='Y-tick thickness:', bg='#f0f0f0').pack(side='left', padx=(8, 0))
        self.effect_ytick_width_var = tk.DoubleVar(value=1.5)
        tk.Spinbox(row4, from_=0.5, to=4.0, increment=0.25, textvariable=self.effect_ytick_width_var, width=4).pack(side='left', padx=3)

        select_frame = tk.LabelFrame(container, text='Optional: Plot Selected Rows Only', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        select_frame.pack(fill='both', expand=True, pady=5)
        self.effect_use_selected_rows_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            select_frame,
            text='Use only selected Grouping-variable rows',
            variable=self.effect_use_selected_rows_var,
            bg='#f0f0f0',
        ).pack(anchor='w', padx=6, pady=(4, 2))

        list_holder = tk.Frame(select_frame, bg='#f0f0f0')
        list_holder.pack(fill='both', expand=True, padx=6, pady=5)
        self.effect_row_listbox = tk.Listbox(list_holder, selectmode='extended', height=10)
        self.effect_row_listbox.pack(side='left', fill='both', expand=True)
        list_scroll = ttk.Scrollbar(list_holder, orient='vertical', command=self.effect_row_listbox.yview)
        list_scroll.pack(side='right', fill='y')
        self.effect_row_listbox.configure(yscrollcommand=list_scroll.set)

        btn_row = tk.Frame(select_frame, bg='#f0f0f0')
        btn_row.pack(fill='x', padx=6, pady=(0, 6))
        tk.Button(
            btn_row,
            text='Select All',
            command=lambda: self.effect_row_listbox.selection_set(0, tk.END),
            bg='#95a5a6',
            fg='white',
        ).pack(side='left', padx=2)
        tk.Button(
            btn_row,
            text='Clear Selection',
            command=lambda: self.effect_row_listbox.selection_clear(0, tk.END),
            bg='#7f8c8d',
            fg='white',
        ).pack(side='left', padx=2)

        tk.Button(
            container,
            text='Generate Effect Plot',
            command=self.generate_effect_size_plot,
            bg='#2196F3',
            fg='white',
            font=('Arial', 12, 'bold'),
            height=2,
        ).pack(pady=10)

        self.effect_status_label = tk.Label(container, text='', bg='#f0f0f0', font=('Arial', 9), fg='blue')
        self.effect_status_label.pack(pady=5)

        self.effect_data = None

    def upload_effect_size_data(self):
        """Upload data for effect-size plotting and apply default column suggestions."""
        file_path = filedialog.askopenfilename(
            title='Select Data File',
            filetypes=[('Excel/CSV', '*.xlsx *.xls *.csv'), ('All Files', '*.*')],
        )
        if not file_path:
            return

        try:
            if file_path.lower().endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                df = pd.read_excel(file_path)

            self.effect_data = df
            self.effect_file_label.config(text=f'OK {os.path.basename(file_path)} ({len(df)} rows)', fg='green')

            cols = [''] + list(df.columns)
            for combo in [
                self.effect_group_combo,
                self.effect_estimate_combo,
                self.effect_sig_combo,
                self.effect_ci_low_combo,
                self.effect_ci_high_combo,
            ]:
                combo['values'] = cols

            self._set_effect_plot_defaults(df)
            self._refresh_effect_row_selector()
            self.effect_status_label.config(text='Data loaded. Defaults set from common lipid-class naming.', fg='green')
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load data:\n{str(e)}')
            logger.error(f'Error loading effect plot data: {e}', exc_info=True)

    def _set_effect_plot_defaults(self, df):
        """Set default mappings using expected lipid-class example names when available.
        
        Prioritizes model-based columns (model_effect, model_ci_*) if available,
        falls back to log2FC and exported CI for backward compatibility.
        """

        def pick(candidates):
            for c in candidates:
                if c in df.columns:
                    return c
            return ''
        
        def pick_comparison_columns():
            """Auto-detect comparison columns (e.g., Control_vs_PD_*) from the data."""
            # Find any _model_effect column or _log2FC column to detect pattern
            for col in df.columns:
                if '_model_effect' in col:
                    # Extract comparison name: "Control_vs_PD" from "Control_vs_PD_model_effect"
                    comp_name = col.replace('_model_effect', '')
                    return comp_name
                elif '_log2FC' in col and '_model' not in col:
                    comp_name = col.replace('_log2FC', '')
                    return comp_name
            return None
        
        comp = pick_comparison_columns()
        
        self.effect_group_var.set(pick(['Class', 'Lipid_Class', 'Class_name', 'Name']))
        
        # Prioritize model-based columns, fall back to log2FC
        if comp:
            self.effect_estimate_var.set(pick([f'{comp}_model_effect', f'{comp}_log2FC', 'log2FC', 'Estimate']))
            self.effect_sig_var.set(pick([f'{comp}_adj_p', 'adj_p', 'p_value_adj']))
            self.effect_ci_low_var.set(pick([f'{comp}_ci_lower_95', 'ci_lower_95', 'CI_Lower_95']))
            self.effect_ci_high_var.set(pick([f'{comp}_ci_upper_95', 'ci_upper_95', 'CI_Upper_95']))
        else:
            self.effect_estimate_var.set(pick(['log2FC', 'Estimate']))
            self.effect_sig_var.set(pick(['adj_p', 'p_value_adj']))
            self.effect_ci_low_var.set(pick(['ci_lower_95', 'CI_Lower_95']))
            self.effect_ci_high_var.set(pick(['ci_upper_95', 'CI_Upper_95']))

    def _refresh_effect_row_selector(self):
        """Populate optional row list from the selected grouping variable column."""
        self.effect_row_listbox.delete(0, tk.END)
        if self.effect_data is None:
            return
        gcol = self.effect_group_var.get().strip()
        if not gcol or gcol not in self.effect_data.columns:
            return

        vals = [str(v) for v in self.effect_data[gcol].dropna().astype(str).tolist() if str(v).strip()]
        unique_vals = sorted(set(vals))
        for v in unique_vals:
            self.effect_row_listbox.insert(tk.END, v)

    def _toggle_effect_top_n_state(self):
        """Enable Top-X spinbox only when Top-X filtering is checked."""
        try:
            if self.effect_limit_top_p_var.get():
                self.effect_top_n_spin.config(state='normal')
            else:
                self.effect_top_n_spin.config(state='disabled')
        except Exception:
            pass

    def generate_effect_size_plot(self):
        """Generate effect-size vs grouping-variable plot with CI and 3-state coloring."""
        if self.effect_data is None:
            messagebox.showwarning('No Data', 'Please upload a data file first.')
            return

        group_col = self.effect_group_var.get().strip()
        est_col = self.effect_estimate_var.get().strip()
        sig_col = self.effect_sig_var.get().strip()
        ci_low_col = self.effect_ci_low_var.get().strip()
        ci_high_col = self.effect_ci_high_var.get().strip()
        required = [group_col, est_col, sig_col, ci_low_col, ci_high_col]
        if any(not c for c in required):
            messagebox.showwarning('Missing Columns', 'Please map all required columns first.')
            return
        missing = [c for c in required if c not in self.effect_data.columns]
        if missing:
            messagebox.showerror('Invalid Columns', f'Missing columns in data:\n{", ".join(missing)}')
            return

        try:
            alpha = float(self.effect_alpha_var.get())
        except Exception:
            alpha = 0.05

        df = self.effect_data.copy()
        df = df[[group_col, est_col, sig_col, ci_low_col, ci_high_col]].copy()
        df[est_col] = pd.to_numeric(df[est_col], errors='coerce')
        df[sig_col] = pd.to_numeric(df[sig_col], errors='coerce')
        df[ci_low_col] = pd.to_numeric(df[ci_low_col], errors='coerce')
        df[ci_high_col] = pd.to_numeric(df[ci_high_col], errors='coerce')
        df[group_col] = df[group_col].astype(str)
        df = df.dropna(subset=[group_col, est_col, sig_col, ci_low_col, ci_high_col]).copy()
        if df.empty:
            messagebox.showwarning('No Rows', 'No valid rows available after numeric conversion.')
            return

        if self.effect_use_selected_rows_var.get():
            idxs = self.effect_row_listbox.curselection()
            selected_vals = {self.effect_row_listbox.get(i) for i in idxs}
            if not selected_vals:
                messagebox.showwarning('No Selection', 'Row filter is enabled, but no rows are selected.')
                return
            df = df[df[group_col].isin(selected_vals)].copy()
            if df.empty:
                messagebox.showwarning('No Rows', 'Selected row filter removed all rows.')
                return

        if self.effect_limit_top_p_var.get():
            try:
                top_n = max(1, int(self.effect_top_n_p_var.get()))
            except Exception:
                top_n = 10
            df = df.nsmallest(top_n, sig_col).copy()
            if df.empty:
                messagebox.showwarning('No Rows', 'Top-X p-value filter removed all rows.')
                return

        def classify(row):
            if row[sig_col] <= alpha:
                return 'Upregulated' if row[est_col] > 0 else 'Downregulated'
            return 'Not significant'

        df['__status__'] = df.apply(classify, axis=1)
        color_map = {
            'Upregulated': '#d73027',
            'Downregulated': '#4575b4',
            'Not significant': '#8c8c8c',
        }

        # Sort by selected method
        sort_method = self.effect_sort_var.get()
        if sort_method == 'Name (Alphabetical)':
            df = df.sort_values(group_col, ascending=True).reset_index(drop=True)
        elif sort_method == 'P-value':
            df = df.sort_values(sig_col, ascending=True).reset_index(drop=True)
        else:  # Effect Size (default)
            df = df.sort_values(est_col, ascending=True).reset_index(drop=True)
        
        x_positions = np.arange(len(df))

        fig, ax = plt.subplots(figsize=(self.effect_fig_w_var.get(), self.effect_fig_h_var.get()))
        
        # Add background vertical lines if enabled
        if self.effect_show_vlines_var.get():
            for i in x_positions:
                ax.axvline(x=i, color='#eeeeee', linestyle='-', linewidth=0.8, alpha=0.6, zorder=0)
        
        for i, row in df.iterrows():
            color = color_map.get(row['__status__'], '#8c8c8c')
            y_value = float(row[est_col])
            y_err_lower = max(0.0, y_value - float(row[ci_low_col]))
            y_err_upper = max(0.0, float(row[ci_high_col]) - y_value)

            ax.errorbar(
                i,
                y_value,
                yerr=[[y_err_lower], [y_err_upper]],
                fmt='o',
                color=color,
                ecolor=color,
                elinewidth=2.0,
                capsize=4,
                capthick=1.5,
                markersize=7,
                markeredgecolor=color,
                zorder=3,
            )

            if self.effect_show_pvalues_var.get():
                def _format_p_label(p_val):
                    if pd.isna(p_val):
                        return 'NA'
                    p_val = float(p_val)
                    if p_val < 1e-4:
                        return f"{p_val:.1e}"
                    # Use one significant figure with plain decimal text for p >= 1e-4
                    from decimal import Decimal, ROUND_HALF_UP
                    d = Decimal(str(p_val))
                    rounded = d.quantize(Decimal(f"1e{d.adjusted()}"), rounding=ROUND_HALF_UP)
                    txt = format(rounded, 'f')
                    if '.' in txt:
                        txt = txt.rstrip('0').rstrip('.')
                    return txt

                p_text = _format_p_label(row[sig_col])
                # Position p-value at top of error bar
                p_y_pos = y_value + y_err_upper
                ax.annotate(
                    p_text,
                    (i, p_y_pos),
                    textcoords='offset points',
                    xytext=(0, 3),
                    ha='center',
                    fontsize=8,
                    color=color,
                    zorder=4,
                )

        ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, zorder=1)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(df[group_col].tolist(), fontsize=self.effect_xtick_fontsize_var.get(), fontweight='bold', rotation=90)
        ax.set_xlabel(self.effect_xlabel_var.get().strip() or 'Lipid Class', fontsize=self.effect_xlabel_fontsize_var.get(), fontweight='bold')
        ax.set_ylabel(self.effect_ylabel_var.get().strip() or 'Effect Size', fontsize=self.effect_ylabel_fontsize_var.get(), fontweight='bold')
        ax.tick_params(axis='x', labelsize=self.effect_xtick_fontsize_var.get(), width=self.effect_xtick_width_var.get())
        ax.tick_params(axis='y', labelsize=self.effect_ytick_fontsize_var.get(), width=self.effect_ytick_width_var.get())
        ax.grid(axis='y', alpha=0.25, linestyle='--', zorder=0)

        from matplotlib.lines import Line2D

        legend_handles = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map['Not significant'], markersize=8, label='Not significant'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map['Upregulated'], markersize=8, label='Upregulated'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map['Downregulated'], markersize=8, label='Downregulated'),
        ]
        ax.legend(handles=legend_handles, loc='best', frameon=True, fontsize=9)

        y_min = np.nanmin(df[ci_low_col].to_numpy(dtype=float))
        y_max = np.nanmax(df[ci_high_col].to_numpy(dtype=float))
        if np.isfinite(y_min) and np.isfinite(y_max):
            y_pad = (y_max - y_min) * 0.08 if y_max != y_min else 1.0
            ax.set_ylim(y_min - y_pad, y_max + y_pad)

        plt.tight_layout()
        save_path = filedialog.asksaveasfilename(
            title='Save Effect Plot',
            defaultextension='.png',
            filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All Files', '*.*')],
        )
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            self.effect_status_label.config(text=f'Saved to {os.path.basename(save_path)}', fg='green')
        else:
            self.effect_status_label.config(text='Plot generated (not saved).', fg='green')
        plt.close(fig)

    # ==================== GLYCAN CLASSIFICATION TAB ====================
    def setup_glycan_classification_tab(self):
        """Setup Glycan Classification utility tab."""
        glycan_frame = ttk.Frame(self.notebook)
        self.notebook.add(glycan_frame, text='Glycan Classification')

        # Create scrollable container
        scroll_host = tk.Frame(glycan_frame, bg='#f0f0f0')
        scroll_host.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(scroll_host, bg='#f0f0f0', highlightthickness=0)
        v_scroll = ttk.Scrollbar(scroll_host, orient='vertical', command=canvas.yview)
        container = tk.Frame(canvas, bg='#f0f0f0')
        container.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas_window = canvas.create_window((0, 0), window=container, anchor='nw')
        canvas.configure(yscrollcommand=v_scroll.set)

        def _fit_glycan_container_width(event):
            canvas.itemconfigure(canvas_window, width=event.width)

        canvas.bind('<Configure>', _fit_glycan_container_width)
        canvas.pack(side='left', fill='both', expand=True)
        v_scroll.pack(side='right', fill='y')

        tk.Label(container, text='Glycan Classification Tool', font=('Arial', 14, 'bold'), bg='#f0f0f0').pack(pady=8)

        # Step 1: Upload
        upload_frame = tk.LabelFrame(container, text='Step 1: Upload Excel File', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        upload_frame.pack(fill='x', pady=5, padx=5)
        row = tk.Frame(upload_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=6, pady=6)
        tk.Button(
            row,
            text='Upload Excel',
            command=self.upload_glycan_data,
            bg='#4CAF50',
            fg='white',
            font=('Arial', 10, 'bold'),
        ).pack(side='left', padx=5)
        self.glycan_file_label = tk.Label(row, text='No file uploaded', bg='#f0f0f0', font=('Arial', 9))
        self.glycan_file_label.pack(side='left', padx=6)

        # Step 2: Column Selection
        col_frame = tk.LabelFrame(container, text='Step 2: Select Columns', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        col_frame.pack(fill='x', pady=5, padx=5)

        # Feature ID column
        row1 = tk.Frame(col_frame, bg='#f0f0f0')
        row1.pack(fill='x', padx=6, pady=5)
        tk.Label(row1, text='Glycan Feature ID Column:', bg='#f0f0f0', width=25, anchor='w').pack(side='left')
        self.glycan_feature_id_var = tk.StringVar()
        self.glycan_feature_combo = ttk.Combobox(row1, textvariable=self.glycan_feature_id_var, state='readonly', width=40)
        self.glycan_feature_combo.pack(side='left', padx=5)

        # Sample columns
        row2 = tk.Frame(col_frame, bg='#f0f0f0')
        row2.pack(fill='x', padx=6, pady=5)
        tk.Button(
            row2,
            text='Select Sample Columns',
            command=self.glycan_select_sample_columns,
            bg='#3498db',
            fg='white',
            font=('Arial', 9, 'bold'),
        ).pack(side='left', padx=5)
        self.glycan_sample_label = tk.Label(row2, text='No columns selected', bg='#f0f0f0', font=('Arial', 9), fg='#555')
        self.glycan_sample_label.pack(side='left', padx=6)

        # Step 3: Classification Summary
        summary_frame = tk.LabelFrame(container, text='Step 3: Preview Classification', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        summary_frame.pack(fill='both', expand=True, pady=5, padx=5)

        summary_label = tk.Label(summary_frame, text='Classification Summary:', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        summary_label.pack(anchor='w', padx=6, pady=(4, 2))

        list_frame = tk.Frame(summary_frame, bg='#f0f0f0')
        list_frame.pack(fill='both', expand=True, padx=6, pady=(2, 6))

        self.glycan_summary_text = tk.Text(list_frame, height=10, wrap='word', font=('Consolas', 9), bg='white')
        self.glycan_summary_text.pack(side='left', fill='both', expand=True)

        summary_scroll = ttk.Scrollbar(list_frame, orient='vertical', command=self.glycan_summary_text.yview)
        summary_scroll.pack(side='right', fill='y')
        self.glycan_summary_text.configure(yscrollcommand=summary_scroll.set)

        # Step 4: Export
        export_frame = tk.LabelFrame(container, text='Step 4: Export Results', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        export_frame.pack(fill='x', pady=5, padx=5)

        row = tk.Frame(export_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=6, pady=6)
        tk.Label(row, text='Output Folder:', bg='#f0f0f0', width=25, anchor='w').pack(side='left')
        self.glycan_output_dir_var = tk.StringVar(value='')
        tk.Entry(row, textvariable=self.glycan_output_dir_var, width=50).pack(side='left', padx=5)
        tk.Button(
            row,
            text='Browse',
            command=self._browse_glycan_output_dir,
            bg='#7f8c8d',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=4)

        # Generate button
        tk.Button(
            container,
            text='Classify Glycans and Export',
            command=self.process_and_export_glycan_data,
            bg='#2196F3',
            fg='white',
            font=('Arial', 12, 'bold'),
            height=2,
        ).pack(fill='x', pady=10, padx=5)

        # Status
        self.glycan_status_label = tk.Label(container, text='', bg='#f0f0f0', font=('Arial', 9), fg='blue')
        self.glycan_status_label.pack(pady=5)

        # Processing Log
        log_frame = tk.LabelFrame(container, text='Processing Log', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        log_frame.pack(fill='both', expand=True, pady=5, padx=5)

        log_inner = tk.Frame(log_frame, bg='#f0f0f0')
        log_inner.pack(fill='both', expand=True, padx=6, pady=6)

        self.glycan_log_text = tk.Text(log_inner, height=14, wrap='word', font=('Consolas', 8), bg='#1e1e1e', fg='#00ff00')
        self.glycan_log_text.pack(side='left', fill='both', expand=True)

        log_scroll = ttk.Scrollbar(log_inner, orient='vertical', command=self.glycan_log_text.yview)
        log_scroll.pack(side='right', fill='y')
        self.glycan_log_text.configure(yscrollcommand=log_scroll.set)

        # Initialize state
        self.glycan_data = None
        self.glycan_sample_columns = []

    def upload_glycan_data(self):
        """Upload Excel file containing glycan data."""
        file_path = filedialog.askopenfilename(
            title='Select Excel File',
            filetypes=[('Excel Files', '*.xlsx *.xls'), ('All Files', '*.*')]
        )

        if not file_path:
            return

        try:
            df = pd.read_excel(file_path)
            self.glycan_data = df

            filename = os.path.basename(file_path)
            self.glycan_file_label.config(text=f'✓ {filename} ({len(df)} rows)', fg='green')

            # Populate feature ID column dropdown
            columns = [''] + list(df.columns)
            self.glycan_feature_combo['values'] = columns

            # Auto-detect feature ID column
            cols_lower = {c.lower(): c for c in df.columns}
            for candidate in ['feature_id', 'glycan_id', 'id', 'glycan', 'feature']:
                if candidate in cols_lower:
                    self.glycan_feature_id_var.set(cols_lower[candidate])
                    break

            self.glycan_status_label.config(text='✓ File uploaded successfully!', fg='green')
            self.glycan_summary_text.delete('1.0', tk.END)
            self.glycan_sample_columns = []
            self.glycan_sample_label.config(text='No columns selected', fg='#555')

        except Exception as e:
            messagebox.showerror('Error', f'Failed to load file: {str(e)}')
            logger.error(f'Error loading glycan data: {e}')
            self.glycan_status_label.config(text='✗ Upload failed', fg='red')

    def glycan_select_sample_columns(self):
        """Select sample columns using statistics verify columns dialog."""
        if self.glycan_data is None:
            messagebox.showwarning('No Data', 'Please upload a file first.')
            return

        selected_columns = self._verify_sample_columns_with_statistics_dialog(
            self.glycan_data,
            'Glycan Classification - Verify Sample Columns'
        )
        if selected_columns is None:
            return

        self.glycan_sample_columns = selected_columns
        if self.glycan_sample_columns:
            self.glycan_sample_label.config(text=f'✓ {len(self.glycan_sample_columns)} columns selected', fg='green')
        else:
            self.glycan_sample_label.config(text='No columns selected', fg='red')

    def _browse_glycan_output_dir(self):
        """Browse for output directory."""
        folder = filedialog.askdirectory(title='Select Output Folder')
        if folder:
            self.glycan_output_dir_var.set(folder)
    def _log_glycan_message(self, message, tag='INFO'):
        """Add a message to the glycan processing log."""
        if not hasattr(self, 'glycan_log_text'):
            return
        timestamp = datetime.now().strftime('%H:%M:%S')
        log_entry = f"[{timestamp}] {tag}: {message}\n"
        self.glycan_log_text.insert(tk.END, log_entry)
        self.glycan_log_text.see(tk.END)
        self.frame.update()


    def process_and_export_glycan_data(self):
        """Process glycan data and export to Excel with detailed logging."""
        # Clear log
        self.glycan_log_text.delete('1.0', tk.END)
        self._log_glycan_message('Starting glycan classification process...', 'START')

        # Validate inputs
        if self.glycan_data is None:
            self._log_glycan_message('No data uploaded - aborting', 'ERROR')
            messagebox.showwarning('No Data', 'Please upload an Excel file first.')
            return

        feature_id_col = self.glycan_feature_id_var.get().strip()
        if not feature_id_col:
            self._log_glycan_message('No feature ID column selected - aborting', 'ERROR')
            messagebox.showwarning('No Column', 'Please select a Feature ID column.')
            return

        if not self.glycan_sample_columns:
            self._log_glycan_message('No sample columns selected - aborting', 'ERROR')
            messagebox.showwarning('No Columns', 'Please select sample columns first.')
            return

        output_dir = self.glycan_output_dir_var.get().strip()
        if not output_dir or not os.path.isdir(output_dir):
            self._log_glycan_message('Invalid output folder - aborting', 'ERROR')
            messagebox.showwarning('Invalid Output', 'Please select a valid output folder.')
            return

        try:
            self.glycan_status_label.config(text='Processing glycan data...', fg='blue')
            self.frame.update()

            self._log_glycan_message(f'Data loaded: {len(self.glycan_data)} rows, {len(self.glycan_data.columns)} columns', 'INFO')
            self._log_glycan_message(f'Feature ID column: "{feature_id_col}"', 'INFO')
            self._log_glycan_message(f'Sample columns: {len(self.glycan_sample_columns)} selected', 'INFO')

            # Process data
            from main_script.glycan_classification import classify_glycan
            
            work_df = self.glycan_data.copy()
            self._log_glycan_message('Applying glycan classification to all rows...', 'INFO')
            
            work_df['Glycan_Class'] = work_df[feature_id_col].astype(str).apply(classify_glycan)
            
            # Show debug info about classifications
            classification_counts = work_df['Glycan_Class'].value_counts(dropna=False)
            self._log_glycan_message('Classification breakdown:', 'INFO')
            
            total_classified = 0
            total_unclassified = 0
            for cls, count in classification_counts.items():
                if cls is None or (isinstance(cls, float) and pd.isna(cls)):
                    self._log_glycan_message(f'  Unclassified (None): {count} rows', 'WARN')
                    total_unclassified = count
                else:
                    self._log_glycan_message(f'  {cls}: {count} rows', 'INFO')
                    total_classified += count
            
            self._log_glycan_message(f'Total classified: {total_classified}, Unclassified: {total_unclassified}', 'INFO')
            
            # Filter to classified only
            result_df = work_df[work_df['Glycan_Class'].notna()].copy()
            
            if result_df.empty:
                self._log_glycan_message('ERROR: No glycans could be classified!', 'ERROR')
                summary_text = "No glycans could be classified.\n\nFirst 10 Feature IDs from data:\n"
                for i, val in enumerate(self.glycan_data[feature_id_col].head(10)):
                    summary_text += f"  {i+1}. {val}\n"
                self.glycan_summary_text.delete('1.0', tk.END)
                self.glycan_summary_text.insert(tk.END, summary_text)
                self.glycan_status_label.config(text='✗ No glycans could be classified - check feature ID format', fg='red')
                self._log_glycan_message('Aborting export - no classified data', 'ERROR')
                return
            
            self._log_glycan_message(f'Grouping {len(result_df)} classified rows by glycan class...', 'INFO')
            
            # Group by classification and sum sample columns
            grouped = result_df.groupby('Glycan_Class', as_index=False)[self.glycan_sample_columns].sum()
            
            self._log_glycan_message(f'Grouped into {len(grouped)} glycan classes', 'INFO')

            # Build QC tables so the grouping decisions are traceable
            qc_grouping_df = work_df.copy()
            qc_grouping_df.insert(0, 'Source_Row', qc_grouping_df.index + 2)
            qc_grouping_df.insert(1, 'Feature_ID', qc_grouping_df[feature_id_col])
            qc_grouping_df['Grouped_Into'] = qc_grouping_df['Glycan_Class'].fillna('Unclassified')
            qc_grouping_df['QC_Status'] = qc_grouping_df['Glycan_Class'].apply(
                lambda value: 'Grouped' if pd.notna(value) else 'Unclassified'
            )

            qc_summary_df = (
                qc_grouping_df.groupby('Grouped_Into', dropna=False)
                .agg(
                    Row_Count=('Source_Row', 'count'),
                    Source_Rows=('Source_Row', lambda values: ', '.join(map(str, values))),
                    Feature_IDs=('Feature_ID', lambda values: ', '.join(map(str, values[:10]))),
                )
                .reset_index()
                .rename(columns={'Grouped_Into': 'Glycan_Class'})
            )

            # Display summary with both counts and abundance to avoid confusion
            class_counts = result_df['Glycan_Class'].value_counts().to_dict()
            summary_text = "Classification Summary:\n" + "="*86 + "\n"
            summary_text += f"{'Class':30s} {'Rows':>10s} {'Total Abundance':>20s}\n"
            summary_text += "-"*86 + "\n"
            for _, row in grouped.iterrows():
                glycan_class = row['Glycan_Class']
                row_count = int(class_counts.get(glycan_class, 0))
                total_abundance = row[self.glycan_sample_columns].sum()
                summary_text += f"{glycan_class:30s} {row_count:10,d} {total_abundance:20,.0f}\n"

            summary_text += "-"*86 + "\n"
            summary_text += f"{'TOTAL':30s} {len(result_df):10,d}"

            self.glycan_summary_text.delete('1.0', tk.END)
            self.glycan_summary_text.insert(tk.END, summary_text)

            # Export to Excel
            self._log_glycan_message('Generating Excel file...', 'INFO')
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_file = os.path.join(output_dir, f'glycan_classification_{timestamp}.xlsx')

            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                grouped.to_excel(writer, sheet_name='Classification', index=False)
                qc_grouping_df.to_excel(writer, sheet_name='QC_Grouping', index=False)
                qc_summary_df.to_excel(writer, sheet_name='QC_Summary', index=False)

            self.glycan_status_label.config(
                text=f'✓ Exported to {os.path.basename(output_file)}',
                fg='green'
            )
            self._log_glycan_message(f'Successfully exported to: {os.path.basename(output_file)}', 'SUCCESS')
            self._log_glycan_message(f'Process complete! {len(grouped)} glycan classes exported.', 'SUCCESS')
            messagebox.showinfo('Success', f'Glycan data classified and exported to:\n{output_file}\n\n{len(grouped)} classifications found')

        except Exception as e:
            self._log_glycan_message(f'FATAL ERROR: {str(e)}', 'ERROR')
            messagebox.showerror('Error', f'Failed to process glycan data: {str(e)}')
            logger.error(f'Error processing glycan data: {e}', exc_info=True)
            self.glycan_status_label.config(text='✗ Processing failed', fg='red')