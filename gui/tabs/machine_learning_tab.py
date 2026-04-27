"""
Machine Learning Tab - ML-based analysis for metabolomics/lipid data

This tab provides machine learning capabilities following the same data loading workflow as Statistics tab:
1. Select data mode (Metabolite/Lipid/Custom)
2. Import Excel (Pos/Neg sheets)
3. Verify columns
4. Configure groups
5. Run ML analysis with automatic filtering and merging

Author: Metabolomics Statistics Tool
Date: 2025
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import pandas as pd
import numpy as np
import threading
import json
import os
import logging
import re
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple
from collections import Counter
from itertools import combinations

# Import shared components
from gui.shared.base_tab import BaseTab
from gui.shared.utils import resolve_runtime_config_path, is_statistics_metadata_col
from gui.shared.column_assignment import show_column_assignment_dialog

logger = logging.getLogger(__name__)


class MachineLearningTab(BaseTab):
    """
    Machine Learning Tab for classification and dimensionality reduction.
    
    Follows the exact same workflow as Statistics tab:
    - Step 1: Import Excel and verify columns
    - Step 2: Configure sample groups
    - Step 3: Configure ML settings
    - Step 4: Run ML analysis (automatic filtering + merging + ML)
    """
    
    def __init__(self, parent, data_manager):
        """Initialize Machine Learning tab"""
        super().__init__(parent, data_manager)
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Initialize data storage (same pattern as Statistics tab)
        self.ml_data_mode = tk.StringVar(value='custom')
        self.pos_df = None
        self.neg_df = None
        self.detected_pos_sample_cols = []
        self.detected_neg_sample_cols = []
        self.detected_sample_cols = []
        self.verified_pos_assignments = {}
        self.verified_neg_assignments = {}
        self.verified_pos_sample_cols = []
        self.verified_neg_sample_cols = []
        self.verified_pos_feature_id_col = None  # Store feature ID column from Pos verification
        self.verified_neg_feature_id_col = None  # Store feature ID column from Neg verification
        self.sample_group_vars = {}
        self.group_definitions = {'Group1': 'Control', 'Group2': 'Disease', 'Group3': 'Treatment', 'Group4': 'Other'}
        self.group_count = 4
        self.auto_assign_patterns = {}  # Store patterns for each group
        self.pairwise_pvalue_map = {}  # Manual mapping: pair key -> selected p-value column
        self.working_folder_var = tk.StringVar(value='')
        self.figure_settings = self._default_figure_settings()
        self.top_n_values = self._default_top_n_values()
        
        # ML results storage
        self.ml_results = None
        self.pairwise_ml_results = None
        self.merged_data = None
        
        # Filtering configuration (same as Statistics tab)
        self.min_samples_per_group_var = tk.StringVar(value='2')
        self.min_samples_type_var = tk.StringVar(value='absolute')
        self.min_samples_percent_var = tk.StringVar(value='50.0')
        
        # Setup UI
        self.setup_ui()
        
        # Setup save/load config functions
        self._setup_ml_config_functions()
        
        # Load saved config on startup
        self._load_ml_config()
        
        logger.info("Machine Learning tab initialized")

    def _default_figure_settings(self) -> Dict[str, float]:
        """Default publication figure settings (user-adjustable)."""
        return {
            # ROC figure
            'roc_width': 8.2,
            'roc_height': 6.8,
            'roc_title_fs': 16.0,
            'roc_label_fs': 14.0,
            'roc_tick_fs': 12.0,
            'roc_legend_fs': 12.0,
            'roc_line_w': 3.0,
            'roc_axis_w': 1.6,

            # Model comparison bar figure
            'comparison_height': 6.0,
            'comparison_base_width': 3.0,
            'comparison_width_per_bar': 1.65,
            'comparison_bar_w': 0.62,
            'comparison_gap_w': 0.34,
            'comparison_error_w': 2.4,
            'comparison_title_fs': 16.0,
            'comparison_label_fs': 14.0,
            'comparison_tick_fs': 12.0,
            'comparison_value_fs': 12.0,
            'comparison_axis_w': 1.6,

            # Horizontal top-metabolite bar figures
            'hbar_width': 10.5,
            'hbar_base_height': 2.4,
            'hbar_height_per_feature': 0.52,
            'hbar_bar_h': 0.55,
            'hbar_gap_h': 0.28,
            'hbar_error_w': 2.2,
            'hbar_title_fs': 16.0,
            'hbar_label_fs': 14.0,
            'hbar_tick_fs': 12.0,
            'hbar_axis_w': 1.6,
        }

    def _default_top_n_values(self) -> Tuple[int, int]:
        """Default top-N ranks used for auto-generated metabolite and Venn figures."""
        return (10, 15)

    def _normalize_top_n_values(self, values: Any) -> Tuple[int, int]:
        """Normalize a pair of top-N values into a sorted unique tuple."""
        default_values = self._default_top_n_values()
        if not isinstance(values, (list, tuple)) or len(values) != 2:
            return default_values

        cleaned = []
        for value in values:
            try:
                cleaned.append(int(value))
            except (TypeError, ValueError):
                return default_values

        if any(v <= 0 for v in cleaned):
            return default_values

        unique_sorted = tuple(sorted(set(cleaned)))
        if len(unique_sorted) != 2:
            return default_values
        return unique_sorted[0], unique_sorted[1]
    
    def setup_ui(self):
        """Create the Machine Learning tab interface"""
        # Main container
        main_frame = tk.Frame(self.frame, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Title
        tk.Label(
            main_frame,
            text='🤖 Machine Learning Analysis',
            font=('Arial', 14, 'bold'),
            bg='#f0f0f0'
        ).pack(pady=(0, 10))
        
        # Create 3-column layout
        body = tk.Frame(main_frame, bg='#f0f0f0')
        body.pack(fill='both', expand=True)
        body.grid_columnconfigure(0, weight=1)  # Left: Steps 1-2
        body.grid_columnconfigure(1, weight=1)  # Middle: Step 3
        body.grid_columnconfigure(2, weight=2)  # Right: Log + Step 4
        body.grid_rowconfigure(0, weight=1)
        
        # ===== LEFT COLUMN: Steps 1-2 =====
        left_frame = tk.LabelFrame(body, text='📋 Steps 1-2: Data & Groups', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 3))
        
        # Add scrollbars to left panel
        left_canvas = tk.Canvas(left_frame, bg='#f0f0f0', highlightthickness=0)
        left_scrollbar_y = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        left_scrollbar_x = ttk.Scrollbar(left_frame, orient="horizontal", command=left_canvas.xview)
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
        
        # Mouse wheel scrolling
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_scrollable.bind("<MouseWheel>", _on_left_mousewheel)
        
        # Create Steps 1-2 in first column
        self._create_step2_import_verify(left_scrollable)
        self._create_step3_group_config(left_scrollable)
        
        # ===== MIDDLE COLUMN: Step 3 =====
        middle_frame = tk.LabelFrame(body, text='⚙️ Step 3: ML Config', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        middle_frame.grid(row=0, column=1, sticky='nsew', padx=3)
        
        # Add scrollbars to middle panel
        middle_canvas = tk.Canvas(middle_frame, bg='#f0f0f0', highlightthickness=0)
        middle_scrollbar_y = ttk.Scrollbar(middle_frame, orient="vertical", command=middle_canvas.yview)
        middle_scrollbar_x = ttk.Scrollbar(middle_frame, orient="horizontal", command=middle_canvas.xview)
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
        
        # Mouse wheel scrolling
        def _on_middle_mousewheel(event):
            middle_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        middle_canvas.bind("<MouseWheel>", _on_middle_mousewheel)
        middle_scrollable.bind("<MouseWheel>", _on_middle_mousewheel)
        
        # Create Step 3 in middle column
        self._create_step4_ml_config(middle_scrollable)
        
        # ===== RIGHT COLUMN: Step 4 + Log/Results =====
        right_frame = tk.LabelFrame(body, text='📈 Step 4: Results', bg='#f0f0f0', font=('Arial', 11, 'bold'))
        right_frame.grid(row=0, column=2, sticky='nsew', padx=(3, 0))
        
        # Add scrollbars to right panel
        right_canvas = tk.Canvas(right_frame, bg='#f0f0f0', highlightthickness=0)
        right_scrollbar_y = ttk.Scrollbar(right_frame, orient="vertical", command=right_canvas.yview)
        right_scrollbar_x = ttk.Scrollbar(right_frame, orient="horizontal", command=right_canvas.xview)
        right_scrollable = tk.Frame(right_canvas, bg='#f0f0f0')
        
        right_scrollable.bind(
            "<Configure>",
            lambda e: right_canvas.configure(scrollregion=right_canvas.bbox("all"))
        )
        
        right_canvas_window = right_canvas.create_window((0, 0), window=right_scrollable, anchor="nw")
        right_canvas.configure(yscrollcommand=right_scrollbar_y.set, xscrollcommand=right_scrollbar_x.set)
        
        def configure_right_scroll(event):
            right_canvas.configure(scrollregion=right_canvas.bbox("all"))
            right_canvas.itemconfig(right_canvas_window, width=event.width)
        
        right_canvas.bind('<Configure>', configure_right_scroll)
        
        right_scrollbar_y.pack(side="right", fill="y", padx=(2, 0))
        right_scrollbar_x.pack(side="bottom", fill="x", pady=(2, 0))
        right_canvas.pack(side="left", fill="both", expand=True)
        
        # Mouse wheel scrolling
        def _on_right_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        right_canvas.bind("<MouseWheel>", _on_right_mousewheel)
        right_scrollable.bind("<MouseWheel>", _on_right_mousewheel)
        
        # Action buttons at top
        action_container = tk.Frame(right_scrollable, bg='#f0f0f0')
        action_container.pack(fill='x', padx=5, pady=5)
        
        self._create_step5_actions(action_container)
        
        # Results text area with scrollbar
        self.ml_log = scrolledtext.ScrolledText(
            right_scrollable,
            wrap=tk.WORD,
            font=("Consolas", 9),
            bg="white",
            fg="black",
            height=40
        )
        self.ml_log.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Initial welcome message
        self.ml_log.insert(tk.END, "="*60 + "\n")
        self.ml_log.insert(tk.END, "🤖 Machine Learning Analysis Tool\n")
        self.ml_log.insert(tk.END, "="*60 + "\n\n")
        self.ml_log.insert(tk.END, "Workflow:\n")
        self.ml_log.insert(tk.END, "1. Import Excel file\n")
        self.ml_log.insert(tk.END, "2. Verify columns\n")
        self.ml_log.insert(tk.END, "3. Configure sample groups\n")
        self.ml_log.insert(tk.END, "4. Run classification analysis\n\n")
        self.ml_log.insert(tk.END, "Ready to begin!\n\n")
        
        logger.info("Machine Learning tab UI created")
    
    def _create_step1_mode_selection(self, parent):
        """Step 1: Data mode (fixed to Custom)."""
        frame = tk.LabelFrame(parent, text='Step 1: Select Data Mode', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        frame.pack(fill='x', padx=5, pady=(5, 10))
        
        mode_frame = tk.Frame(frame, bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        tk.Label(mode_frame, text='Data type:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(side='left', padx=(0, 10))
        tk.Label(
            mode_frame,
            text='Custom (fixed)',
            bg='#f0f0f0',
            fg='#2c3e50',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=5)
        
        # Mode description
        self.mode_desc_label = tk.Label(
            frame,
            text='Custom mode for preprocessed/combined data',
            bg='#e3f2fd',
            fg='#1565c0',
            font=('Arial', 8, 'italic'),
            wraplength=280,
            justify='left',
            padx=5,
            pady=5
        )
        self.mode_desc_label.pack(fill='x', padx=5, pady=(0, 5))
    
    def _create_step2_import_verify(self, parent):
        """Step 1: Import and verify data"""
        frame = tk.LabelFrame(parent, text='Step 1: Import & Verify Data', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        frame.pack(fill='x', padx=5, pady=(0, 10))
        
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}

        top_btn_row = tk.Frame(frame, bg='#f0f0f0')
        top_btn_row.pack(fill='x', padx=5, pady=5)
        top_btn_row.grid_columnconfigure(0, weight=1)
        top_btn_row.grid_columnconfigure(1, weight=1)

        tk.Button(
            top_btn_row,
            text='📂 Import Excel',
            command=self._import_excel,
            bg='#8e44ad',
            fg='white',
            **btn_style
        ).grid(row=0, column=0, sticky='ew', padx=(0, 3))

        tk.Button(
            top_btn_row,
            text='📁 Working Folder',
            command=self._select_working_folder,
            bg='#16a085',
            fg='white',
            **btn_style
        ).grid(row=0, column=1, sticky='ew', padx=(3, 0))

        self.working_folder_label = tk.Label(
            frame,
            text='Working folder: Not set',
            bg='#f0f0f0',
            fg='#666666',
            font=('Arial', 8),
            anchor='w',
            justify='left'
        )
        self.working_folder_label.pack(fill='x', padx=7, pady=(0, 4))
        
        tk.Button(
            frame,
            text='🔍 Verify Columns',
            command=self._verify_columns,
            bg='#2980b9',
            fg='white',
            **btn_style
        ).pack(fill='x', padx=5, pady=5)
        
        # Status label
        self.data_status_label = tk.Label(
            frame,
            text="⚠️ No data loaded",
            fg="red",
            bg='#f0f0f0',
            font=('Arial', 9)
        )
        self.data_status_label.pack(pady=5)

    def _select_working_folder(self):
        """Select working folder used for automated figure exports."""
        selected_dir = filedialog.askdirectory(title='Select Working Folder for Figure Outputs')
        if not selected_dir:
            return

        self.working_folder_var.set(selected_dir)
        if hasattr(self, 'working_folder_label'):
            self.working_folder_label.config(text=f"Working folder: {selected_dir}", fg='#2c3e50')

        self._save_ml_config()
        self._log(f"📁 Working folder set: {selected_dir}\n")

    def _open_figure_settings_dialog(self):
        """Popup to configure figure size, fonts, and line/bar thickness."""
        dialog = tk.Toplevel(self.root)
        dialog.title('Configure Figure Settings')
        dialog.geometry('760x780')
        dialog.transient(self.root)
        dialog.grab_set()

        canvas = tk.Canvas(dialog, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(dialog, orient='vertical', command=canvas.yview)
        content = tk.Frame(canvas, bg='#f0f0f0')
        content.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=content, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        tk.Label(content, text='Manual Figure Controls', bg='#f0f0f0', fg='#2c3e50', font=('Arial', 12, 'bold')).pack(anchor='w', padx=10, pady=(10, 4))
        tk.Label(
            content,
            text='These settings are applied to ROC, model comparison bars, top-metabolite plots, and Venn diagrams.',
            bg='#f0f0f0',
            fg='#555555',
            font=('Arial', 9)
        ).pack(anchor='w', padx=10, pady=(0, 8))

        field_defs = {
            'ROC': [
                ('roc_width', 'Width'), ('roc_height', 'Height'),
                ('roc_title_fs', 'Title Font'), ('roc_label_fs', 'Axis Label Font'),
                ('roc_tick_fs', 'Tick Font'), ('roc_legend_fs', 'Legend Font'),
                ('roc_line_w', 'Curve Line Width'), ('roc_axis_w', 'Axis Line Width')
            ],
            'Model Comparison Bar': [
                ('comparison_height', 'Height'), ('comparison_base_width', 'Base Width'),
                ('comparison_width_per_bar', 'Width per Bar'), ('comparison_bar_w', 'Bar Width'),
                ('comparison_gap_w', 'Gap Between Bars'), ('comparison_error_w', 'Error Bar Width'),
                ('comparison_title_fs', 'Title Font'), ('comparison_label_fs', 'Axis Label Font'),
                ('comparison_tick_fs', 'Tick Font'), ('comparison_value_fs', 'Value Label Font'),
                ('comparison_axis_w', 'Axis Line Width')
            ],
            'Top Metabolites Horizontal Bar': [
                ('hbar_width', 'Width'), ('hbar_base_height', 'Base Height'),
                ('hbar_height_per_feature', 'Height per Feature'), ('hbar_bar_h', 'Bar Height'),
                ('hbar_gap_h', 'Gap Between Bars'), ('hbar_error_w', 'Error Bar Width'),
                ('hbar_title_fs', 'Title Font'), ('hbar_label_fs', 'Axis Label Font'),
                ('hbar_tick_fs', 'Tick Font'), ('hbar_axis_w', 'Axis Line Width')
            ]
        }

        setting_vars = {}
        defaults = self._default_figure_settings()

        for section, fields in field_defs.items():
            section_frame = tk.LabelFrame(content, text=section, bg='#f0f0f0', font=('Arial', 10, 'bold'))
            section_frame.pack(fill='x', padx=10, pady=6)

            for key, label in fields:
                row = tk.Frame(section_frame, bg='#f0f0f0')
                row.pack(fill='x', padx=8, pady=2)
                tk.Label(row, text=f'{label}:', width=22, anchor='w', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
                var = tk.StringVar(value=str(self.figure_settings.get(key, defaults.get(key, 1.0))))
                setting_vars[key] = var
                ttk.Entry(row, textvariable=var, width=12).pack(side='left')

        top_n_frame = tk.LabelFrame(content, text='Top-N Figure Selection', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        top_n_frame.pack(fill='x', padx=10, pady=6)
        tk.Label(
            top_n_frame,
            text='Controls the top-metabolite bar charts and the Venn output for 2-3 model comparisons.',
            bg='#f0f0f0',
            fg='#555555',
            font=('Arial', 9)
        ).pack(anchor='w', padx=8, pady=(6, 4))

        top_n_defaults = self._default_top_n_values()
        top_n_vars = []
        for idx, label in enumerate(['First top-N', 'Second top-N']):
            row = tk.Frame(top_n_frame, bg='#f0f0f0')
            row.pack(fill='x', padx=8, pady=2)
            tk.Label(row, text=f'{label}:', width=22, anchor='w', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
            var = tk.StringVar(value=str(self.top_n_values[idx] if idx < len(self.top_n_values) else top_n_defaults[idx]))
            top_n_vars.append(var)
            ttk.Spinbox(row, from_=1, to=500, increment=1, textvariable=var, width=12).pack(side='left')

        tk.Label(
            top_n_frame,
            text='Example pairs: 10 and 15, 15 and 20, 20 and 25.',
            bg='#f0f0f0',
            fg='#666666',
            font=('Arial', 8)
        ).pack(anchor='w', padx=8, pady=(0, 6))

        btn_frame = tk.Frame(content, bg='#f0f0f0')
        btn_frame.pack(fill='x', padx=10, pady=(10, 12))

        def _apply_settings(close_after=False):
            updated = self.figure_settings.copy()
            for key, var in setting_vars.items():
                raw = var.get().strip()
                try:
                    val = float(raw)
                except ValueError:
                    messagebox.showerror('Invalid Value', f'Invalid numeric value for {key}: {raw}')
                    return

                if val <= 0:
                    messagebox.showerror('Invalid Value', f'{key} must be > 0')
                    return
                updated[key] = val

            try:
                top_n_raw = [int(var.get().strip()) for var in top_n_vars]
            except ValueError:
                messagebox.showerror('Invalid Value', 'Top-N values must be whole numbers.')
                return

            if any(v <= 0 for v in top_n_raw):
                messagebox.showerror('Invalid Value', 'Top-N values must be greater than 0.')
                return

            top_n_values = tuple(sorted(set(top_n_raw)))
            if len(top_n_values) != 2:
                messagebox.showerror('Invalid Value', 'Top-N values must be two different numbers.')
                return

            self.figure_settings = updated
            self.top_n_values = top_n_values
            self._save_ml_config()
            self._log('Figure settings updated. Future auto-generated figures will use the new values.\n')
            if close_after:
                dialog.destroy()

        def _reset_defaults():
            for key, var in setting_vars.items():
                var.set(str(defaults.get(key, 1.0)))
            for idx, var in enumerate(top_n_vars):
                var.set(str(top_n_defaults[idx]))

        tk.Button(btn_frame, text='Apply', command=lambda: _apply_settings(False), bg='#27ae60', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=(0, 6))
        tk.Button(btn_frame, text='Apply & Close', command=lambda: _apply_settings(True), bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=6)
        tk.Button(btn_frame, text='Reset Defaults', command=_reset_defaults, bg='#f39c12', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=6)
        tk.Button(btn_frame, text='Cancel', command=dialog.destroy, bg='#95a5a6', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=6)
    
    def _create_step3_group_config(self, parent):
        """Step 2: Configure groups"""
        frame = tk.LabelFrame(parent, text='Step 2: Configure Groups', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        frame.pack(fill='x', padx=5, pady=(0, 10))
        
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        
        # Configure Groups button (disabled until verification)
        self.configure_groups_btn = tk.Button(
            frame,
            text='⚙️ Config Groups',
            command=self._configure_groups,
            bg='#9b59b6',
            fg='white',
            state='disabled',
            **btn_style
        )
        self.configure_groups_btn.pack(fill='x', padx=5, pady=5)
        
        # Group IDs & Labels section (same as Statistics tab)
        group_ids_frame = tk.LabelFrame(frame, text='Group IDs & Labels', bg='#f0f0f0')
        group_ids_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        # Buttons for add/remove groups
        buttons_frame = tk.Frame(group_ids_frame, bg='#f0f0f0')
        buttons_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        tk.Button(
            buttons_frame,
            text='+ Add Group',
            command=self.add_group,
            bg='#27ae60',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            bd=2,
            pady=2
        ).pack(side='left', padx=(0, 5))
        
        tk.Button(
            buttons_frame,
            text='- Remove',
            command=self.remove_group,
            bg='#e74c3c',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            bd=2,
            pady=2
        ).pack(side='left')
        
        # Scrollable frame for group list
        groups_canvas = tk.Canvas(group_ids_frame, bg='#f0f0f0', height=120, highlightthickness=0)
        groups_scrollbar = ttk.Scrollbar(group_ids_frame, orient="vertical", command=groups_canvas.yview)
        self.groups_scrollable_frame = tk.Frame(groups_canvas, bg='#f0f0f0')
        
        self.groups_scrollable_frame.bind(
            "<Configure>",
            lambda e: groups_canvas.configure(scrollregion=groups_canvas.bbox("all"))
        )
        
        groups_canvas.create_window((0, 0), window=self.groups_scrollable_frame, anchor="nw")
        groups_canvas.configure(yscrollcommand=groups_scrollbar.set)
        
        groups_canvas.pack(side='left', fill='both', expand=True, padx=5, pady=(0, 5))
        groups_scrollbar.pack(side='right', fill='y', pady=(0, 5))
        
        # Initialize group ID variables
        self.group_id_vars = {}
        
        # Populate initial groups
        self.refresh_group_ui()
        
        # ========== Feature Filters (right under groups) ==========
        feature_filter_frame = tk.LabelFrame(frame, text='Feature Filters (combinable)', bg='#f0f0f0', font=('Arial', 9, 'bold'))
        feature_filter_frame.pack(fill='x', padx=5, pady=(5, 5))
        
        tk.Label(feature_filter_frame, text='Select filters to apply before ML analysis:', bg='#f0f0f0', font=('Arial', 7, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2,2))

        # --- Replicate Filtering (first in feature filters) ---
        replicate_filter_frame = tk.LabelFrame(feature_filter_frame, text='Replicate Filtering', bg='#f0f0f0')
        replicate_filter_frame.pack(fill='x', padx=5, pady=(5, 5))

        tk.Label(replicate_filter_frame, text='Min samples/group:', bg='#f0f0f0', font=('Arial', 8)).pack(anchor='w', padx=5, pady=(5, 2))

        type_frame = tk.Frame(replicate_filter_frame, bg='#f0f0f0')
        type_frame.pack(fill='x', padx=10)

        tk.Radiobutton(
            type_frame,
            text='Absolute',
            variable=self.min_samples_type_var,
            value='absolute',
            bg='#f0f0f0',
            font=('Arial', 8)
        ).pack(side='left')

        tk.Radiobutton(
            type_frame,
            text='Percentage',
            variable=self.min_samples_type_var,
            value='percentage',
            bg='#f0f0f0',
            font=('Arial', 8)
        ).pack(side='left')

        value_frame = tk.Frame(replicate_filter_frame, bg='#f0f0f0')
        value_frame.pack(fill='x', padx=10, pady=(2, 5))

        tk.Label(value_frame, text='Count:', bg='#f0f0f0', font=('Arial', 8)).pack(side='left')
        tk.Entry(value_frame, textvariable=self.min_samples_per_group_var, width=6).pack(side='left', padx=5)

        tk.Label(value_frame, text='%:', bg='#f0f0f0', font=('Arial', 8)).pack(side='left', padx=(5, 0))
        tk.Entry(value_frame, textvariable=self.min_samples_percent_var, width=6).pack(side='left', padx=5)
        
        # --- Endogenous Filter ---
        self.filter_endogenous_var = tk.BooleanVar(value=False)
        endogenous_check = tk.Checkbutton(feature_filter_frame, text='Endogenous="Yes" only', 
                                         variable=self.filter_endogenous_var, bg='#f0f0f0', font=('Arial', 8))
        endogenous_check.pack(anchor='w', padx=5, pady=(2, 0))
        
        # --- Has HMDB ID Filter ---
        self.filter_has_hmdb_var = tk.BooleanVar(value=False)
        hmdb_check = tk.Checkbutton(feature_filter_frame, text='Has HMDB ID (non-empty)', 
                                   variable=self.filter_has_hmdb_var, bg='#f0f0f0', font=('Arial', 8))
        hmdb_check.pack(anchor='w', padx=5, pady=(2, 0))
        
        # --- P-Value Filter ---
        pvalue_filter_row = tk.Frame(feature_filter_frame, bg='#f0f0f0')
        pvalue_filter_row.pack(fill='x', padx=5, pady=(2, 5))
        
        self.filter_pvalue_var = tk.BooleanVar(value=False)
        pvalue_check = tk.Checkbutton(pvalue_filter_row, text='P-Value <', 
                                     variable=self.filter_pvalue_var, bg='#f0f0f0', font=('Arial', 8))
        pvalue_check.pack(side='left')
        
        self.filter_pvalue_threshold_var = tk.StringVar(value='0.05')
        pvalue_entry = tk.Entry(pvalue_filter_row, textvariable=self.filter_pvalue_threshold_var, width=6, font=('Arial', 8))
        pvalue_entry.pack(side='left', padx=(2, 5))
        
        tk.Label(pvalue_filter_row, text='(default: 0.05)', bg='#f0f0f0', font=('Arial', 7), fg='#7f8c8d').pack(side='left')

        pairwise_verify_row = tk.Frame(feature_filter_frame, bg='#f0f0f0')
        pairwise_verify_row.pack(fill='x', padx=5, pady=(0, 5))
        tk.Button(
            pairwise_verify_row,
            text='Verify Pairwise P-Value Columns',
            command=self._open_pairwise_pvalue_verification_dialog,
            bg='#2c3e50',
            fg='white',
            font=('Arial', 8, 'bold'),
            relief='raised',
            bd=2,
            padx=6,
            pady=2
        ).pack(anchor='w')
        
        # Info label
        tk.Label(feature_filter_frame, text='ⓘ Columns must be assigned in Verify Columns.', 
                bg='#fff3cd', fg='#856404', font=('Arial', 7), justify='left', padx=3, pady=3).pack(fill='x', padx=5, pady=(0, 5))

    
    def _create_step4_ml_config(self, parent):
        """Step 3: ML configuration - content goes directly in parent (no extra wrapper)"""
        # Analysis type (fixed)
        tk.Label(parent, text='Analysis Type:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 2))

        self.analysis_type_var = tk.StringVar(value='classification')
        tk.Label(
            parent,
            text='🎯 Classification (fixed)',
            bg='#f0f0f0',
            fg='#2c3e50',
            font=('Arial', 9, 'bold')
        ).pack(anchor='w', padx=10, pady=1)

        self.generate_pairwise_addon_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            parent,
            text='Pairwise add-on (multiclass only)',
            variable=self.generate_pairwise_addon_var,
            command=self._save_ml_config,
            bg='#f0f0f0',
            font=('Arial', 9, 'bold'),
            anchor='w',
            justify='left'
        ).pack(anchor='w', padx=10, pady=(1, 4))
        
        # Model selection (for classification/feature importance)
        self.model_frame = tk.LabelFrame(parent, text='Model (Select one or multiple)', bg='#f0f0f0')
        self.model_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Label(self.model_frame, text='Single or multiple selection (auto-detected):', 
                bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 2))
        
        self.model_checkboxes = {}
        models = ['Random Forest', 'Gradient Boosting', 'SVM (RBF)', 'Logistic Regression', 'Linear Discriminant Analysis']
        for model in models:
            var = tk.BooleanVar(value=(model == 'Random Forest'))  # Default select RF
            self.model_checkboxes[model] = var
            tk.Checkbutton(
                self.model_frame,
                text=model,
                variable=var,
                bg='#f0f0f0'
            ).pack(anchor='w', padx=15, pady=1)
        
        # Parameters
        param_frame = tk.LabelFrame(parent, text='Parameters', bg='#f0f0f0')
        param_frame.pack(fill='x', padx=5, pady=5)
        
        # Test size
        row = tk.Frame(param_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Test Size:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.test_size_var = tk.StringVar(value='0.3')
        ttk.Spinbox(row, from_=0.0, to=0.9, increment=0.05, textvariable=self.test_size_var, width=10).pack(side='left')
        tk.Label(row, text='(0 = CV only)', bg='#f0f0f0', font=('Arial', 7), fg='#7f8c8d').pack(side='left', padx=(5, 0))
        
        # CV folds
        self.cv_frame = tk.Frame(param_frame, bg='#f0f0f0')
        self.cv_frame.pack(fill='x', padx=5, pady=2)
        tk.Label(self.cv_frame, text='CV Folds:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.cv_folds_var = tk.StringVar(value='5')
        ttk.Spinbox(self.cv_frame, from_=2, to=10, increment=1, textvariable=self.cv_folds_var, width=10).pack(side='left')
        
        # Scaling
        row = tk.Frame(param_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Scaling:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.scaling_var = tk.StringVar(value='standard')
        ttk.Combobox(row, textvariable=self.scaling_var, values=['standard', 'robust', 'none'], state='readonly', width=12).pack(side='left')

        # Class weighting for imbalanced classes
        row = tk.Frame(param_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Class Weight:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.class_weight_var = tk.StringVar(value='none')
        ttk.Combobox(row, textvariable=self.class_weight_var, values=['none', 'balanced'], state='readonly', width=12).pack(side='left')
        tk.Label(row, text="  'balanced' for minority groups", font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)
        
        # Regularization settings (for linear models)
        reg_frame = tk.LabelFrame(parent, text='⚙️ Regularization (Linear Models)', bg='#f0f0f0')
        reg_frame.pack(fill='x', padx=5, pady=5)
        
        # Regularization type
        row = tk.Frame(reg_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Type:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.regularization_type_var = tk.StringVar(value='l2')
        ttk.Combobox(row, textvariable=self.regularization_type_var, 
                    values=['l2', 'l1', 'elasticnet'], state='readonly', width=12).pack(side='left')
        tk.Label(row, text='  L2=Ridge, L1=Lasso', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)
        
        # Regularization strength
        row = tk.Frame(reg_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Strength (C):', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.regularization_strength_var = tk.StringVar(value='medium')
        strength_combo = ttk.Combobox(row, textvariable=self.regularization_strength_var, 
                                     values=['strong', 'medium', 'weak'], state='readonly', width=12)
        strength_combo.pack(side='left')
        tk.Label(row, text='  Lower=stronger', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)
        
        # Max iterations
        row = tk.Frame(reg_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Max Iterations:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.max_iter_var = tk.StringVar(value='1000')
        ttk.Combobox(row, textvariable=self.max_iter_var, 
                    values=['500', '1000', '2000', '5000'], state='readonly', width=12).pack(side='left')
        
        # Robustness settings
        robust_frame = tk.LabelFrame(parent, text='🔁 Robustness Testing', bg='#f0f0f0')
        robust_frame.pack(fill='x', padx=5, pady=5)
        
        # Repeated runs
        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Repeated Runs:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.repeated_runs_var = tk.StringVar(value='1')
        ttk.Combobox(row, textvariable=self.repeated_runs_var, 
                    values=['1', '5', '10', '20', '30'], state='readonly', width=12).pack(side='left')
        tk.Label(row, text='  10-20 recommended', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)
        
        # Base random seed for reproducibility
        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Base Seed:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.base_seed_var = tk.StringVar(value='42')
        seed_spin = ttk.Spinbox(row, from_=0, to=1000, increment=1, textvariable=self.base_seed_var, width=8)
        seed_spin.pack(side='left', padx=5)
        tk.Label(row, text='(Runs: seed, seed+1,...)', 
            font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)
        
        # Stability tracking
        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.stability_tracking_var = tk.BooleanVar(value=True)
        tk.Checkbutton(row, text='Enable feature stability tracking', 
                      variable=self.stability_tracking_var, bg='#f0f0f0',
                      font=('Arial', 9)).pack(anchor='w', padx=5)
        
        # Stability threshold
        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Stability Threshold:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.stability_threshold_var = tk.StringVar(value='70')
        ttk.Spinbox(row, from_=50, to=100, increment=5, textvariable=self.stability_threshold_var, width=10).pack(side='left')
        tk.Label(row, text='  % (keep features ≥ this)', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)

        # Auto plot generation
        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.auto_generate_plots_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            row,
            text='Auto-generate figures (ROC, model bars, top metabolites)',
            variable=self.auto_generate_plots_var,
            command=self._save_ml_config,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(anchor='w', padx=5)

        row = tk.Frame(robust_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=(2, 4))
        tk.Button(
            row,
            text='Configure Figures',
            command=self._open_figure_settings_dialog,
            bg='#2c3e50',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            bd=2,
            padx=6,
            pady=2
        ).pack(anchor='w', padx=5)

        # Advanced validation/tuning controls
        adv_frame = tk.LabelFrame(parent, text='🧪 Advanced Validation & Tuning', bg='#f0f0f0')
        adv_frame.pack(fill='x', padx=5, pady=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.tune_hyperparameters_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            row,
            text='Enable hyperparameter tuning (inner CV when nested CV is ON)',
            variable=self.tune_hyperparameters_var,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(anchor='w', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Tuning Strategy:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.tuning_strategy_var = tk.StringVar(value='grid')
        ttk.Combobox(
            row,
            textvariable=self.tuning_strategy_var,
            values=['grid', 'random'],
            state='readonly',
            width=12
        ).pack(side='left')
        tk.Label(row, text='  grid=full, random=faster', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Random Iterations:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.tuning_iter_var = tk.StringVar(value='20')
        ttk.Spinbox(row, from_=5, to=200, increment=5, textvariable=self.tuning_iter_var, width=10).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.use_repeated_cv_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row,
            text='Use Repeated Stratified CV',
            variable=self.use_repeated_cv_var,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(anchor='w', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='CV Repeats:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.cv_repeats_var = tk.StringVar(value='3')
        ttk.Spinbox(row, from_=1, to=20, increment=1, textvariable=self.cv_repeats_var, width=10).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.nested_cv_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row,
            text='Enable nested CV (outer estimate + inner tuning)',
            variable=self.nested_cv_var,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(anchor='w', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='SVM Calibration:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.calibration_method_var = tk.StringVar(value='none')
        ttk.Combobox(
            row,
            textvariable=self.calibration_method_var,
            values=['none', 'sigmoid', 'isotonic'],
            state='readonly',
            width=12
        ).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Permutation Runs:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.permutation_test_runs_var = tk.StringVar(value='0')
        ttk.Spinbox(row, from_=0, to=500, increment=10, textvariable=self.permutation_test_runs_var, width=10).pack(side='left')
        tk.Label(row, text='  0 disables test', font=('Arial', 7), fg='#666', bg='#f0f0f0').pack(side='left', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Imputation:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.imputation_method_var = tk.StringVar(value='half_min')
        ttk.Combobox(
            row,
            textvariable=self.imputation_method_var,
            values=['none', 'half_min', 'median_per_group', 'median_global', 'knn'],
            state='readonly',
            width=16
        ).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='KNN Neighbors:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.imputation_knn_neighbors_var = tk.StringVar(value='5')
        ttk.Spinbox(row, from_=1, to=20, increment=1, textvariable=self.imputation_knn_neighbors_var, width=10).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        self.auto_skip_scaling_tree_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            row,
            text='Auto skip scaling for tree models (RF/GB)',
            variable=self.auto_skip_scaling_tree_var,
            bg='#f0f0f0',
            font=('Arial', 9)
        ).pack(anchor='w', padx=5)

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Feature Select:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.feature_selection_method_var = tk.StringVar(value='none')
        ttk.Combobox(
            row,
            textvariable=self.feature_selection_method_var,
            values=['none', 'variance', 'univariate', 'lasso', 'rf_rfe'],
            state='readonly',
            width=16
        ).pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Variance %ile:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.variance_percentile_var = tk.StringVar(value='10')
        self.variance_percentile_spin = ttk.Spinbox(row, from_=0, to=100, increment=5, textvariable=self.variance_percentile_var, width=10)
        self.variance_percentile_spin.pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='Select K:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.univariate_k_var = tk.StringVar(value='50')
        self.univariate_k_spin = ttk.Spinbox(row, from_=5, to=5000, increment=5, textvariable=self.univariate_k_var, width=10)
        self.univariate_k_spin.pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='LASSO C:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.lasso_c_var = tk.StringVar(value='0.1')
        self.lasso_c_combo = ttk.Combobox(row, textvariable=self.lasso_c_var, values=['0.01', '0.1', '1.0', '10.0'], state='readonly', width=10)
        self.lasso_c_combo.pack(side='left')

        row = tk.Frame(adv_frame, bg='#f0f0f0')
        row.pack(fill='x', padx=5, pady=2)
        tk.Label(row, text='RFE Features:', width=15, anchor='w', bg='#f0f0f0').pack(side='left')
        self.rfe_n_features_var = tk.StringVar(value='50')
        self.rfe_n_features_spin = ttk.Spinbox(row, from_=5, to=5000, increment=5, textvariable=self.rfe_n_features_var, width=10)
        self.rfe_n_features_spin.pack(side='left')

        # Restore defaults button
        button_row = tk.Frame(adv_frame, bg='#f0f0f0')
        button_row.pack(fill='x', padx=5, pady=(8, 4))
        tk.Button(
            button_row,
            text='Restore Defaults',
            command=self._restore_ml_defaults,
            bg='#34495e',
            fg='white',
            font=('Arial', 8, 'bold'),
            relief='raised',
            bd=2,
            padx=6,
            pady=2
        ).pack(anchor='w')

        # Keep feature-selection controls in sync with selected method.
        self.feature_selection_method_var.trace_add('write', lambda *_: self._update_feature_selection_controls())
        self._update_feature_selection_controls()
    
    def _create_step5_actions(self, parent):
        """Step 4: Action buttons"""
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 5}
        
        self.run_button = tk.Button(
            parent,
            text='▶️ Run Analysis',
            command=self._run_analysis,
            bg='#27ae60',
            fg='white',
            state='disabled',
            **btn_style
        )
        self.run_button.pack(fill='x', padx=0, pady=(0, 5))
        
        self.test_models_button = tk.Button(
            parent,
            text='🔬 Test Models',
            command=self._test_models,
            bg='#3498db',
            fg='white',
            state='disabled',
            **btn_style
        )
        self.test_models_button.pack(fill='x', padx=0, pady=(0, 5))
        
        self.clear_log_button = tk.Button(
            parent,
            text='🗑️ Clear Log',
            command=self._clear_log,
            bg='#e74c3c',
            fg='white',
            **btn_style
        )
        self.clear_log_button.pack(fill='x', padx=0, pady=0)
    
    def _on_mode_change(self):
        """Handle data mode change"""
        # Machine Learning tab is intentionally fixed to custom mode.
        if self.ml_data_mode.get() != 'custom':
            self.ml_data_mode.set('custom')
        if hasattr(self, 'mode_desc_label'):
            self.mode_desc_label.config(text='Custom mode for preprocessed/combined data')
        self._log("\n📌 Data mode fixed to: CUSTOM\n")
    
    def _on_analysis_type_change(self):
        """Handle analysis type change"""
        # Analysis type is fixed to classification in this tab.
        if self.analysis_type_var.get() != 'classification':
            self.analysis_type_var.set('classification')
        self.model_frame.pack(fill='x', padx=5, pady=5)
        self.cv_frame.pack(fill='x', padx=5, pady=2)

    def _pairwise_key(self, g1: str, g2: str) -> str:
        """Stable key for a pair of group names."""
        a = str(g1).strip()
        b = str(g2).strip()
        return '||'.join(sorted([a, b], key=lambda x: x.lower()))

    def _get_available_pvalue_columns_for_pairwise(self) -> List[str]:
        """Collect likely p-value columns from loaded/verified data."""
        candidates = []

        for df in (getattr(self, 'pos_df', None), getattr(self, 'neg_df', None), getattr(self, 'merged_data', None)):
            if isinstance(df, pd.DataFrame):
                candidates.extend([str(c) for c in df.columns])

        # Include verified p-value assignments (both mapping directions).
        for mapping in (getattr(self, 'verified_pos_assignments', {}) or {}, getattr(self, 'verified_neg_assignments', {}) or {}):
            for k, v in mapping.items():
                ks = str(k).strip().lower()
                vs = str(v).strip().lower()
                if ks in {'pvalue', 'p-value', 'p value'} and isinstance(v, str) and v.strip():
                    candidates.append(v)
                if vs in {'pvalue', 'p-value', 'p value'} and isinstance(k, str) and k.strip():
                    candidates.append(k)

        seen = set()
        unique = []
        for c in candidates:
            if c not in seen:
                seen.add(c)
                unique.append(c)

        # Prefer p-value-like names.
        def _is_pv(name: str) -> bool:
            n = re.sub(r'[^a-z0-9]+', '', str(name).lower())
            for k in ('pvalue', 'pval', 'padj', 'adjp', 'adjustedp', 'fdr', 'qvalue'):
                if k in n:
                    return True
            return False

        filtered = [c for c in unique if _is_pv(c)]
        return filtered if filtered else unique

    def _open_pairwise_pvalue_verification_dialog(self):
        """Popup to let users map each pairwise comparison to a p-value column."""
        if not getattr(self, 'generate_pairwise_addon_var', None) or not self.generate_pairwise_addon_var.get():
            messagebox.showinfo(
                'Pairwise Mode Required',
                'Enable Pairwise add-on first, then verify pairwise p-value columns.'
            )
            return

        groups = sorted({str(var.get()).strip() for var in getattr(self, 'sample_group_vars', {}).values() if str(var.get()).strip()})
        if len(groups) < 2:
            groups = sorted({str(v).strip() for v in getattr(self, 'group_definitions', {}).values() if str(v).strip()})

        if len(groups) < 2:
            messagebox.showwarning('No Groups', 'Define at least two groups before verifying pairwise p-values.')
            return

        pairs = list(combinations(groups, 2))
        pvalue_columns = self._get_available_pvalue_columns_for_pairwise()

        dialog = tk.Toplevel(self.root)
        dialog.title('Verify Pairwise P-Value Columns')
        dialog.geometry('760x520')
        dialog.transient(self.root)
        dialog.grab_set()

        header = tk.Label(
            dialog,
            text='Map each pairwise comparison to its p-value column.',
            font=('Arial', 10, 'bold'),
            anchor='w',
            justify='left'
        )
        header.pack(fill='x', padx=10, pady=(10, 6))

        hint = tk.Label(
            dialog,
            text='Required when P-Value filter is enabled for pairwise ML.',
            font=('Arial', 8),
            fg='#555555',
            anchor='w',
            justify='left'
        )
        hint.pack(fill='x', padx=10, pady=(0, 8))

        body = tk.Frame(dialog)
        body.pack(fill='both', expand=True, padx=10, pady=(0, 10))

        canvas = tk.Canvas(body, highlightthickness=0)
        scrollbar = ttk.Scrollbar(body, orient='vertical', command=canvas.yview)
        scroll_frame = tk.Frame(canvas)
        scroll_frame.bind('<Configure>', lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=scroll_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        tk.Label(scroll_frame, text='Pairwise Comparison', font=('Arial', 9, 'bold')).grid(row=0, column=0, sticky='w', padx=(4, 12), pady=(4, 6))
        tk.Label(scroll_frame, text='P-Value Column', font=('Arial', 9, 'bold')).grid(row=0, column=1, sticky='w', padx=4, pady=(4, 6))

        combo_vars = {}
        combo_values = pvalue_columns if pvalue_columns else ['(No p-value columns found)']
        for i, (g1, g2) in enumerate(pairs, start=1):
            pair_label = f'{g1} vs {g2}'
            pair_key = self._pairwise_key(g1, g2)
            tk.Label(scroll_frame, text=pair_label, font=('Arial', 9)).grid(row=i, column=0, sticky='w', padx=(4, 12), pady=3)

            default_value = combo_values[0] if combo_values else ''
            current = self.pairwise_pvalue_map.get(pair_key, default_value)
            if current not in combo_values:
                combo_values_with_current = combo_values + [current]
            else:
                combo_values_with_current = combo_values

            var = tk.StringVar(value=current)
            combo = ttk.Combobox(scroll_frame, textvariable=var, values=combo_values_with_current, state='readonly', width=54)
            combo.grid(row=i, column=1, sticky='w', padx=4, pady=3)
            combo_vars[pair_key] = var

        footer = tk.Frame(dialog)
        footer.pack(fill='x', padx=10, pady=(0, 10))

        def _save_and_close():
            updated = {}
            for k, var in combo_vars.items():
                selected = str(var.get()).strip()
                if selected and selected != '(No p-value columns found)':
                    updated[k] = selected
            self.pairwise_pvalue_map = updated
            self._save_ml_config()
            self._log(f"✅ Saved pairwise p-value verification for {len(updated)} pair(s)\n")
            dialog.destroy()

        tk.Button(footer, text='Save', command=_save_and_close, bg='#27ae60', fg='white', font=('Arial', 9, 'bold')).pack(side='left')
        tk.Button(footer, text='Cancel', command=dialog.destroy, bg='#e74c3c', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=(8, 0))

    def _update_feature_selection_controls(self):
        """Enable only the feature-selection parameters relevant to the selected method."""
        method = str(self.feature_selection_method_var.get()).strip().lower()

        def _set_state(widget, enabled: bool):
            if widget is not None:
                widget.config(state='normal' if enabled else 'disabled')

        _set_state(getattr(self, 'variance_percentile_spin', None), method == 'variance')
        _set_state(getattr(self, 'univariate_k_spin', None), method == 'univariate')
        _set_state(getattr(self, 'lasso_c_combo', None), method == 'lasso')
        _set_state(getattr(self, 'rfe_n_features_spin', None), method == 'rf_rfe')

    def _restore_ml_defaults(self):
        """Restore the built-in default values and save them as the current config."""
        self.test_size_var.set('0.3')
        self.cv_folds_var.set('5')
        self.scaling_var.set('standard')
        self.class_weight_var.set('none')
        self.regularization_type_var.set('l2')
        self.regularization_strength_var.set('medium')
        self.max_iter_var.set('1000')
        self.repeated_runs_var.set('1')
        self.base_seed_var.set('42')
        self.stability_tracking_var.set(True)
        self.stability_threshold_var.set('70')
        self.auto_generate_plots_var.set(True)
        self.generate_pairwise_addon_var.set(False)
        self.pairwise_pvalue_map = {}
        self.top_n_values = self._default_top_n_values()

        self.tune_hyperparameters_var.set(True)
        self.tuning_strategy_var.set('grid')
        self.tuning_iter_var.set('20')
        self.use_repeated_cv_var.set(False)
        self.cv_repeats_var.set('3')
        self.nested_cv_var.set(False)
        self.calibration_method_var.set('none')
        self.permutation_test_runs_var.set('0')
        self.imputation_method_var.set('half_min')
        self.imputation_knn_neighbors_var.set('5')
        self.auto_skip_scaling_tree_var.set(False)
        self.feature_selection_method_var.set('none')
        self.variance_percentile_var.set('10')
        self.univariate_k_var.set('50')
        self.lasso_c_var.set('0.1')
        self.rfe_n_features_var.set('50')

        self._update_feature_selection_controls()
        self._save_ml_config()
        self._log("\n♻️ Restored ML defaults.\n")
    
    def _import_excel(self):
        """Import Excel file with Pos/Neg sheets (same logic as Statistics tab)"""
        mode = self.ml_data_mode.get()
        
        title = {
            'metabolite': 'Select Metabolite Excel (with Pos_id/Neg_id)',
            'lipid': 'Select Lipid Excel (with Positive_Lipids/Negative_Lipids)',
            'custom': 'Select Preprocessed Excel'
        }.get(mode, 'Select Excel File')
        
        file_path = filedialog.askopenfilename(
            title=title,
            filetypes=[('Excel files', '*.xlsx *.xls'), ('All files', '*.*')]
        )
        
        if not file_path:
            return
        
        try:
            self._log(f"\n{'='*60}\n")
            self._log(f"📂 Importing {mode.upper()} Excel file\n")
            self._log(f"{'='*60}\n")
            self._log(f"File: {os.path.basename(file_path)}\n\n")
            
            xl = pd.ExcelFile(file_path)
            sheet_names = [str(name) for name in xl.sheet_names]
            self._log(f"Found {len(sheet_names)} sheets: {', '.join(sheet_names)}\n\n")
            
            # Reset previous data
            self.pos_df = None
            self.neg_df = None
            self.detected_pos_sample_cols = []
            self.detected_neg_sample_cols = []
            self.detected_sample_cols = []
            
            if mode == 'custom':
                # Custom mode: Load single combined sheet
                self._import_custom_sheet(xl, file_path)
                return
            
            # Find Pos/Neg sheets based on mode
            pos_sheet = None
            neg_sheet = None
            
            if mode == 'lipid':
                for candidate in ['Positive_Lipids', 'Positive_Lipid', 'Pos_Lipids', 'Pos_Lipid']:
                    if candidate in xl.sheet_names:
                        pos_sheet = candidate
                        break
                for candidate in ['Negative_Lipids', 'Negative_Lipid', 'Neg_Lipids', 'Neg_Lipid']:
                    if candidate in xl.sheet_names:
                        neg_sheet = candidate
                        break
            else:  # metabolite
                for candidate in ['Pos_id', 'Positive', 'Pos', 'POS']:
                    if candidate in xl.sheet_names:
                        pos_sheet = candidate
                        break
                for candidate in ['Neg_id', 'Negative', 'Neg', 'NEG']:
                    if candidate in xl.sheet_names:
                        neg_sheet = candidate
                        break
            
            if not pos_sheet and not neg_sheet:
                messagebox.showerror(
                    'No Expected Sheets',
                    f'Could not find expected sheets for {mode} mode.\n\n'
                    f'Expected: Pos_id/Neg_id (metabolite) or Positive_Lipids/Negative_Lipids (lipid)\n\n'
                    f'Try using Custom mode if your data is preprocessed.'
                )
                self._log("❌ No expected sheets found\n")
                return
            
            # Load sheets
            if pos_sheet:
                self.pos_df = xl.parse(pos_sheet)
                self._log(f"✅ Loaded {pos_sheet}: {len(self.pos_df)} rows, {len(self.pos_df.columns)} columns\n")
                
                # Detect sample columns (same logic as Statistics tab)
                self.detected_pos_sample_cols = self._detect_sample_columns(self.pos_df, mode)
                self._log(f"   Detected {len(self.detected_pos_sample_cols)} sample columns\n")
            
            if neg_sheet:
                self.neg_df = xl.parse(neg_sheet)
                self._log(f"✅ Loaded {neg_sheet}: {len(self.neg_df)} rows, {len(self.neg_df.columns)} columns\n")
                
                # Detect sample columns
                self.detected_neg_sample_cols = self._detect_sample_columns(self.neg_df, mode)
                self._log(f"   Detected {len(self.detected_neg_sample_cols)} sample columns\n")
            
            # Union of sample columns
            seen = set()
            union_cols = []
            for cols in [self.detected_pos_sample_cols, self.detected_neg_sample_cols]:
                for c in cols:
                    if c not in seen:
                        seen.add(c)
                        union_cols.append(c)
            
            self.detected_sample_cols = union_cols
            self._log(f"\n📊 Total unique sample columns: {len(union_cols)}\n")
            
            # Update status
            self.data_status_label.config(
                text=f"✅ Loaded: {len(union_cols)} samples (not verified yet)",
                fg="orange"
            )
            
            self._log(f"\n{'='*60}\n")
            self._log("✅ Import complete! Next: Click '🔍 Verify Columns'\n")
            self._log(f"{'='*60}\n\n")
            
        except Exception as e:
            logger.error(f"Error importing Excel: {e}", exc_info=True)
            self._log(f"\n❌ Import failed: {str(e)}\n")
            messagebox.showerror("Import Error", f"Failed to import Excel:\n{str(e)}")
    
    def _import_custom_sheet(self, xl: pd.ExcelFile, file_path: str):
        """Import custom preprocessed sheet"""
        if len(xl.sheet_names) == 0:
            messagebox.showerror("No Sheets", "No sheets found in Excel file")
            return
        
        # Use first sheet or ask user
        sheet_name = xl.sheet_names[0]
        if len(xl.sheet_names) > 1:
            from tkinter import simpledialog
            sheet_name = simpledialog.askstring(
                "Select Sheet",
                f"Multiple sheets found:\n" + "\n".join(str(name) for name in xl.sheet_names) + f"\n\nEnter sheet name:",
                initialvalue=str(sheet_name)
            )
            if not sheet_name or sheet_name not in xl.sheet_names:
                self._log("❌ Invalid sheet selected\n")
                return
        
        df = xl.parse(sheet_name)
        self._log(f"✅ Loaded {sheet_name}: {len(df)} rows, {len(df.columns)} columns\n")
        
        # Detect columns
        sample_cols = self._detect_sample_columns(df, 'custom')
        
        # Store as "merged" data
        self.pos_df = df
        self.neg_df = None
        self.detected_sample_cols = sample_cols
        self.detected_pos_sample_cols = sample_cols
        
        self._log(f"📊 Detected {len(sample_cols)} sample columns\n")
        
        # Update status
        self.data_status_label.config(
            text=f"✅ Custom data loaded: {len(sample_cols)} samples",
            fg="orange"
        )
        
        self._log("\n✅ Custom data imported! Next: Click '🔍 Verify Columns'\n\n")
    
    def _detect_sample_columns(self, df: pd.DataFrame, mode: str) -> List[str]:
        """
        Detect sample columns using same logic as Statistics tab.
        
        Sample columns are numeric columns that are NOT feature/metadata columns.
        """
        from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
        
        if mode == 'metabolite' or mode == 'custom':
            # Use metabolite detection
            _, sample_cols = detect_feature_and_sample_columns(df)
            return sample_cols
        else:
            # Lipid mode: manual detection
            sample_cols = []
            for col in df.columns:
                # Skip if it's a lipid feature column
                if self._is_lipid_feature_col(col):
                    continue
                # Skip if it's a metadata column
                if is_statistics_metadata_col(col):
                    continue
                # Include if numeric
                if pd.api.types.is_numeric_dtype(df[col]):
                    sample_cols.append(col)
            return sample_cols
    
    def _is_lipid_feature_col(self, col_name: str) -> bool:
        """Check if column is a lipid feature column"""
        normalized = str(col_name).lower().strip().replace(' ', '').replace('_', '')
        lipid_features = [
            'lipidid', 'class', 'lipidgroup', 'charge', 'calcmz', 'basert', 'subclass',
            'adduction', 'ionformula', 'molstructure', 'obsmz', 'obsrt', 'ppmdiff', 'polarity',
            'classname', 'lipidclass'
        ]
        return normalized in lipid_features
    
    def _get_feature_id_column(self, df: Optional[pd.DataFrame] = None) -> Optional[str]:
        """Get the verified feature ID column name.
        
        Returns the feature ID column name from verified column assignments,
        or tries to find a suitable column in the dataframe.
        """
        # Priority 1: User verified feature ID column
        if hasattr(self, 'verified_pos_feature_id_col') and self.verified_pos_feature_id_col:
            return self.verified_pos_feature_id_col
        if hasattr(self, 'verified_neg_feature_id_col') and self.verified_neg_feature_id_col:
            return self.verified_neg_feature_id_col
        
        # Priority 2: Check common ID column names in dataframe
        if df is not None:
            for col in ['Name', 'Metabolite', 'Compound', 'LipidIon', 'Feature_ID', 'ID', 
                       'Molecule', 'CompoundName', 'Compound_Name', 'metabolite_name']:
                if col in df.columns:
                    return col
        
        return None
    
    def _verify_columns(self):
        """Verify columns using shared column assignment dialog (same as Statistics tab)"""
        if self.pos_df is None and self.neg_df is None:
            messagebox.showwarning("No Data", "Please import Excel data first")
            return
        
        def _verify_thread():
            """Background thread for column verification"""
            try:
                mode = self.ml_data_mode.get()
                tab_type = 'statistics_metabolite' if mode == 'metabolite' else 'statistics_lipid' if mode == 'lipid' else 'statistics_metabolite'
                
                self._log("\n🔍 Verifying columns...\n")
                
                # Verify input sheet
                if self.pos_df is not None:
                    self._log("📂 Verifying input sheet...\n")
                    
                    pos_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=self.pos_df,
                        tab_type=tab_type,
                        auto_calculate=False,
                        dialog_title=f"Input {mode.capitalize()} - Verify Columns",
                        detected_sample_cols=self.detected_pos_sample_cols,
                        allow_skip=True
                    )
                    
                    if pos_result is None:
                        self._log("❌ Verification cancelled\n")
                        return
                    
                    if pos_result.get('skipped'):
                        self._log("⏭ Skipped input sheet\n")
                        self.pos_df = None
                        self.verified_pos_sample_cols = []
                        self.verified_pos_assignments = {}
                    else:
                        self.verified_pos_assignments = pos_result['assignments']
                        self.verified_pos_sample_cols = pos_result.get('sample_cols', [])
                        # Extract feature ID column from assignments
                        if 'feature_id_col' in pos_result:
                            self.verified_pos_feature_id_col = pos_result['feature_id_col']
                        self._log(f"✅ Input sheet verified: {len(self.verified_pos_sample_cols)} samples")
                        if self.verified_pos_feature_id_col:
                            self._log(f" (Feature ID: {self.verified_pos_feature_id_col})")
                        self._log("\n")
                
                # Verify Negative sheet
                if self.neg_df is not None:
                    self._log("📂 Verifying Negative sheet...\n")
                    
                    neg_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=self.neg_df,
                        tab_type=tab_type,
                        auto_calculate=False,
                        dialog_title=f"Negative {mode.capitalize()} - Verify Columns",
                        detected_sample_cols=self.detected_neg_sample_cols,
                        allow_skip=True
                    )
                    
                    if neg_result is None:
                        self._log("❌ Verification cancelled\n")
                        return
                    
                    if neg_result.get('skipped'):
                        self._log("⏭ Skipped Negative sheet\n")
                        self.neg_df = None
                        self.verified_neg_sample_cols = []
                        self.verified_neg_assignments = {}
                    else:
                        self.verified_neg_assignments = neg_result['assignments']
                        self.verified_neg_sample_cols = neg_result.get('sample_cols', [])
                        # Extract feature ID column from assignments
                        if 'feature_id_col' in neg_result:
                            self.verified_neg_feature_id_col = neg_result['feature_id_col']
                        self._log(f"✅ Negative verified: {len(self.verified_neg_sample_cols)} samples")
                        if self.verified_neg_feature_id_col:
                            self._log(f" (Feature ID: {self.verified_neg_feature_id_col})")
                        self._log("\n")
                
                # Union of verified sample columns
                seen = set()
                union_cols = []
                for cols in [self.verified_pos_sample_cols, self.verified_neg_sample_cols]:
                    for c in cols:
                        if c not in seen:
                            seen.add(c)
                            union_cols.append(c)
                
                if len(union_cols) == 0:
                    self.root.after(0, lambda: messagebox.showerror(
                        "No Samples",
                        "No sample columns verified. Please try again."
                    ))
                    return
                
                # Create sample group vars for Configure Groups
                self.sample_group_vars = {}
                for col in union_cols:
                    self.sample_group_vars[col] = tk.StringVar(value='')
                
                # Update status
                def update_ui():
                    self.data_status_label.config(
                        text=f"✅ Verified: {len(union_cols)} samples (groups not assigned)",
                        fg="orange"
                    )
                    self.configure_groups_btn.config(state='normal')
                
                self.root.after(0, update_ui)
                
                self._log(f"\n✅ Column verification complete!\n")
                self._log(f"Total verified samples: {len(union_cols)}\n")
                self._log(f"\nNext: Click '⚙️ Configure Groups'\n\n")
                
            except Exception as e:
                logger.error(f"Error verifying columns: {e}", exc_info=True)
                self._log(f"\n❌ Verification error: {str(e)}\n")
                self.root.after(0, lambda: messagebox.showerror("Verification Error", str(e)))
        
        thread = threading.Thread(target=_verify_thread, daemon=True)
        thread.start()
    
    def _configure_groups(self):
        """Configure sample groups using pattern matching (same as Statistics tab)"""
        if not self.sample_group_vars:
            messagebox.showwarning("No Samples", "Please verify columns first")
            return
        
        # Sync group_definitions from entry fields before opening dialog
        for group_id, var in self.group_id_vars.items():
            label = var.get().strip()
            if label:
                self.group_definitions[group_id] = label
        
        # Open pattern matching dialog (simplified version of Statistics tab's auto_assign_groups)
        dialog = tk.Toplevel(self.root)
        dialog.title('Configure Sample Groups')
        dialog.geometry('600x700')
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Main frame with scrollbar
        main_frame = tk.Frame(dialog, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        canvas = tk.Canvas(main_frame, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable = tk.Frame(canvas, bg='#f0f0f0')
        
        scrollable.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Title
        tk.Label(
            scrollable,
            text='Auto-Assign Groups by Pattern',
            font=('Arial', 12, 'bold'),
            bg='#f0f0f0'
        ).pack(pady=(0, 10))
        
        # Pattern definition
        pattern_frame = tk.LabelFrame(scrollable, text='Define Patterns', bg='#f0f0f0')
        pattern_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        tk.Label(
            pattern_frame,
            text='Enter keywords/patterns for each group (one per line):',
            bg='#f0f0f0'
        ).pack(anchor='w', padx=5, pady=5)
        
        pattern_vars = {}
        for group_id, label in self.group_definitions.items():
            group_frame = tk.LabelFrame(pattern_frame, text=f'{group_id}: {label}', bg='#f0f0f0')
            group_frame.pack(fill='x', padx=5, pady=5)
            
            pattern_text = tk.Text(group_frame, height=3, font=('Arial', 9))
            pattern_text.pack(fill='x', padx=5, pady=5)
            
            # Load saved patterns if available
            if group_id in self.auto_assign_patterns:
                saved_patterns = self.auto_assign_patterns[group_id]
                if isinstance(saved_patterns, list):
                    pattern_text.insert('1.0', '\n'.join(saved_patterns))
                elif isinstance(saved_patterns, str):
                    pattern_text.insert('1.0', saved_patterns)
            
            pattern_vars[group_id] = pattern_text
        
        # Current assignments display
        assign_frame = tk.LabelFrame(scrollable, text='Current Assignments', bg='#f0f0f0')
        assign_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Header
        header = tk.Frame(assign_frame, bg='#e8f4f8')
        header.pack(fill='x', pady=(0, 2))
        tk.Label(header, text='Sample', bg='#e8f4f8', font=('Arial', 9, 'bold'), width=30, anchor='w').pack(side='left', padx=5)
        tk.Label(header, text='Group', bg='#e8f4f8', font=('Arial', 9, 'bold'), width=20, anchor='w').pack(side='left', padx=5)
        
        # Sample rows with dropdowns
        assignment_canvas = tk.Canvas(assign_frame, bg='white', height=200)
        assignment_scrollbar = ttk.Scrollbar(assign_frame, orient='vertical', command=assignment_canvas.yview)
        assignment_scrollable = tk.Frame(assignment_canvas, bg='white')
        
        assignment_scrollable.bind("<Configure>", lambda e: assignment_canvas.configure(scrollregion=assignment_canvas.bbox("all")))
        assignment_canvas.create_window((0, 0), window=assignment_scrollable, anchor="nw")
        assignment_canvas.configure(yscrollcommand=assignment_scrollbar.set)
        
        assignment_scrollbar.pack(side="right", fill="y")
        assignment_canvas.pack(side="left", fill="both", expand=True)
        
        group_labels = list(self.group_definitions.values())
        
        for col_name in sorted(self.sample_group_vars.keys()):
            row = tk.Frame(assignment_scrollable, bg='white')
            row.pack(fill='x', pady=1)
            
            tk.Label(row, text=col_name, bg='white', font=('Arial', 8), width=30, anchor='w').pack(side='left', padx=5)
            
            ttk.Combobox(
                row,
                values=[''] + group_labels,
                textvariable=self.sample_group_vars[col_name],
                state='readonly',
                width=18,
                font=('Arial', 8)
            ).pack(side='left', padx=5)
        
        # Auto-apply patterns button
        def apply_patterns():
            for group_id, pattern_text in pattern_vars.items():
                patterns = pattern_text.get('1.0', tk.END).strip().split('\n')
                group_label = self.group_definitions[group_id]
                
                count = 0
                for pattern in patterns:
                    pattern = pattern.strip()
                    if not pattern:
                        continue
                    
                    for col_name in self.sample_group_vars.keys():
                        if pattern.lower() in col_name.lower():
                            self.sample_group_vars[col_name].set(group_label)
                            count += 1
                
                if count > 0:
                    self._log(f"✓ Assigned {count} samples to {group_label} using patterns\n")
        
        btn_frame = tk.Frame(scrollable, bg='#f0f0f0')
        btn_frame.pack(fill='x', pady=10)
        
        tk.Button(
            btn_frame,
            text='🔄 Apply Patterns',
            command=apply_patterns,
            bg='#3498db',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=5)
        
        def save_and_close():
            # Sync latest group labels from entry fields before saving config
            for group_id, var in self.group_id_vars.items():
                label = var.get().strip()
                if label:
                    self.group_definitions[group_id] = label

            # Save patterns for future use
            for group_id, pattern_text in pattern_vars.items():
                patterns = pattern_text.get('1.0', tk.END).strip()
                if patterns:
                    self.auto_assign_patterns[group_id] = patterns.split('\n')
                else:
                    self.auto_assign_patterns[group_id] = []

            # Persist group definitions + patterns immediately
            self._save_ml_config()
            
            # Collect assignments
            group_map = {}
            for col, var in self.sample_group_vars.items():
                group = var.get().strip()
                if group:
                    group_map[col] = group
            
            if len(group_map) == 0:
                messagebox.showwarning("No Assignments", "Please assign at least one sample to a group")
                return
            
            # Check unassigned
            unassigned = [c for c in self.sample_group_vars.keys() if c not in group_map]
            if unassigned:
                response = messagebox.askyesno(
                    "Unassigned Samples",
                    f"{len(unassigned)} samples not assigned.\n\nContinue anyway?"
                )
                if not response:
                    return
            
            # Update status
            n_groups = len(set(group_map.values()))
            self.data_status_label.config(
                text=f"✅ Ready: {len(group_map)} samples, {n_groups} groups",
                fg="green"
            )
            self.run_button.config(state='normal')
            self.test_models_button.config(state='normal')
            
            self._log(f"\n✅ Group configuration complete!\n")
            self._log(f"Assigned samples: {len(group_map)}\n")
            self._log(f"Groups: {n_groups}\n")
            
            # Log group counts
            counts = Counter(group_map.values())
            for group, count in sorted(counts.items()):
                self._log(f"  • {group}: {count} samples\n")
            
            self._log(f"\nReady for ML analysis!\n\n")
            
            dialog.destroy()
        
        tk.Button(
            btn_frame,
            text='✅ Save & Close',
            command=save_and_close,
            bg='#27ae60',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=5)
        
        tk.Button(
            btn_frame,
            text='❌ Cancel',
            command=dialog.destroy,
            bg='#e74c3c',
            fg='white',
            font=('Arial', 9, 'bold')
        ).pack(side='left', padx=5)
    
    def refresh_group_ui(self):
        """Refresh the group UI with current group definitions"""
        # Clear existing group entries
        for widget in self.groups_scrollable_frame.winfo_children():
            widget.destroy()

        # Helper: ensure each group StringVar has a write-trace to auto-save edits
        def _ensure_group_trace(gid: str, var: tk.StringVar):
            if getattr(var, '_ml_trace_bound', False):
                return

            def update_group_def(*args, group_id=gid, entry_var=var):
                new_label = entry_var.get().strip()
                if new_label:
                    self.group_definitions[group_id] = new_label
                    self._save_ml_config()

            var.trace_add('write', update_group_def)
            setattr(var, '_ml_trace_bound', True)
        
        # Recreate group entries
        for i, (group_id, default_label) in enumerate(self.group_definitions.items()):
            id_frame = tk.Frame(self.groups_scrollable_frame, bg='#f0f0f0')
            id_frame.pack(fill='x', padx=3, pady=2)
            tk.Label(id_frame, text=f'{group_id}:', bg='#f0f0f0', width=8).pack(side='left')
            
            # Create or reuse StringVar
            if group_id not in self.group_id_vars:
                self.group_id_vars[group_id] = tk.StringVar(value=default_label)
            else:
                # Keep UI var in sync with persisted definitions on refresh
                if self.group_id_vars[group_id].get() != default_label:
                    self.group_id_vars[group_id].set(default_label)
            
            # Create entry widget with the StringVar
            entry_var = self.group_id_vars[group_id]
            _ensure_group_trace(group_id, entry_var)
            tk.Entry(id_frame, textvariable=entry_var, font=('Arial', 9), width=20).pack(side='left', padx=(5,5))
    
    def add_group(self):
        """Add a new group to the group definitions"""
        self.group_count += 1
        new_group_id = f'Group{self.group_count}'
        self.group_definitions[new_group_id] = f'Group{self.group_count}'
        self.group_id_vars[new_group_id] = tk.StringVar(value=f'Group{self.group_count}')
        self.refresh_group_ui()
        self._save_ml_config()
        self._log(f"✅ Added {new_group_id}\n")
    
    def remove_group(self):
        """Remove the last group from group definitions"""
        if self.group_count <= 2:  # Don't allow less than 2 groups
            messagebox.showwarning("Minimum Groups", "At least 2 groups are required for ML analysis.")
            return
        
        last_group_id = f'Group{self.group_count}'
        if last_group_id in self.group_definitions:
            del self.group_definitions[last_group_id]
            if last_group_id in self.group_id_vars:
                del self.group_id_vars[last_group_id]
            self.group_count -= 1
            self.refresh_group_ui()
            self._save_ml_config()
            self._log(f"✅ Removed {last_group_id}\n")
    
    def _run_analysis(self):
        """Run ML analysis with automatic filtering and merging"""
        if self.auto_generate_plots_var.get() or self.generate_pairwise_addon_var.get():
            working_folder = self.working_folder_var.get().strip()
            if not working_folder:
                messagebox.showwarning(
                    "Working Folder Required",
                    "Figure generation and/or pairwise add-on output is enabled.\n\nPlease select a Working Folder first (Step 1)."
                )
                return
            if not os.path.isdir(working_folder):
                messagebox.showerror(
                    "Invalid Working Folder",
                    "The selected Working Folder does not exist.\n\nPlease choose a valid folder before running analysis."
                )
                return

        self.run_button.config(state='disabled', text='⏳ Running...')
        
        thread = threading.Thread(target=self._run_analysis_thread, daemon=True)
        thread.start()
    
    def _run_analysis_thread(self):
        """Background thread for ML analysis"""
        try:
            from time import perf_counter
            analysis_start = perf_counter()

            self._log(f"\n{'='*60}\n")
            self._log("🚀 Starting ML Analysis Pipeline\n")
            self._log(f"{'='*60}\n\n")
            
            # Step 1: Filter data
            self._log("Step 1: Filtering data...\n")
            merged_df, group_map, feature_cols = self._filter_and_merge_data()
            
            if merged_df is None or len(merged_df) == 0:
                self._log("❌ No data available after filtering\n")
                self.root.after(0, lambda: messagebox.showerror(
                    "No Data",
                    "No data available after filtering"
                ))
                return
            
            self._log(f"✅ Filtered data ready: {len(merged_df)} features, {len(feature_cols)} samples\n\n")

            # Determine pairwise-only mode early so we can avoid global p-value prefilter.
            pre_unique_groups = sorted({str(v) for v in group_map.values() if str(v).strip()})
            pairwise_only_precalc = bool(self.generate_pairwise_addon_var.get()) and len(pre_unique_groups) >= 3
            
            # ========== Feature Filters (can be combined) ==========
            filters_applied = []
            
            # Helper: Get column names from verified assignments
            def get_assigned_columns(col_type):
                """Find actual column names assigned to a type from verified assignments.

                Note: the column assignment dialog returns a mapping of {Column Type -> Column Name}.
                Older code paths may provide {Column Name -> Column Type}, so support both.
                """

                def _extract(mapping):
                    results = []
                    if not mapping:
                        return results
                    target_norm = str(col_type).strip().lower()

                    # Case A: mapping is {type -> column}
                    for key, value in mapping.items():
                        key_norm = str(key).strip().lower()
                        if key_norm == target_norm and isinstance(value, str) and value.strip():
                            results.append(value)

                    # Case B: mapping is {column -> type}
                    for key, value in mapping.items():
                        value_norm = str(value).strip().lower()
                        if value_norm == target_norm and isinstance(key, str) and key.strip():
                            results.append(key)

                    # De-duplicate while preserving order
                    return list(dict.fromkeys(results))

                matches = []
                matches.extend(_extract(getattr(self, 'verified_pos_assignments', {}) or {}))
                matches.extend(_extract(getattr(self, 'verified_neg_assignments', {}) or {}))
                return list(dict.fromkeys(matches))

            def _debug_log_assignment_state(col_type, df):
                """Log assignment + column presence info to help debug resolution issues."""
                try:
                    pos_map = getattr(self, 'verified_pos_assignments', {}) or {}
                    neg_map = getattr(self, 'verified_neg_assignments', {}) or {}
                    assigned = get_assigned_columns(col_type)

                    self._log(f"  🔎 DEBUG: '{col_type}' assigned columns: {assigned if assigned else 'None'}\n")
                    self._log(f"  🔎 DEBUG: verified_pos_assignments keys example: {list(pos_map.keys())[:8]}\n")
                    self._log(f"  🔎 DEBUG: verified_neg_assignments keys example: {list(neg_map.keys())[:8]}\n")

                    if df is not None:
                        # show any df columns that look related
                        target = str(col_type).strip().lower().replace('_', '').replace(' ', '')
                        related = []
                        for c in df.columns:
                            c_norm = str(c).strip().lower().replace('_', '').replace(' ', '')
                            if target in c_norm or c_norm in [target, f"is{target}"]:
                                related.append(c)
                        if related:
                            self._log(f"  🔎 DEBUG: merged dataframe related columns: {related[:10]}\n")
                except Exception:
                    # Never let debug logging break the analysis
                    pass

            def resolve_column_name(assigned_name, df_columns):
                """Resolve assigned column name against dataframe columns (case-insensitive)."""
                if not assigned_name:
                    return None
                if assigned_name in df_columns:
                    return assigned_name
                assigned_norm = str(assigned_name).strip().lower().replace(' ', '').replace('_', '')
                for col in df_columns:
                    col_norm = str(col).strip().lower().replace(' ', '').replace('_', '')
                    if col_norm == assigned_norm:
                        return col
                return None

            def choose_best_pvalue_column(candidates):
                """Choose best p-value column among candidates (prefer adjusted p-values)."""
                if not candidates:
                    return None
                def _priority(col_name: str) -> tuple:
                    normalized = re.sub(r'[^a-z0-9]', '', str(col_name).lower())
                    if any(keyword in normalized for keyword in ['adjp', 'padj', 'fdr', 'qvalue', 'adjustedp']):
                        return (0, len(normalized))
                    if normalized.endswith('pvalue') or normalized in ['p', 'pval']:
                        return (1, len(normalized))
                    return (2, len(normalized))
                return min(candidates, key=_priority)
            
            # --- Endogenous Filter ---
            if hasattr(self, 'filter_endogenous_var') and self.filter_endogenous_var.get():
                self._log("🔬 Applying Endogenous filter (from Endogenous_Source)...\n")
                
                # Primary: use verified Endogenous_Source (this project uses it as the endogenous signal)
                source_cols = get_assigned_columns('Endogenous_Source')
                if source_cols:
                    source_cols = [resolve_column_name(c, merged_df.columns) for c in source_cols]
                    source_cols = [c for c in source_cols if c is not None]
                    source_cols = list(dict.fromkeys(source_cols))
                if source_cols:
                    self._log(f"  ℹ️ Using Endogenous_Source column(s): {source_cols}\n")

                # Back-compat: if a dataset truly has an Endogenous yes/no column, still support it
                endogenous_cols = get_assigned_columns('Endogenous')
                resolved_endogenous_cols = []
                if endogenous_cols:
                    resolved_endogenous_cols = [resolve_column_name(c, merged_df.columns) for c in endogenous_cols]
                    resolved_endogenous_cols = [c for c in resolved_endogenous_cols if c is not None]
                    resolved_endogenous_cols = list(dict.fromkeys(resolved_endogenous_cols))
                endogenous_col = resolved_endogenous_cols[0] if resolved_endogenous_cols else None
                
                # Fallback: search for column name (case-insensitive)
                if not endogenous_col:
                    for col in merged_df.columns:
                        if col.lower().replace('_', '').replace(' ', '') in ['endogenous', 'isendogenous']:
                            endogenous_col = col
                            break

                # Fallback: also accept endogenous_source column naming
                if not endogenous_col and not source_cols:
                    for col in merged_df.columns:
                        if col.lower().replace('_', '').replace(' ', '') in ['endogenoussource', 'endogenous_source']:
                            source_cols = [col]
                            self._log(f"  ℹ️ Using '{col}' for endogenous filtering (fallback match).\n")
                            break
                
                if (resolved_endogenous_cols or source_cols or (endogenous_col and endogenous_col in merged_df.columns)):
                    before_count = len(merged_df)
                    truthy = {'yes', 'y', 'true', '1', 'endogenous'}
                    if resolved_endogenous_cols:
                        endogenous_mask = pd.Series(False, index=merged_df.index)
                        for col in resolved_endogenous_cols:
                            endogenous_mask |= merged_df[col].astype(str).str.strip().str.lower().isin(truthy)
                    elif source_cols:
                        # DEBUG: Show unique values in Endogenous_Source column
                        for col in source_cols:
                            unique_vals = merged_df[col].astype(str).str.strip().str.lower().unique()
                            unique_preview = list(unique_vals[:10])
                            self._log(f"  🔎 DEBUG: Unique values in '{col}': {unique_preview}\n")
                        
                        # Match endogenous: value contains 'endo' OR equals common truthy values
                        # Exclude if value explicitly contains 'exo' (exogenous)
                        endogenous_mask = pd.Series(False, index=merged_df.index)
                        for col in source_cols:
                            series = merged_df[col].astype(str).str.strip().str.lower()
                            # Option 1: contains 'endo' (but not 'exo')
                            contains_endo = series.str.contains('endo', na=False) & ~series.str.contains('exo', na=False)
                            # Option 2: common truthy values (yes, y, true, 1)
                            is_truthy = series.isin(truthy)
                            col_mask = contains_endo | is_truthy
                            endogenous_mask |= col_mask
                    else:
                        endogenous_mask = merged_df[endogenous_col].astype(str).str.strip().str.lower().isin(truthy)
                    merged_df = merged_df[endogenous_mask].copy()
                    after_count = len(merged_df)
                    self._log(f"  ✅ Endogenous only: {before_count} → {after_count} features ({before_count - after_count} excluded)\n")
                    filters_applied.append('Endogenous (Source)')
                else:
                    # Warning - column not assigned
                    self._log("  ⚠️ Endogenous_Source column not found! Please assign it in Verify Columns.\n")
                    _debug_log_assignment_state('Endogenous_Source', merged_df)
                    self.root.after(0, lambda: messagebox.showwarning(
                        "Column Not Assigned",
                        "Endogenous filter is enabled but 'Endogenous_Source' column was not found.\n\n"
                        "Please assign it in 'Verify Columns' (type: Endogenous_Source) or disable the filter."
                    ))
            
            # --- Has HMDB ID Filter ---
            if hasattr(self, 'filter_has_hmdb_var') and self.filter_has_hmdb_var.get():
                self._log("🔬 Applying HMDB ID filter...\n")
                
                # First try to get from verified assignments
                hmdb_cols = get_assigned_columns('HMDB ID')
                hmdb_col = resolve_column_name(hmdb_cols[0], merged_df.columns) if hmdb_cols else None
                
                # Fallback: search for column name (case-insensitive)
                if not hmdb_col:
                    for col in merged_df.columns:
                        col_lower = col.lower().replace('_', '').replace(' ', '')
                        if col_lower in ['hmdbid', 'hmdb', 'hmdb_id']:
                            hmdb_col = col
                            break
                
                if hmdb_col and hmdb_col in merged_df.columns:
                    before_count = len(merged_df)
                    # Keep rows where HMDB ID is not empty/null/nan
                    hmdb_mask = merged_df[hmdb_col].astype(str).str.strip().replace('', pd.NA).replace('nan', pd.NA).replace('None', pd.NA).notna()
                    merged_df = merged_df[hmdb_mask].copy()
                    after_count = len(merged_df)
                    self._log(f"  ✅ Has HMDB ID: {before_count} → {after_count} features ({before_count - after_count} excluded)\n")
                    filters_applied.append('Has HMDB ID')
                else:
                    # Warning - column not assigned
                    self._log(f"  ⚠️ HMDB_ID column not found! Please assign it in Verify Columns.\n")
                    _debug_log_assignment_state('HMDB ID', merged_df)
                    self.root.after(0, lambda: messagebox.showwarning(
                        "Column Not Assigned",
                        "HMDB ID filter is enabled but 'HMDB_ID' column was not found.\n\n"
                        "Please assign it in 'Verify Columns' or disable the filter."
                    ))
            
            # --- P-Value Filter ---
            if hasattr(self, 'filter_pvalue_var') and self.filter_pvalue_var.get():
                if pairwise_only_precalc:
                    self._log("🔬 Pairwise-only mode detected: skipping global P-Value prefilter.\n")
                    self._log("  ℹ️ Pair-specific mapped P-Value columns will be enforced in pairwise stage.\n")
                else:
                    try:
                        pvalue_threshold = float(self.filter_pvalue_threshold_var.get().strip())
                    except (ValueError, TypeError):
                        pvalue_threshold = 0.05
                        self._log(f"  ⚠️ Invalid p-value threshold, using default 0.05\n")

                    self._log(f"🔬 Applying P-Value filter (< {pvalue_threshold})...\n")
                    pvalue_columns = []

                    # Prefer verified assignment(s): use ALL assigned p-value columns (row-wise OR)
                    pvalue_cols = get_assigned_columns('P-Value')
                    if pvalue_cols:
                        resolved = [resolve_column_name(c, merged_df.columns) for c in pvalue_cols]
                        resolved = [c for c in resolved if c is not None]
                        # De-duplicate while preserving order
                        resolved = list(dict.fromkeys(resolved))
                        pvalue_columns = resolved
                        if pvalue_columns:
                            preview = ', '.join(pvalue_columns[:5])
                            more = ' ...' if len(pvalue_columns) > 5 else ''
                            self._log(f"  ℹ️ Using {len(pvalue_columns)} assigned P-Value column(s) (OR across): {preview}{more}\n")

                    if pvalue_columns:
                        before_count = len(merged_df)
                        # Keep rows where ANY p-value column < threshold
                        combined_mask = pd.Series(False, index=merged_df.index)
                        for col in pvalue_columns:
                            pvals = pd.to_numeric(merged_df[col], errors='coerce')
                            combined_mask |= (pvals < pvalue_threshold)
                        merged_df = merged_df[combined_mask].copy()
                        after_count = len(merged_df)
                        self._log(
                            f"  ✅ Any P-Value < {pvalue_threshold}: {before_count} → {after_count} features ({before_count - after_count} excluded)\n"
                        )
                        filters_applied.append(f'Any P-Value<{pvalue_threshold}')
                    else:
                        self._log("  ❌ P-Value filter requires user-assigned P-Value column(s).\n")
                        self.root.after(0, lambda: messagebox.showerror(
                            "Column Not Assigned",
                            "P-Value filter is enabled but no assigned P-Value column is available.\n\n"
                            "Please assign P-Value column(s) in 'Verify Columns' or disable the filter."
                        ))
                        return
            
            # Summary of filters
            if filters_applied:
                self._log(f"\n✅ Filters applied: {', '.join(filters_applied)}\n")
                self._log(f"   Final feature count: {len(merged_df)}\n\n")
            
            # Check if any data remains
            if len(merged_df) == 0:
                self._log("❌ No features remain after applying filters\n")
                self.root.after(0, lambda: messagebox.showerror(
                    "No Data",
                    "No features remain after applying filters.\n"
                    "Try disabling some filters or adjusting thresholds."
                ))
                return
            
            # Store for export
            self.merged_data = merged_df
            
            # Step 2: Run ML analysis
            self._log("Step 2: Running ML analysis...\n")
            
            from main_script.ml_models import MetabolomicsMLAnalysis, format_classification_summary, format_pca_summary
            
            # Get feature ID column from verified columns
            feature_id_col = self._get_feature_id_column(merged_df)
            
            ml_analysis = MetabolomicsMLAnalysis(
                data_df=merged_df,
                group_assignments=group_map,
                feature_columns=feature_cols,
                feature_id_col=feature_id_col
            )
            
            analysis_type = 'classification'
            self.analysis_type_var.set('classification')
            
            if analysis_type == 'classification' or analysis_type == 'feature_importance':
                # Get selected models from checkboxes
                selected_models = [model for model, var in self.model_checkboxes.items() if var.get()]
                
                if not selected_models:
                    messagebox.showwarning("No Model Selected", "Please select at least one model.")
                    return
                
                # Auto-detect: single vs multi-model
                is_multi_model = len(selected_models) > 1
                
                cv_folds = int(self.cv_folds_var.get())
                test_size_str = self.test_size_var.get().strip()
                
                # Validate and fix test_size (0 = CV-only mode, no held-out test set)
                try:
                    test_size = float(test_size_str)
                    if test_size < 0.0 or test_size >= 1.0:
                        test_size = 0.3  # Default if invalid
                        self._log(f"⚠️ Invalid test size {test_size_str}, using default 0.3\n")
                    elif test_size == 0.0:
                        if is_multi_model:
                            messagebox.showerror("Invalid Test Size", 
                                "Multi-model comparison requires test_size > 0 for held-out test set.")
                            return
                        self._log("ℹ️ Test size = 0: Using CV-only mode (no held-out test set)\n")
                        self._log("   All samples used for cross-validation for more reliable estimates\n")
                except (ValueError, TypeError):
                    test_size = 0.3
                    self._log(f"⚠️ Invalid test size format, using default 0.3\n")
                
                scaling = self.scaling_var.get()
                class_weight_mode = self.class_weight_var.get().strip().lower()
                class_weight = 'balanced' if class_weight_mode == 'balanced' else None
                
                # Get regularization settings
                reg_type = self.regularization_type_var.get()
                reg_strength = self.regularization_strength_var.get()
                max_iter = int(self.max_iter_var.get())
                
                # Map regularization strength to C value
                c_value_map = {'strong': 0.1, 'medium': 1.0, 'weak': 10.0}
                c_value = c_value_map.get(reg_strength, 1.0)
                
                # Get robustness settings
                repeated_runs = int(self.repeated_runs_var.get())
                base_seed = int(self.base_seed_var.get())
                stability_tracking = self.stability_tracking_var.get()
                stability_threshold = float(self.stability_threshold_var.get())
                tune_hyperparameters = self.tune_hyperparameters_var.get()
                tuning_strategy = self.tuning_strategy_var.get().strip().lower()
                tuning_iter = int(self.tuning_iter_var.get())
                use_repeated_cv = self.use_repeated_cv_var.get()
                cv_repeats = int(self.cv_repeats_var.get())
                nested_cv = self.nested_cv_var.get()
                calibration_method = self.calibration_method_var.get().strip().lower()
                if calibration_method == 'none':
                    calibration_method = None
                permutation_test_runs = int(self.permutation_test_runs_var.get())
                imputation_method = self.imputation_method_var.get().strip().lower()
                imputation_knn_neighbors = int(self.imputation_knn_neighbors_var.get())
                auto_skip_scaling_for_trees = self.auto_skip_scaling_tree_var.get()
                feature_selection_method = self.feature_selection_method_var.get().strip().lower()
                variance_percentile = float(self.variance_percentile_var.get())
                univariate_k = int(self.univariate_k_var.get())
                lasso_c = float(self.lasso_c_var.get())
                rfe_n_features = int(self.rfe_n_features_var.get())

                pairwise_requested = bool(self.generate_pairwise_addon_var.get())
                unique_groups = sorted({str(v) for v in group_map.values() if str(v).strip()})
                pairwise_only_mode = pairwise_requested and len(unique_groups) >= 3

                if pairwise_only_mode and hasattr(self, 'filter_pvalue_var') and self.filter_pvalue_var.get():
                    pairwise_map = getattr(self, 'pairwise_pvalue_map', {}) or {}
                    missing_pairs = []
                    invalid_pairs = []
                    for g1, g2 in combinations(unique_groups, 2):
                        pair_key = self._pairwise_key(g1, g2)
                        selected_col = str(pairwise_map.get(pair_key, '')).strip()
                        if not selected_col:
                            missing_pairs.append(f"{g1} vs {g2}")
                        elif selected_col not in merged_df.columns:
                            invalid_pairs.append(f"{g1} vs {g2} -> {selected_col}")

                    if missing_pairs or invalid_pairs:
                        if missing_pairs:
                            self._log("❌ Pairwise p-value mapping is incomplete.\n")
                            self._log(f"   Missing mappings: {', '.join(missing_pairs)}\n")
                        if invalid_pairs:
                            self._log("❌ Some mapped pairwise p-value columns are not in the current data.\n")
                            self._log(f"   Invalid mappings: {', '.join(invalid_pairs)}\n")
                        self.root.after(0, lambda: messagebox.showerror(
                            "Pairwise P-Value Mapping Required",
                            "Pairwise mode + P-Value filter requires an explicit p-value column for every pair.\n\n"
                            "Open 'Verify Pairwise P-Value Columns' and map each pair to a valid column, then rerun."
                        ))
                        return
                
                if is_multi_model:
                    self._log(f"\n{'='*60}\n")
                    self._log("🔬 MULTI-MODEL COMPARISON (Equal Training Conditions)\n")
                    self._log(f"{'='*60}\n\n")
                    self._log(f"Models: {', '.join(selected_models)}\n")
                    self._log(f"Repeated runs: {repeated_runs} (reporting mean ± std across runs)\n")
                    self._log(f"Base seed: {base_seed} (seed, seed+1, seed+2, ...)\n")
                else:
                    model_name = selected_models[0]
                    self._log(f"Model: {model_name}\n")
                
                self._log(f"Test size: {test_size}\n")
                self._log(f"CV folds: {cv_folds}\n")
                self._log(f"Scaling: {scaling}\n")
                self._log(f"Class weight: {class_weight_mode}\n")
                self._log(f"Regularization: {reg_type.upper()} | C={c_value} ({reg_strength})\n")
                self._log(f"Max iterations: {max_iter}\n")
                if pairwise_only_mode:
                    self._log("Pairwise mode: ON (pairwise-only; main multiclass run skipped)\n")
                else:
                    self._log(f"Pairwise mode: {'OFF' if not pairwise_requested else 'ON (not applicable: <3 groups)'}\n")
                if not is_multi_model:
                    self._log(f"Repeated runs: {repeated_runs}\n")
                    self._log(f"Hyperparameter tuning: {'ON' if tune_hyperparameters else 'OFF'} ({tuning_strategy})\n")
                    self._log(f"Repeated CV: {'ON' if use_repeated_cv else 'OFF'} | repeats={cv_repeats}\n")
                    self._log(f"Nested CV: {'ON' if nested_cv else 'OFF'}\n")
                    self._log(f"Calibration: {calibration_method or 'none'}\n")
                    self._log(f"Permutation runs: {permutation_test_runs}\n")
                    self._log(f"Imputation: {imputation_method} (knn={imputation_knn_neighbors})\n")
                    self._log(f"Auto-skip scaling (trees): {'ON' if auto_skip_scaling_for_trees else 'OFF'}\n")
                    self._log(f"Feature selection: {feature_selection_method}\n")
                    if feature_selection_method == 'variance':
                        self._log(f"Variance percentile: {variance_percentile}\n")
                    elif feature_selection_method == 'univariate':
                        self._log(f"Univariate K: {univariate_k}\n")
                    elif feature_selection_method == 'lasso':
                        self._log(f"LASSO C: {lasso_c}\n")
                    elif feature_selection_method == 'rf_rfe':
                        self._log(f"RFE features: {rfe_n_features}\n")
                if class_weight == 'balanced':
                    class_counts = Counter(group_map.values())
                    total_n = sum(class_counts.values())
                    if total_n > 0:
                        self._log("Estimated balanced class weights (N / n_i):\n")
                        for group_name, count in sorted(class_counts.items()):
                            weight = total_n / max(count, 1)
                            self._log(f"  {group_name}: n={count}, weight~{weight:.2f}\n")
                if repeated_runs > 1 and not is_multi_model:
                    self._log(f"Base seed: {base_seed} (runs use: {base_seed}, {base_seed+1}, {base_seed+2}, ...)\n")
                    self._log(f"Stability tracking: {'ON' if stability_tracking else 'OFF'}\n")
                    if stability_tracking:
                        self._log(f"Stability threshold: {stability_threshold}%\n")
                self._log("\n")

                if pairwise_only_mode:
                    self.ml_results = None
                    self._log("\n🔀 Running pairwise-only analysis (main multiclass run is skipped)...\n")
                    pairwise_info = self._generate_pairwise_ml_addon_outputs(
                        merged_df=merged_df,
                        group_map=group_map,
                        feature_cols=feature_cols,
                        selected_models=selected_models,
                        is_multi_model=is_multi_model,
                        cv_folds=cv_folds,
                        test_size=test_size,
                        scaling=scaling,
                        reg_type=reg_type,
                        c_value=c_value,
                        max_iter=max_iter,
                        class_weight=class_weight,
                        repeated_runs=repeated_runs,
                        base_seed=base_seed,
                        tune_hyperparameters=tune_hyperparameters,
                        tuning_strategy=tuning_strategy,
                        tuning_iter=tuning_iter,
                        use_repeated_cv=use_repeated_cv,
                        cv_repeats=cv_repeats,
                        nested_cv=nested_cv,
                        calibration_method=calibration_method,
                        permutation_test_runs=permutation_test_runs,
                        imputation_method=imputation_method,
                        imputation_knn_neighbors=imputation_knn_neighbors,
                        auto_skip_scaling_for_trees=auto_skip_scaling_for_trees,
                        feature_selection_method=feature_selection_method,
                        variance_percentile=variance_percentile,
                        univariate_k=univariate_k,
                        lasso_c=lasso_c,
                        rfe_n_features=rfe_n_features,
                        stability_tracking=stability_tracking,
                        stability_threshold=stability_threshold,
                    )
                    self.pairwise_ml_results = pairwise_info
                    self._log(f"   Output folder: {pairwise_info.get('output_root', 'N/A')}\n")
                    self._log(f"   Excel file: {pairwise_info.get('excel_path', 'N/A')}\n")
                    if pairwise_info.get('auc_plot_path'):
                        self._log(f"   AUC comparison plot: {pairwise_info.get('auc_plot_path')}\n")
                else:
                    # Route to single or multi-model flow
                    if is_multi_model:
                        # Multi-model comparison
                        results = ml_analysis.run_multi_model_comparison(
                            model_names=selected_models,
                            test_size=test_size,
                            scaling_method=scaling,
                            regularization_type=reg_type,
                            C=c_value,
                            max_iter=max_iter,
                            class_weight=class_weight,
                            repeated_runs=repeated_runs,
                            random_state=base_seed
                        )
                    else:
                        # Single model with optional repeats
                        model_name = selected_models[0]
                        results = ml_analysis.run_classification(
                            model_name=model_name,
                            test_size=test_size,
                            cv_folds=cv_folds,
                            scaling_method=scaling,
                            regularization_type=reg_type,
                            C=c_value,
                            max_iter=max_iter,
                            class_weight=class_weight,
                            random_state=base_seed,
                            repeated_runs=repeated_runs,
                            tune_hyperparameters=tune_hyperparameters,
                            tuning_strategy=tuning_strategy,
                            tuning_iter=tuning_iter,
                            use_repeated_cv=use_repeated_cv,
                            cv_repeats=cv_repeats,
                            nested_cv=nested_cv,
                            calibration_method=calibration_method,
                            permutation_test_runs=permutation_test_runs,
                            imputation_method=imputation_method,
                            imputation_knn_neighbors=imputation_knn_neighbors,
                            auto_skip_scaling_for_trees=auto_skip_scaling_for_trees,
                            feature_selection_method=feature_selection_method,
                            variance_percentile=variance_percentile,
                            univariate_k=univariate_k,
                            lasso_C=lasso_c,
                            rfe_n_features=rfe_n_features,
                            stability_tracking=stability_tracking,
                            stability_threshold=stability_threshold,
                        )

                    self.ml_results = results

                    if results is None:
                        self._log("\n❌ All runs failed. No results to summarize.\n")
                        self.root.after(0, lambda: messagebox.showerror(
                            "Analysis Error",
                            "All runs failed. Please check your filters and data columns."
                        ))
                        return

                    # Check for overfitting and display warnings
                    self._check_overfitting_warnings(results, len(merged_df), len(feature_cols))

                    # Display results
                    if analysis_type == 'feature_importance':
                        # Show feature importance
                        if results.get('feature_importances'):
                            self._log(f"\n{'='*60}\n")
                            self._log("⭐ Feature Importance Analysis\n")
                            self._log(f"{'='*60}\n\n")

                            top_features = results['feature_importances']['top_features']
                            self._log(f"Top {len(top_features)} Features:\n\n")

                            for rank, (feature, importance) in enumerate(top_features, 1):
                                bar_length = int(importance * 40)
                                bar = "█" * bar_length + "░" * (40 - bar_length)
                                self._log(f"{rank:2d}. {feature[:35]:35s} {importance:.6f}\n")
                                self._log(f"    {bar}\n\n")
                        else:
                            summary = format_classification_summary(results)
                            self._log(summary)
                    else:
                        summary = format_classification_summary(results)
                        self._log(summary)

                    # Optional publication figure generation for classification outputs
                    if self.auto_generate_plots_var.get():
                        try:
                            base_dir = self.working_folder_var.get().strip()
                            if not base_dir:
                                raise ValueError("Working folder is not set")
                            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                            figures_root = os.path.join(base_dir, f"ML_Figures_{timestamp}")
                            figure_info = ml_analysis.generate_publication_figures(
                                results,
                                figures_root,
                                top_n_values=self.top_n_values,
                                figure_settings=self.figure_settings
                            )

                            if figure_info and isinstance(figure_info, dict):
                                self._log("\n📊 Auto-generated figures:\n")
                                self._log(f"   Output folder: {figure_info.get('output_root', figures_root)}\n")
                                combined_roc_path = figure_info.get('combined_roc_path')
                                if combined_roc_path:
                                    self._log(f"   Combined ROC: {combined_roc_path}\n")
                                model_dirs = figure_info.get('model_dirs', {})
                                if model_dirs:
                                    for model_name, model_dir in model_dirs.items():
                                        self._log(f"   - {model_name}: {model_dir}\n")
                            else:
                                self._log("\nℹ️ Figure generation skipped (insufficient data for plots).\n")
                        except Exception as fig_err:
                            logger.warning(f"Figure generation failed: {fig_err}", exc_info=True)
                            self._log(f"\n⚠️ Figure generation failed: {fig_err}\n")
                    
            elif analysis_type == 'pca':
                try:
                    n_components = 10  # Default
                    scaling = self.scaling_var.get()
                    
                    self._log(f"Components: {n_components}\n")
                    self._log(f"Scaling: {scaling}\n\n")
                    
                    results = ml_analysis.run_pca(
                        n_components=n_components,
                        scaling_method=scaling
                    )
                    
                    self.ml_results = results
                    
                    # Safe extraction of PCA results
                    if isinstance(results, dict):
                        self._log(f"\n✅ PCA Complete!\n")
                        self._log(f"Components extracted: {results.get('n_components', 'N/A')}\n")
                        if 'explained_variance_ratio' in results:
                            self._log(f"Explained Variance Ratio:\n")
                            for i, var in enumerate(results['explained_variance_ratio'][:10]):
                                self._log(f"  PC{i+1}: {var:.4f}\n")
                    else:
                        self._log(f"PCA completed successfully\n")
                except Exception as pca_err:
                    self._log(f"\n⚠️ PCA error: {str(pca_err)}\n")
                    self._log(f"Check data dimensions and ensure sufficient samples per group\n")
                    
            elif analysis_type == 'lda':
                try:
                    scaling = self.scaling_var.get()
                    
                    self._log(f"Scaling: {scaling}\n\n")
                    
                    results = ml_analysis.run_lda(scaling_method=scaling)
                    self.ml_results = results
                    
                    # Safe extraction of LDA results
                    if isinstance(results, dict):
                        self._log(f"\n✅ LDA Complete!\n")
                        self._log(f"Discriminants: {results.get('n_components', 'N/A')}\n")
                        if 'explained_variance' in results and isinstance(results['explained_variance'], (list, tuple)):
                            self._log(f"Explained Variance:\n")
                            for i, var in enumerate(results['explained_variance']):
                                self._log(f"  LD{i+1}: {var:.4f}\n")
                    else:
                        self._log(f"LDA completed successfully\n")
                except Exception as lda_err:
                    self._log(f"\n⚠️ LDA error: {str(lda_err)}\n")
                    self._log(f"Check that you have at least 2 groups and sufficient samples\n")
            
            self._log(f"\n{'='*60}\n")
            self._log("✅ ML Analysis Complete!\n")
            self._log(f"{'='*60}\n\n")
            self._log(f"⏱️ Total analysis time: {perf_counter() - analysis_start:.1f}s\n\n")

            auto_export_path = None
            if self.ml_results is not None:
                auto_export_path = self._auto_export_main_results(ml_analysis, self.ml_results)
                if auto_export_path:
                    self._log(f"📥 Main ML results auto-exported: {auto_export_path}\n")

            # Explicit completion notice for the UI and log
            pairwise_only_done = self.ml_results is None and self.pairwise_ml_results is not None
            if pairwise_only_done:
                self._log("🎉 Pairwise-only analysis finished successfully. Outputs were saved to the pairwise folder/Excel.\n\n")
            else:
                self._log("🎉 Analysis finished successfully. Results were auto-exported to your Working Folder.\n\n")
            if hasattr(self, 'data_status_label'):
                self.root.after(0, lambda: self.data_status_label.config(text='✅ ML analysis complete', fg='green'))

            if self.ml_results is not None:
                self.root.after(0, lambda: messagebox.showinfo(
                    "Analysis Complete",
                    "ML analysis finished successfully.\n\n"
                    f"Main results were auto-exported to:\n{auto_export_path or 'Working Folder'}"
                ))
            else:
                self.root.after(0, lambda: messagebox.showinfo(
                    "Pairwise Analysis Complete",
                    "Pairwise-only ML analysis finished successfully.\n\n"
                    "Outputs were auto-exported to the pairwise folder in your Working Folder."
                ))
            
        except Exception as e:
            logger.error(f"ML analysis error: {e}", exc_info=True)
            error_msg = str(e)
            self._log(f"\n❌ Analysis error: {error_msg}\n")
            self.root.after(0, lambda msg=error_msg: messagebox.showerror("Analysis Error", msg))
            
        finally:
            self.root.after(0, lambda: self.run_button.config(state='normal', text='▶️ Run ML Analysis'))

    def _safe_path_token(self, value: Any) -> str:
        """Convert free text into a filesystem-safe token."""
        token = re.sub(r'[^A-Za-z0-9._-]+', '_', str(value).strip())
        token = re.sub(r'_+', '_', token).strip('._')
        return token or 'group'

    def _safe_sheet_name(self, value: str) -> str:
        """Convert arbitrary text into an Excel-safe sheet name."""
        cleaned = re.sub(r'[\\/*?:\[\]]+', '_', str(value))
        cleaned = cleaned.strip()
        return cleaned[:31] or 'Sheet'

    def _auto_export_main_results(self, ml_analysis, results: Dict[str, Any]) -> Optional[str]:
        """Automatically export main ML results to the working folder."""
        try:
            if results is None:
                return None

            base_dir = self.working_folder_var.get().strip()
            if not base_dir:
                raise ValueError("Working folder is not set")

            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            export_root = os.path.join(base_dir, f"ML_Results_{timestamp}")
            os.makedirs(export_root, exist_ok=True)

            final_path = os.path.join(export_root, 'ml_results.xlsx')
            temp_path = os.path.join(export_root, 'ml_results_core.xlsx')

            ml_analysis.results = results
            ml_analysis.export_results_to_excel(temp_path, 'classification')

            temp_wb = pd.ExcelFile(temp_path)
            with pd.ExcelWriter(final_path, engine='openpyxl') as writer:
                for sheet_name in temp_wb.sheet_names:
                    sheet_name_str = str(sheet_name)
                    df = pd.read_excel(temp_wb, sheet_name=sheet_name_str)
                    df.to_excel(writer, sheet_name=sheet_name_str, index=False)
                temp_wb.close()

                if self.merged_data is not None:
                    self.merged_data.to_excel(writer, sheet_name='Filtered_Dataset', index=False)

                settings_data = {
                    'Parameter': [
                        'Analysis Date',
                        'Model',
                        'Test Size',
                        'CV Folds',
                        'Scaling Method',
                        'Class Weight',
                        'Regularization Type',
                        'Regularization Strength',
                        'Max Iterations',
                        'Repeated Runs',
                        'Base Seed',
                        'Seed Strategy',
                        'Imputation Method',
                        'Imputation KNN Neighbors',
                        'Auto Skip Scaling (Trees)',
                        'Feature Selection Method',
                        'Features Before Selection',
                        'Features After Selection',
                        'Stability Tracking',
                        'Stability Threshold',
                        'Feature Filter',
                    ],
                    'Value': [
                        pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
                        (", ".join(results.get('models_trained', [])) if results.get('comparison_type') == 'multi_model' else results.get('model_name', 'N/A')),
                        results.get('test_size', 'N/A'),
                        results.get('cv_folds', 'N/A'),
                        results.get('scaling_method', 'N/A'),
                        results.get('class_weight', 'none'),
                        self.analysis_type_var.get(),
                        self.regularization_strength_var.get(),
                        self.max_iter_var.get(),
                        results.get('n_repeats', '1'),
                        self.base_seed_var.get(),
                        "Fixed base seed with increments (seed, seed+1, seed+2, ...)",
                        results.get('imputation_method', 'N/A'),
                        results.get('imputation_knn_neighbors', 'N/A'),
                        'ON' if results.get('auto_skip_scaling_for_trees') else 'OFF',
                        results.get('feature_selection_method', 'none'),
                        results.get('n_features_before_selection', results.get('n_features', 'N/A')),
                        results.get('n_features_after_selection', results.get('n_features', 'N/A')),
                        'ON' if self.stability_tracking_var.get() else 'OFF',
                        f"{self.stability_threshold_var.get()}%",
                        results.get('feature_filter', 'None'),
                    ]
                }
                pd.DataFrame(settings_data).to_excel(writer, sheet_name='Analysis_Settings', index=False)

            if os.path.exists(temp_path):
                os.remove(temp_path)
            return final_path

        except Exception as e:
            logger.error(f"Auto export error: {e}", exc_info=True)
            self._log(f"⚠️ Auto-export failed: {e}\n")
            return None

    def _generate_pairwise_ml_addon_outputs(
        self,
        merged_df,
        group_map,
        feature_cols,
        selected_models,
        is_multi_model,
        cv_folds,
        test_size,
        scaling,
        reg_type,
        c_value,
        max_iter,
        class_weight,
        repeated_runs,
        base_seed,
        tune_hyperparameters,
        tuning_strategy,
        tuning_iter,
        use_repeated_cv,
        cv_repeats,
        nested_cv,
        calibration_method,
        permutation_test_runs,
        imputation_method,
        imputation_knn_neighbors,
        auto_skip_scaling_for_trees,
        feature_selection_method,
        variance_percentile,
        univariate_k,
        lasso_c,
        rfe_n_features,
        stability_tracking,
        stability_threshold,
    ):
        """Generate pairwise (A vs B) ML outputs for multiclass data."""
        from main_script.ml_models import MetabolomicsMLAnalysis, format_classification_summary
        from time import perf_counter

        base_dir = self.working_folder_var.get().strip()
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        addon_root = os.path.join(base_dir, f"Pairwise_ML_{timestamp}")
        os.makedirs(addon_root, exist_ok=True)
        addon_start = perf_counter()

        def _norm_token(text: Any) -> str:
            return re.sub(r'[^a-z0-9]+', '', str(text).lower())

        def _is_pvalue_like(col_name: str) -> bool:
            n = _norm_token(col_name)
            keys = ('pvalue', 'pval', 'padj', 'adjp', 'adjustedp', 'fdr', 'qvalue')
            return any(k in n for k in keys)

        def _extract_verified_pvalue_columns() -> List[str]:
            cols = []
            target_names = {'pvalue', 'p-value', 'p value'}

            def _collect(mapping: Dict[str, Any]):
                if not mapping:
                    return
                for key, value in mapping.items():
                    key_norm = str(key).strip().lower()
                    value_norm = str(value).strip().lower()
                    # Support both {type: column} and {column: type} formats.
                    if key_norm in target_names and isinstance(value, str) and value.strip():
                        cols.append(value)
                    if value_norm in target_names and isinstance(key, str) and key.strip():
                        cols.append(key)

            _collect(getattr(self, 'verified_pos_assignments', {}) or {})
            _collect(getattr(self, 'verified_neg_assignments', {}) or {})

            unique = []
            seen = set()
            for c in cols:
                if c not in seen:
                    unique.append(c)
                    seen.add(c)
            return unique

        def _resolve_pairwise_pvalue_columns(g1: str, g2: str, all_columns: List[str], verified_cols: List[str]) -> List[str]:
            manual_map = getattr(self, 'pairwise_pvalue_map', {}) or {}
            pair_key = self._pairwise_key(g1, g2)

            manual_col = manual_map.get(pair_key)
            if manual_col and manual_col in all_columns:
                return [manual_col]
            return []

        unique_groups = sorted({str(v) for v in group_map.values() if str(v).strip()})
        feature_id_col = self._get_feature_id_column(merged_df)
        verified_pvalue_cols = _extract_verified_pvalue_columns()

        pvalue_enabled = bool(getattr(self, 'filter_pvalue_var', None) and self.filter_pvalue_var.get())
        try:
            pvalue_threshold = float(self.filter_pvalue_threshold_var.get().strip())
        except Exception:
            pvalue_threshold = 0.05

        pairwise_runs = []
        summary_rows = []

        for idx, (g1, g2) in enumerate(combinations(unique_groups, 2)):
            pair_label = f"{g1} vs {g2}"
            pair_token = f"{self._safe_path_token(g1)}_vs_{self._safe_path_token(g2)}"
            pair_seed = int(base_seed) + idx
            pair_start = perf_counter()
            pvalue_cols_used = []

            pair_feature_cols = [c for c in feature_cols if str(group_map.get(c, '')).strip() in {g1, g2}]
            pair_group_map = {c: group_map[c] for c in pair_feature_cols if c in group_map}

            if len(pair_feature_cols) < 4 or len(set(pair_group_map.values())) < 2:
                msg = "Insufficient samples for pairwise modeling"
                pairwise_runs.append({'pair': pair_label, 'status': 'failed', 'error': msg})
                summary_rows.append({'Pair': pair_label, 'Status': 'Failed', 'Reason': msg})
                self._log(f"   ⚠️ {pair_label}: {msg}\n")
                continue

            pair_df = merged_df
            if pvalue_enabled:
                pvalue_cols_used = _resolve_pairwise_pvalue_columns(g1, g2, list(merged_df.columns), verified_pvalue_cols)
                if pvalue_cols_used:
                    before_n = len(pair_df)
                    mask = pd.Series(False, index=pair_df.index)
                    for col in pvalue_cols_used:
                        mask |= (pd.to_numeric(pair_df[col], errors='coerce') < pvalue_threshold)
                    pair_df = pair_df[mask].copy()
                    after_n = len(pair_df)
                    self._log(
                        f"   🔎 {pair_label}: pairwise p-value filter (<{pvalue_threshold}) "
                        f"using {len(pvalue_cols_used)} column(s) -> {before_n} to {after_n} features\n"
                    )
                else:
                    msg = "Missing mapped pairwise p-value column"
                    pairwise_runs.append({'pair': pair_label, 'status': 'failed', 'error': msg})
                    summary_rows.append({'Pair': pair_label, 'Status': 'Failed', 'Reason': msg})
                    self._log(f"   ❌ {pair_label}: {msg}\n")
                    continue

            if len(pair_df) == 0:
                msg = "No features remain after pairwise p-value filtering"
                pairwise_runs.append({'pair': pair_label, 'status': 'failed', 'error': msg})
                summary_rows.append({'Pair': pair_label, 'Status': 'Failed', 'Reason': msg})
                self._log(f"   ⚠️ {pair_label}: {msg}\n")
                continue

            pair_analysis = MetabolomicsMLAnalysis(
                data_df=pair_df,
                group_assignments=pair_group_map,
                feature_columns=pair_feature_cols,
                feature_id_col=feature_id_col,
            )

            try:
                model_start = perf_counter()
                if is_multi_model:
                    pair_results = pair_analysis.run_multi_model_comparison(
                        model_names=selected_models,
                        test_size=test_size,
                        scaling_method=scaling,
                        regularization_type=reg_type,
                        C=c_value,
                        max_iter=max_iter,
                        class_weight=class_weight,
                        repeated_runs=repeated_runs,
                        random_state=pair_seed,
                    )
                else:
                    pair_results = pair_analysis.run_classification(
                        model_name=selected_models[0],
                        test_size=test_size,
                        cv_folds=cv_folds,
                        scaling_method=scaling,
                        regularization_type=reg_type,
                        C=c_value,
                        max_iter=max_iter,
                        class_weight=class_weight,
                        random_state=pair_seed,
                        repeated_runs=repeated_runs,
                        tune_hyperparameters=tune_hyperparameters,
                        tuning_strategy=tuning_strategy,
                        tuning_iter=tuning_iter,
                        use_repeated_cv=use_repeated_cv,
                        cv_repeats=cv_repeats,
                        nested_cv=nested_cv,
                        calibration_method=calibration_method,
                        permutation_test_runs=permutation_test_runs,
                        imputation_method=imputation_method,
                        imputation_knn_neighbors=imputation_knn_neighbors,
                        auto_skip_scaling_for_trees=auto_skip_scaling_for_trees,
                        feature_selection_method=feature_selection_method,
                        variance_percentile=variance_percentile,
                        univariate_k=univariate_k,
                        lasso_C=lasso_c,
                        rfe_n_features=rfe_n_features,
                        stability_tracking=stability_tracking,
                        stability_threshold=stability_threshold,
                    )
                model_elapsed = perf_counter() - model_start
                self._log(f"   ⏱️ {pair_label}: model analysis completed in {model_elapsed:.1f}s\n")

                figure_start = perf_counter()
                pair_folder = os.path.join(addon_root, f"Pair_{pair_token}")
                os.makedirs(pair_folder, exist_ok=True)

                figure_info = pair_analysis.generate_publication_figures(
                    pair_results,
                    pair_folder,
                    top_n_values=self.top_n_values,
                    figure_settings=self.figure_settings,
                )
                figure_elapsed = perf_counter() - figure_start
                self._log(f"   ⏱️ {pair_label}: figure generation completed in {figure_elapsed:.1f}s\n")
                if isinstance(figure_info, dict) and figure_info.get('combined_roc_path'):
                    self._log(f"   📈 {pair_label}: combined ROC -> {figure_info.get('combined_roc_path')}\n")

                excel_start = perf_counter()
                summary_rows.append({
                    'Pair': pair_label,
                    'Status': 'OK',
                    'Reason': '',
                    'Features After Pair Filter': len(pair_df),
                    'Pairwise P-Value Columns': ', '.join(pvalue_cols_used),
                    'Classes': ', '.join(pair_results.get('class_labels', [])) if isinstance(pair_results, dict) else '',
                    'Model': ', '.join(pair_results.get('models_trained', [])) if pair_results.get('comparison_type') == 'multi_model' else pair_results.get('model_name', ''),
                    'Accuracy': pair_results.get('test_accuracy', None),
                    'AUC': pair_results.get('auc', None),
                    'Output Folder': figure_info.get('output_root', pair_folder) if isinstance(figure_info, dict) else pair_folder,
                })

                pairwise_runs.append({
                    'pair': pair_label,
                    'status': 'ok',
                    'results': pair_results,
                    'summary': format_classification_summary(pair_results),
                    'output_folder': figure_info.get('output_root', pair_folder) if isinstance(figure_info, dict) else pair_folder,
                    'pvalue_columns': pvalue_cols_used,
                    'n_features_after_pair_filter': len(pair_df),
                })

                pair_elapsed = perf_counter() - pair_start
                self._log(f"   ✅ {pair_label}: generated in {pair_folder} (total {pair_elapsed:.1f}s)\n")

            except Exception as e:
                pairwise_runs.append({'pair': pair_label, 'status': 'failed', 'error': str(e)})
                summary_rows.append({'Pair': pair_label, 'Status': 'Failed', 'Reason': str(e)})
                self._log(f"   ❌ {pair_label}: {e}\n")

        # Create cross-pair metric comparison bar graphs.
        auc_rows = []
        acc_rows = []
        for entry in pairwise_runs:
            if entry.get('status') != 'ok':
                continue
            pair_name = str(entry.get('pair', 'pair'))
            res = entry.get('results', {}) or {}
            if res.get('comparison_type') == 'multi_model':
                for model_name in res.get('models_trained', []):
                    model_res = (res.get('model_results', {}) or {}).get(model_name, {}) or {}
                    auc_val = model_res.get('auc_mean')
                    if auc_val is None:
                        auc_val = model_res.get('auc')
                    acc_val = model_res.get('test_accuracy_mean')
                    if acc_val is None:
                        acc_val = model_res.get('test_accuracy')
                    if auc_val is not None:
                        auc_rows.append({
                            'Pair': pair_name,
                            'Model': model_name,
                            'Mean': float(auc_val),
                            'Std': float(model_res.get('auc_std') or 0.0),
                        })
                    if acc_val is not None:
                        acc_rows.append({
                            'Pair': pair_name,
                            'Model': model_name,
                            'Mean': float(acc_val),
                            'Std': float(model_res.get('test_accuracy_std') or 0.0),
                        })
            else:
                auc_val = res.get('auc')
                acc_val = res.get('test_accuracy')
                if auc_val is not None:
                    auc_rows.append({
                        'Pair': pair_name,
                        'Model': res.get('model_name', 'Model'),
                        'Mean': float(auc_val),
                        'Std': float(res.get('auc_std') or 0.0),
                    })
                if acc_val is not None:
                    acc_rows.append({
                        'Pair': pair_name,
                        'Model': res.get('model_name', 'Model'),
                        'Mean': float(acc_val),
                        'Std': float(res.get('test_accuracy_std') or 0.0),
                    })

        auc_plot_path = None
        acc_plot_path = None

        def _plot_pair_metric(metric_df: pd.DataFrame, metric_label: str, filename: str) -> Optional[str]:
            if metric_df.empty:
                return None
            try:
                import matplotlib.pyplot as plt

                comparison_value_fs = float(self.figure_settings.get('comparison_value_fs', 12.0))

                def _safe_float(value: Any) -> Optional[float]:
                    try:
                        return float(value)
                    except (TypeError, ValueError):
                        return None

                metric_pivot = metric_df.pivot_table(index='Pair', columns='Model', values='Mean', aggfunc='mean')
                std_pivot = metric_df.pivot_table(index='Pair', columns='Model', values='Std', aggfunc='mean').reindex_like(metric_pivot)
                fig_w = max(8.5, 2.4 + 1.7 * len(metric_pivot.index))
                fig_h = 6.0
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                metric_pivot.plot(kind='bar', ax=ax, width=0.8)
                ax.set_ylim(0.0, 1.0)
                ax.set_ylabel(metric_label, fontweight='bold')
                ax.set_xlabel('Pairwise Comparison', fontweight='bold')
                ax.set_title(f'Pairwise {metric_label} Comparison', fontweight='bold')
                ax.grid(axis='y', alpha=0.25)
                ax.tick_params(axis='x', rotation=20)
                for lbl in ax.get_xticklabels() + ax.get_yticklabels():
                    lbl.set_fontweight('bold')
                leg = ax.legend(title='Model', framealpha=0.95)
                if leg:
                    for txt in leg.get_texts():
                        txt.set_fontweight('bold')

                pair_index = list(metric_pivot.index)
                top_candidates = []
                for pair_idx, pair_name in enumerate(pair_index):
                    if pair_name not in std_pivot.index:
                        continue
                    for model_idx, model_name in enumerate(metric_pivot.columns):
                        value = _safe_float(metric_pivot.loc[pair_name, model_name])
                        if value is None:
                            continue
                        std_val = _safe_float(std_pivot.loc[pair_name, model_name]) or 0.0
                        bar = None
                        for container in ax.containers:
                            if model_idx < len(container):
                                bar = container[model_idx]
                                break
                        if bar is None:
                            continue
                        x = bar.get_x() + bar.get_width() / 2
                        top = value + std_val
                        top_candidates.append(top)
                        ax.errorbar(x, value, yerr=std_val, fmt='none', ecolor='black', elinewidth=1.1, capsize=3, capthick=1.1, zorder=4)
                        ax.text(
                            x,
                            top + 0.03,
                            f"{value:.3f}±{std_val:.3f}",
                            ha='center',
                            va='bottom',
                            fontsize=comparison_value_fs,
                            fontweight='bold',
                            rotation=0,
                            clip_on=False,
                        )

                if top_candidates:
                    ax.set_ylim(0.0, max(1.0, max(top_candidates) + 0.12))

                ax.margins(y=0.12)
                fig.tight_layout()
                out_path = os.path.join(addon_root, filename)
                fig.savefig(out_path, dpi=300, bbox_inches='tight')
                plt.close(fig)
                self._log(f"\n📊 Pairwise {metric_label.lower()} comparison plot saved: {out_path}\n")
                return out_path
            except Exception as plot_err:
                self._log(f"\n⚠️ Could not generate pairwise {metric_label.lower()} comparison plot: {plot_err}\n")
                return None

        auc_plot_path = _plot_pair_metric(pd.DataFrame(auc_rows), 'AUC', 'pairwise_auc_comparison.png')
        acc_plot_path = _plot_pair_metric(pd.DataFrame(acc_rows), 'Accuracy', 'pairwise_accuracy_comparison.png')

        excel_write_start = perf_counter()
        excel_path = os.path.join(addon_root, 'pairwise_ml_results.xlsx')
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            summary_export_rows = []
            for entry in pairwise_runs:
                if entry.get('status') != 'ok':
                    summary_export_rows.append({
                        'Pair': entry.get('pair', 'pair'),
                        'Status': 'Failed',
                        'Reason': entry.get('error', 'Unknown error'),
                        'Features After Pair Filter': entry.get('n_features_after_pair_filter', ''),
                        'Accuracy': '',
                        'AUC': '',
                        'Output Folder': entry.get('output_folder', ''),
                    })
                    continue

                res = entry.get('results', {}) or {}
                summary_export_rows.append({
                    'Pair': entry.get('pair', 'pair'),
                    'Status': 'OK',
                    'Reason': '',
                    'Features After Pair Filter': entry.get('n_features_after_pair_filter', ''),
                    'Accuracy': res.get('test_accuracy', None),
                    'Accuracy Std': res.get('test_accuracy_std', None),
                    'AUC': res.get('auc', None),
                    'AUC Std': res.get('auc_std', None),
                    'Output Folder': entry.get('output_folder', ''),
                })

            pd.DataFrame(summary_export_rows).to_excel(writer, sheet_name='Summary', index=False)

            for entry in pairwise_runs:
                if entry.get('status') != 'ok':
                    continue

                pair_name = str(entry.get('pair', 'pair'))

                summary_text = str(entry.get('summary', '')).strip()
                if summary_text:
                    text_df = pd.DataFrame({'Summary': summary_text.splitlines()})
                    text_df.to_excel(writer, sheet_name=self._safe_sheet_name(f"T_{pair_name}"), index=False)

        excel_elapsed = perf_counter() - excel_write_start
        total_elapsed = perf_counter() - addon_start

        self._log(f"⏱️ Pairwise add-on Excel export time: {excel_elapsed:.1f}s\n")
        self._log(f"\n⏱️ Pairwise add-on total time: {total_elapsed:.1f}s\n")

        return {
            'output_root': addon_root,
            'excel_path': excel_path,
            'auc_plot_path': auc_plot_path,
            'accuracy_plot_path': acc_plot_path,
            'pairs': pairwise_runs,
            'elapsed_seconds': total_elapsed,
        }
    
    def _run_repeated_classification(self, ml_analysis, model_name, test_size, cv_folds, scaling,
                                     reg_type, c_value, max_iter, n_repeats, base_seed, 
                                     stability_tracking, stability_threshold,
                                     class_weight):
        """Run classification multiple times and aggregate results"""
        from collections import Counter
        
        self._log(f"🔁 Running {n_repeats} repeated analyses for robustness...\n\n")
        
        all_results = []
        feature_selections = []  # Track which features are selected in each run
        for run_idx in range(n_repeats):
            # Determine random seed for this run
            # Always increment seed to ensure different data splits for true robustness testing
            random_seed = base_seed + run_idx
            
            self._log(f"[Run {run_idx + 1}/{n_repeats}] seed={random_seed}...\n")
            
            try:
                result = ml_analysis.run_classification(
                    model_name=model_name,
                    test_size=test_size,
                    cv_folds=cv_folds,
                    scaling_method=scaling,
                    regularization_type=reg_type,
                    C=c_value,
                    max_iter=max_iter,
                    random_state=random_seed,
                    class_weight=class_weight
                )
                
                all_results.append(result)
                
                # Track feature importance/selection if available
                if stability_tracking and result.get('feature_importances'):
                    top_features = [f[0] for f in result['feature_importances']['top_features'][:20]]
                    feature_selections.append(top_features)
                
                # Log results (handle CV-only mode where test_accuracy is None)
                if result.get('cv_only_mode') or result['test_accuracy'] is None:
                    self._log(f"  CV: {result['cv_mean_accuracy']:.4f} ± {result['cv_std_accuracy']:.4f}\n")
                else:
                    self._log(f"  CV: {result['cv_mean_accuracy']:.4f} | Test: {result['test_accuracy']:.4f}\n")
                
            except Exception as e:
                self._log(f"  ⚠️ Run failed: {str(e)}\n")
        
        # Aggregate results
        self._log(f"\n{'='*60}\n")
        self._log("📊 ROBUSTNESS ANALYSIS RESULTS\n")
        self._log(f"{'='*60}\n\n")
        
        if not all_results:
            self._log("❌ All runs failed\n")
            return None
        
        # Calculate mean and std for metrics
        cv_means = [r['cv_mean_accuracy'] for r in all_results]
        cv_stds = [r['cv_std_accuracy'] for r in all_results]
        train_accs = [r['train_accuracy'] for r in all_results]
        
        # Check if CV-only mode (test_accuracy is None)
        cv_only_mode = all_results[0].get('cv_only_mode', False) or all_results[0]['test_accuracy'] is None
        test_accs = [r['test_accuracy'] for r in all_results if r['test_accuracy'] is not None]
        
        self._log(f"Performance Summary (n={len(all_results)} runs):\n")
        self._log(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
        self._log(f"  CV Accuracy:      {np.mean(cv_means):.4f} ± {np.std(cv_means):.4f}\n")
        if cv_only_mode:
            self._log(f"  Test Accuracy:    N/A (CV-only mode - all samples in CV)\n")
        else:
            self._log(f"  Test Accuracy:    {np.mean(test_accs):.4f} ± {np.std(test_accs):.4f}\n")
        self._log(f"  Train Accuracy:   {np.mean(train_accs):.4f} ± {np.std(train_accs):.4f}\n\n")
        
        # Feature stability analysis
        if stability_tracking and feature_selections:
            self._log(f"\n⭐ Feature Stability Analysis\n")
            self._log(f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
            
            # Count feature occurrences
            all_features = [f for run_features in feature_selections for f in run_features]
            feature_counts = Counter(all_features)
            
            # Calculate frequency percentage
            total_runs = len(feature_selections)
            stable_features = []
            
            self._log(f"Feature selection frequency (appears in ≥{stability_threshold}% of runs):\n\n")
            
            for feature, count in feature_counts.most_common(30):
                frequency_pct = (count / total_runs) * 100
                if frequency_pct >= stability_threshold:
                    stable_features.append((feature, frequency_pct))
                    self._log(f"  {feature:50s} {count:2d}/{total_runs} ({frequency_pct:5.1f}%)\n")
            
            if stable_features:
                self._log(f"\n✅ {len(stable_features)} stable features identified (≥{stability_threshold}%)\n")
            else:
                self._log(f"\n⚠️ No features met {stability_threshold}% threshold\n")
                self._log(f"   Consider: lower threshold, more runs, or stronger regularization\n")
        
        # Return the best performing run result with aggregated stats
        if cv_only_mode or not test_accs:
            best_idx = int(np.argmax(cv_means))
        else:
            best_idx = int(np.argmax(test_accs))
        best_result = all_results[best_idx].copy()
        
        # Add robustness statistics
        best_result['robustness_stats'] = {
            'n_runs': len(all_results),
            'cv_mean': np.mean(cv_means),
            'cv_std': np.std(cv_means),
            'test_mean': np.mean(test_accs) if test_accs else None,
            'test_std': np.std(test_accs) if test_accs else None,
            'train_mean': np.mean(train_accs),
            'train_std': np.std(train_accs)
        }
        
        if stability_tracking and feature_selections:
            best_result['stable_features'] = stable_features
        
        return best_result
    
    def _test_models(self):
        """Open dialog to test multiple models with different parameters"""
        dialog = tk.Toplevel(self.root)
        dialog.title("🔬 Test Multiple Models")
        dialog.geometry("650x750")
        dialog.resizable(True, True)
        dialog.minsize(600, 600)
        
        # Create a canvas with scrollbar for the entire dialog
        canvas = tk.Canvas(dialog, highlightthickness=0)
        scrollbar = ttk.Scrollbar(dialog, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pack canvas and scrollbar
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Mouse wheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        # Clean up mouse wheel binding when dialog closes
        def on_close():
            canvas.unbind_all("<MouseWheel>")
            dialog.destroy()
        dialog.protocol("WM_DELETE_WINDOW", on_close)
        
        # Main frame inside scrollable area
        main_frame = tk.Frame(scrollable_frame)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Title
        title_label = tk.Label(main_frame, text="Model Comparison & Hyperparameter Testing", 
                              font=('Arial', 12, 'bold'))
        title_label.pack(pady=(0, 10))
        
        # Instructions
        info_text = "Select models and parameter ranges to test. The system will evaluate all combinations and recommend the best setup."
        info_label = tk.Label(main_frame, text=info_text, wraplength=550, justify='left', fg='#555')
        info_label.pack(pady=(0, 10))
        
        # Model selection
        model_frame = tk.LabelFrame(main_frame, text="Select Models to Test", font=('Arial', 10, 'bold'))
        model_frame.pack(fill='x', pady=5)
        
        model_vars = {}
        models = ['Random Forest', 'Gradient Boosting', 'SVM (RBF)', 'Logistic Regression', 'Linear Discriminant Analysis']
        for model in models:
            var = tk.BooleanVar(value=True)
            model_vars[model] = var
            tk.Checkbutton(model_frame, text=model, variable=var, font=('Arial', 9)).pack(anchor='w', padx=10, pady=2)
        
        # Parameter options
        param_frame = tk.LabelFrame(main_frame, text="Parameter Ranges to Test", font=('Arial', 10, 'bold'))
        param_frame.pack(fill='x', pady=5)
        
        # Test size options
        tk.Label(param_frame, text="Test Size:", font=('Arial', 9, 'bold')).pack(anchor='w', padx=10, pady=(5,2))
        test_size_vars = {}
        for size in ['0.2', '0.25', '0.3', '0.35']:
            var = tk.BooleanVar(value=(size == '0.3'))
            test_size_vars[size] = var
            tk.Checkbutton(param_frame, text=size, variable=var, font=('Arial', 9)).pack(anchor='w', padx=20, pady=1)
        
        # Scaling options
        tk.Label(param_frame, text="Scaling Method:", font=('Arial', 9, 'bold')).pack(anchor='w', padx=10, pady=(10,2))
        scaling_vars = {}
        for method in ['standard', 'robust', 'none']:
            var = tk.BooleanVar(value=(method == 'standard'))
            scaling_vars[method] = var
            tk.Checkbutton(param_frame, text=method.title(), variable=var, font=('Arial', 9)).pack(anchor='w', padx=20, pady=1)
        
        # CV folds
        cv_frame = tk.Frame(param_frame)
        cv_frame.pack(fill='x', padx=10, pady=(10,5))
        tk.Label(cv_frame, text="CV Folds:", font=('Arial', 9, 'bold')).pack(side='left')
        cv_folds_var = tk.StringVar(value='5')
        ttk.Spinbox(cv_frame, from_=3, to=10, increment=1, textvariable=cv_folds_var, width=8).pack(side='left', padx=10)
        
        # Results display options
        options_frame = tk.LabelFrame(main_frame, text="Output Options", font=('Arial', 10, 'bold'))
        options_frame.pack(fill='x', pady=5)
        
        show_all_var = tk.BooleanVar(value=False)
        tk.Checkbutton(options_frame, text="Show all results (not just top 5)", 
                      variable=show_all_var, font=('Arial', 9)).pack(anchor='w', padx=10, pady=5)
        
        # Action buttons
        button_frame = tk.Frame(main_frame)
        button_frame.pack(fill='x', pady=10)
        
        def run_test():
            # Collect selected models
            selected_models = [model for model, var in model_vars.items() if var.get()]
            if not selected_models:
                messagebox.showwarning("No Models", "Please select at least one model to test.")
                return
            
            # Collect parameter options
            selected_test_sizes = [size for size, var in test_size_vars.items() if var.get()]
            selected_scalings = [method for method, var in scaling_vars.items() if var.get()]
            
            if not selected_test_sizes:
                messagebox.showwarning("No Test Sizes", "Please select at least one test size.")
                return
            
            if not selected_scalings:
                messagebox.showwarning("No Scaling", "Please select at least one scaling method.")
                return
            
            # Close dialog and run tests
            canvas.unbind_all("<MouseWheel>")
            dialog.destroy()
            class_weight_mode = self.class_weight_var.get().strip().lower()
            class_weight = 'balanced' if class_weight_mode == 'balanced' else None

            self._run_model_comparison(selected_models, selected_test_sizes, selected_scalings,
                                      int(cv_folds_var.get()), show_all_var.get(), class_weight)
        
        tk.Button(button_frame, text="▶️ Run Tests", command=run_test, 
                 bg='#27ae60', fg='white', font=('Arial', 10, 'bold'), pady=8).pack(side='left', expand=True, fill='x', padx=(0,5))
        tk.Button(button_frame, text="✖️ Cancel", command=on_close,
                 bg='#e74c3c', fg='white', font=('Arial', 10, 'bold'), pady=8).pack(side='left', expand=True, fill='x', padx=(5,0))
        
        # Center dialog
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.update_idletasks()
        x = (dialog.winfo_screenwidth() // 2) - (dialog.winfo_width() // 2)
        y = (dialog.winfo_screenheight() // 2) - (dialog.winfo_height() // 2)
        dialog.geometry(f"+{x}+{y}")
    
    def _run_model_comparison(self, models, test_sizes, scalings, cv_folds, show_all, class_weight=None):
        """Run model comparison in background thread"""
        self.test_models_button.config(state='disabled', text='⏳ Testing...')
        self.run_button.config(state='disabled')
        
        thread = threading.Thread(target=self._model_comparison_thread, 
                                 args=(models, test_sizes, scalings, cv_folds, show_all, class_weight), 
                                 daemon=True)
        thread.start()
    
    def _model_comparison_thread(self, models, test_sizes, scalings, cv_folds, show_all, class_weight):
        """Background thread for model comparison"""
        try:
            self._log(f"\n{'='*60}\n")
            self._log("🔬 MODEL COMPARISON & TESTING\n")
            self._log(f"{'='*60}\n\n")
            
            # Prepare data
            merged_df, group_map, feature_cols = self._filter_and_merge_data()
            
            if merged_df is None or len(merged_df) == 0:
                self._log("❌ No data available\n")
                return
            
            from main_script.ml_models import MetabolomicsMLAnalysis
            
            # Get feature ID column from verified columns
            feature_id_col = self._get_feature_id_column(merged_df)
            
            ml_analysis = MetabolomicsMLAnalysis(
                data_df=merged_df,
                group_assignments=group_map,
                feature_columns=feature_cols,
                feature_id_col=feature_id_col
            )
            
            # Test all combinations
            results = []
            total_tests = len(models) * len(test_sizes) * len(scalings)
            current_test = 0
            
            self._log(f"Testing {total_tests} combinations...\n")
            self._log(f"Models: {', '.join(models)}\n")
            self._log(f"Test Sizes: {', '.join(test_sizes)}\n")
            self._log(f"Scaling: {', '.join(scalings)}\n")
            self._log(f"CV Folds: {cv_folds}\n\n")
            self._log(f"Class Weight: {'balanced' if class_weight == 'balanced' else 'none'}\n\n")
            
            for model_name in models:
                for test_size in test_sizes:
                    for scaling in scalings:
                        current_test += 1
                        self._log(f"[{current_test}/{total_tests}] Testing {model_name} | test_size={test_size} | scaling={scaling}...\n")
                        
                        try:
                            # Run classification
                            result = ml_analysis.run_classification(
                                model_name=model_name,
                                test_size=float(test_size),
                                cv_folds=cv_folds,
                                scaling_method=scaling,
                                class_weight=class_weight
                            )
                            
                            # Store key metrics (handle CV-only mode)
                            results.append({
                                'model': model_name,
                                'test_size': test_size,
                                'scaling': scaling,
                                'class_weight': 'balanced' if class_weight == 'balanced' else 'none',
                                'cv_mean': result['cv_mean_accuracy'],
                                'cv_std': result['cv_std_accuracy'],
                                'test_acc': result['test_accuracy'] if result['test_accuracy'] is not None else 'N/A',
                                'train_acc': result['train_accuracy']
                            })
                            
                        except Exception as e:
                            self._log(f"  ⚠️ Failed: {str(e)}\n")
            
            # Sort by CV mean accuracy (primary) and test accuracy (secondary)
            # Handle N/A test_acc by treating it as 0 for sorting
            results.sort(key=lambda x: (x['cv_mean'], x['test_acc'] if isinstance(x['test_acc'], (int, float)) else 0), reverse=True)
            
            # Display results
            self._log(f"\n{'='*60}\n")
            self._log("📊 MODEL COMPARISON RESULTS\n")
            self._log(f"{'='*60}\n\n")
            
            # Show top results or all
            display_count = len(results) if show_all else min(5, len(results))
            
            for i, r in enumerate(results[:display_count], 1):
                rank_emoji = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}."
                self._log(f"{rank_emoji} {r['model']}\n")
                self._log(f"   Test Size: {r['test_size']} | Scaling: {r['scaling']} | Class Weight: {r['class_weight']}\n")
                self._log(f"   CV Accuracy: {r['cv_mean']:.4f} ± {r['cv_std']:.4f}\n")
                self._log(f"   Test Accuracy: {r['test_acc']:.4f}\n")
                self._log(f"   Train Accuracy: {r['train_acc']:.4f}\n\n")
            
            if not show_all and len(results) > 5:
                self._log(f"... and {len(results) - 5} more combinations tested\n\n")
            
            # Recommendation
            best = results[0]
            self._log(f"\n{'='*60}\n")
            self._log("✅ RECOMMENDED CONFIGURATION\n")
            self._log(f"{'='*60}\n")
            self._log(f"Model: {best['model']}\n")
            self._log(f"Test Size: {best['test_size']}\n")
            self._log(f"Scaling: {best['scaling']}\n")
            self._log(f"Class Weight: {best['class_weight']}\n")
            self._log(f"Expected CV Accuracy: {best['cv_mean']:.4f} ± {best['cv_std']:.4f}\n")
            self._log(f"\n💡 Use these settings for your final analysis!\n")
            
            # Auto-apply best settings
            self.root.after(0, lambda: self._apply_best_settings(best))
            
        except Exception as e:
            logger.error(f"Model comparison error: {e}", exc_info=True)
            self._log(f"\n❌ Error during model comparison: {str(e)}\n")
            
        finally:
            self.root.after(0, lambda: self.test_models_button.config(state='normal', text='🔬 Test Models'))
            self.root.after(0, lambda: self.run_button.config(state='normal'))
    
    def _apply_best_settings(self, best_config):
        """Apply the best configuration to the UI"""
        # Set best model checkbox
        best_model = best_config.get('model', 'Random Forest')
        for model, var in self.model_checkboxes.items():
            var.set(model == best_model)
        
        self.test_size_var.set(best_config['test_size'])
        self.scaling_var.set(best_config['scaling'])
        if 'class_weight' in best_config:
            self.class_weight_var.set(best_config['class_weight'])
        self._log("\n✅ Best settings applied to UI!\n")
        messagebox.showinfo("Settings Applied", 
                           f"Best configuration applied:\n\n"
                           f"Model: {best_model}\n"
                           f"Test Size: {best_config['test_size']}\n"
                           f"Scaling: {best_config['scaling']}\n"
                           f"Class Weight: {best_config.get('class_weight', 'none')}\n\n"
                           f"You can now run your final analysis!")
    
    def _check_overfitting_warnings(self, results, n_features, n_samples):
        """Check for signs of overfitting and display warnings"""
        
        # Calculate feature-to-sample ratio
        feature_ratio = n_features / n_samples if n_samples > 0 else 0
        
        # Get metrics (handle CV-only mode)
        test_acc = results.get('test_accuracy')  # Can be None in CV-only mode
        train_acc = results.get('train_accuracy', 0)
        cv_mean = results.get('cv_mean_accuracy', 0)
        cv_only_mode = results.get('cv_only_mode', False) or test_acc is None
        
        # Check for robustness stats if available
        robustness = results.get('robustness_stats', {})
        test_std = robustness.get('test_std', 0) if robustness else 0
        
        # Flags for overfitting
        warnings = []
        critical_warnings = []
        
        # RED FLAG 1: Perfect or near-perfect test accuracy (skip if CV-only)
        if not cv_only_mode and test_acc is not None and test_acc >= 0.95:
            if test_std == 0:
                critical_warnings.append(
                    "🚨 CRITICAL: Test accuracy = 1.00 with zero variance across all runs\n"
                    "   → This indicates the model is exploiting high-dimensional separability\n"
                    "   → NOT true generalization - data is trivially separable\n"
                )
            else:
                warnings.append(
                    "⚠️  Test accuracy very high (≥0.95) - possible overfitting\n"
                )
        
        # RED FLAG 2: High feature-to-sample ratio
        if feature_ratio > 5:
            severity = "CRITICAL" if feature_ratio > 20 else "WARNING"
            critical_warnings.append(
                f"🚨 {severity}: {n_features} features vs {n_samples} samples (ratio: {feature_ratio:.1f}:1)\n"
                f"   → Curse of dimensionality - model can memorize, not learn\n"
                f"   → Recommendation: Use feature selection or dimensionality reduction\n"
            )
        
        # RED FLAG 3: Near-perfect CV with small sample size
        if cv_mean >= 0.95 and n_samples < 30:
            warnings.append(
                f"⚠️  CV accuracy very high ({cv_mean:.3f}) with only {n_samples} samples\n"
                "   → Likely residual overfitting even with cross-validation\n"
            )
        
        # RED FLAG 4: Training accuracy = 1.0 (perfect fit)
        if train_acc >= 0.999:
            warnings.append(
                "⚠️  Perfect training accuracy (1.00)\n"
                "   → Model is memorizing training data\n"
            )
        
        # Display warnings if any
        if critical_warnings or warnings:
            self._log(f"\n{'='*60}\n")
            self._log("⚠️  OVERFITTING RISK ASSESSMENT\n")
            self._log(f"{'='*60}\n\n")
            
            if critical_warnings:
                self._log("CRITICAL ISSUES:\n")
                for warning in critical_warnings:
                    self._log(warning + "\n")
            
            if warnings:
                self._log("WARNINGS:\n")
                for warning in warnings:
                    self._log(warning + "\n")
            
            self._log("\n💡 RECOMMENDATIONS TO FIX OVERFITTING:\n")
            self._log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
            
            if feature_ratio > 10:
                self._log("1. REDUCE DIMENSIONALITY (Critical):\n")
                self._log("   • Use PCA/LDA to reduce to 5-10 components\n")
                self._log("   • Apply univariate feature selection (t-test, fold-change)\n")
                self._log("   • Use stability threshold ≥80% to keep only robust features\n\n")
            
            self._log("2. INCREASE REGULARIZATION:\n")
            self._log("   • Change to 'Strong' regularization (C=0.1)\n")
            self._log("   • Try L1 (Lasso) for automatic feature selection\n")
            self._log("   • Increase filtering threshold (e.g., 80% samples per group)\n\n")
            
            self._log("3. MODEL SELECTION:\n")
            self._log("   • Avoid complex models (Random Forest, Gradient Boosting) with n<30\n")
            self._log("   • Stick to linear models (Logistic Regression, LDA) with strong regularization\n\n")
            
            self._log("4. VALIDATION:\n")
            self._log("   • Acquire more samples (aim for n≥50 per group minimum)\n")
            self._log("   • Use external validation set (different cohort/timepoint)\n")
            self._log("   • Increase repeated runs to 20-30 to assess true stability\n\n")
            
            self._log("⚠️  Current results are NOT publishable without addressing these issues.\n")
            self._log("   Reviewers will reject claims based on memorization, not biology.\n\n")
    
    def _filter_and_merge_data(self) -> Tuple[Optional[pd.DataFrame], Dict[str, str], List[str]]:
        """
        Filter input data with proper feature filtering.
        
        Returns:
            (merged_df, group_map, feature_columns)
        """
        mode = self.ml_data_mode.get()
        
        # Get group assignments
        group_map = {}
        for col, var in self.sample_group_vars.items():
            group = var.get().strip()
            if group:
                group_map[col] = group
        
        if len(group_map) == 0:
            self._log("❌ No group assignments\n")
            return None, {}, []
        
        # Get per-group feature detection thresholds.
        group_thresholds = self._get_group_detection_thresholds(group_map)
        threshold_type = self.min_samples_type_var.get()
        if threshold_type == 'percentage':
            try:
                percent = float(self.min_samples_percent_var.get())
            except Exception:
                percent = 50.0
            parts = [f"{grp}>={group_thresholds[grp]}" for grp in sorted(group_thresholds.keys())]
            self._log(f"Feature filtering thresholds ({percent:g}% per group): {', '.join(parts)}\n\n")
        else:
            min_samples = min(group_thresholds.values()) if group_thresholds else 1
            self._log(f"Feature filtering threshold: {min_samples} samples per group\n\n")
        
        # Apply sample group filtering first
        filtered_map, excluded = self._filter_groups_by_threshold(group_map, group_thresholds)
        
        if excluded:
            self._log(f"⚠️  Excluded sample groups (below threshold):\n")
            for grp, cnt in excluded.items():
                self._log(f"  • {grp}: {cnt} samples\n")
        
        if len(filtered_map) == 0:
            self._log("❌ No samples remain after filtering\n")
            return None, {}, []
        
        self._log(f"✅ Retained: {len(filtered_map)} samples across {len(set(filtered_map.values()))} groups\n\n")
        
        # Filter features on input sheet(s)
        frames = []
        
        if self.pos_df is not None and len(self.verified_pos_sample_cols) > 0:
            valid_cols = [c for c in self.verified_pos_sample_cols if c in filtered_map]
            if valid_cols:
                pos_subset = self.pos_df.copy()
                features_before = len(pos_subset)
                
                # Filter features: must meet each group's required detection count.
                pos_subset = self._filter_features_by_group(pos_subset, valid_cols, filtered_map, group_thresholds)
                features_after = len(pos_subset)
                
                if len(pos_subset) > 0:
                    # Keep all non-sample columns plus valid sample columns
                    keep_cols = [c for c in pos_subset.columns if c not in self.verified_pos_sample_cols or c in valid_cols]
                    pos_subset = pos_subset[keep_cols]
                    frames.append(pos_subset)
                    self._log(f"✅ Input sheet: {features_before} → {features_after} features ({features_before - features_after} filtered), {len(valid_cols)} samples\n")
        
        if self.neg_df is not None and len(self.verified_neg_sample_cols) > 0:
            valid_cols = [c for c in self.verified_neg_sample_cols if c in filtered_map]
            if valid_cols:
                neg_subset = self.neg_df.copy()
                features_before = len(neg_subset)
                
                # Filter features: must meet each group's required detection count.
                neg_subset = self._filter_features_by_group(neg_subset, valid_cols, filtered_map, group_thresholds)
                features_after = len(neg_subset)
                
                if len(neg_subset) > 0:
                    keep_cols = [c for c in neg_subset.columns if c not in self.verified_neg_sample_cols or c in valid_cols]
                    neg_subset = neg_subset[keep_cols]
                    neg_subset['Polarity_Source'] = 'Negative'
                    frames.append(neg_subset)
                    self._log(f"✅ Negative: {features_before} → {features_after} features ({features_before - features_after} filtered), {len(valid_cols)} samples\n")
        
        if len(frames) == 0:
            self._log("❌ No valid data frames after filtering\n")
            return None, {}, []

        # Combine frames if needed, but do not apply area-based sorting.
        if len(frames) == 1:
            combined = frames[0].copy()
        else:
            combined = pd.concat(frames, ignore_index=True)

        sample_cols = list(filtered_map.keys())
        combined = combined.reset_index(drop=True)
        
        # Remove duplicates (keep first = highest Area)
        before = len(combined)
        
        if mode == 'lipid':
            id_col = self._get_id_column(combined)
            if id_col and id_col in combined.columns:
                combined = combined.drop_duplicates(subset=[id_col], keep='first')
                self._log(f"Deduplicated by {id_col}: removed {before - len(combined)} rows\n")
        else:
            if 'Name' in combined.columns:
                combined = combined.drop_duplicates(subset=['Name'], keep='first')
                self._log(f"Deduplicated by Name: removed {before - len(combined)} rows\n")
        
        # Remove legacy helper columns if present.
        combined = combined.drop(columns=['Polarity_Source', 'Area (Max.)'], errors='ignore')
        
        self._log(f"\n✅ Final filtered data: {len(combined)} features × {len(sample_cols)} samples\n")
        
        return combined, filtered_map, sample_cols
    
    def _filter_features_by_group(self, df: pd.DataFrame, sample_cols: List[str], 
                                   group_map: Dict[str, str], group_thresholds: Dict[str, int]) -> pd.DataFrame:
        """
        Filter features to keep only those with >= threshold detections in EACH group.
        
        Example: If threshold=8, and feature has Control=7, PD=10, it's removed
        because Control didn't meet the threshold.
        """
        # Group samples by their group assignment
        groups_samples = {}
        for sample, group in group_map.items():
            if sample in sample_cols:
                if group not in groups_samples:
                    groups_samples[group] = []
                groups_samples[group].append(sample)
        
        # For each feature (row), check if it meets threshold in ALL groups
        keep_mask = []
        
        for idx, row in df.iterrows():
            meets_threshold_all_groups = True
            
            for group, group_sample_cols in groups_samples.items():
                # Count non-zero, non-NaN values in this group's samples
                values = row[group_sample_cols]
                detected_count = ((values > 0) & values.notna()).sum()
                
                required = int(group_thresholds.get(group, 1))
                if detected_count < required:
                    meets_threshold_all_groups = False
                    break
            
            keep_mask.append(meets_threshold_all_groups)
        
        return df[keep_mask].reset_index(drop=True)
    
    def _get_group_detection_thresholds(self, group_map: Dict[str, str]) -> Dict[str, int]:
        """Calculate required detected-sample counts per group for feature filtering."""
        counts = Counter(group_map.values())
        if not counts:
            return {}

        threshold_type = self.min_samples_type_var.get()

        if threshold_type == 'percentage':
            percent = float(self.min_samples_percent_var.get())
            return {
                grp: max(1, int(np.ceil(cnt * percent / 100.0)))
                for grp, cnt in counts.items()
            }
        else:
            threshold = max(1, int(self.min_samples_per_group_var.get()))
            return {grp: threshold for grp in counts.keys()}

    def _get_min_samples_threshold(self, group_map: Dict[str, str]) -> int:
        """Backward-compatible helper returning the minimum per-group threshold."""
        thresholds = self._get_group_detection_thresholds(group_map)
        return min(thresholds.values()) if thresholds else 1
    
    def _filter_groups_by_threshold(self, group_map: Dict[str, str], group_thresholds: Dict[str, int]) -> Tuple[Dict[str, str], Dict[str, int]]:
        """Filter groups below threshold"""
        counts = Counter(group_map.values())
        excluded = {
            grp: cnt
            for grp, cnt in counts.items()
            if cnt < int(group_thresholds.get(grp, 1))
        }
        filtered = {
            sample: grp
            for sample, grp in group_map.items()
            if counts[grp] >= int(group_thresholds.get(grp, 1))
        }
        return filtered, excluded
    
    def _get_id_column(self, df: pd.DataFrame) -> Optional[str]:
        """Get the ID column for deduplication"""
        mode = self.ml_data_mode.get()
        
        if mode == 'lipid':
            for col in ['LipidID', 'Lipid_ID', 'lipidid']:
                if col in df.columns:
                    return col
        else:
            for col in ['Name', 'Feature', 'Metabolite']:
                if col in df.columns:
                    return col
        
        return None
    
    def _export_results(self):
        """Export ML results to Excel"""
        if self.ml_results is None:
            messagebox.showwarning("No Results", "No results to export")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save ML Results",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
        )
        
        if not file_path:
            return
        
        try:
            # Export main results to temporary file first.
            # This avoids creating an empty destination workbook if export fails upstream.
            from main_script.ml_models import MetabolomicsMLAnalysis

            ml_analysis = MetabolomicsMLAnalysis(
                data_df=self.merged_data if self.merged_data is not None else pd.DataFrame(),
                group_assignments={},
                feature_columns=[]
            )
            ml_analysis.results = self.ml_results

            analysis_type = self.analysis_type_var.get()
            if analysis_type == 'feature_importance':
                analysis_type = 'classification'

            temp_file = file_path.replace('.xlsx', '_temp.xlsx')
            ml_analysis.export_results_to_excel(temp_file, analysis_type)

            if not os.path.exists(temp_file):
                raise RuntimeError("Temporary export file was not created.")

            temp_wb = pd.ExcelFile(temp_file)
            if not temp_wb.sheet_names:
                temp_wb.close()
                raise RuntimeError("Temporary export produced no visible sheets.")

            with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                # Copy all generated sheets from temp export
                for sheet_name in temp_wb.sheet_names:
                    sheet_name_str = str(sheet_name)
                    df = pd.read_excel(temp_wb, sheet_name_str)
                    df.to_excel(writer, sheet_name=sheet_name_str, index=False)
                temp_wb.close()
                
                # Add filtered dataset sheet
                if self.merged_data is not None:
                    self.merged_data.to_excel(writer, sheet_name='Filtered_Dataset', index=False)
                    self._log(f"✅ Added 'Filtered_Dataset' sheet with {len(self.merged_data)} features\n")
                
                # Add reproducibility settings sheet
                settings_data = {
                    'Parameter': [
                        'Analysis Date',
                        'Model',
                        'Test Size',
                        'CV Folds',
                        'Scaling Method',
                        'Class Weight',
                        'Regularization Type',
                        'Regularization Strength',
                        'Max Iterations',
                        'Repeated Runs',
                        'Base Seed',
                        'Seed Strategy',
                        'Imputation Method',
                        'Imputation KNN Neighbors',
                        'Auto Skip Scaling (Trees)',
                        'Feature Selection Method',
                        'Features Before Selection',
                        'Features After Selection',
                        'Stability Tracking',
                        'Stability Threshold',
                        'Feature Filter',
                    ],
                    'Value': [
                        pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
                        (", ".join(self.ml_results.get('models_trained', [])) if self.ml_results.get('comparison_type') == 'multi_model' else self.ml_results.get('model_name', 'N/A')),
                        self.ml_results.get('test_size', 'N/A'),
                        self.ml_results.get('cv_folds', 'N/A'),
                        self.ml_results.get('scaling_method', 'N/A'),
                        self.ml_results.get('class_weight', 'none'),
                        self.analysis_type_var.get(),
                        self.regularization_strength_var.get(),
                        self.max_iter_var.get(),
                        self.ml_results.get('n_repeats', '1'),
                        self.base_seed_var.get(),
                        "Fixed base seed with increments (seed, seed+1, seed+2, ...)",
                        self.ml_results.get('imputation_method', 'N/A'),
                        self.ml_results.get('imputation_knn_neighbors', 'N/A'),
                        'ON' if self.ml_results.get('auto_skip_scaling_for_trees') else 'OFF',
                        self.ml_results.get('feature_selection_method', 'none'),
                        self.ml_results.get('n_features_before_selection', self.ml_results.get('n_features', 'N/A')),
                        self.ml_results.get('n_features_after_selection', self.ml_results.get('n_features', 'N/A')),
                        'ON' if self.stability_tracking_var.get() else 'OFF',
                        f"{self.stability_threshold_var.get()}%",
                        self.ml_results.get('feature_filter', 'None'),
                    ]
                }
                pd.DataFrame(settings_data).to_excel(writer, sheet_name='Analysis_Settings', index=False)
                self._log(f"✅ Added 'Analysis_Settings' sheet for reproducibility\n")

                if os.path.exists(temp_file):
                    os.remove(temp_file)
            
            self._log(f"\n✅ Results exported to: {os.path.basename(file_path)}\n\n")
            messagebox.showinfo("Export Complete", f"Results exported successfully!\n\nFile: {os.path.basename(file_path)}")
            
        except Exception as e:
            logger.error(f"Export error: {e}", exc_info=True)
            self._log(f"\n❌ Export failed: {str(e)}\n")
            messagebox.showerror("Export Error", f"Failed to export results:\n{str(e)}")
    
    def _log(self, message: str):
        """Thread-safe logging"""
        def append():
            self.ml_log.insert(tk.END, message)
            self.ml_log.see(tk.END)
            self.ml_log.update_idletasks()
        
        try:
            self.root.after(0, append)
        except Exception:
            # Fallback for cases where root is not available
            print(message, end='')
    
    def _clear_log(self):
        """Clear the log text area and show welcome message"""
        self.ml_log.delete('1.0', tk.END)
        
        # Re-insert welcome message
        self.ml_log.insert(tk.END, "="*60 + "\n")
        self.ml_log.insert(tk.END, "🤖 Machine Learning Analysis Tool\n")
        self.ml_log.insert(tk.END, "="*60 + "\n\n")
        self.ml_log.insert(tk.END, "Workflow:\n")
        self.ml_log.insert(tk.END, "1. Import Excel file\n")
        self.ml_log.insert(tk.END, "2. Verify columns\n")
        self.ml_log.insert(tk.END, "3. Configure sample groups\n")
        self.ml_log.insert(tk.END, "4. Run classification analysis\n\n")
        self.ml_log.insert(tk.END, "Ready to begin!\n\n")
    
    def _setup_ml_config_functions(self):
        """Setup config save/load infrastructure"""
        self._ml_config_loaded = False
    
    def _ml_config_file(self):
        """Get path to ml config file"""
        try:
            return resolve_runtime_config_path('ml_config.json')
        except Exception:
            return None
    
    def _gather_ml_config(self) -> dict:
        """Gather all ML tab settings into a config dict"""
        # Get selected models from checkboxes
        selected_models = [model for model, var in self.model_checkboxes.items() if var.get()]
        primary_model = selected_models[0] if selected_models else 'Random Forest'
        
        config = {
            'group_definitions': self.group_definitions.copy(),
            'group_count': len(self.group_definitions),
            'auto_assign_patterns': self.auto_assign_patterns.copy(),
            'ml_data_mode': 'custom',
            'analysis_type': 'classification',
            'model_type': primary_model,
            'selected_models': selected_models,
            'test_size': self.test_size_var.get(),
            'cv_folds': self.cv_folds_var.get(),
            'scaling': self.scaling_var.get(),
            'class_weight': self.class_weight_var.get(),
            'regularization_type': self.regularization_type_var.get(),
            'regularization_strength': self.regularization_strength_var.get(),
            'max_iter': self.max_iter_var.get(),
            'repeated_runs': self.repeated_runs_var.get(),
            'base_seed': self.base_seed_var.get(),
            'stability_tracking': self.stability_tracking_var.get(),
            'stability_threshold': self.stability_threshold_var.get(),
            'tune_hyperparameters': self.tune_hyperparameters_var.get(),
            'tuning_strategy': self.tuning_strategy_var.get(),
            'tuning_iter': self.tuning_iter_var.get(),
            'use_repeated_cv': self.use_repeated_cv_var.get(),
            'cv_repeats': self.cv_repeats_var.get(),
            'nested_cv': self.nested_cv_var.get(),
            'calibration_method': self.calibration_method_var.get(),
            'permutation_test_runs': self.permutation_test_runs_var.get(),
            'imputation_method': self.imputation_method_var.get(),
            'imputation_knn_neighbors': self.imputation_knn_neighbors_var.get(),
            'auto_skip_scaling_for_trees': self.auto_skip_scaling_tree_var.get(),
            'feature_selection_method': self.feature_selection_method_var.get(),
            'variance_percentile': self.variance_percentile_var.get(),
            'univariate_k': self.univariate_k_var.get(),
            'lasso_c': self.lasso_c_var.get(),
            'rfe_n_features': self.rfe_n_features_var.get(),
            'auto_generate_plots': self.auto_generate_plots_var.get(),
            'generate_pairwise_addon': self.generate_pairwise_addon_var.get(),
            'pairwise_pvalue_map': self.pairwise_pvalue_map.copy(),
            'working_folder': self.working_folder_var.get(),
            'figure_settings': self.figure_settings.copy(),
            'top_n_values': list(self.top_n_values),
        }
        return config
    
    def _save_ml_config(self):
        """Save ML config to JSON file"""
        try:
            if not self._ml_config_loaded:
                return  # Don't save if we're still loading
            
            config_file = self._ml_config_file()
            if not config_file:
                return
            
            config = self._gather_ml_config()
            
            # Ensure config directory exists
            config_dir = os.path.dirname(config_file)
            os.makedirs(config_dir, exist_ok=True)
            
            # Write to JSON
            with open(config_file, 'w') as f:
                json.dump(config, f, indent=2)
            
            logger.info(f"ML config saved to {config_file}")
        
        except Exception as e:
            logger.error(f"Error saving ML config: {e}")
    
    def _load_ml_config(self):
        """Load ML config from JSON file"""
        try:
            config_file = self._ml_config_file()
            if not config_file or not os.path.exists(config_file):
                self._ml_config_loaded = True
                return
            
            # Check if file is empty
            if os.path.getsize(config_file) == 0:
                logger.warning(f"ML config file is empty: {config_file}")
                self._ml_config_loaded = True
                return
            
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            # Apply group definitions (critical: must load group_count and reconstruct group_id_vars)
            gd = config.get('group_definitions') or {}
            if gd:
                self.group_definitions = gd.copy()
                self.group_count = config.get('group_count', len(gd))
                # Reconstruct group_id_vars with loaded labels
                self.group_id_vars = {gid: tk.StringVar(value=label) for gid, label in self.group_definitions.items()}
                self.refresh_group_ui()
                logger.info(f"Loaded {self.group_count} groups: {list(self.group_definitions.values())}")
            
            # Load auto-assign patterns
            if 'auto_assign_patterns' in config:
                self.auto_assign_patterns = config.get('auto_assign_patterns', {}).copy()
                logger.info(f"Loaded auto-assign patterns for {len(self.auto_assign_patterns)} groups")
            
            # Apply ML settings (only for UI variables that exist)
            if 'test_size' in config:
                self.test_size_var.set(config.get('test_size', '0.3'))
            if 'cv_folds' in config:
                self.cv_folds_var.set(config.get('cv_folds', '5'))
            if 'scaling' in config:
                self.scaling_var.set(config.get('scaling', 'standard'))
            if 'class_weight' in config:
                self.class_weight_var.set(config.get('class_weight', 'none'))
            # Analysis type is fixed for ML tab.
            self.analysis_type_var.set('classification')
            # Restore model checkboxes from config
            if 'selected_models' in config:
                selected_models = config.get('selected_models', [])
                for model, var in self.model_checkboxes.items():
                    var.set(model in selected_models)
            elif 'model_type' in config:
                # Fallback: if old config format, restore just the primary model
                primary_model = config.get('model_type', 'Random Forest')
                for model, var in self.model_checkboxes.items():
                    var.set(model == primary_model)
            if 'regularization_type' in config:
                self.regularization_type_var.set(config.get('regularization_type', 'l2'))
            if 'regularization_strength' in config:
                self.regularization_strength_var.set(config.get('regularization_strength', 'medium'))
            if 'max_iter' in config:
                self.max_iter_var.set(config.get('max_iter', '1000'))
            if 'repeated_runs' in config:
                self.repeated_runs_var.set(config.get('repeated_runs', '1'))
            if 'base_seed' in config:
                self.base_seed_var.set(config.get('base_seed', '42'))
            if 'stability_tracking' in config:
                self.stability_tracking_var.set(config.get('stability_tracking', True))
            if 'stability_threshold' in config:
                self.stability_threshold_var.set(config.get('stability_threshold', '70'))
            if 'tune_hyperparameters' in config:
                self.tune_hyperparameters_var.set(config.get('tune_hyperparameters', True))
            if 'tuning_strategy' in config:
                self.tuning_strategy_var.set(config.get('tuning_strategy', 'grid'))
            if 'tuning_iter' in config:
                self.tuning_iter_var.set(config.get('tuning_iter', '20'))
            if 'use_repeated_cv' in config:
                self.use_repeated_cv_var.set(config.get('use_repeated_cv', False))
            if 'cv_repeats' in config:
                self.cv_repeats_var.set(config.get('cv_repeats', '3'))
            if 'nested_cv' in config:
                self.nested_cv_var.set(config.get('nested_cv', False))
            if 'calibration_method' in config:
                self.calibration_method_var.set(config.get('calibration_method', 'none'))
            if 'permutation_test_runs' in config:
                self.permutation_test_runs_var.set(config.get('permutation_test_runs', '0'))
            if 'imputation_method' in config:
                self.imputation_method_var.set(config.get('imputation_method', 'half_min'))
            if 'imputation_knn_neighbors' in config:
                self.imputation_knn_neighbors_var.set(config.get('imputation_knn_neighbors', '5'))
            if 'auto_skip_scaling_for_trees' in config:
                self.auto_skip_scaling_tree_var.set(config.get('auto_skip_scaling_for_trees', False))
            if 'feature_selection_method' in config:
                self.feature_selection_method_var.set(config.get('feature_selection_method', 'none'))
            if 'variance_percentile' in config:
                self.variance_percentile_var.set(config.get('variance_percentile', '10'))
            if 'univariate_k' in config:
                self.univariate_k_var.set(config.get('univariate_k', '50'))
            if 'lasso_c' in config:
                self.lasso_c_var.set(config.get('lasso_c', '0.1'))
            if 'rfe_n_features' in config:
                self.rfe_n_features_var.set(config.get('rfe_n_features', '50'))
            if 'auto_generate_plots' in config:
                self.auto_generate_plots_var.set(config.get('auto_generate_plots', True))
            if 'generate_pairwise_addon' in config:
                self.generate_pairwise_addon_var.set(config.get('generate_pairwise_addon', False))
            if 'pairwise_pvalue_map' in config and isinstance(config.get('pairwise_pvalue_map'), dict):
                self.pairwise_pvalue_map = config.get('pairwise_pvalue_map', {}).copy()
            if 'working_folder' in config:
                working_folder = config.get('working_folder', '')
                self.working_folder_var.set(working_folder)
                if hasattr(self, 'working_folder_label'):
                    if working_folder:
                        self.working_folder_label.config(text=f"Working folder: {working_folder}", fg='#2c3e50')
                    else:
                        self.working_folder_label.config(text='Working folder: Not set', fg='#666666')
            if 'figure_settings' in config and isinstance(config.get('figure_settings'), dict):
                defaults = self._default_figure_settings()
                loaded = config.get('figure_settings', {})
                merged = defaults.copy()
                for k, v in loaded.items():
                    try:
                        val = float(v)
                        if val > 0:
                            merged[k] = val
                    except (ValueError, TypeError):
                        continue
                self.figure_settings = merged
            if 'top_n_values' in config:
                self.top_n_values = self._normalize_top_n_values(config.get('top_n_values'))

            self._update_feature_selection_controls()

            # Enforce fixed mode regardless of legacy saved values.
            self.ml_data_mode.set('custom')
            if hasattr(self, 'mode_desc_label'):
                self.mode_desc_label.config(text='Custom mode for preprocessed/combined data')
            
            logger.info(f"ML config loaded from {config_file}")
            self._ml_config_loaded = True
        
        except json.JSONDecodeError as e:
            logger.warning(f"ML config file has invalid JSON (will use defaults): {e}")
            self._ml_config_loaded = True
        except Exception as e:
            logger.error(f"Error loading ML config: {e}")
            self._ml_config_loaded = True

