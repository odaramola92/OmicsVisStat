import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext, simpledialog
import logging
import datetime
import threading
import os
import json
import pandas as pd
import time
import seaborn as sns
from typing import Optional, Dict, List, Any
import traceback

# Import shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.shared.utils import resolve_runtime_config_path
# Column verification removed - automated pairwise detection handles all cases

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class VisualizationTab(BaseTab):
    """Visualization Tab - Placeholder for data visualization
    
    This tab will handle visualization of metabolite and pathway data.
    """
    
    def __init__(self, parent, data_manager):
        """Initialize Visualization tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles (runs only once)
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Setup memory_store as reference to data_manager's memory store
        self.memory_store = self.data_manager.memory_store
        
        # Create UI
        self.setup_ui()
        print("[OK] Visualization Tab initialized")
    
    def _enable_mousewheel_scroll(self, widget):
        """Enable mousewheel scrolling on a canvas widget."""
        def _on_mousewheel(event):
            try:
                if widget.winfo_exists():
                    widget.yview_scroll(int(-1 * (event.delta / 120)), "units")
            except tk.TclError:
                pass
        
        widget.bind("<Enter>", lambda e: widget.bind("<MouseWheel>", _on_mousewheel))
        widget.bind("<Leave>", lambda e: widget.unbind("<MouseWheel>"))
    
    def setup_ui(self):
        """Setup the UI for the Visualization tab with full user control."""
        # Create internal notebook for Visualization and Statistics sub-tabs
        self.notebook = ttk.Notebook(self.frame)
        self.notebook.pack(fill='both', expand=True, padx=0, pady=0)
        
        self.setup_visualization_tab()
    
    def setup_visualization_tab(self):
        from main_script.metabolites_visualization import (
            PCAParams, VolcanoParams, BoxplotParams, BargraphParams, HeatmapParams, ROCParams,
            CommonVizContext, VizResults
        )
        
        # Create Visualization sub-tab within this tab
        self.viz_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.viz_tab, text='')
        
        # Initialize group definitions if not already present
        if not hasattr(self, 'viz_group_definitions'):
            self.viz_group_definitions = {}
        # Auto-import group definitions & assignments from Statistics config if available
        # (Priority before any visualization-specific defaults). This gives a seamless
        # cross-tab experience so the user does not need to redefine groups.
        if not self.viz_group_definitions:
            try:
                stats_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'statistics_config.json')
                if os.path.exists(stats_cfg_path):
                    with open(stats_cfg_path, 'r', encoding='utf-8') as f:
                        stats_cfg = json.load(f)
                    loaded_gd = stats_cfg.get('group_definitions') or {}
                    # group_definitions is stored as {GroupID: Label}
                    if loaded_gd:
                        self.viz_group_definitions = loaded_gd.copy()
                        try:
                            if hasattr(self, 'viz_preferred_group_order'):
                                self.viz_preferred_group_order.set(','.join(list(self.viz_group_definitions.values())))
                        except Exception:
                            pass
                        # Also import sample_group_assignments mapping if present
                        sample_map = stats_cfg.get('sample_group_assignments') or {}
                        if sample_map:
                            # Normalize mapping values to labels (already labels in stats config)
                            self.viz_group_mapping = sample_map.copy()
                        self.log_viz_message(f"🔄 Imported {len(self.viz_group_definitions)} groups from Statistics config: {list(self.viz_group_definitions.values())}")
            except Exception as e:
                try:
                    self.log_viz_message(f"⚠️ Failed importing Statistics groups: {e}")
                except Exception:
                    pass
        
        # Define visualization config auto-save handler
        if not hasattr(self, '_viz_config_changed'):
            def _viz_config_changed():
                """Auto-save visualization configuration when changes occur."""
                try:
                    if hasattr(self, '_save_viz_config'):
                        self._save_viz_config()
                except Exception:
                    pass
            self._viz_config_changed = _viz_config_changed
        
        # Initialize visualization parameters
        self.viz_params = {
            'pca': PCAParams(),
            'volcano': VolcanoParams(),
            'venn': None,
            'boxplot': BoxplotParams(),
            'bargraph': BargraphParams(),
            'heatmap': HeatmapParams(),
            'roc': ROCParams()
        }
        
        self.viz_selected = {
            'pca': tk.BooleanVar(value=True),
            'volcano': tk.BooleanVar(value=False),
            'venn': tk.BooleanVar(value=False),
            'boxplot': tk.BooleanVar(value=False),
            'bargraph': tk.BooleanVar(value=False),
            'heatmap': tk.BooleanVar(value=False),
            'roc': tk.BooleanVar(value=False)
        }
        
        # Initialize verified column assignments (populated by column dialog)
        self.verified_assignments = None
        
        # === NEW: State tracking for workflow enforcement ===
        self.groups_configured = False
        self.stat_cols_configured = False
        
        # Stat column assignments: {(g1, g2): {'pvalue': col, 'log2fc': col, ...}, ...}
        self.stat_column_assignments = {
            'id_column': None,  # Name, metabolite_id, LipidID
            'comparisons': {}   # Populated by Configure Stat Columns dialog
        }
        
        self.viz_color_map = {}
        self.viz_cancel_flag = threading.Event()
        
        # Initialize metabolite list dictionaries for per-comparison lists
        self.heatmap_metabolite_lists = {}
        self.roc_metabolite_lists = {}
        
        # Initialize auto-size variables
        self.heatmap_auto_size = tk.BooleanVar(value=True)
        
        # Main container with 3 panels (mimic Statistics tab layout so the
        # visualization tab expands to use available space instead of leaving
        # large blank areas). We use a canvas with an inner frame and a grid
        # of three columns (left config, center settings, right log/progress).
        self.viz_canvas = tk.Canvas(self.viz_tab, bg='#f0f0f0', highlightthickness=0)
        viz_scroll = ttk.Scrollbar(self.viz_tab, orient='vertical', command=self.viz_canvas.yview)
        viz_inner = tk.Frame(self.viz_canvas, bg='#f0f0f0')
        viz_inner.bind('<Configure>', lambda e: self.viz_canvas.configure(scrollregion=self.viz_canvas.bbox('all')))
        canvas_window = self.viz_canvas.create_window((0, 0), window=viz_inner, anchor='nw')
        self.viz_canvas.configure(yscrollcommand=viz_scroll.set)

        def _configure_viz_canvas(event):
            # Ensure the inner window width tracks the canvas width so
            # child frames expand horizontally. Do NOT force inner height.
            try:
                self.viz_canvas.itemconfig(canvas_window, width=event.width)
            except Exception:
                try:
                    self.viz_canvas.itemconfig(canvas_window, width=event.width)
                except Exception:
                    pass
            try:
                self.viz_canvas.configure(scrollregion=self.viz_canvas.bbox('all'))
            except Exception:
                pass

        self.viz_canvas.bind('<Configure>', _configure_viz_canvas)
        self.viz_canvas.pack(side='left', fill='both', expand=True, padx=0, pady=0)
        viz_scroll.pack(side='right', fill='y')

        # Container inside canvas where we place the three column layout
        container = tk.Frame(viz_inner, bg='#f0f0f0')
        container.pack(fill='both', expand=True, padx=5, pady=0)

        # Body frame with 3 columns: left (config), center (settings), right (log)
        body = tk.Frame(container, bg='#f0f0f0')
        body.pack(fill='both', expand=True)
        
        # CRITICAL: Ensure container expands to fill available canvas height
        # This allows body.grid_rowconfigure(0, weight=1) to work properly
        def _on_canvas_configure(event):
            # Get the height of the canvas viewport
            canvas_height = self.viz_canvas.winfo_height()
            # Get current scroll region height
            scroll_height = self.viz_canvas.bbox('all')[3] if self.viz_canvas.bbox('all') else 0
            # If content is smaller than canvas, expand the container to fill
            if scroll_height < canvas_height:
                viz_inner.configure(height=canvas_height)
                self.viz_canvas.itemconfig(canvas_window, height=canvas_height)
        
        self.viz_canvas.bind('<Configure>', _on_canvas_configure, add=True)
        # keep a reference for potential runtime layout adjustments
        self.viz_body = body

        # Column minimums and weights to distribute available width similar to statistics tab
        viz_column_mins = getattr(self, 'viz_column_mins', (320, 600, 300))
        try:
            c0, c1, c2 = viz_column_mins
        except Exception:
            c0, c1, c2 = 320, 600, 300
        body.grid_columnconfigure(0, weight=1, minsize=c0)
        body.grid_columnconfigure(1, weight=3, minsize=c1)
        body.grid_columnconfigure(2, weight=1, minsize=c2)
        body.grid_rowconfigure(0, weight=1, minsize=500)  # Increased minimum height for better visibility of columns 1 & 2
        # Ensure any additional rows don't take space
        for i in range(1, 10):
            body.grid_rowconfigure(i, weight=0, minsize=0)

        # LEFT: scrollable config column
        left_frame = tk.Frame(body, bg='#f0f0f0')
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 5))
        
        # Add scrollbar to left column
        left_canvas = tk.Canvas(left_frame, bg='#f0f0f0', highlightthickness=0)
        left_scrollbar = ttk.Scrollbar(left_frame, orient='vertical', command=left_canvas.yview)

        left_inner = tk.Frame(left_canvas, bg='#f0f0f0')
        left_inner.bind('<Configure>', lambda e: left_canvas.configure(scrollregion=left_canvas.bbox('all')))
        left_canvas_window = left_canvas.create_window((0, 0), window=left_inner, anchor='nw')
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        
        def _configure_left_canvas(event):
            try:
                left_canvas.itemconfig(left_canvas_window, width=event.width)
            except Exception:
                pass
            try:
                left_canvas.configure(scrollregion=left_canvas.bbox('all'))
            except Exception:
                pass
        
        left_canvas.bind('<Configure>', _configure_left_canvas)
        left_scrollbar.pack(side='right', fill='y')
        left_canvas.pack(side='left', fill='both', expand=True)
        
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), 'units')
        
        left_canvas.bind('<MouseWheel>', _on_left_mousewheel)
        left_inner.bind('<MouseWheel>', _on_left_mousewheel)

        # CENTER: scrollable settings column
        center_frame = tk.Frame(body, bg='#f0f0f0')
        center_frame.grid(row=0, column=1, sticky='nsew', padx=5)
        
        # Add scrollbar to center column
        center_canvas = tk.Canvas(center_frame, bg='#f0f0f0', highlightthickness=0)
        center_scrollbar = ttk.Scrollbar(center_frame, orient='vertical', command=center_canvas.yview)
        center_inner = tk.Frame(center_canvas, bg='#f0f0f0')
        center_inner.bind('<Configure>', lambda e: center_canvas.configure(scrollregion=center_canvas.bbox('all')))
        center_canvas_window = center_canvas.create_window((0, 0), window=center_inner, anchor='nw')
        center_canvas.configure(yscrollcommand=center_scrollbar.set)
        
        def _configure_center_canvas(event):
            try:
                center_canvas.itemconfig(center_canvas_window, width=event.width)
            except Exception:
                pass
            try:
                center_canvas.configure(scrollregion=center_canvas.bbox('all'))
            except Exception:
                pass
        
        center_canvas.bind('<Configure>', _configure_center_canvas)
        center_scrollbar.pack(side='right', fill='y')
        center_canvas.pack(side='left', fill='both', expand=True)
        
        def _on_center_mousewheel(event):
            center_canvas.yview_scroll(int(-1*(event.delta/120)), 'units')
        
        center_canvas.bind('<MouseWheel>', _on_center_mousewheel)
        center_inner.bind('<MouseWheel>', _on_center_mousewheel)

        # RIGHT: progress + log column
        right_frame = tk.Frame(body, bg='#f0f0f0')
        right_frame.grid(row=0, column=2, sticky='nsew', padx=(5, 0))

        # Setup panels (they will create their own internal scroll areas where needed)
        self.setup_viz_left_panel(left_inner)
        self.setup_viz_center_panel(center_inner)
        self.setup_viz_right_panel(right_frame)
        
        # Load saved visualization config (group definitions, mappings, etc.) at startup
        try:
            if hasattr(self, '_load_viz_config'):
                self._load_viz_config()
        except Exception as e:
            try:
                self.log_viz_message(f"⚠️ Failed loading saved visualization config: {e}")
            except Exception:
                pass

    def setup_viz_left_panel(self, parent):
        """Setup Panel 1: Data source, status, colors, output."""
        # Data source selection
        data_frame = ttk.LabelFrame(parent, text="Data Source", padding=10)
        data_frame.pack(fill='x', padx=5, pady=5)

        # Visualization mode radio (Metabolite vs Lipid vs Custom) at top
        mode_frame = ttk.Frame(data_frame)
        mode_frame.pack(fill='x', pady=(0,4))
        mode_label = ttk.Label(mode_frame, text="Mode:")
        mode_label.pack(side='left')
        self._create_tooltip(mode_label, "Metabolite: metabolite analysis (Pos_id/Neg_id sheets)\nLipid: lipid-specific analysis (Pos_lipid/Neg_lipid sheets)\nCustom: custom sheet names")
        self.viz_mode = tk.StringVar(value='metabolite')  # metabolite | lipid | custom
        for txt, val in [("Metabolite", 'metabolite'), ("Lipid", 'lipid'), ("Custom", 'custom')]:
            rb = ttk.Radiobutton(mode_frame, text=txt, value=val, variable=self.viz_mode, command=self.update_viz_data_status)
            rb.pack(side='left', padx=2)
        
        self.viz_data_source = tk.StringVar(value='memory')
        ttk.Radiobutton(data_frame, text="Use statistics from current session", 
                       variable=self.viz_data_source, value='memory',
                       command=self.update_viz_import_state).pack(anchor='w')
        ttk.Radiobutton(data_frame, text="Import statistical results Excel file", 
                       variable=self.viz_data_source, value='file',
                       command=self.update_viz_import_state).pack(anchor='w')
        
        # File import row
        file_row = ttk.Frame(data_frame)
        file_row.pack(fill='x', pady=2)
        
        self.viz_import_file = tk.StringVar()
        self.viz_file_entry = ttk.Entry(file_row, textvariable=self.viz_import_file, width=25)
        self.viz_file_entry.pack(side='left', expand=True, fill='x')
        self.viz_browse_btn = tk.Button(file_row, text="Browse", command=self.browse_viz_import_file,
                                        bg='#95a5a6', fg='white', font=('Arial', 9), relief='raised')
        self.viz_browse_btn.pack(side='right', padx=(2, 0))
        
        tk.Button(data_frame, text="Load & Analyze", command=self.load_viz_data,
                  bg='#16a085', fg='white', font=('Arial', 9, 'bold'), relief='raised').pack(pady=2)
        
        # Verify Columns button (for file import mode)
        self.verify_cols_btn = tk.Button(data_frame, text="🔍 Verify Columns", command=self.verify_viz_columns,
                  bg='#e67e22', fg='white', font=('Arial', 9, 'bold'), relief='raised', state='disabled')
        self.verify_cols_btn.pack(pady=2, fill='x')

        # Group & Stat Column Configuration
        config_frame = ttk.LabelFrame(parent, text="Configuration (Required)", padding=10)
        config_frame.pack(fill='x', padx=5, pady=5)
        
        # Step 1: Configure Groups (disabled until columns verified)
        group_row = ttk.Frame(config_frame)
        group_row.pack(fill='x', pady=2)
        self.configure_groups_btn = tk.Button(group_row, text="⚙️ Configure Groups", 
                  command=self._configure_groups_wrapper,
                  bg='#9b59b6', fg='white', font=('Arial', 9, 'bold'), relief='raised', state='disabled')
        self.configure_groups_btn.pack(side='left', fill='x', expand=True)
        self.groups_status_label = ttk.Label(group_row, text="⚠️", font=('Arial', 12))
        self.groups_status_label.pack(side='right', padx=5)
        
        # Step 2: Configure Stat Columns
        stat_row = ttk.Frame(config_frame)
        stat_row.pack(fill='x', pady=2)
        self.configure_stat_cols_btn = tk.Button(stat_row, text="Configure Stat Columns", 
                  command=self.configure_stat_columns,
                  bg='#2980b9', fg='white', font=('Arial', 9, 'bold'), relief='raised')
        self.configure_stat_cols_btn.pack(side='left', fill='x', expand=True)
        self.stat_cols_status_label = ttk.Label(stat_row, text="⚠️", font=('Arial', 12))
        self.stat_cols_status_label.pack(side='right', padx=5)
        
        # Configuration status message
        self.config_status = tk.StringVar(value="⚠️ Complete steps 1 & 2 before generating plots")
        ttk.Label(config_frame, textvariable=self.config_status, font=('Arial', 9, 'italic')).pack(anchor='w', pady=2)
        
        # Group color customization
        color_frame = ttk.LabelFrame(parent, text="Group Colors", padding=10)
        color_frame.pack(fill='x', padx=5, pady=5)
        
        # Group status
        self.viz_group_status = tk.StringVar(value="No groups configured")
        ttk.Label(color_frame, textvariable=self.viz_group_status, font=('Arial', 9)).pack(anchor='w', pady=2)
        
        # Color controls container
        self.viz_color_frame = ttk.Frame(color_frame)
        self.viz_color_frame.pack(fill='x', pady=2)
        
        # Output directory selection
        output_frame = ttk.LabelFrame(parent, text="Output Settings", padding=10)
        output_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(output_frame, text="Output Directory:").pack(anchor='w')
        output_row = ttk.Frame(output_frame)
        output_row.pack(fill='x', pady=2)
        
        self.viz_output_dir = tk.StringVar(value=os.path.join(os.getcwd(), 'visualizations'))
        ttk.Entry(output_row, textvariable=self.viz_output_dir, width=20).pack(side='left', expand=True, fill='x')
        tk.Button(output_row, text="Browse", command=self.choose_viz_output_dir,
                  bg='#95a5a6', fg='white', font=('Arial', 9), relief='raised').pack(side='right', padx=(2, 0))

        # Preferred group order entry
        order_frame = ttk.LabelFrame(parent, text="Preferred Group Order (comma separated)", padding=8)
        order_frame.pack(fill='x', padx=5, pady=5)
        self.viz_preferred_group_order = tk.StringVar(value="Group1,Group2,Group3")
        ttk.Entry(order_frame, textvariable=self.viz_preferred_group_order).pack(fill='x')
        ttk.Label(order_frame, text="Leave blank to use detected order.", font=('Arial',8,'italic')).pack(anchor='w', pady=(2,0))

        # Reset settings button (removed duplicate save/load - now only in global header)
        settings_vis_frame = ttk.Frame(parent)
        settings_vis_frame.pack(fill='x', padx=5, pady=(6,0))
        tk.Button(settings_vis_frame, text='Reset Visualization Defaults', command=self.reset_visualization_defaults,
                  bg='#e74c3c', fg='white', font=('Arial', 8), relief='raised').pack(side='right', padx=2)
        
        # Initialize the import state (must be after all buttons are created)
        self.update_viz_import_state()

    def setup_viz_center_panel(self, parent):
        """Setup Panel 2: Plot settings with internal tabs."""
        settings_frame = ttk.LabelFrame(parent, text="Plot Settings", padding=10)
        settings_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Action buttons at TOP of center panel for easy accessibility
        action_frame = ttk.Frame(settings_frame)
        action_frame.pack(fill='x', pady=(0, 10))
        
        # Main action buttons - store references for enabling/disabling
        btn_frame = ttk.Frame(action_frame)
        btn_frame.pack(fill='x')
        
        self.generate_selected_btn = tk.Button(btn_frame, text="🎯 Generate Selected Plots", 
                  command=self.generate_selected_viz_plots, state='disabled',
                  bg='#27ae60', fg='white', font=('Arial', 9, 'bold'), relief='raised')
        self.generate_selected_btn.pack(side='left', padx=5)
        self.generate_all_btn = tk.Button(btn_frame, text="📊 Generate All Plots", 
                  command=self.generate_all_viz_plots, state='disabled',
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold'), relief='raised')
        self.generate_all_btn.pack(side='left', padx=5)
        tk.Button(btn_frame, text="⏹️ Stop", 
              command=self.stop_viz_plots,
              bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'), relief='raised').pack(side='left', padx=5)
        tk.Button(btn_frame, text="📁 Open Folder", 
                  command=self.open_viz_output_folder,
                  bg='#3498db', fg='white', font=('Arial', 9, 'bold'), relief='raised').pack(side='right', padx=5)
        
        # Create notebook for plot type tabs
        self.viz_notebook = ttk.Notebook(settings_frame)
        self.viz_notebook.pack(fill='both', expand=True, pady=(0, 0))
        
        # Create tabs for each plot type
        self.create_pca_tab()
        self.create_volcano_tab() 
        self.create_venn_tab()
        self.create_boxplot_tab()
        self.create_bargraph_tab()
        self.create_heatmap_tab()
        self.create_roc_tab()
        self.create_metabolite_tab()

    def setup_viz_right_panel(self, parent):
        """Setup Panel 3: Progress and log."""
        # Progress section
        progress_frame = ttk.LabelFrame(parent, text="Progress", padding=10)
        progress_frame.pack(fill='x', padx=5, pady=5)
        
        self.viz_progress = ttk.Progressbar(progress_frame, mode='determinate')
        self.viz_progress.pack(fill='x', pady=2)
        
        self.viz_progress_label = ttk.Label(progress_frame, text="Ready")
        self.viz_progress_label.pack(anchor='w', pady=2)
        
        # Log section
        log_frame = ttk.LabelFrame(parent, text="Visualization Log", padding=10)
        log_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Slightly increased log height for better visibility
        self.viz_log = scrolledtext.ScrolledText(log_frame, height=18, width=35)
        self.viz_log.pack(fill='both', expand=True)
        
        # Clear log button
        tk.Button(log_frame, text="Clear Log", command=self.clear_viz_log,
              bg='#95a5a6', fg='white', font=('Arial', 8), relief='raised').pack(anchor='e', pady=2)

    def clear_viz_log(self):
        """Clear the visualization log."""
        self.viz_log.delete(1.0, tk.END)

    def create_pca_tab(self):
        """Create PCA settings tab."""
        pca_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(pca_tab, text="PCA")
        self.create_pca_panel(pca_tab)

    def create_volcano_tab(self):
        """Create Volcano plot settings tab."""
        volcano_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(volcano_tab, text="Volcano")
        self.create_volcano_panel(volcano_tab)

    def create_venn_tab(self):
        """Create Venn settings tab (placed after Volcano)."""
        venn_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(venn_tab, text="Venn")
        self.create_venn_panel(venn_tab)

    def create_boxplot_tab(self):
        """Create Boxplot settings tab."""
        boxplot_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(boxplot_tab, text="Boxplots")
        self.create_boxplot_panel(boxplot_tab)

    def create_bargraph_tab(self):
        """Create Bar graph settings tab."""
        bargraph_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(bargraph_tab, text="Bar Graphs")
        self.create_bargraph_panel(bargraph_tab)

    def create_heatmap_tab(self):
        """Create Heatmap settings tab."""
        heatmap_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(heatmap_tab, text="Heatmaps")
        self.create_heatmap_panel(heatmap_tab)

    def create_roc_tab(self):
        """Create ROC settings tab."""
        roc_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(roc_tab, text="ROC")
        self.create_roc_panel(roc_tab)

    def create_metabolite_tab(self):
        """Create Metabolite selection tab."""
        metabolite_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(metabolite_tab, text="Custom List")
        self.create_metabolite_selection_panel(metabolite_tab)

    def create_pca_panel(self, parent):
        """Create PCA parameter panel."""
        panel = ttk.LabelFrame(parent, text="PCA Analysis", padding=10)
        panel.pack(fill='x', padx=5, pady=2)
        
        # Enable checkbox
        pca_enable_cb = ttk.Checkbutton(panel, text="Generate PCA plots", 
                       variable=self.viz_selected['pca'])
        pca_enable_cb.pack(anchor='w')
        self._create_tooltip(pca_enable_cb, "Enable to generate Principal Component Analysis plots showing sample grouping and variance")
        
        # Comparison selection control (for pairwise PCA only) - MOVED TO TOP
        comp_select_frame = ttk.LabelFrame(panel, text="🔍 Pairwise Comparison Selection", padding=5)
        comp_select_frame.pack(fill='x', pady=4)
        
        info_label = ttk.Label(comp_select_frame, text="Control which pairwise comparisons to generate (only affects pairwise mode)", 
                              font=('Arial', 9, 'italic'), foreground='#666')
        info_label.pack(anchor='w', pady=(0, 5))
        
        self.pca_selected_comparisons = None  # Will store list of (g1, g2) tuples
        self.pca_comp_status_label = ttk.Label(comp_select_frame, text="All pairwise comparisons will be plotted")
        self.pca_comp_status_label.pack(anchor='w', pady=(0, 5))
        
        ttk.Button(comp_select_frame, text="Configure Pairwise Comparisons...", 
                   command=self._configure_pca_comparisons).pack(anchor='w')
        
        # Parameters frame
        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        
        # Components
        ttk.Label(params_frame, text="Max Components:").grid(row=0, column=0, sticky='w', padx=2)
        self.pca_components = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=2, to=50, textvariable=self.pca_components, width=10).grid(row=0, column=1, padx=2)
        
        # Enable 3D by default since checkboxes are removed
        self.pca_3d = tk.BooleanVar(value=True)
        
        self.pca_interactive = tk.BooleanVar(value=False)
        pca_interactive_cb = ttk.Checkbutton(params_frame, text="Interactive 3D (HTML)", variable=self.pca_interactive)
        pca_interactive_cb.grid(row=0, column=2, padx=10, sticky='w')
        self._create_tooltip(pca_interactive_cb, "Create interactive HTML 3D plots (can rotate/zoom in browser)")
        
        # Enable scree by default since checkbox is removed
        self.pca_scree = tk.BooleanVar(value=True)
        
        # Loadings
        self.pca_loadings = tk.BooleanVar(value=False)
        pca_loadings_cb = ttk.Checkbutton(params_frame, text="Export loadings (top", variable=self.pca_loadings)
        pca_loadings_cb.grid(row=1, column=0, sticky='w')
        self._create_tooltip(pca_loadings_cb, "Export top metabolites contributing most to principal components (feature importance)")
        self.pca_loadings_k = tk.IntVar(value=30)
        ttk.Spinbox(params_frame, from_=10, to=100, textvariable=self.pca_loadings_k, width=8).grid(row=1, column=1, padx=2)
        ttk.Label(params_frame, text=")").grid(row=1, column=2, sticky='w')

        # Figure size controls
        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        self.pca_fig_width = tk.DoubleVar(value=10.0)
        self.pca_fig_height = tk.DoubleVar(value=5.0)
        self.pca_fig_dpi = tk.IntVar(value=250)
        ttk.Label(size_frame, text="Width:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=30.0, increment=0.5, textvariable=self.pca_fig_width, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=30.0, increment=0.5, textvariable=self.pca_fig_height, width=6).grid(row=0, column=3, padx=2)
        
        # Title options
        title_frame = ttk.LabelFrame(panel, text="Title Options", padding=5)
        title_frame.pack(fill='x', pady=4)
        self.pca_skip_title = tk.BooleanVar(value=False)
        self.pca_show_legend = tk.BooleanVar(value=True)  # Show legend by default for all PCA plots
        ttk.Checkbutton(title_frame, text="Skip writing titles on plots", variable=self.pca_skip_title).grid(row=0, column=0, sticky='w', padx=5)
        ttk.Checkbutton(title_frame, text="Show legend on plots", variable=self.pca_show_legend).grid(row=0, column=1, sticky='w', padx=15)
        ttk.Label(size_frame, text="DPI:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.pca_fig_dpi, width=6).grid(row=0, column=5, padx=2)
        
        # Point size controls
        point_frame = ttk.LabelFrame(panel, text="Point Sizes")
        point_frame.pack(fill='x', pady=4)
        self.pca_point_size_2d = tk.DoubleVar(value=70)
        self.pca_point_size_3d = tk.DoubleVar(value=50)
        ttk.Label(point_frame, text="2D Point Size:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(point_frame, from_=5, to=100, increment=5, textvariable=self.pca_point_size_2d, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(point_frame, text="3D Point Size:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(point_frame, from_=5, to=100, increment=5, textvariable=self.pca_point_size_3d, width=6).grid(row=0, column=3, padx=2)
        
        # 3D viewing angles controls
        view_frame = ttk.LabelFrame(panel, text="3D View Angles (degrees)")
        view_frame.pack(fill='x', pady=4)
        self.pca_view_azim = tk.DoubleVar(value=-60)
        self.pca_view_elev = tk.DoubleVar(value=30)
        
        azim_label = ttk.Label(view_frame, text="Azimuth (rotation):")
        azim_label.grid(row=0, column=0, padx=2, sticky='w')
        azim_spin = ttk.Spinbox(view_frame, from_=-180, to=180, increment=5, textvariable=self.pca_view_azim, width=6)
        azim_spin.grid(row=0, column=1, padx=2)
        self._create_tooltip(azim_label, "Horizontal rotation angle for 3D plots (left/right perspective). Try -60 for good group separation.")
        self._create_tooltip(azim_spin, "Horizontal rotation angle for 3D plots (left/right perspective). Try -60 for good group separation.")
        
        elev_label = ttk.Label(view_frame, text="Elevation (tilt):")
        elev_label.grid(row=0, column=2, padx=8, sticky='w')
        elev_spin = ttk.Spinbox(view_frame, from_=0, to=90, increment=5, textvariable=self.pca_view_elev, width=6)
        elev_spin.grid(row=0, column=3, padx=2)
        self._create_tooltip(elev_label, "Vertical tilt angle for 3D plots (up/down view). Try 30 for good perspective.")
        self._create_tooltip(elev_spin, "Vertical tilt angle for 3D plots (up/down view). Try 30 for good perspective.")
        
        # Save options controls
        save_frame = ttk.LabelFrame(panel, text="💾 Save Options")
        save_frame.pack(fill='x', pady=4)
        self.pca_save_2d = tk.BooleanVar(value=True)
        self.pca_save_3d = tk.BooleanVar(value=False)
        self.pca_save_excel = tk.BooleanVar(value=False)
        ttk.Checkbutton(save_frame, text="Save 2D Plots", variable=self.pca_save_2d).grid(row=0, column=0, padx=5, sticky='w')
        ttk.Checkbutton(save_frame, text="Save 3D Plots", variable=self.pca_save_3d).grid(row=0, column=1, padx=5, sticky='w')
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.pca_save_excel).grid(row=0, column=2, padx=5, sticky='w')
        
        # Output format selection
        format_frame = ttk.LabelFrame(panel, text="Output Format", padding=5)
        format_frame.pack(fill='x', pady=4)
        self.pca_output_format = tk.StringVar(value='png')
        ttk.Radiobutton(format_frame, text="PNG (Raster)", variable=self.pca_output_format, value='png').pack(anchor='w')
        ttk.Radiobutton(format_frame, text="SVG (Vector)", variable=self.pca_output_format, value='svg').pack(anchor='w')
        
        # Font size controls
        font_frame = ttk.LabelFrame(panel, text="🔤 Font Sizes (pt)")
        font_frame.pack(fill='x', pady=4)
        self.pca_xlabel_fontsize = tk.IntVar(value=14)
        self.pca_ylabel_fontsize = tk.IntVar(value=14)
        self.pca_title_fontsize = tk.IntVar(value=14)
        self.pca_tick_fontsize = tk.IntVar(value=12)
        self.pca_legend_fontsize = tk.IntVar(value=12)
        
        # Add trace to sync ylabel with xlabel for combined control
        def _sync_pca_axis_fonts(*args):
            self.pca_ylabel_fontsize.set(self.pca_xlabel_fontsize.get())
        self.pca_xlabel_fontsize.trace_add('write', _sync_pca_axis_fonts)
        
        ttk.Label(font_frame, text="Axis Labels (PC1, PC2, PC3):").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.pca_xlabel_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(font_frame, text="Title:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.pca_title_fontsize, width=5).grid(row=0, column=3, padx=2)
        ttk.Label(font_frame, text="Tick Labels:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.pca_tick_fontsize, width=5).grid(row=0, column=5, padx=2)
        
        ttk.Label(font_frame, text="Legend:").grid(row=1, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=24, textvariable=self.pca_legend_fontsize, width=5).grid(row=1, column=1, padx=2)

    def create_volcano_panel(self, parent):
        """Create Volcano plot parameter panel."""
        panel = ttk.LabelFrame(parent, text="Volcano Plots", padding=10)
        panel.pack(fill='x', padx=5, pady=2)
        
        # Enable checkbox
        ttk.Checkbutton(panel, text="Generate volcano plots", 
                       variable=self.viz_selected['volcano']).pack(anchor='w')
        
        # Parameters
        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        
        ttk.Label(params_frame, text="P-value threshold:").grid(row=0, column=0, sticky='w', padx=2)
        self.volcano_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(params_frame, textvariable=self.volcano_p_thresh, width=10).grid(row=0, column=1, padx=2)
        
        ttk.Label(params_frame, text="Fold-change threshold:").grid(row=0, column=2, sticky='w', padx=10)
        self.volcano_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(params_frame, textvariable=self.volcano_fc_thresh, width=10).grid(row=0, column=3, padx=2)
        # Option to skip fold-change cutoff and use p-value only
        self.volcano_skip_fc = tk.BooleanVar(value=True)
        def _on_skip_fc_change():
            try:
                state = 'disabled' if self.volcano_skip_fc.get() else 'normal'
                # Find the entry widget (it's the grid child at row=0, column=3 in params_frame)
                for w in params_frame.grid_slaves(row=0, column=3):
                    try:
                        w.configure(state=state)
                    except Exception:
                        pass
            except Exception:
                pass
        ttk.Checkbutton(params_frame, text="Skip FC cutoff", variable=self.volcano_skip_fc, command=_on_skip_fc_change).grid(row=1, column=2, columnspan=2, sticky='w', padx=4)
        _on_skip_fc_change()
        
        # Annotation options
        self.volcano_annotate = tk.BooleanVar(value=False)
        volcano_annotate_cb = ttk.Checkbutton(params_frame, text="Annotate top", variable=self.volcano_annotate)
        volcano_annotate_cb.grid(row=1, column=0, sticky='w')
        self._create_tooltip(volcano_annotate_cb, "Add metabolite name labels for the top N most significant points on volcano plot")
        self.volcano_top_n = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=0, to=50, textvariable=self.volcano_top_n, width=8).grid(row=1, column=1, padx=2)
        ttk.Label(params_frame, text="Label size:").grid(row=2, column=0, padx=2, sticky='e')
        self.volcano_annot_fontsize = tk.IntVar(value=8)
        ttk.Spinbox(params_frame, from_=4, to=20, textvariable=self.volcano_annot_fontsize, width=8).grid(row=2, column=1, padx=2)

        # Figure size controls
        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        # Match metabolites_visualization.VolcanoParams defaults (5 x 4 inches)
        self.volcano_fig_width = tk.DoubleVar(value=6.0)
        self.volcano_fig_height = tk.DoubleVar(value=5.0)
        self.volcano_fig_dpi = tk.IntVar(value=250)
        ttk.Label(size_frame, text="Width:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=30.0, increment=0.5, textvariable=self.volcano_fig_width, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=30.0, increment=0.5, textvariable=self.volcano_fig_height, width=6).grid(row=0, column=3, padx=2)
        ttk.Label(size_frame, text="DPI:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.volcano_fig_dpi, width=6).grid(row=0, column=5, padx=2)
        
        # Point size controls
        point_frame = ttk.LabelFrame(panel, text="Point Sizes")
        point_frame.pack(fill='x', pady=4)
        self.volcano_point_size_sig = tk.DoubleVar(value=35)
        self.volcano_point_size_nonsig = tk.DoubleVar(value=35)
        ttk.Label(point_frame, text="Significant Points:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(point_frame, from_=5, to=100, increment=5, textvariable=self.volcano_point_size_sig, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(point_frame, text="Non-Significant:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(point_frame, from_=5, to=100, increment=5, textvariable=self.volcano_point_size_nonsig, width=6).grid(row=0, column=3, padx=2)
        
        # Font size controls
        font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        font_frame.pack(fill='x', pady=4)
        self.volcano_xlabel_fontsize = tk.IntVar(value=14)
        self.volcano_ylabel_fontsize = tk.IntVar(value=14)
        self.volcano_title_fontsize = tk.IntVar(value=14)
        self.volcano_tick_fontsize = tk.IntVar(value=14)
        self.volcano_count_fontsize = tk.IntVar(value=16)
        self.volcano_total_fontsize = tk.IntVar(value=16)
        self.volcano_legend_fontsize = tk.IntVar(value=12)
        
        ttk.Label(font_frame, text="Axis Labels:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.volcano_xlabel_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(font_frame, text="Title:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.volcano_title_fontsize, width=5).grid(row=0, column=3, padx=2)
        ttk.Label(font_frame, text="Tick Labels:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=30, textvariable=self.volcano_tick_fontsize, width=5).grid(row=0, column=5, padx=2)
        
        ttk.Label(font_frame, text="Count Text:").grid(row=1, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=24, textvariable=self.volcano_count_fontsize, width=5).grid(row=1, column=1, padx=2)
        ttk.Label(font_frame, text="Total Text:").grid(row=1, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=24, textvariable=self.volcano_total_fontsize, width=5).grid(row=1, column=3, padx=2)
        ttk.Label(font_frame, text="Legend:").grid(row=1, column=4, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=6, to=24, textvariable=self.volcano_legend_fontsize, width=5).grid(row=1, column=5, padx=2)
        
        # Count background style
        count_bg_frame = ttk.LabelFrame(panel, text="Count Background Style", padding=5)
        count_bg_frame.pack(fill='x', pady=4)
        self.volcano_count_background = tk.StringVar(value='colored')
        ttk.Radiobutton(count_bg_frame, text="Colored Background", variable=self.volcano_count_background, value='colored').pack(anchor='w')
        ttk.Radiobutton(count_bg_frame, text="Transparent Background", variable=self.volcano_count_background, value='transparent').pack(anchor='w')
        
        # Output format selection
        format_frame = ttk.LabelFrame(panel, text="Output Format", padding=5)
        format_frame.pack(fill='x', pady=4)
        self.volcano_output_format = tk.StringVar(value='png')
        ttk.Radiobutton(format_frame, text="PNG (Raster)", variable=self.volcano_output_format, value='png').pack(anchor='w')
        ttk.Radiobutton(format_frame, text="SVG (Vector)", variable=self.volcano_output_format, value='svg').pack(anchor='w')
        
        # Title options
        title_frame = ttk.LabelFrame(panel, text="Title Options", padding=5)
        title_frame.pack(fill='x', pady=4)
        self.volcano_skip_title = tk.BooleanVar(value=False)
        ttk.Checkbutton(title_frame, text="Skip writing titles on plots", variable=self.volcano_skip_title).pack(anchor='w')
        
        # Excel export control
        save_frame = ttk.LabelFrame(panel, text="💾 Save Options", padding=5)
        save_frame.pack(fill='x', pady=4)
        self.volcano_save_excel = tk.BooleanVar(value=False)
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.volcano_save_excel).pack(anchor='w')
        
        # Comparison selection control
        comp_select_frame = ttk.LabelFrame(panel, text="🔍 Comparison Selection", padding=5)
        comp_select_frame.pack(fill='x', pady=4)
        
        self.volcano_selected_comparisons = None  # Will store list of (g1, g2) tuples
        self.volcano_comp_status_label = ttk.Label(comp_select_frame, text="All comparisons will be plotted")
        self.volcano_comp_status_label.pack(anchor='w', pady=(0, 5))
        
        ttk.Button(comp_select_frame, text="Configure Comparisons...", 
                   command=self._configure_volcano_comparisons).pack(anchor='w')
        
        # end volcano panel

    def create_venn_panel(self, parent):
        """Create Venn configuration panel."""
        panel = ttk.LabelFrame(parent, text="Venn Diagrams", padding=10)
        panel.pack(fill='both', expand=True, padx=5, pady=2)

        # Display canonicalization helper: keep internal pair tuples as-is but
        # always display the base group first (user expects 'Control' first).
        BASE_VENN_GROUP = 'control'
        def _canon_pair(a, b):
            """Return a tuple (first, second) for display where the BASE_VENN_GROUP
            appears first if either a or b matches it (case-insensitive). Otherwise
            preserve original order."""
            try:
                if str(a).strip().lower() == BASE_VENN_GROUP:
                    return a, b
                if str(b).strip().lower() == BASE_VENN_GROUP:
                    return b, a
            except Exception:
                pass
            return a, b

        # Enable checkbox
        if 'venn' not in self.viz_selected:
            self.viz_selected['venn'] = tk.BooleanVar(value=False)
        ttk.Checkbutton(panel, text="Generate Venn diagrams", variable=self.viz_selected['venn']).pack(anchor='w')

        # Add quick setup section at the top
        quick_frame = ttk.LabelFrame(panel, text="Quick Setup (Recommended)", padding=8)
        quick_frame.pack(fill='x', pady=(5,10))
        
        # info_label = ttk.Label(quick_frame, text="Generate Venn diagrams comparing all your groups (like in R):", 
        #                        foreground='blue')
        # info_label.pack(anchor='w', pady=(0,4))
        
        def _auto_generate_venn():
            """Automatically create Venn diagrams with pairwise comparisons. Supports multiple Venns for 4+ groups."""
            try:
                groups = self.get_viz_groups()
                if not groups or len(groups) < 2:
                    messagebox.showwarning("Auto-Generate Venn", 
                                         "Need at least 2 groups configured.\n\n"
                                         "Please configure groups in the Groups tab first.")
                    return
                
                from itertools import combinations
                # Generate all pairwise comparisons in canonical order (baseline first)
                all_pairs = [_canon_pair(a, b) for a, b in combinations(groups, 2)]
                
                # If >3 groups, use multi-Venn dialog
                if len(groups) > 3:
                    venn_configs = self._show_comparison_selector(all_pairs, groups)
                    if not venn_configs:  # User cancelled or no selection
                        return
                    
                    # Extract all unique pairs from all venns
                    all_selected_pairs = []
                    pair_set = set()
                    for config in venn_configs:
                        for pair in config['pairs']:
                            if pair not in pair_set:
                                all_selected_pairs.append(pair)
                                pair_set.add(pair)
                    
                    self.venn_pairs_list = all_selected_pairs
                    
                    # Map each venn's pairs to indices in the consolidated pairs list
                    self.venn_specs = []
                    for config in venn_configs:
                        indices = []
                        for pair in config['pairs']:
                            try:
                                idx = self.venn_pairs_list.index(pair) + 1
                                indices.append(idx)
                            except ValueError:
                                pass
                        self.venn_specs.append({
                            'name': config['name'],
                            'indices': indices
                        })
                else:
                    # Simple case: 2-3 groups = single Venn
                    self.venn_pairs_list = list(all_pairs)
                    venn_name = f"All_Groups_Venn"
                    self.venn_specs = [{'name': venn_name, 'indices': list(range(1, len(all_pairs)+1))}]
                
                # Update displays (display-only canonical order: Control first)
                self.venn_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_pairs_list]))
                self.venn_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_specs]))
                
                # ALSO auto-generate All Molecules config with same structure
                self.venn_allmol_pairs_list = list(self.venn_pairs_list)  # Copy same pairs
                self.venn_allmol_specs = [{'name': spec['name'], 'indices': list(spec['indices'])} for spec in self.venn_specs]  # Copy specs
                
                # Update All Molecules displays
                self.venn_allmol_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_allmol_pairs_list]))
                self.venn_allmol_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_allmol_specs]))
                
                # Update the pairs index display (show canonical display order)
                if hasattr(self, 'venn_pairs_index_label'):
                    lines = [f"{i+1}:{_canon_pair(a,b)[0]}|{_canon_pair(a,b)[1]}" for i,(a,b) in enumerate(self.venn_pairs_list)]
                    self.venn_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
                
                # Update All Molecules pairs index display
                if hasattr(self, 'venn_allmol_pairs_index_label'):
                    lines = [f"{i+1}:{_canon_pair(a,b)[0]}|{_canon_pair(a,b)[1]}" for i,(a,b) in enumerate(self.venn_allmol_pairs_list)]
                    self.venn_allmol_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
                
                # Update the summary display
                if hasattr(self, '_update_venn_summary'):
                    self._update_venn_summary()
                
                msg = f"✅ Auto-configured {len(self.venn_specs)} Venn diagram(s)!\n\n"
                msg += f"Groups detected: {', '.join(groups)}\n"
                msg += f"Total comparisons: {len(self.venn_pairs_list)}\n"
                msg += f"Venn diagrams created: {len(self.venn_specs)}\n"
                msg += f"All Molecules Venns: {len(self.venn_allmol_specs)} (same structure)\n\n"
                msg += "Click 'Generate Plots' to visualize."
                
                self.log_viz_message(f"Auto-generated {len(self.venn_specs)} Venn(s) + {len(self.venn_allmol_specs)} All Molecules: {len(self.venn_pairs_list)} comparisons from {len(groups)} groups")
                messagebox.showinfo("Auto-Generate Complete", msg)
                
            except Exception as e:
                messagebox.showerror("Auto-Generate Error", f"Failed to auto-generate Venn:\n{e}")
        
        auto_btn = ttk.Button(quick_frame, text="🔄 Auto-Generate Venn from All Groups", 
                             command=_auto_generate_venn, width=40)
        auto_btn.pack(pady=4)
        
        # Add edit button for modifying existing comparisons
        def _edit_comparisons():
            """Edit existing comparisons - remove unwanted ones."""
            groups = self.get_viz_groups()
            if not groups or len(groups) < 2:
                messagebox.showwarning("Edit Comparisons", 
                                      "Need at least 2 groups configured.\n\n"
                                      "Please configure groups in the Groups tab first.")
                return
            
            # Generate fresh pairs from current groups
            from itertools import combinations
            all_pairs = [_canon_pair(a, b) for a, b in combinations(groups, 2)]
            
            if not all_pairs:
                messagebox.showinfo("Edit Comparisons", 
                                  "No comparisons available with current groups.")
                return
            
            # Pass current selection to pre-check existing pairs
            venn_configs = self._show_comparison_selector(all_pairs, groups)
            
            if venn_configs is not None:  # User clicked OK (could be empty list)
                # Extract all unique pairs from the configs
                all_selected_pairs = []
                pair_set = set()
                for config in venn_configs:
                    for pair in config['pairs']:
                        if pair not in pair_set:
                            all_selected_pairs.append(pair)
                            pair_set.add(pair)
                
                if not all_selected_pairs:
                    response = messagebox.askyesno(
                        "Remove All Comparisons?",
                        "This will remove all comparisons.\n\nContinue?",
                        icon='warning'
                    )
                    if not response:
                        return
                
                # Update pairs list
                self.venn_pairs_list = all_selected_pairs
                self.venn_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_pairs_list]))
                
                # Update Venn specs from the configs
                self.venn_specs = []
                for config in venn_configs:
                    indices = []
                    for pair in config['pairs']:
                        try:
                            idx = self.venn_pairs_list.index(pair) + 1
                            indices.append(idx)
                        except ValueError:
                            pass
                    if indices:
                        self.venn_specs.append({
                            'name': config['name'],
                            'indices': indices
                        })
                
                if not self.venn_specs:
                    self.venn_specs_var.set("")
                else:
                    self.venn_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_specs]))
                
                # Update displays
                if hasattr(self, 'venn_pairs_index_label') and all_selected_pairs:
                    lines = [f"{i+1}:{a}|{b}" for i,(a,b) in enumerate(self.venn_pairs_list)]
                    self.venn_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
                elif hasattr(self, 'venn_pairs_index_label'):
                    self.venn_pairs_index_label.configure(text='Pairs Index: (none)')
                
                if hasattr(self, '_update_venn_summary'):
                    self._update_venn_summary()
                
                self.log_viz_message(f"Edited comparisons: {len(all_selected_pairs)} remaining")
                messagebox.showinfo("Edit Complete", 
                                  f"Updated to {len(all_selected_pairs)} comparison(s)")
        
        edit_btn = ttk.Button(quick_frame, text="✏️ Edit/Remove Comparisons", 
                             command=_edit_comparisons, width=40)
        edit_btn.pack(pady=(2,4))
        
        hint_label = ttk.Label(quick_frame, 
                              text="(Auto-Generate creates all pairs, then Edit to select which ones to keep)",
                              font=('TkDefaultFont', 8), foreground='gray')
        hint_label.pack(anchor='w')

        # Thresholds aligned with volcano defaults
        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        ttk.Label(params_frame, text="P-value threshold:").grid(row=0, column=0, sticky='w', padx=2)
        self.venn_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(params_frame, textvariable=self.venn_p_thresh, width=10).grid(row=0, column=1, padx=2)
        ttk.Label(params_frame, text="Fold-change threshold:").grid(row=0, column=2, sticky='w', padx=10)
        self.venn_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(params_frame, textvariable=self.venn_fc_thresh, width=10).grid(row=0, column=3, padx=2)
        self.venn_skip_fc = tk.BooleanVar(value=True)
        def _on_venn_skip_fc_change():
            try:
                state = 'disabled' if self.venn_skip_fc.get() else 'normal'
                for w in params_frame.grid_slaves(row=0, column=3):
                    try: w.configure(state=state)
                    except Exception: pass
            except Exception:
                pass
        ttk.Checkbutton(params_frame, text="Skip FC cutoff", variable=self.venn_skip_fc, command=_on_venn_skip_fc_change).grid(row=1, column=2, columnspan=2, sticky='w', padx=4)
        _on_venn_skip_fc_change()
        
        # ===================================================================
        # All Molecule Comparison Across Sample Groups
        # ===================================================================
        all_mol_frame = ttk.LabelFrame(panel, text="All Molecule Comparison Across Sample Groups", padding=8)
        all_mol_frame.pack(fill='x', pady=(10,4))
        
        ttk.Label(all_mol_frame, text="Generate Venn diagrams without statistical filters (p-value/FC).",
                 font=('TkDefaultFont', 9), foreground='#444').pack(anchor='w', pady=(0,6))
        
        # Enable/disable checkbox
        self.venn_generate_all_molecules = tk.BooleanVar(value=True)
        ttk.Checkbutton(all_mol_frame, text="Generate All Molecules Venn diagrams", 
                       variable=self.venn_generate_all_molecules).pack(anchor='w', pady=2)
        
        ttk.Label(all_mol_frame, text="Minimum non-zero/non-NaN samples required per group to consider metabolite 'present':",
                 font=('TkDefaultFont', 8), foreground='#555').pack(anchor='w', pady=(4,2))
        
        min_presence_frame = ttk.Frame(all_mol_frame)
        min_presence_frame.pack(fill='x', pady=2)
        
        self.venn_allmol_min_presence_type = tk.StringVar(value='count')
        ttk.Radiobutton(min_presence_frame, text="Count:", variable=self.venn_allmol_min_presence_type, 
                       value='count').grid(row=0, column=0, sticky='w', padx=2)
        self.venn_allmol_min_presence_count = tk.IntVar(value=3)
        ttk.Spinbox(min_presence_frame, from_=1, to=50, textvariable=self.venn_allmol_min_presence_count, 
                   width=8).grid(row=0, column=1, padx=2)
        ttk.Label(min_presence_frame, text="samples").grid(row=0, column=2, sticky='w', padx=2)
        
        ttk.Radiobutton(min_presence_frame, text="Percentage:", variable=self.venn_allmol_min_presence_type, 
                       value='percentage').grid(row=1, column=0, sticky='w', padx=2)
        self.venn_allmol_min_presence_percent = tk.DoubleVar(value=50.0)
        ttk.Spinbox(min_presence_frame, from_=0, to=100, increment=5, textvariable=self.venn_allmol_min_presence_percent, 
                   width=8).grid(row=1, column=1, padx=2)
        ttk.Label(min_presence_frame, text="%").grid(row=1, column=2, sticky='w', padx=2)
        
        # All Molecules Comparison Configuration
        allmol_config_frame = ttk.LabelFrame(all_mol_frame, text="Configure Comparisons for All Molecules Venn", padding=6)
        allmol_config_frame.pack(fill='x', pady=(8,4))
        
        ttk.Label(allmol_config_frame, text="Add group pairs for All Molecules Venn analysis:").pack(anchor='w', pady=(0,4))
        
        allmol_pairs_row = ttk.Frame(allmol_config_frame)
        allmol_pairs_row.pack(fill='x', pady=2)
        ttk.Label(allmol_pairs_row, text="Comparison (A|B):").pack(side='left')
        self.venn_allmol_pair_entry = tk.StringVar()
        ttk.Entry(allmol_pairs_row, textvariable=self.venn_allmol_pair_entry, width=30).pack(side='left', padx=4)
        
        if not hasattr(self, 'venn_allmol_pairs_list'):
            self.venn_allmol_pairs_list = []  # list[tuple[str,str]]
        if not hasattr(self, 'venn_allmol_specs'):
            self.venn_allmol_specs = []  # list[dict{name, indices:list[int]}]
        
        def _add_allmol_pair():
            txt = self.venn_allmol_pair_entry.get().strip()
            if '|' in txt:
                a, b = [t.strip() for t in txt.split('|', 1)]
                if a and b:
                    self.venn_allmol_pairs_list.append((a, b))
                    self.venn_allmol_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_allmol_pairs_list]))
                    self.venn_allmol_pair_entry.set('')
                    _refresh_allmol_pairs_index_label()
        
        def _remove_allmol_pair():
            try:
                idx_str = simpledialog.askstring("Remove Pair", "Enter pair index to remove:")
                if idx_str:
                    idx = int(idx_str) - 1
                    if 0 <= idx < len(self.venn_allmol_pairs_list):
                        removed = self.venn_allmol_pairs_list.pop(idx)
                        self.venn_allmol_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_allmol_pairs_list]))
                        _refresh_allmol_pairs_index_label()
                        self.log_viz_message(f"Removed All Molecules pair {idx+1}: {removed[0]}|{removed[1]}")
                    else:
                        messagebox.showwarning("Invalid Index", f"Index {idx+1} is out of range.")
            except ValueError:
                messagebox.showerror("Invalid Input", "Please enter a valid number.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to remove pair: {e}")
        
        def _quick_import_allmol_pairs():
            """Auto-generate all pairwise combinations for All Molecules Venn."""
            try:
                groups = self.get_viz_groups()
                if not groups or len(groups) < 2:
                    messagebox.showwarning("Quick Import", "Need at least 2 groups configured.")
                    return
                from itertools import combinations
                new_pairs = list(combinations(groups, 2))
                if new_pairs:
                    self.venn_allmol_pairs_list.extend(new_pairs)
                    self.venn_allmol_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_allmol_pairs_list]))
                    _refresh_allmol_pairs_index_label()
                    self.log_viz_message(f"All Molecules Quick import: added {len(new_pairs)} pairs from {len(groups)} groups")
                    messagebox.showinfo("Quick Import", f"Added {len(new_pairs)} pairwise pairs.")
            except Exception as e:
                messagebox.showerror("Quick Import Error", f"Failed: {e}")
        
        allmol_buttons_line = ttk.Frame(allmol_config_frame)
        allmol_buttons_line.pack(fill='x', pady=(4,2))
        ttk.Button(allmol_buttons_line, text="Add", command=_add_allmol_pair).pack(side='left', padx=2)
        ttk.Button(allmol_buttons_line, text="Remove", command=_remove_allmol_pair).pack(side='left', padx=2)
        ttk.Button(allmol_buttons_line, text="Quick Import", command=_quick_import_allmol_pairs).pack(side='left', padx=2)
        
        def _edit_allmol_comparisons():
            """Launch multi-selection dialog for All Molecules comparisons."""
            try:
                groups = self.get_viz_groups()
                if not groups or len(groups) < 2:
                    messagebox.showwarning("Edit Comparisons", "Need at least 2 groups configured.")
                    return
                
                from itertools import combinations
                all_pairs = [_canon_pair(a, b) for a, b in combinations(groups, 2)]
                
                if len(groups) > 3:
                    venn_configs = self._show_comparison_selector(all_pairs, groups, title="All Molecules: Create Venn Diagrams")
                    if not venn_configs:
                        return
                    
                    # Extract all unique pairs
                    all_selected_pairs = []
                    pair_set = set()
                    for config in venn_configs:
                        for pair in config['pairs']:
                            if pair not in pair_set:
                                all_selected_pairs.append(pair)
                                pair_set.add(pair)
                    
                    self.venn_allmol_pairs_list = all_selected_pairs
                    
                    # Map pairs to indices
                    self.venn_allmol_specs = []
                    for config in venn_configs:
                        indices = []
                        for pair in config['pairs']:
                            try:
                                idx = self.venn_allmol_pairs_list.index(pair) + 1
                                indices.append(idx)
                            except ValueError:
                                pass
                        self.venn_allmol_specs.append({'name': config['name'], 'indices': indices})
                    
                    # Update displays
                    self.venn_allmol_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_allmol_pairs_list]))
                    self.venn_allmol_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_allmol_specs]))
                    
                    if hasattr(self, 'venn_allmol_pairs_index_label'):
                        lines = [f"{i+1}:{_canon_pair(a,b)[0]}|{_canon_pair(a,b)[1]}" for i,(a,b) in enumerate(self.venn_allmol_pairs_list)]
                        self.venn_allmol_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
                    
                    self.log_viz_message(f"All Molecules: Configured {len(self.venn_allmol_specs)} Venn(s) with {len(self.venn_allmol_pairs_list)} comparisons")
                else:
                    messagebox.showinfo("Edit Comparisons", "Use this for 4+ groups. For 2-3 groups, use 'Quick Import'.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to edit comparisons: {e}")
        
        ttk.Button(allmol_buttons_line, text="✏️ Edit Comparisons (Multi-Select)", command=_edit_allmol_comparisons).pack(side='left', padx=6)
        
        self.venn_allmol_pairs_var = tk.StringVar(value="")
        ttk.Label(allmol_config_frame, textvariable=self.venn_allmol_pairs_var, wraplength=500).pack(anchor='w', pady=(2,6))
        
        # All Molecules Venn specs builder
        allmol_venn_row = ttk.Frame(allmol_config_frame)
        allmol_venn_row.pack(fill='x', pady=2)
        ttk.Label(allmol_venn_row, text="Venn Name:").pack(side='left')
        self.venn_allmol_name_entry = tk.StringVar()
        ttk.Entry(allmol_venn_row, textvariable=self.venn_allmol_name_entry, width=24).pack(side='left', padx=4)
        ttk.Label(allmol_venn_row, text="Comparisons to include (indices):").pack(side='left', padx=6)
        self.venn_allmol_indices_entry = tk.StringVar()
        ttk.Entry(allmol_venn_row, textvariable=self.venn_allmol_indices_entry, width=24).pack(side='left')
        
        def _add_allmol_venn_spec():
            name = self.venn_allmol_name_entry.get().strip() or f"AllMolVenn{len(self.venn_allmol_specs)+1}"
            idxs = []
            raw = self.venn_allmol_indices_entry.get().strip()
            if raw:
                for x in raw.split(','):
                    try:
                        k = int(x.strip())
                        if k >= 1: idxs.append(k)
                    except Exception:
                        pass
            self.venn_allmol_specs.append({'name': name, 'indices': idxs})
            self.venn_allmol_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_allmol_specs]))
            self.venn_allmol_name_entry.set('')
            self.venn_allmol_indices_entry.set('')
        
        def _remove_allmol_venn_spec():
            try:
                venn_name = simpledialog.askstring("Remove Venn", "Enter All Molecules Venn name to remove:")
                if venn_name:
                    found = False
                    for i, spec in enumerate(self.venn_allmol_specs):
                        if spec['name'] == venn_name:
                            removed = self.venn_allmol_specs.pop(i)
                            self.venn_allmol_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_allmol_specs]))
                            self.log_viz_message(f"Removed All Molecules Venn: {removed['name']}")
                            found = True
                            break
                    if not found:
                        messagebox.showwarning("Not Found", f"Venn '{venn_name}' not found.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to remove Venn: {e}")
        
        allmol_venn_buttons_row = ttk.Frame(allmol_config_frame)
        allmol_venn_buttons_row.pack(fill='x', pady=(2,2))
        ttk.Button(allmol_venn_buttons_row, text="Add Venn", command=_add_allmol_venn_spec).pack(side='left', padx=2)
        ttk.Button(allmol_venn_buttons_row, text="Remove Venn", command=_remove_allmol_venn_spec).pack(side='left', padx=2)
        
        self.venn_allmol_specs_var = tk.StringVar(value="")
        ttk.Label(allmol_config_frame, text="Defined All Molecules Venns:").pack(anchor='w')
        ttk.Label(allmol_config_frame, textvariable=self.venn_allmol_specs_var, wraplength=500).pack(anchor='w')
        
        def _refresh_allmol_pairs_index_label():
            if self.venn_allmol_pairs_list:
                lines = [f"{i+1}:{_canon_pair(a,b)[0]}|{_canon_pair(a,b)[1]}" for i,(a,b) in enumerate(self.venn_allmol_pairs_list)]
                self.venn_allmol_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
            else:
                self.venn_allmol_pairs_index_label.configure(text='Pairs Index: (none)')
        
        self.venn_allmol_pairs_index_label = ttk.Label(allmol_config_frame, text='Pairs Index: (none)')
        self.venn_allmol_pairs_index_label.pack(anchor='w', pady=(4,6))

        # Font size controls
        venn_font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        venn_font_frame.pack(fill='x', pady=4)
        self.venn_number_fontsize = tk.IntVar(value=16)  # Bigger default as requested
        self.venn_label_fontsize = tk.IntVar(value=11)
        
        ttk.Label(venn_font_frame, text="Venn Numbers (bold):").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(venn_font_frame, from_=8, to=36, textvariable=self.venn_number_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(venn_font_frame, text="Venn Labels:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(venn_font_frame, from_=6, to=30, textvariable=self.venn_label_fontsize, width=5).grid(row=0, column=3, padx=2)

        # Output format selection
        format_frame = ttk.LabelFrame(panel, text="Output Format", padding=5)
        format_frame.pack(fill='x', pady=4)
        self.venn_output_format = tk.StringVar(value='png')
        ttk.Radiobutton(format_frame, text="PNG (Raster)", variable=self.venn_output_format, value='png').pack(anchor='w')
        ttk.Radiobutton(format_frame, text="SVG (Vector)", variable=self.venn_output_format, value='svg').pack(anchor='w')

        # Advanced/Manual configuration (collapsible)
        advanced_frame = ttk.LabelFrame(panel, text="Advanced: Manual Configuration (Optional)", padding=8)
        advanced_frame.pack(fill='both', expand=True, pady=6)
        
        adv_info = ttk.Label(advanced_frame, 
                            text="For custom Venn configurations, manually define comparisons below.\n"
                                 "Most users can skip this and use Auto-Generate above.",
                            foreground='gray', font=('TkDefaultFont', 8))
        adv_info.pack(anchor='w', pady=(0,4))

        # Comparisons builder
        ttk.Label(advanced_frame, text="Available group pairs come from current groups.").pack(anchor='w')

        pairs_row = ttk.Frame(advanced_frame)
        pairs_row.pack(fill='x', pady=2)
        ttk.Label(pairs_row, text="Comparison (A|B):").pack(side='left')
        self.venn_pair_entry = tk.StringVar()
        ttk.Entry(pairs_row, textvariable=self.venn_pair_entry, width=30).pack(side='left', padx=4)
        if not hasattr(self, 'venn_pairs_list'):
            self.venn_pairs_list = []  # list[tuple[str,str]]
        def _add_pair():
            txt = self.venn_pair_entry.get().strip()
            if '|' in txt:
                a, b = [t.strip() for t in txt.split('|', 1)]
                if a and b:
                    self.venn_pairs_list.append((a, b))
                    self.venn_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_pairs_list]))
                    self.venn_pair_entry.set('')
                    _refresh_pairs_index_label()
        def _remove_pair():
            try:
                idx_str = simpledialog.askstring("Remove Pair", "Enter pair index to remove:")
                if idx_str:
                    idx = int(idx_str) - 1
                    if 0 <= idx < len(self.venn_pairs_list):
                        removed = self.venn_pairs_list.pop(idx)
                        self.venn_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_pairs_list]))
                        _refresh_pairs_index_label()
                        self.log_viz_message(f"Removed pair {idx+1}: {removed[0]}|{removed[1]}")
                    else:
                        messagebox.showwarning("Invalid Index", f"Index {idx+1} is out of range.")
            except ValueError:
                messagebox.showerror("Invalid Input", "Please enter a valid number.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to remove pair: {e}")
        def _quick_import_pairs():
            """Auto-generate all pairwise combinations from current groups."""
            try:
                groups = self.get_viz_groups()
                if not groups or len(groups) < 2:
                    messagebox.showwarning("Quick Import", "Need at least 2 groups configured. Please configure groups first.")
                    return
                from itertools import combinations
                # Preserve configured group order (do not sort alphabetically)
                new_pairs = list(combinations(groups, 2))
                if new_pairs:
                    self.venn_pairs_list.extend(new_pairs)
                    self.venn_pairs_var.set(', '.join([f"{_canon_pair(x,y)[0]}|{_canon_pair(x,y)[1]}" for x,y in self.venn_pairs_list]))
                    _refresh_pairs_index_label()
                    self.log_viz_message(f"Quick import: added {len(new_pairs)} pairwise combinations from {len(groups)} groups")
                    messagebox.showinfo("Quick Import", f"Added {len(new_pairs)} pairwise pairs from detected groups.")
            except Exception as e:
                messagebox.showerror("Quick Import Error", f"Failed to import pairs: {e}")
        
        # Move the action buttons to the next line under the comparison entry
        buttons_line = ttk.Frame(advanced_frame)
        buttons_line.pack(fill='x', pady=(4,2))
        ttk.Button(buttons_line, text="Add", command=_add_pair).pack(side='left', padx=2)
        ttk.Button(buttons_line, text="Remove", command=_remove_pair).pack(side='left', padx=2)
        ttk.Button(buttons_line, text="Quick Import", command=_quick_import_pairs).pack(side='left', padx=2)

        self.venn_pairs_var = tk.StringVar(value="")
        ttk.Label(advanced_frame, textvariable=self.venn_pairs_var, wraplength=500).pack(anchor='w', pady=(2,6))

        # Venn specs builder
        venn_row = ttk.Frame(advanced_frame)
        venn_row.pack(fill='x', pady=2)
        ttk.Label(venn_row, text="Venn Name:").pack(side='left')
        self.venn_name_entry = tk.StringVar()
        ttk.Entry(venn_row, textvariable=self.venn_name_entry, width=24).pack(side='left', padx=4)
        ttk.Label(venn_row, text="Comparisons to include (A,B)").pack(side='left', padx=6)
        self.venn_indices_entry = tk.StringVar()
        ttk.Entry(venn_row, textvariable=self.venn_indices_entry, width=24).pack(side='left')

        if not hasattr(self, 'venn_specs'):
            self.venn_specs = []  # list[dict{name, indices:list[int]}]
        def _add_venn_spec():
            name = self.venn_name_entry.get().strip() or f"Venn{len(self.venn_specs)+1}"
            idxs = []
            raw = self.venn_indices_entry.get().strip()
            if raw:
                for x in raw.split(','):
                    try:
                        k = int(x.strip())
                        if k >= 1: idxs.append(k)
                    except Exception:
                        pass
            self.venn_specs.append({'name': name, 'indices': idxs})
            self.venn_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_specs]))
            self.venn_name_entry.set('')
            self.venn_indices_entry.set('')
        def _remove_venn_spec():
            try:
                venn_name = simpledialog.askstring("Remove Venn", "Enter Venn name to remove:")
                if venn_name:
                    found = False
                    for i, spec in enumerate(self.venn_specs):
                        if spec['name'] == venn_name:
                            removed = self.venn_specs.pop(i)
                            self.venn_specs_var.set('; '.join([f"{v['name']}:[{','.join(map(str,v['indices']))}]" for v in self.venn_specs]))
                            self.log_viz_message(f"Removed Venn: {removed['name']}")
                            found = True
                            break
                    if not found:
                        messagebox.showwarning("Not Found", f"Venn '{venn_name}' not found.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to remove Venn: {e}")
        venn_buttons_row = ttk.Frame(advanced_frame)
        venn_buttons_row.pack(fill='x', pady=(2,2))
        ttk.Button(venn_buttons_row, text="Add Venn", command=_add_venn_spec).pack(side='left', padx=2)
        ttk.Button(venn_buttons_row, text="Remove Venn", command=_remove_venn_spec).pack(side='left', padx=2)

        self.venn_specs_var = tk.StringVar(value="")
        ttk.Label(advanced_frame, text="Defined Venns:").pack(anchor='w')
        ttk.Label(advanced_frame, textvariable=self.venn_specs_var, wraplength=500).pack(anchor='w')

        # Helper: list pairs with indices
        def _refresh_pairs_index_label():
            if self.venn_pairs_list:
                lines = [f"{i+1}:{_canon_pair(a,b)[0]}|{_canon_pair(a,b)[1]}" for i,(a,b) in enumerate(self.venn_pairs_list)]
                self.venn_pairs_index_label.configure(text='Pairs Index: ' + '; '.join(lines))
            else:
                self.venn_pairs_index_label.configure(text='Pairs Index: (none)')
        self.venn_pairs_index_label = ttk.Label(advanced_frame, text='Pairs Index: (none)')
        self.venn_pairs_index_label.pack(anchor='w', pady=(4,6))
        # Update indices label whenever pairs list changes via Add
        pairs_row.bind_all('<<VennPairsChanged>>', lambda e: _refresh_pairs_index_label())
        
        # Add visual summary at bottom
        summary_frame = ttk.LabelFrame(panel, text="Current Configuration Summary", padding=8)
        summary_frame.pack(fill='x', pady=(10,5))
        
        self.venn_summary_label = ttk.Label(summary_frame, 
                                            text="No Venn configured yet. Use 'Auto-Generate' button above.",
                                            foreground='gray')
        self.venn_summary_label.pack(anchor='w', pady=4)
        
        def _update_venn_summary():
            """Update the visual summary of current Venn configuration."""
            try:
                if not self.venn_specs or not self.venn_pairs_list:
                    self.venn_summary_label.configure(
                        text="No Venn configured yet. Use 'Auto-Generate' button above.",
                        foreground='gray')
                    return
                
                summary_lines = []
                for spec in self.venn_specs:
                    name = spec.get('name', 'Unnamed')
                    indices = spec.get('indices', [])
                    comparisons = []
                    for idx in indices:
                        if 1 <= idx <= len(self.venn_pairs_list):
                            g1, g2 = self.venn_pairs_list[idx-1]
                            ca, cb = _canon_pair(g1, g2)
                            comparisons.append(f"{ca} vs {cb}")
                    
                    if comparisons:
                        summary_lines.append(f"📊 {name}:")
                        summary_lines.append(f"   Comparisons: {', '.join(comparisons)}")
                        summary_lines.append(f"   (Total: {len(comparisons)} sets)")
                
                if summary_lines:
                    self.venn_summary_label.configure(
                        text='\n'.join(summary_lines),
                        foreground='darkgreen')
                else:
                    self.venn_summary_label.configure(
                        text="Venn defined but no valid comparisons selected.",
                        foreground='orange')
            except Exception:
                pass
        
        # Store the update function so auto-generate can call it
        self._update_venn_summary = _update_venn_summary
        
        # Bind to update summary when specs change
        self.venn_specs_var.trace_add('write', lambda *args: _update_venn_summary())
    
    
        # --------------- Tooltip Helper ---------------
    
    def _create_tooltip(self, widget, text):
        tooltip = {'win': None}
        def enter(_):
            if tooltip['win'] is not None:
                return
            try:
                x, y, cx, cy = widget.bbox("insert") if hasattr(widget, 'bbox') else (0, 0, 0, 0)
            except Exception: 
                x = y = 0
            x += widget.winfo_rootx() + 25
            y += widget.winfo_rooty() + 20
            win = tk.Toplevel(widget)
            win.wm_overrideredirect(True)
            win.wm_geometry(f"+{x}+{y}")
            lbl = tk.Label(win, text=text, justify='left', background='#ffffe0', relief='solid', borderwidth=1, font=('Arial', 8))
            lbl.pack(ipadx=4, ipady=2)
            tooltip['win'] = win
        def leave(_):
            if tooltip['win'] is not None:
                try:
                    tooltip['win'].destroy()
                except Exception:
                    pass
                tooltip['win'] = None
        widget.bind('<Enter>', enter)
        widget.bind('<Leave>', leave)

    # --------------- Comparison Selector Dialog (Reusable) ---------------
    def _show_comparison_selector_dialog(self, all_pairs, current_selection=None, title="Select Comparisons"):
        """Show a dialog to select which comparisons to include.
        
        Parameters
        ----------
        all_pairs : list of tuple
            All possible comparison pairs (g1, g2)
        current_selection : list of tuple, optional
            Currently selected pairs
        title : str
            Dialog title
            
        Returns
        -------
        list of tuple or None
            Selected pairs, or None if cancelled
        """
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("500x600")
        dialog.transient(self.root)
        dialog.grab_set()
        
        result = {'selection': None}
        
        # Header
        header_frame = ttk.Frame(dialog, padding=10)
        header_frame.pack(fill='x')
        ttk.Label(header_frame, text=title, font=('Arial', 12, 'bold')).pack(anchor='w')
        ttk.Label(header_frame, text="Select which comparisons to include in the analysis:",
                 font=('Arial', 9)).pack(anchor='w', pady=(5, 0))
        
        # Checkbox list
        list_frame = ttk.Frame(dialog, padding=10)
        list_frame.pack(fill='both', expand=True)
        
        # Add scrollbar
        canvas = tk.Canvas(list_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(list_frame, orient='vertical', command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mouse wheel / two-finger scrolling
        # Ensure canvas can be scrolled with mouse wheel
        if hasattr(self, '_enable_mousewheel_scroll'):
            self._enable_mousewheel_scroll(canvas)
        
        # Create checkboxes for each pair
        check_vars = {}
        current_set = set(current_selection) if current_selection else set()
        
        for g1, g2 in all_pairs:
            pair_key = (g1, g2)
            var = tk.BooleanVar(value=(pair_key in current_set or (g2, g1) in current_set))
            check_vars[pair_key] = var
            
            cb = ttk.Checkbutton(scrollable_frame, text=f"{g1} vs {g2}", variable=var)
            cb.pack(anchor='w', padx=5, pady=2)
        
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        # Selection controls
        control_frame = ttk.Frame(dialog, padding=10)
        control_frame.pack(fill='x')
        
        def select_all():
            for var in check_vars.values():
                var.set(True)
        
        def deselect_all():
            for var in check_vars.values():
                var.set(False)
        
        ttk.Button(control_frame, text="Select All", command=select_all).pack(side='left', padx=5)
        ttk.Button(control_frame, text="Deselect All", command=deselect_all).pack(side='left', padx=5)
        
        # Count label
        count_label = ttk.Label(control_frame, text="", font=('Arial', 9))
        count_label.pack(side='left', padx=20)
        
        def update_count():
            count = sum(1 for var in check_vars.values() if var.get())
            count_label.configure(text=f"Selected: {count}/{len(all_pairs)}")
        
        # Update count when checkboxes change
        for var in check_vars.values():
            var.trace_add('write', lambda *args: update_count())
        
        update_count()  # Initial count
        
        # Buttons
        button_frame = ttk.Frame(dialog, padding=10)
        button_frame.pack(fill='x')
        
        def on_ok():
            selected = [pair for pair, var in check_vars.items() if var.get()]
            result['selection'] = selected
            dialog.destroy()
        
        def on_cancel():
            result['selection'] = None
            dialog.destroy()
        
        ttk.Button(button_frame, text="OK", command=on_ok, width=12).pack(side='right', padx=5)
        ttk.Button(button_frame, text="Cancel", command=on_cancel, width=12).pack(side='right', padx=5)
        
        # Wait for dialog to close
        dialog.wait_window()
        
        return result['selection']

    def _configure_volcano_comparisons(self):
        """Show dialog to select which comparisons to plot in volcano plots."""
        # Get all possible comparisons from the data
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        result = self._show_comparison_selector_dialog(
            all_pairs, 
            self.volcano_selected_comparisons,
            title="Select Volcano Comparisons"
        )
        
        if result is not None:
            self.volcano_selected_comparisons = result
            # Update status label
            if not result:  # Empty list
                self.volcano_comp_status_label.config(text="⚠️ No comparisons selected - no plots will be generated")
            elif len(result) == len(all_pairs):
                self.volcano_comp_status_label.config(text="All comparisons will be plotted")
            else:
                self.volcano_comp_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons")

    def _configure_boxplot_comparisons(self):
        """Show dialog to select which comparisons to plot in boxplots."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        result = self._show_comparison_selector_dialog(
            all_pairs, 
            self.boxplot_selected_comparisons,
            title="Select Boxplot Comparisons"
        )
        
        if result is not None:
            self.boxplot_selected_comparisons = result
            if not result:
                self.boxplot_comp_status_label.config(text="⚠️ No comparisons selected - no plots will be generated")
            elif len(result) == len(all_pairs):
                self.boxplot_comp_status_label.config(text="All comparisons will be plotted")
            else:
                self.boxplot_comp_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons")
    
    def _configure_boxplot_annotations(self):
        """Show dialog to select which comparisons to annotate in boxplots."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        result = self._show_comparison_selector_dialog(
            all_pairs, 
            self.boxplot_annotate_comparisons,
            title="Select Comparisons to Annotate"
        )
        
        if result is not None:
            self.boxplot_annotate_comparisons = result
            if not result:
                self.boxplot_annotate_status_label.config(text="⚠️ No comparisons will be annotated")
            elif len(result) == len(all_pairs):
                self.boxplot_annotate_status_label.config(text="All comparisons will be annotated (if 'Significance stars' enabled)")
            else:
                self.boxplot_annotate_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons to annotate")

    def _configure_boxplot_groups(self):
        """Show dialog to select which groups to display in boxplots."""
        groups = self.get_viz_groups()
        if not groups:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        if len(groups) < 2:
            messagebox.showinfo("Info", "You have only one group. At least 2 groups are needed for meaningful boxplots.")
            return
        
        # Create dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("Select Groups to Display in Boxplots")
        dialog.geometry("450x500")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="Select which groups to include in boxplots:",
                 font=('Arial', 10, 'bold')).pack(pady=10, padx=10)
        ttk.Label(dialog, text="Only selected groups will appear in the boxplot visualizations.",
                 font=('Arial', 9), foreground='#666').pack(pady=5, padx=10)
        
        # Frame for checkboxes
        checkbox_frame = ttk.Frame(dialog)
        checkbox_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Add scrollbar for many groups
        canvas = tk.Canvas(checkbox_frame)
        scrollbar = ttk.Scrollbar(checkbox_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mouse wheel scrolling
        def _on_mousewheel(event):
            try:
                if canvas.winfo_exists():
                    canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except tk.TclError:
                pass
        canvas.bind("<MouseWheel>", _on_mousewheel)
        
        # Create checkbox variables
        group_vars = {}
        for group in groups:
            # Default: all selected (either None means all, or all are in the list)
            if self.boxplot_selected_groups is None:
                default_value = True  # None means all groups selected
            else:
                default_value = group in self.boxplot_selected_groups
            var = tk.BooleanVar(value=default_value)
            group_vars[group] = var
            ttk.Checkbutton(scrollable_frame, text=group, variable=var).pack(anchor='w', pady=2)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill='x', padx=10, pady=10)
        
        def select_all():
            for var in group_vars.values():
                var.set(True)
        
        def deselect_all():
            for var in group_vars.values():
                var.set(False)
        
        def apply_selection():
            selected = [group for group, var in group_vars.items() if var.get()]
            
            if len(selected) < 2:
                messagebox.showwarning("Invalid Selection", 
                                     "Please select at least 2 groups for boxplots.\n"
                                     "Boxplots require multiple groups for comparison.")
                return
            
            # If all groups are selected, set to None (meaning "all groups")
            if len(selected) == len(groups):
                self.boxplot_selected_groups = None
            else:
                self.boxplot_selected_groups = selected
            
            # Update status label
            if len(selected) == len(groups):
                self.boxplot_groups_status_label.config(text="All groups will be displayed")
            else:
                group_list = ", ".join(selected)
                self.boxplot_groups_status_label.config(text=f"Selected: {', '.join(selected[:3])}{'...' if len(selected) > 3 else ''} ({len(selected)}/{len(groups)} groups)")
            
            canvas.unbind("<MouseWheel>")
            dialog.destroy()
        
        def cancel():
            canvas.unbind("<MouseWheel>")
            dialog.destroy()
        
        ttk.Button(btn_frame, text="Select All", command=select_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Deselect All", command=deselect_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Apply", command=apply_selection).pack(side='right', padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel).pack(side='right', padx=5)

    def _configure_bargraph_comparisons(self):
        """Show dialog to select which comparisons to plot in bar graphs."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return

        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))

        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return

        result = self._show_comparison_selector_dialog(
            all_pairs,
            self.bargraph_selected_comparisons,
            title="Select Bar Graph Comparisons"
        )

        if result is not None:
            self.bargraph_selected_comparisons = result
            if not result:
                self.bargraph_comp_status_label.config(text="⚠️ No comparisons selected - no plots will be generated")
            elif len(result) == len(all_pairs):
                self.bargraph_comp_status_label.config(text="All comparisons will be plotted")
            else:
                self.bargraph_comp_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons")

    def _configure_bargraph_annotations(self):
        """Show dialog to select which comparisons to annotate in bar graphs."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return

        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))

        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return

        result = self._show_comparison_selector_dialog(
            all_pairs,
            self.bargraph_annotate_comparisons,
            title="Select Bar Graph Comparisons to Annotate"
        )

        if result is not None:
            self.bargraph_annotate_comparisons = result
            if not result:
                self.bargraph_annotate_status_label.config(text="⚠️ No comparisons will be annotated")
            elif len(result) == len(all_pairs):
                self.bargraph_annotate_status_label.config(text="All comparisons will be annotated (if 'Significance stars' enabled)")
            else:
                self.bargraph_annotate_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons to annotate")

    def _configure_bargraph_groups(self):
        """Show dialog to select which groups to display in bar graphs."""
        groups = self.get_viz_groups()
        if not groups:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return

        if len(groups) < 2:
            messagebox.showinfo("Info", "You have only one group. At least 2 groups are needed for meaningful bar graphs.")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Select Groups to Display in Bar Graphs")
        dialog.geometry("450x500")
        dialog.transient(self.root)
        dialog.grab_set()

        ttk.Label(dialog, text="Select which groups to include in bar graphs:",
                 font=('Arial', 10, 'bold')).pack(pady=10, padx=10)
        ttk.Label(dialog, text="Only selected groups will appear in the bar graph visualizations.",
                 font=('Arial', 9), foreground='#666').pack(pady=5, padx=10)

        checkbox_frame = ttk.Frame(dialog)
        checkbox_frame.pack(fill='both', expand=True, padx=10, pady=10)

        canvas = tk.Canvas(checkbox_frame)
        scrollbar = ttk.Scrollbar(checkbox_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        def _on_mousewheel(event):
            try:
                if canvas.winfo_exists():
                    canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except tk.TclError:
                pass
        canvas.bind("<MouseWheel>", _on_mousewheel)

        group_vars = {}
        for group in groups:
            if self.bargraph_selected_groups is None:
                default_value = True
            else:
                default_value = group in self.bargraph_selected_groups
            var = tk.BooleanVar(value=default_value)
            group_vars[group] = var
            ttk.Checkbutton(scrollable_frame, text=group, variable=var).pack(anchor='w', pady=2)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill='x', padx=10, pady=10)

        def select_all():
            for var in group_vars.values():
                var.set(True)

        def deselect_all():
            for var in group_vars.values():
                var.set(False)

        def apply_selection():
            selected = [group for group, var in group_vars.items() if var.get()]

            if len(selected) < 2:
                messagebox.showwarning("Invalid Selection",
                                     "Please select at least 2 groups for bar graphs.\n"
                                     "Bar graphs require multiple groups for comparison.")
                return

            if len(selected) == len(groups):
                self.bargraph_selected_groups = None
            else:
                self.bargraph_selected_groups = selected

            if len(selected) == len(groups):
                self.bargraph_groups_status_label.config(text="All groups will be displayed")
            else:
                self.bargraph_groups_status_label.config(text=f"Selected: {', '.join(selected[:3])}{'...' if len(selected) > 3 else ''} ({len(selected)}/{len(groups)} groups)")

            canvas.unbind("<MouseWheel>")
            dialog.destroy()

        def cancel():
            canvas.unbind("<MouseWheel>")
            dialog.destroy()

        ttk.Button(btn_frame, text="Select All", command=select_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Deselect All", command=deselect_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Apply", command=apply_selection).pack(side='right', padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel).pack(side='right', padx=5)

    def _configure_roc_groups(self):
        """Show dialog to select which groups to display in ROC curves."""
        groups = self.get_viz_groups()
        if not groups:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        if len(groups) < 2:
            messagebox.showinfo("Info", "You have only one group. At least 2 groups are needed for ROC curves.")
            return
        
        # Create dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("Select Groups to Display in ROC Curves")
        dialog.geometry("450x500")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="Select which groups to include in ROC curves:",
                 font=('Arial', 10, 'bold')).pack(pady=10, padx=10)
        ttk.Label(dialog, text="Only selected groups will appear in the ROC visualizations.",
                 font=('Arial', 9), foreground='#666').pack(pady=5, padx=10)
        
        # Frame for checkboxes
        checkbox_frame = ttk.Frame(dialog)
        checkbox_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Add scrollbar for many groups
        canvas = tk.Canvas(checkbox_frame)
        scrollbar = ttk.Scrollbar(checkbox_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mouse wheel scrolling
        def _on_mousewheel(event):
            try:
                if canvas.winfo_exists():
                    canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except tk.TclError:
                pass
        canvas.bind("<MouseWheel>", _on_mousewheel)
        
        # Create checkbox variables
        group_vars = {}
        for group in groups:
            # Default: all selected (either None means all, or all are in the list)
            if not hasattr(self, 'roc_selected_groups') or self.roc_selected_groups is None:
                default_value = True  # None means all groups selected
            else:
                default_value = group in self.roc_selected_groups
            var = tk.BooleanVar(value=default_value)
            group_vars[group] = var
            ttk.Checkbutton(scrollable_frame, text=group, variable=var).pack(anchor='w', pady=2)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill='x', padx=10, pady=10)
        
        def select_all():
            for var in group_vars.values():
                var.set(True)
        
        def deselect_all():
            for var in group_vars.values():
                var.set(False)
        
        def apply_selection():
            selected = [group for group, var in group_vars.items() if var.get()]
            
            if len(selected) < 2:
                messagebox.showwarning("Invalid Selection", 
                                     "Please select at least 2 groups for ROC curves.\n"
                                     "ROC curves require multiple groups for comparison.")
                return
            
            # If all groups are selected, set to None (meaning "all groups")
            if len(selected) == len(groups):
                self.roc_selected_groups = None
            else:
                self.roc_selected_groups = selected
            
            # Update status label
            if len(selected) == len(groups):
                self.roc_groups_status_label.config(text="All groups will be displayed")
            else:
                self.roc_groups_status_label.config(text=f"Selected: {', '.join(selected[:3])}{'...' if len(selected) > 3 else ''} ({len(selected)}/{len(groups)} groups)")
            
            canvas.unbind("<MouseWheel>")
            dialog.destroy()
        
        def cancel():
            canvas.unbind("<MouseWheel>")
            dialog.destroy()
        
        ttk.Button(btn_frame, text="Select All", command=select_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Deselect All", command=deselect_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Apply", command=apply_selection).pack(side='right', padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel).pack(side='right', padx=5)

    def _configure_boxplot_filter_comparison(self):
        """Show dialog to select a SINGLE comparison for boxplot significance filter (All/Specific mode)."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        # Create a dialog to select ONE comparison
        dialog = tk.Toplevel(self.root)
        dialog.title("Select Reference Comparison for Filter")
        dialog.geometry("400x450")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="Select ONE comparison to use for 'All' or 'Specific' filter mode:",
                 font=('Arial', 10, 'bold')).pack(pady=10, padx=10)
        ttk.Label(dialog, text="This determines which p-value/FC columns are used for filtering.",
                 font=('Arial', 9), foreground='#666').pack(pady=5, padx=10)
        
        # Frame for radio buttons
        radio_frame = ttk.Frame(dialog)
        radio_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        selected_var = tk.StringVar()
        if self.boxplot_filter_comparison:
            g1, g2 = self.boxplot_filter_comparison
            selected_var.set(f"{g1}_vs_{g2}")
        
        # Add scrollbar for many comparisons
        canvas = tk.Canvas(radio_frame)
        scrollbar = ttk.Scrollbar(radio_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mouse wheel / two-finger scrolling
        if hasattr(self, '_enable_mousewheel_scroll'):
            self._enable_mousewheel_scroll(canvas)
        
        for g1, g2 in all_pairs:
            ttk.Radiobutton(scrollable_frame, text=f"{g1} vs {g2}", 
                          variable=selected_var, value=f"{g1}_vs_{g2}").pack(anchor='w', pady=2)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Buttons
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill='x', padx=10, pady=10)
        
        result_holder = [None]
        
        def on_ok():
            selection = selected_var.get()
            if selection:
                # Parse back to tuple
                parts = selection.split('_vs_')
                if len(parts) == 2:
                    result_holder[0] = (parts[0], parts[1])
            dialog.destroy()
        
        def on_cancel():
            dialog.destroy()
        
        ttk.Button(btn_frame, text="OK", command=on_ok).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side='left', padx=5)
        
        dialog.wait_window()
        
        if result_holder[0]:
            self.boxplot_filter_comparison = result_holder[0]
            g1, g2 = result_holder[0]
            self.boxplot_filter_comp_status_label.config(text=f"Selected: {g1} vs {g2}")

    def _configure_heatmap_comparisons(self):
        """Configure which comparisons to include in heatmap generation with combined options."""
        if not hasattr(self, 'viz_group_definitions') or not self.viz_group_definitions:
            messagebox.showwarning("No Groups", "Please define visualization groups first.")
            return
        
        # Generate all possible pairwise comparisons
        groups = list(self.viz_group_definitions.values())
        from itertools import combinations
        all_pairs = list(combinations(groups, 2))
        
        if not all_pairs:
            messagebox.showwarning("Not Enough Groups", "Need at least 2 groups for comparisons.")
            return
        
        # Create dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("Configure Heatmap Comparisons")
        dialog.geometry("600x550")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Main frame
        main_frame = ttk.Frame(dialog, padding=10)
        main_frame.pack(fill='both', expand=True)
        
        ttk.Label(main_frame, text="Select Comparisons for Heatmap Generation", 
                font=('Arial', 12, 'bold')).pack(pady=(0,10))

        # Keep primary actions at the top for faster access.
        top_btn_frame = ttk.Frame(main_frame)
        top_btn_frame.pack(fill='x', pady=(0, 8))
        
        # Comparison selection frame with scrollbar
        comp_frame = ttk.LabelFrame(main_frame, text="Individual Comparisons", padding=5)
        comp_frame.pack(fill='both', expand=True, pady=5)
        
        canvas = tk.Canvas(comp_frame, height=200)
        scrollbar = ttk.Scrollbar(comp_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Enable mouse wheel / two-finger scrolling
        if hasattr(self, '_enable_mousewheel_scroll'):
            self._enable_mousewheel_scroll(canvas)
        
        # Create checkboxes for each comparison
        comp_vars = {}
        for g1, g2 in all_pairs:
            var = tk.BooleanVar(value=True)  # AUTO-CHECK ALL by default
            comp_vars[(g1, g2)] = var
            ttk.Checkbutton(scrollable_frame, text=f"{g1} vs {g2}", 
                        variable=var).pack(anchor='w', padx=5, pady=2)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Combined heatmap options frame
        combined_frame = ttk.LabelFrame(main_frame, text="Combined Heatmap Options", padding=10)
        combined_frame.pack(fill='x', pady=5)
        
        # Combined heatmap checkbox
        combined_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(combined_frame, text="Generate Combined Heatmap", 
                    variable=combined_var).pack(anchor='w', pady=(0,10))
        
        # Combined mode selection
        mode_frame = ttk.Frame(combined_frame)
        mode_frame.pack(fill='x', padx=20)
        
        # Metabolite scope selection (new)
        ttk.Label(mode_frame, text="Metabolite Scope:", 
                font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        
        metabolite_scope_var = tk.StringVar(value='selected')
        
        scope_info = {
            'selected': "Selected: Use only metabolites from selected comparisons",
            'all': "All: Use all significant metabolites across all comparisons"
        }
        
        for scope, description in scope_info.items():
            scope_frame = ttk.Frame(mode_frame)
            scope_frame.pack(fill='x', pady=2)
            ttk.Radiobutton(scope_frame, text=scope.capitalize(), value=scope, 
                        variable=metabolite_scope_var).pack(side='left')
            ttk.Label(scope_frame, text=description, font=('Arial', 8, 'italic'), 
                    foreground='#666').pack(side='left', padx=(10,0))
        
        # Combination Mode selection
        ttk.Label(mode_frame, text="Combination Mode:", 
                font=('Arial', 9, 'bold')).pack(anchor='w', pady=(10,5))
        
        combined_mode_var = tk.StringVar(value='union')
        
        mode_info = {
            'union': "Union: All significant metabolites (respects scope setting above)",
            'intersection': "Intersection: Metabolites significant in ALL selected comparisons (respects scope)",
            'concatenate': "Concatenate: Lists appended in comparison order (first occurrence kept)"
        }
        
        for mode, description in mode_info.items():
            rb_frame = ttk.Frame(mode_frame)
            rb_frame.pack(fill='x', pady=2)
            ttk.Radiobutton(rb_frame, text=mode.capitalize(), value=mode, 
                        variable=combined_mode_var).pack(side='left')
            ttk.Label(rb_frame, text=description, font=('Arial', 8, 'italic'), 
                    foreground='#666').pack(side='left', padx=(10,0))
        
        # Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill='x', pady=(10,0))
        
        def apply_selection():
            selected = [(g1, g2) for (g1, g2), var in comp_vars.items() if var.get()]
            
            if not selected:
                messagebox.showwarning("No Selection", "Please select at least one comparison.")
                return
            
            # Store selected comparisons
            self.heatmap_selected_comparisons = selected
            
            # Update status label
            if len(selected) == len(all_pairs):
                self.heatmap_comp_status_label.config(text="All comparisons will be plotted")
            else:
                self.heatmap_comp_status_label.config(
                    text=f"{len(selected)} comparison(s) selected: " + 
                        ", ".join([f"{g1} vs {g2}" for g1, g2 in selected[:3]]) +
                        (f" and {len(selected)-3} more" if len(selected) > 3 else "")
                )
            
            # Store combined heatmap settings
            if hasattr(self, 'viz_params') and 'heatmap' in self.viz_params:
                self.viz_params['heatmap'].combined = combined_var.get()
                self.viz_params['heatmap'].combined_mode = combined_mode_var.get()
                self.viz_params['heatmap'].metabolite_scope = metabolite_scope_var.get()
            
            # Also store in attributes for settings persistence
            if not hasattr(self, 'heatmap_combined'):
                self.heatmap_combined = tk.BooleanVar()
            if not hasattr(self, 'heatmap_combined_mode'):
                self.heatmap_combined_mode = tk.StringVar()
            if not hasattr(self, 'heatmap_metabolite_scope'):
                self.heatmap_metabolite_scope = tk.StringVar()
            self.heatmap_combined.set(combined_var.get())
            self.heatmap_combined_mode.set(combined_mode_var.get())
            self.heatmap_metabolite_scope.set(metabolite_scope_var.get())
            
            logger.info(f"Heatmap comparisons configured: {len(selected)} selected, "
                    f"Combined: {combined_var.get()}, Mode: {combined_mode_var.get()}")
            
            dialog.destroy()
        
        def select_all():
            for var in comp_vars.values():
                var.set(True)
        
        def deselect_all():
            for var in comp_vars.values():
                var.set(False)
        
        ttk.Button(btn_frame, text="Select All", command=select_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Deselect All", command=deselect_all).pack(side='left', padx=5)
        ttk.Button(top_btn_frame, text="Cancel", command=dialog.destroy).pack(side='right', padx=5)
        ttk.Button(top_btn_frame, text="Apply", command=apply_selection).pack(side='right', padx=5)
        
        # Load existing selections if any
        if hasattr(self, 'heatmap_selected_comparisons') and self.heatmap_selected_comparisons:
            for g1, g2 in all_pairs:
                if (g1, g2) not in self.heatmap_selected_comparisons and (g2, g1) not in self.heatmap_selected_comparisons:
                    comp_vars[(g1, g2)].set(False)
        
        # Load existing combined settings
        if hasattr(self, 'heatmap_combined'):
            try:
                combined_var.set(self.heatmap_combined.get())
            except:
                pass
        if hasattr(self, 'heatmap_combined_mode'):
            try:
                combined_mode_var.set(self.heatmap_combined_mode.get())
            except:
                pass
        if hasattr(self, 'heatmap_metabolite_scope'):
            try:
                metabolite_scope_var.set(self.heatmap_metabolite_scope.get())
            except:
                pass

    def _configure_heatmap_metabolite_lists(self):
        """Show dialog to manage per-comparison metabolite lists for heatmaps."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        # Create dialog for managing per-comparison lists
        dialog = tk.Toplevel(self.root)
        dialog.title("Manage Heatmap Metabolite Lists")
        dialog.geometry("600x500")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Instructions
        info_frame = ttk.Frame(dialog, padding=10)
        info_frame.pack(fill='x')
        ttk.Label(info_frame, text="Upload different metabolite lists for different comparisons.", 
                 font=('Arial', 10, 'bold')).pack(anchor='w')
        ttk.Label(info_frame, text="Lists uploaded here take priority over the general metabolite list.", 
                 font=('Arial', 9, 'italic'), foreground='#666').pack(anchor='w')
        
        # Scrollable frame for comparison list entries
        canvas = tk.Canvas(dialog)
        scrollbar = ttk.Scrollbar(dialog, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Store entry widgets
        list_entries = {}
        
        for g1, g2 in all_pairs:
            pair_frame = ttk.LabelFrame(scrollable_frame, text=f"{g1} vs {g2}", padding=5)
            pair_frame.pack(fill='x', padx=5, pady=3)
            
            current_file = self.heatmap_metabolite_lists.get((g1, g2), "")
            
            entry_var = tk.StringVar(value=current_file)
            list_entries[(g1, g2)] = entry_var
            
            entry = ttk.Entry(pair_frame, textvariable=entry_var, width=50)
            entry.pack(side='left', expand=True, fill='x', padx=(0, 5))
            
            def browse_func(pair=( g1, g2), var=entry_var):
                filename = filedialog.askopenfilename(
                    title=f"Select metabolite list for {pair[0]} vs {pair[1]}",
                    filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
                )
                if filename:
                    var.set(filename)
            
            ttk.Button(pair_frame, text="Browse", command=browse_func).pack(side='left')
            
            def clear_func(pair=(g1, g2), var=entry_var):
                var.set("")
            
            ttk.Button(pair_frame, text="Clear", command=clear_func).pack(side='left', padx=(2, 0))
        
        canvas.pack(side='left', fill='both', expand=True, padx=(10, 0), pady=10)
        scrollbar.pack(side='right', fill='y', pady=10)
        
        # Buttons
        button_frame = ttk.Frame(dialog, padding=10)
        button_frame.pack(fill='x')
        
        def save_and_close():
            # Update the metabolite lists dictionary
            for pair, var in list_entries.items():
                filepath = var.get().strip()
                if filepath:
                    self.heatmap_metabolite_lists[pair] = filepath
                elif pair in self.heatmap_metabolite_lists:
                    del self.heatmap_metabolite_lists[pair]
            
            dialog.destroy()
        
        ttk.Button(button_frame, text="Save", command=save_and_close).pack(side='right', padx=5)
        ttk.Button(button_frame, text="Cancel", command=dialog.destroy).pack(side='right')
        
        dialog.wait_window()
    
    def _configure_roc_comparisons(self):
        """Show dialog to select which comparisons to plot in ROC curves."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        result = self._show_comparison_selector_dialog(
            all_pairs, 
            self.roc_selected_comparisons,
            title="Select ROC Comparisons"
        )
        
        if result is not None:
            self.roc_selected_comparisons = result
            if not result:
                self.roc_comp_status_label.config(text="⚠️ No comparisons selected - no plots will be generated")
            elif len(result) == len(all_pairs):
                self.roc_comp_status_label.config(text="All comparisons will be plotted")
            else:
                self.roc_comp_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} comparisons")
    
    def _configure_roc_metabolite_lists(self):
        """Show dialog to manage per-comparison metabolite lists for ROC curves."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        # Create dialog for managing per-comparison lists
        dialog = tk.Toplevel(self.root)
        dialog.title("Manage ROC Metabolite Lists")
        dialog.geometry("600x500")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Instructions
        info_frame = ttk.Frame(dialog, padding=10)
        info_frame.pack(fill='x')
        ttk.Label(info_frame, text="Upload different metabolite lists for different comparisons.", 
                 font=('Arial', 10, 'bold')).pack(anchor='w')
        ttk.Label(info_frame, text="Lists uploaded here take priority over the general metabolite list.", 
                 font=('Arial', 9, 'italic'), foreground='#666').pack(anchor='w')
        
        # Scrollable frame for comparison list entries
        canvas = tk.Canvas(dialog)
        scrollbar = ttk.Scrollbar(dialog, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Store entry widgets
        list_entries = {}
        
        for g1, g2 in all_pairs:
            pair_frame = ttk.LabelFrame(scrollable_frame, text=f"{g1} vs {g2}", padding=5)
            pair_frame.pack(fill='x', padx=5, pady=3)
            
            current_file = self.roc_metabolite_lists.get((g1, g2), "")
            
            entry_var = tk.StringVar(value=current_file)
            list_entries[(g1, g2)] = entry_var
            
            entry = ttk.Entry(pair_frame, textvariable=entry_var, width=50)
            entry.pack(side='left', expand=True, fill='x', padx=(0, 5))
            
            def browse_func(pair=(g1, g2), var=entry_var):
                filename = filedialog.askopenfilename(
                    title=f"Select metabolite list for {pair[0]} vs {pair[1]}",
                    filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
                )
                if filename:
                    var.set(filename)
            
            ttk.Button(pair_frame, text="Browse", command=browse_func).pack(side='left')
            
            def clear_func(pair=(g1, g2), var=entry_var):
                var.set("")
            
            ttk.Button(pair_frame, text="Clear", command=clear_func).pack(side='left', padx=(2, 0))
        
        canvas.pack(side='left', fill='both', expand=True, padx=(10, 0), pady=10)
        scrollbar.pack(side='right', fill='y', pady=10)
        
        # Buttons
        button_frame = ttk.Frame(dialog, padding=10)
        button_frame.pack(fill='x')
        
        def save_and_close():
            # Update the metabolite lists dictionary
            for pair, var in list_entries.items():
                filepath = var.get().strip()
                if filepath:
                    self.roc_metabolite_lists[pair] = filepath
                elif pair in self.roc_metabolite_lists:
                    del self.roc_metabolite_lists[pair]
            
            dialog.destroy()
        
        ttk.Button(button_frame, text="Save", command=save_and_close).pack(side='right', padx=5)
        ttk.Button(button_frame, text="Cancel", command=dialog.destroy).pack(side='right')
        
        dialog.wait_window()
    
    def _configure_pca_comparisons(self):
        """Show dialog to select which pairwise comparisons to generate in PCA."""
        groups = self.get_viz_groups()
        if not groups or len(groups) < 2:
            messagebox.showwarning("No Data", "Please load data and select groups first.")
            return
        
        all_pairs = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                all_pairs.append((g1, g2))
        
        if not all_pairs:
            messagebox.showwarning("Insufficient Groups", "Need at least 2 groups for comparisons.")
            return
        
        result = self._show_comparison_selector_dialog(
            all_pairs, 
            self.pca_selected_comparisons,
            title="Select PCA Pairwise Comparisons"
        )
        
        if result is not None:
            self.pca_selected_comparisons = result
            if not result:
                self.pca_comp_status_label.config(text="⚠️ No pairwise comparisons selected - only all-groups PCA will be generated")
            elif len(result) == len(all_pairs):
                self.pca_comp_status_label.config(text="All pairwise comparisons will be plotted")
            else:
                self.pca_comp_status_label.config(text=f"Selected: {len(result)}/{len(all_pairs)} pairwise comparisons")

    def _show_comparison_selector(self, all_pairs, groups, title="Create Venn Diagrams - Multi-Selection Mode"):
        """Show dialog to create multiple Venn configurations (for 4+ groups)."""
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("850x650")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Header
        header = ttk.Label(dialog, 
                          text=f"You have {len(groups)} groups → {len(all_pairs)} possible comparisons.\n"
                               f"Create one or more Venn diagrams (2-4 comparisons each).\n"
                               f"Give each Venn a name, select comparisons, click 'Add This Venn', then create another or 'Done'.",
                          font=('TkDefaultFont', 10, 'bold'),
                          foreground='darkblue',
                          justify='left')
        header.pack(pady=10, padx=10, anchor='w')
        
        # Main top row: Quick Filters on left, Action buttons on right
        top_row = ttk.Frame(dialog)
        top_row.pack(fill='x', padx=10, pady=(0,10))
        
        # Quick filters section (left side)
        filter_frame = ttk.LabelFrame(top_row, text="⚡ Quick Filters", padding=10)
        filter_frame.pack(side='left', fill='both', expand=True)
        
        # Quick filter functions (defined early for button commands)
        def filter_by_group(group_name):
            """Select only comparisons involving a specific group."""
            for (g1, g2), var in pair_vars.items():
                if g1 == group_name or g2 == group_name:
                    var.set(True)
                else:
                    var.set(False)
            update_count_label()
        
        def select_all():
            for var in pair_vars.values():
                var.set(True)
            update_count_label()
        
        def deselect_all():
            for var in pair_vars.values():
                var.set(False)
            update_count_label()
        
        def select_first_n(n):
            """Select first n comparisons."""
            for i, var in enumerate(pair_vars.values()):
                var.set(i < n)
            update_count_label()
        
        # Add filter buttons
        ttk.Button(filter_frame, text="✓ Select All", 
                  command=select_all).grid(row=0, column=0, padx=2, pady=2, sticky='ew')
        ttk.Button(filter_frame, text="✗ Deselect All", 
                  command=deselect_all).grid(row=0, column=1, padx=2, pady=2, sticky='ew')
        ttk.Button(filter_frame, text="First 3", 
                  command=lambda: select_first_n(3)).grid(row=0, column=2, padx=2, pady=2, sticky='ew')
        ttk.Button(filter_frame, text="First 2", 
                  command=lambda: select_first_n(2)).grid(row=0, column=3, padx=2, pady=2, sticky='ew')
        
        # Add group filter buttons
        ttk.Label(filter_frame, text="Include only comparisons with:").grid(row=1, column=0, columnspan=4, sticky='w', pady=(5,2))
        for i, group in enumerate(groups):
            btn = ttk.Button(filter_frame, 
                           text=f"{group}",
                           command=lambda g=group: filter_by_group(g))
            btn.grid(row=2, column=i % 4, padx=2, pady=2, sticky='ew')
        
        # Action buttons section (right side)
        button_container = ttk.Frame(top_row, padding=10)
        button_container.pack(side='right', fill='y')
        
        ttk.Label(button_container, text="Actions:", font=('TkDefaultFont', 10, 'bold')).pack(anchor='w', pady=(0,8))
        
        # Venn name entry
        name_frame = ttk.Frame(button_container)
        name_frame.pack(fill='x', pady=(0,8))
        ttk.Label(name_frame, text="Venn Name:").pack(anchor='w')
        venn_name_var = tk.StringVar(value="")
        ttk.Entry(name_frame, textvariable=venn_name_var, width=18).pack(fill='x', pady=2)
        
        # Accumulated Venns list
        accumulated_venns = []
        
        # Display canonicalization helper: keep internal pair tuples as-is but
        # always display the base group first (user expects 'Control' first).
        BASE_VENN_GROUP = 'control'
        def _canon_pair(a, b):
            if a.lower() == BASE_VENN_GROUP or b.lower() == BASE_VENN_GROUP:
                if a.lower() == BASE_VENN_GROUP:
                    return a, b
                else:
                    return b, a
            else:
                # Alphabetical order for consistency
                return (a, b) if a < b else (b, a)
        
        def on_add_venn():
            """Add current selection as a new Venn spec."""
            selected = [(g1, g2) for (g1, g2), var in pair_vars.items() if var.get()]
            if not selected:
                messagebox.showwarning("No Selection", "Please select at least one comparison.")
                return
            
            # Validate count (now supports up to 4 comparisons!)
            if len(selected) > 4:
                messagebox.showwarning("Too Many Comparisons", 
                    f"You selected {len(selected)} comparisons.\n"
                    f"Venn diagrams support 2-4 comparisons.\n"
                    f"Please deselect some comparisons.")
                return
            
            # Get venn name
            venn_name = venn_name_var.get().strip()
            if not venn_name:
                venn_name = f"Venn{len(accumulated_venns)+1}"
            
            # Check for duplicate names
            if any(v['name'] == venn_name for v in accumulated_venns):
                messagebox.showwarning("Duplicate Name", 
                    f"A Venn named '{venn_name}' already exists.\nPlease use a different name.")
                return
            
            # Get indices of selected pairs in all_pairs
            indices = []
            for i, pair in enumerate(all_pairs):
                if pair_vars.get(pair, tk.BooleanVar(value=False)).get():
                    indices.append(i+1)
            
            accumulated_venns.append({
                'name': venn_name,
                'pairs': selected.copy(),
                'indices': indices.copy()
            })
            
            # Update display
            update_accumulated_display()
            
            # Clear for next
            venn_name_var.set("")
            for var in pair_vars.values():
                var.set(False)
            update_count_label()
            
            self.log_viz_message(f"Added Venn '{venn_name}' with {len(selected)} comparisons")
        
        def on_done():
            """Finish and apply all accumulated Venns."""
            if not accumulated_venns:
                response = messagebox.askyesno("No Venns Created", 
                    "You haven't added any Venn configurations yet.\nClose anyway?")
                if response:
                    dialog.destroy()
                return
            dialog.destroy()
        
        def on_cancel():
            """Cancel all changes."""
            accumulated_venns.clear()
            dialog.destroy()
        
        ttk.Button(button_container, text="➕ Add This Venn", 
                  command=on_add_venn, width=18).pack(pady=3, fill='x')
        ttk.Button(button_container, text="✓ Done", 
                  command=on_done, width=18).pack(pady=3, fill='x')
        ttk.Button(button_container, text="✗ Cancel All", 
                  command=on_cancel, width=18).pack(pady=3, fill='x')
        
        # Selection area with scrollbar
        selection_frame = ttk.LabelFrame(dialog, text="Select Comparisons for Current Venn", padding=10)
        selection_frame.pack(fill='both', expand=True, padx=10, pady=(0,10))
        
        canvas = tk.Canvas(selection_frame)
        scrollbar = ttk.Scrollbar(selection_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Checkbox variables for each pair
        pair_vars = {}
        for i, (g1, g2) in enumerate(all_pairs):
            var = tk.BooleanVar(value=False)  # Start unchecked
            pair_vars[(g1, g2)] = var
            display_a, display_b = _canon_pair(g1, g2)
            cb = ttk.Checkbutton(scrollable_frame, 
                               text=f"{display_a} vs {display_b}",
                               variable=var)
            cb.grid(row=i, column=0, sticky='w', pady=2, padx=5)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Count label
        count_label = ttk.Label(dialog, text="", foreground='darkgreen', font=('TkDefaultFont', 9, 'bold'))
        count_label.pack(pady=5)
        
        def update_count_label():
            selected_count = sum(1 for var in pair_vars.values() if var.get())
            if selected_count == 0:
                count_label.configure(text="No comparisons selected", foreground='gray')
            elif selected_count <= 4:
                count_label.configure(text=f"✓ {selected_count} comparison(s) selected - Perfect for Venn!", foreground='darkgreen')
            else:
                count_label.configure(text=f"⚠️ {selected_count} selected - Too many! Venn supports 2-4 only", foreground='red')
        
        # Bind checkbox changes
        for var in pair_vars.values():
            var.trace_add('write', lambda *args: update_count_label())
        
        update_count_label()
        
        # Accumulated Venns display
        accumulated_frame = ttk.LabelFrame(dialog, text="📊 Accumulated Venn Diagrams", padding=10)
        accumulated_frame.pack(fill='x', padx=10, pady=(0,10))
        
        accumulated_text = tk.Text(accumulated_frame, height=5, wrap='word', state='disabled',
                                   font=('TkDefaultFont', 9))
        accumulated_text.pack(fill='x')
        
        def update_accumulated_display():
            """Update the display of accumulated Venns."""
            accumulated_text.config(state='normal')
            accumulated_text.delete('1.0', 'end')
            if accumulated_venns:
                for venn in accumulated_venns:
                    comparisons = ', '.join([f"{_canon_pair(a,b)[0]} vs {_canon_pair(a,b)[1]}" for a, b in venn['pairs']])
                    accumulated_text.insert('end', f"• {venn['name']}: {comparisons}\n")
            else:
                accumulated_text.insert('end', "No Venns added yet. Select comparisons above and click 'Add This Venn'.")
            accumulated_text.config(state='disabled')
        
        update_accumulated_display()
        
        # Wait for dialog
        dialog.wait_window()
        
        # Return accumulated venns
        return accumulated_venns if accumulated_venns else None

    def _run_venn_from_gui(self):
        """Legacy method - now redirects to _generate_venn_plots."""
        try:
            if not self._confirm_venn_filter_alignment():
                self.log_viz_message("⚠️ Venn generation cancelled: filter confirmation not accepted.")
                return
            if not getattr(self, 'viz_group_mapping', None):
                self.log_viz_message("❌ Error: No group configuration found. Please configure groups first.")
                return
            complete_df = self.get_current_complete_df()
            if complete_df is None or complete_df.empty:
                self.log_viz_message("❌ Error: No data available. Please load or compute statistics first.")
                return
            groups = self.get_viz_groups()
            sample_to_group = self.get_viz_sample_to_group_mapping()
            sample_cols = list(sample_to_group.keys())
            output_dir = self.viz_output_dir.get()
            os.makedirs(output_dir, exist_ok=True)
            self._generate_venn_plots(complete_df, groups, sample_cols, sample_to_group, output_dir)
        except Exception as e:
            self.log_viz_message(f"Venn error: {e}")

    def _generate_venn_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir):
        """Generate Venn diagrams based on configured specs.
        
        Generates filtered Venns (p-value/FC cutoffs) and optionally All Molecules Venns
        if the user has checked the "Generate All Molecules" option.
        """
        try:
            # Debug: Log sample mapping to verify FC columns are included
            logger.info(f"🔍 Venn sample mapping debug:")
            logger.info(f"   └─ Total sample columns: {len(sample_cols)}")
            for grp in groups:
                grp_cols = [c for c in sample_cols if sample_to_group.get(c) == grp]
                logger.info(f"   └─ {grp}: {len(grp_cols)} samples = {grp_cols}")
            
            from main_script.metabolites_visualization import run_venn_analysis, VennParams, VennSpec
            ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'venn'))
            os.makedirs(ctx.output_dir, exist_ok=True)
            
            # Gather params
            generate_all_molecules = self.venn_generate_all_molecules.get() if hasattr(self, 'venn_generate_all_molecules') else False
            
            # Standard cutoff parameters
            pval = self.venn_p_thresh.get() if hasattr(self, 'venn_p_thresh') else 0.05
            fc = 0.0 if (hasattr(self, 'venn_skip_fc') and self.venn_skip_fc.get()) else (self.venn_fc_thresh.get() if hasattr(self, 'venn_fc_thresh') else 2.0)
            
            # All Molecules parameters
            allmol_min_presence_type = self.venn_allmol_min_presence_type.get() if hasattr(self, 'venn_allmol_min_presence_type') else 'count'
            allmol_min_presence_count = self.venn_allmol_min_presence_count.get() if hasattr(self, 'venn_allmol_min_presence_count') else 3
            allmol_min_presence_percent = self.venn_allmol_min_presence_percent.get() if hasattr(self, 'venn_allmol_min_presence_percent') else 50.0
            
            # Debug: log what parameters were read from GUI
            self.log_viz_message(f"DEBUG: Venn All Molecules params from GUI:")
            self.log_viz_message(f"  └─ min_presence_type = '{allmol_min_presence_type}'")
            self.log_viz_message(f"  └─ min_presence_count = {allmol_min_presence_count}")
            self.log_viz_message(f"  └─ min_presence_percent = {allmol_min_presence_percent}")
            
            # Font sizes
            venn_number_fontsize = self.venn_number_fontsize.get() if hasattr(self, 'venn_number_fontsize') else 16
            venn_label_fontsize = self.venn_label_fontsize.get() if hasattr(self, 'venn_label_fontsize') else 11
            
            # Build specs from indices (all user-configured Venn specs)
            specs = []
            for spec in getattr(self, 'venn_specs', []):
                comps = []
                for idx in spec.get('indices', []):
                    if 1 <= idx <= len(self.venn_pairs_list):
                        comps.append(self.venn_pairs_list[idx-1])
                if comps:
                    specs.append(VennSpec(name=spec.get('name', f"Venn{len(specs)+1}"), comparisons=comps))
            
            # Check if user has configured Venn specs using the setup button
            if not specs:
                self.log_viz_message('❌ Error: No Venn diagrams configured.')
                self.log_viz_message('   Please use the "Setup Venns" button to configure Venn diagrams before generating.')
                self.log_viz_message('   Or uncheck "Generate Venn diagrams" if you don\'t need them.')
                return
            
            all_files_created = []
            all_errors = []
            
            # Count unique groups involved in comparisons to detect single comparison scenario
            unique_groups = set()
            for spec in specs:
                for g1, g2 in spec.comparisons:
                    unique_groups.add(g1)
                    unique_groups.add(g2)
            
            # For single comparison (2 groups), skip filtered Venn and generate only All Molecules
            skip_filtered_venn = len(unique_groups) <= 2
            
            # ============================================================
            # 1. Generate FILTERED Venns (p-value/FC cutoffs) - Skip if only 2 groups
            # ============================================================
            if not skip_filtered_venn:
                # Get output format from GUI
                venn_output_format = getattr(self, 'venn_output_format', tk.StringVar(value='png')).get()
                
                params_filtered = VennParams(
                    p_threshold=pval, 
                    fc_threshold=fc, 
                    venn_specs=specs,
                    venn_number_fontsize=venn_number_fontsize,
                    venn_label_fontsize=venn_label_fontsize,
                    skip_all_cutoffs=False,
                    min_presence_type='count',
                    min_presence_count=1,
                    min_presence_percent=0.0,
                    output_format=venn_output_format
                )
                
                self.log_viz_message(f"Venn: Generating {len(specs)} filtered Venn diagram(s) with p<{pval}, FC≥{fc if fc > 0 else 'none'}")
                
                result_filtered = run_venn_analysis(ctx, params_filtered)
                
                # Log summary
                self.log_viz_message(result_filtered.summary)
                if hasattr(result_filtered, 'venn_summaries') and result_filtered.venn_summaries:
                    for venn_summary in result_filtered.venn_summaries:
                        self.log_viz_message(f"  └─ {venn_summary}")
                
                all_files_created.extend(result_filtered.files_created)
                all_errors.extend(result_filtered.errors)
            else:
                self.log_viz_message('ℹ️ Venn: Skipping filtered Venn - only 1 comparison (2 groups) detected.')
                self.log_viz_message('   Filtered Venn would create empty groups for a single comparison.')
            
            # ============================================================
            # 2. Generate ALL MOLECULES Venns - ONLY IF CHECKBOX IS CHECKED
            # ============================================================
            if generate_all_molecules:
                # Build All Molecules specs from separate configuration
                specs_allmol = []
                allmol_pairs = getattr(self, 'venn_allmol_pairs_list', [])
                allmol_specs_config = getattr(self, 'venn_allmol_specs', [])
                
                for spec in allmol_specs_config:
                    comps = []
                    for idx in spec.get('indices', []):
                        if 1 <= idx <= len(allmol_pairs):
                            comps.append(allmol_pairs[idx-1])
                    if comps:
                        specs_allmol.append(VennSpec(name=spec.get('name', f"AllMolVenn{len(specs_allmol)+1}"), comparisons=comps))
                
                if not specs_allmol:
                    self.log_viz_message('⚠️ All Molecules Venn: No configurations defined. Skipping.')
                else:
                    # Get output format from GUI
                    venn_output_format = getattr(self, 'venn_output_format', tk.StringVar(value='png')).get()
                    
                    params_allmol = VennParams(
                        p_threshold=1.0,  # Disable p-value filtering
                        fc_threshold=0.0,  # Disable FC filtering
                        venn_specs=specs_allmol,
                        venn_number_fontsize=venn_number_fontsize,
                        venn_label_fontsize=venn_label_fontsize,
                        skip_all_cutoffs=True,
                        min_presence_type=allmol_min_presence_type,
                        min_presence_count=allmol_min_presence_count,
                        min_presence_percent=allmol_min_presence_percent,
                        output_format=venn_output_format
                    )
                    
                    allmol_desc = f"{allmol_min_presence_count} samples" if allmol_min_presence_type == 'count' else f"{allmol_min_presence_percent}%"
                    self.log_viz_message(f"Venn: Generating {len(specs_allmol)} All Molecules Venn diagram(s) (min {allmol_desc})")
                    
                    result_allmol = run_venn_analysis(ctx, params_allmol)
                    
                    # Log summary
                    self.log_viz_message(result_allmol.summary)
                    if hasattr(result_allmol, 'venn_summaries') and result_allmol.venn_summaries:
                        for venn_summary in result_allmol.venn_summaries:
                            self.log_viz_message(f"  └─ {venn_summary}")
                    
                    all_files_created.extend(result_allmol.files_created)
                    all_errors.extend(result_allmol.errors)
            
            # ============================================================
            # Final Summary
            # ============================================================
            if all_files_created:
                self.log_viz_message(f"Venn: Total {len(all_files_created)} files saved to {ctx.output_dir}")
            
            if all_errors:
                for e in all_errors:
                    self.log_viz_message(f"Venn warning: {e}")
        except Exception as e:
            import traceback
            self.log_viz_message(f"❌ Venn generation failed: {e}")
            self.log_viz_message(traceback.format_exc())

    def create_heatmap_panel(self, parent):
        """Create Heatmap parameter panel - REORGANIZED VERSION."""
        panel = ttk.LabelFrame(parent, text="Heatmaps", padding=10)
        panel.pack(fill='both', expand=True, padx=5, pady=2)

        ttk.Checkbutton(panel, text="Generate heatmaps",
                        variable=self.viz_selected['heatmap']).pack(anchor='w')

        # Comparison and Metabolite List controls side-by-side (MOVED UP for visibility)
        comp_metab_frame = ttk.LabelFrame(panel, text="🔍 Comparison & Metabolite Configuration", padding=5)
        comp_metab_frame.pack(fill='x', pady=4)
        
        # Create two columns
        left_col = ttk.Frame(comp_metab_frame)
        left_col.pack(side='left', fill='both', expand=True, padx=(0,5))
        right_col = ttk.Frame(comp_metab_frame)
        right_col.pack(side='right', fill='both', expand=True, padx=(5,0))
        
        # Left: Comparison Selection
        ttk.Label(left_col, text="Comparison Selection", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.heatmap_selected_comparisons = None  # Will store list of (g1, g2) tuples
        self.heatmap_comp_status_label = ttk.Label(left_col, text="All comparisons will be plotted", 
                                                font=('Arial', 8), foreground='#666')
        self.heatmap_comp_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(left_col, text="Configure Comparisons...", 
                command=self._configure_heatmap_comparisons).pack(anchor='w')
        
        # Right: Per-Comparison Metabolite Lists
        ttk.Label(right_col, text="Per-Comparison Lists", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.heatmap_metabolite_lists = {}  # Will store {(g1, g2): filepath}
        self.heatmap_skip_unlisted = tk.BooleanVar(value=False)
        info_label = ttk.Label(right_col, text="Upload different lists for comparisons", 
                            font=('Arial', 8, 'italic'), foreground='#666')
        info_label.pack(anchor='w', pady=(0, 5))
        ttk.Checkbutton(right_col, text="Skip comparisons without lists", 
                    variable=self.heatmap_skip_unlisted).pack(anchor='w', pady=(0, 5))
        ttk.Button(right_col, text="Manage Metabolite Lists...", 
                command=self._configure_heatmap_metabolite_lists).pack(anchor='w')
        
        # Custom Metabolite List - Top Right
        custom_list_frame = ttk.LabelFrame(panel, text="Custom Metabolite List", padding=5)
        custom_list_frame.pack(fill='x', pady=4)
        self.heatmap_custom_list = tk.StringVar()
        list_row = ttk.Frame(custom_list_frame)
        list_row.pack(fill='x', pady=2)
        ttk.Entry(list_row, textvariable=self.heatmap_custom_list, width=40).pack(side='left', expand=True, fill='x', padx=(0,5))
        ttk.Button(list_row, text="Browse", command=lambda: self.browse_metabolite_list('heatmap')).pack(side='right')
        self.heatmap_use_custom_only = tk.BooleanVar(value=False)
        ttk.Checkbutton(custom_list_frame, text='Use custom list ONLY (ignore p-value/FC filters)', 
                    variable=self.heatmap_use_custom_only).pack(anchor='w', pady=(2,0))

        # Basic parameters
        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        ttk.Label(params_frame, text="Max metabolites:").grid(row=0, column=0, sticky='w', padx=2)
        self.heatmap_max = tk.IntVar(value=0)
        ttk.Spinbox(params_frame, from_=0, to=500, textvariable=self.heatmap_max, width=10).grid(row=0, column=1, padx=2)
        ttk.Label(params_frame, text="(0=all)").grid(row=0, column=2, sticky='w')
        
        # Divider line options
        divider_frame = ttk.LabelFrame(panel, text="Divider Lines", padding=5)
        divider_frame.pack(fill='x', pady=4)
        self.heatmap_show_fc_divider = tk.BooleanVar(value=True)
        ttk.Checkbutton(divider_frame, text="Show dividing line between up/down regulated", 
                    variable=self.heatmap_show_fc_divider).pack(anchor='w', padx=2, pady=2)
        self.heatmap_show_sample_divider = tk.BooleanVar(value=True)
        ttk.Checkbutton(divider_frame, text="Show dividing line between sample groups", 
                    variable=self.heatmap_show_sample_divider).pack(anchor='w', padx=2, pady=2)

        # Significance thresholds frame
        sig_frame = ttk.LabelFrame(panel, text="Significance Filters")
        sig_frame.pack(fill='x', pady=4)
        ttk.Label(sig_frame, text="P-value <").grid(row=0, column=0, sticky='w', padx=2)
        self.heatmap_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(sig_frame, textvariable=self.heatmap_p_thresh, width=8).grid(row=0, column=1, padx=2)
        ttk.Label(sig_frame, text="|Fold Change| ≥").grid(row=0, column=2, sticky='w', padx=8)
        self.heatmap_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(sig_frame, textvariable=self.heatmap_fc_thresh, width=8).grid(row=0, column=3, padx=2)
        # Option to skip fold-change cutoff for heatmaps
        self.heatmap_skip_fc = tk.BooleanVar(value=True)
        def _on_heatmap_skip_fc_change():
            try:
                state = 'disabled' if self.heatmap_skip_fc.get() else 'normal'
                for w in sig_frame.grid_slaves(row=0, column=3):
                    try:
                        w.configure(state=state)
                    except Exception:
                        pass
            except Exception:
                pass
        ttk.Checkbutton(sig_frame, text="Skip FC cutoff", variable=self.heatmap_skip_fc, command=_on_heatmap_skip_fc_change).grid(row=0, column=4, padx=4, sticky='w')
        _on_heatmap_skip_fc_change()
        ttk.Label(sig_frame, text="Mode:").grid(row=1, column=0, sticky='w', padx=2)
        self.heatmap_filter_mode = tk.StringVar(value='any')
        ttk.Radiobutton(sig_frame, text='Any', value='any', variable=self.heatmap_filter_mode).grid(row=1, column=1, sticky='w')
        ttk.Radiobutton(sig_frame, text='All', value='all', variable=self.heatmap_filter_mode).grid(row=1, column=2, sticky='w')
        ttk.Radiobutton(sig_frame, text='Specific', value='specific', variable=self.heatmap_filter_mode).grid(row=1, column=3, sticky='w', padx=0)
        ttk.Label(sig_frame, text="Comparison:", font=('Arial', 8)).grid(row=1, column=4, sticky='w', padx=0)
        self.heatmap_specific_comparison = tk.StringVar(value="")
        specific_entry = ttk.Entry(sig_frame, textvariable=self.heatmap_specific_comparison, width=15)
        specific_entry.grid(row=1, column=5, sticky='ew', padx=0)
        self._create_tooltip(specific_entry, "Enter comparison (e.g., 'GroupA|GroupB' or 'AD|Control')\nUsed when Mode='Specific'")

        # Color scaling frame
        scale_frame = ttk.LabelFrame(panel, text="Color Scaling")
        scale_frame.pack(fill='x', pady=4)
        
        # Default scale (-3 to 3) - checked by default
        self.heatmap_use_fixed_scale = tk.BooleanVar(value=True)
        ttk.Checkbutton(scale_frame, text="Use default scale (-3 to 3)", 
                       variable=self.heatmap_use_fixed_scale).grid(row=0, column=0, columnspan=4, sticky='w', padx=2, pady=2)
        
        # Auto scale option
        self.heatmap_auto_scale = tk.BooleanVar(value=False)
        ttk.Checkbutton(scale_frame, text="Auto scale (5th-95th percentile)", 
                       variable=self.heatmap_auto_scale).grid(row=1, column=0, columnspan=4, sticky='w', padx=2, pady=2)
        
        # Manual scale option
        ttk.Label(scale_frame, text="Manual scale:", font=('', 9, 'bold')).grid(row=2, column=0, columnspan=4, sticky='w', padx=2, pady=(8,2))
        ttk.Label(scale_frame, text="vmin:").grid(row=3, column=0, sticky='w', padx=2)
        self.heatmap_vmin = tk.DoubleVar(value=-3.0)
        ttk.Spinbox(scale_frame, from_=-20, to=0, increment=0.5, textvariable=self.heatmap_vmin, width=6).grid(row=3, column=1, padx=2)
        ttk.Label(scale_frame, text="vmax:").grid(row=3, column=2, sticky='w', padx=(12,2))
        self.heatmap_vmax = tk.DoubleVar(value=3.0)
        ttk.Spinbox(scale_frame, from_=0, to=20, increment=0.5, textvariable=self.heatmap_vmax, width=6).grid(row=3, column=3, padx=2)
        ttk.Label(scale_frame, text="(Used when both above are unchecked)", font=('', 8, 'italic')).grid(row=4, column=0, columnspan=4, sticky='w', padx=2, pady=(0,4))

        # Figure size controls
        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        self.heatmap_auto_size = tk.BooleanVar(value=True)
        ttk.Checkbutton(size_frame, text="Auto-size (recommended)", variable=self.heatmap_auto_size).grid(row=0, column=0, columnspan=2, sticky='w', pady=(0,4))
        self.heatmap_fig_width = tk.DoubleVar(value=8.0)
        self.heatmap_fig_height = tk.DoubleVar(value=6.0)
        self.heatmap_fig_dpi = tk.IntVar(value=190)
        ttk.Label(size_frame, text="Width:").grid(row=1, column=0, padx=2, sticky='w')
        width_spin = ttk.Spinbox(size_frame, from_=2.0, to=40.0, increment=0.5, textvariable=self.heatmap_fig_width, width=6)
        width_spin.grid(row=1, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=1, column=2, padx=8, sticky='w')
        height_spin = ttk.Spinbox(size_frame, from_=2.0, to=40.0, increment=0.5, textvariable=self.heatmap_fig_height, width=6)
        height_spin.grid(row=1, column=3, padx=2)
        ttk.Label(size_frame, text="DPI:").grid(row=1, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.heatmap_fig_dpi, width=6).grid(row=1, column=5, padx=2)
        
        # Function to enable/disable manual size controls
        def toggle_size_controls():
            state = 'disabled' if self.heatmap_auto_size.get() else 'normal'
            width_spin.config(state=state)
            height_spin.config(state=state)
        
        self.heatmap_auto_size.trace_add('write', lambda *args: toggle_size_controls())
        toggle_size_controls()  # Initial state

        # Layout controls for dendrogram and colorbar
        layout_frame = ttk.LabelFrame(panel, text="Layout (Dendrogram & Colorbar)", padding=5)
        layout_frame.pack(fill='x', pady=4)
        
        # Clustering control - OFF by default so split lines work correctly
        self.heatmap_cluster = tk.BooleanVar(value=False)
        ttk.Checkbutton(layout_frame, text="Cluster rows (hierarchical clustering)", 
                       variable=self.heatmap_cluster).grid(row=0, column=0, columnspan=2, sticky='w', pady=(0,5))
        
        # Show colorbar by default
        self.heatmap_show_colorbar = tk.BooleanVar(value=True)
        show_cbar_cb = ttk.Checkbutton(layout_frame, text="Show colorbar (top)", variable=self.heatmap_show_colorbar)
        show_cbar_cb.grid(row=1, column=0, columnspan=2, sticky='w')
        # Dendrogram width (percent of heatmap width column in GridSpec ratios)
        ttk.Label(layout_frame, text="Dendrogram width (%):").grid(row=2, column=0, padx=2, sticky='w')
        self.heatmap_dendro_width_pct = tk.DoubleVar(value=10.0)
        dendro_spin = ttk.Spinbox(layout_frame, from_=0.0, to=40.0, increment=1.0, textvariable=self.heatmap_dendro_width_pct, width=6)
        dendro_spin.grid(row=2, column=1, padx=2)
        # Colorbar height (in inches - fixed distance from heatmap)
        ttk.Label(layout_frame, text="Colorbar height (inches):").grid(row=2, column=2, padx=8, sticky='w')
        self.heatmap_cbar_height_inches = tk.DoubleVar(value=0.3)
        cbar_spin = ttk.Spinbox(layout_frame, from_=0.3, to=3.0, increment=0.1, textvariable=self.heatmap_cbar_height_inches, width=6)
        cbar_spin.grid(row=2, column=3, padx=2)
        # Disable colorbar height control when colorbar is hidden
        def toggle_cbar_controls(*_):
            state = 'normal' if self.heatmap_show_colorbar.get() else 'disabled'
            try:
                cbar_spin.config(state=state)
            except Exception:
                pass
        self.heatmap_show_colorbar.trace_add('write', toggle_cbar_controls)
        toggle_cbar_controls()

        # Font size controls - simplified for heatmaps
        font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        font_frame.pack(fill='x', pady=4)
        self.heatmap_feature_fontsize = tk.IntVar(value=10)  # For rows (metabolites/features)
        self.heatmap_sample_fontsize = tk.IntVar(value=10)   # For columns (samples)
        self.heatmap_title_fontsize = tk.IntVar(value=14)    # For plot title
        
        ttk.Label(font_frame, text="Feature labels (rows):").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=20, textvariable=self.heatmap_feature_fontsize, width=5).grid(row=0, column=1, padx=2)
        
        ttk.Label(font_frame, text="Sample labels (columns):").grid(row=0, column=2, padx=12, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=20, textvariable=self.heatmap_sample_fontsize, width=5).grid(row=0, column=3, padx=2)
        
        ttk.Label(font_frame, text="Title:").grid(row=1, column=0, padx=2, sticky='w', pady=(4,0))
        ttk.Spinbox(font_frame, from_=8, to=24, textvariable=self.heatmap_title_fontsize, width=5).grid(row=1, column=1, padx=2, pady=(4,0))

        # Excel export control 
        save_frame = ttk.LabelFrame(panel, text="💾 Save Options", padding=5)
        save_frame.pack(fill='x', pady=4)
        self.heatmap_save_excel = tk.BooleanVar(value=False)
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.heatmap_save_excel).pack(anchor='w')
        
        # Output format selection
        format_frame = ttk.LabelFrame(panel, text="Output Format", padding=5)
        format_frame.pack(fill='x', pady=4)
        self.heatmap_output_format = tk.StringVar(value='png')
        ttk.Radiobutton(format_frame, text="PNG (Raster)", variable=self.heatmap_output_format, value='png').pack(anchor='w')
        ttk.Radiobutton(format_frame, text="SVG (Vector)", variable=self.heatmap_output_format, value='svg').pack(anchor='w')

    def create_boxplot_panel(self, parent):
        """Create Boxplot parameter panel - REORGANIZED VERSION."""
        panel = ttk.LabelFrame(parent, text="Boxplots", padding=10)
        panel.pack(fill='both', expand=True, padx=5, pady=2)

        ttk.Checkbutton(panel, text="Generate boxplots",
                        variable=self.viz_selected['boxplot']).pack(anchor='w')

        # Group Selection and Annotation Configuration - Side by Side (MOVED UP for visibility)
        comp_metab_frame = ttk.LabelFrame(panel, text="🔍 Group & Annotation Configuration", padding=5)
        comp_metab_frame.pack(fill='x', pady=4)
        
        # Create two columns
        left_column = ttk.Frame(comp_metab_frame)
        left_column.pack(side='left', fill='both', expand=True, padx=(0, 5))
        
        right_column = ttk.Frame(comp_metab_frame)
        right_column.pack(side='left', fill='both', expand=True, padx=(5, 0))
        
        # LEFT: Group selection for boxplots
        ttk.Label(left_column, text="Groups to Display", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.boxplot_selected_groups = None  # Will store list of group names to display
        self.boxplot_groups_status_label = ttk.Label(left_column, text="All groups will be displayed", 
                                font=('Arial', 8), foreground='#666')
        self.boxplot_groups_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(left_column, text="Set Groups...", 
                   command=self._configure_boxplot_groups).pack(anchor='w')
        
        # RIGHT: Annotation comparison selection
        ttk.Label(right_column, text="Annotation Comparisons", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.boxplot_annotate_comparisons = None
        self.boxplot_annotate_status_label = ttk.Label(right_column, text="All will be annotated", 
                                font=('Arial', 8), foreground='#666')
        self.boxplot_annotate_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(right_column, text="Configure Annotation...", 
                   command=self._configure_boxplot_annotations).pack(anchor='w')
        
        # Custom Metabolite List - Top Right
        custom_list_frame = ttk.LabelFrame(panel, text="Custom Metabolite List", padding=5)
        custom_list_frame.pack(fill='x', pady=4)
        self.boxplot_custom_list = tk.StringVar()
        list_row = ttk.Frame(custom_list_frame)
        list_row.pack(fill='x', pady=2)
        ttk.Entry(list_row, textvariable=self.boxplot_custom_list, width=40).pack(side='left', expand=True, fill='x', padx=(0,5))
        ttk.Button(list_row, text="Browse", command=lambda: self.browse_metabolite_list('boxplot')).pack(side='right')
        self.boxplot_use_custom_only = tk.BooleanVar(value=False)
        ttk.Checkbutton(custom_list_frame, text='Use custom list ONLY (ignore p-value/FC filters)', 
                    variable=self.boxplot_use_custom_only).pack(anchor='w', pady=(2,0))

        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        ttk.Label(params_frame, text="Top N metabolites:").grid(row=0, column=0, sticky='w', padx=2)
        self.boxplot_top_n = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=1, to=100, textvariable=self.boxplot_top_n, width=8).grid(row=0, column=1, padx=2)
        self.boxplot_no_limit = tk.BooleanVar(value=False)
        ttk.Checkbutton(params_frame, text="No limit (plot all)", variable=self.boxplot_no_limit).grid(row=0, column=2, padx=10)
        self.boxplot_annotate = tk.BooleanVar(value=True)
        ttk.Checkbutton(params_frame, text="Show significance stars", variable=self.boxplot_annotate).grid(row=1, column=0, columnspan=3, sticky='w', pady=5)

        # Significance filters
        sig_frame = ttk.LabelFrame(panel, text="Significance Filters")
        sig_frame.pack(fill='x', pady=4)
        ttk.Label(sig_frame, text="P-value <").grid(row=0, column=0, sticky='w', padx=2)
        self.boxplot_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(sig_frame, textvariable=self.boxplot_p_thresh, width=8).grid(row=0, column=1, padx=2)
        ttk.Label(sig_frame, text="|Fold Change| ≥").grid(row=0, column=2, sticky='w', padx=8)
        self.boxplot_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(sig_frame, textvariable=self.boxplot_fc_thresh, width=8).grid(row=0, column=3, padx=2)
        self.boxplot_skip_fc = tk.BooleanVar(value=True)
        def _on_boxplot_skip_fc_change():
            try:
                state = 'disabled' if self.boxplot_skip_fc.get() else 'normal'
                for w in sig_frame.grid_slaves(row=0, column=3):
                    try:
                        w.configure(state=state)
                    except Exception:
                        pass
            except Exception:
                pass
        ttk.Checkbutton(sig_frame, text="Skip FC cutoff", variable=self.boxplot_skip_fc, command=_on_boxplot_skip_fc_change).grid(row=0, column=4, padx=4, sticky='w')
        _on_boxplot_skip_fc_change()
        ttk.Label(sig_frame, text="Mode:").grid(row=1, column=0, sticky='w', padx=2)
        self.boxplot_filter_mode = tk.StringVar(value='any')
        ttk.Radiobutton(sig_frame, text='Any', value='any', variable=self.boxplot_filter_mode).grid(row=1, column=1, sticky='w')
        ttk.Radiobutton(sig_frame, text='All', value='all', variable=self.boxplot_filter_mode).grid(row=1, column=2, sticky='w')
        ttk.Radiobutton(sig_frame, text='Specific', value='specific', variable=self.boxplot_filter_mode).grid(row=1, column=3, sticky='w', padx=0)
        ttk.Label(sig_frame, text="Comparison:", font=('Arial', 8)).grid(row=1, column=4, sticky='w', padx=0)
        self.boxplot_specific_comparison = tk.StringVar(value="")
        specific_entry = ttk.Entry(sig_frame, textvariable=self.boxplot_specific_comparison, width=15)
        specific_entry.grid(row=1, column=5, sticky='ew', padx=0)
        self._create_tooltip(specific_entry, "Enter comparison (e.g., 'GroupA|GroupB' or 'AD|Control')\nUsed when Mode='Specific'")

        # Figure size controls
        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        self.boxplot_fig_width = tk.DoubleVar(value=3.0)
        self.boxplot_fig_height = tk.DoubleVar(value=3.0)
        self.boxplot_fig_dpi = tk.IntVar(value=240)
        ttk.Label(size_frame, text="Width:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.boxplot_fig_width, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.boxplot_fig_height, width=6).grid(row=0, column=3, padx=2)
        ttk.Label(size_frame, text="DPI:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.boxplot_fig_dpi, width=6).grid(row=0, column=5, padx=2)

        # Font size controls
        font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        font_frame.pack(fill='x', pady=4)
        self.boxplot_xlabel_fontsize = tk.IntVar(value=12)
        self.boxplot_ylabel_fontsize = tk.IntVar(value=10)
        self.boxplot_title_fontsize = tk.IntVar(value=14)
        self.boxplot_tick_fontsize = tk.IntVar(value=10)
        
        ttk.Label(font_frame, text="X-label:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.boxplot_xlabel_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(font_frame, text="Y-label:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.boxplot_ylabel_fontsize, width=5).grid(row=0, column=3, padx=2)
        ttk.Label(font_frame, text="Title:").grid(row=1, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.boxplot_title_fontsize, width=5).grid(row=1, column=1, padx=2)
        ttk.Label(font_frame, text="Tick:").grid(row=1, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.boxplot_tick_fontsize, width=5).grid(row=1, column=3, padx=2)
        
        # Title wrap width control
        ttk.Label(font_frame, text="Title wrap:").grid(row=2, column=0, padx=2, sticky='w')
        self.boxplot_title_wrap_width = tk.IntVar(value=25)
        ttk.Spinbox(font_frame, from_=20, to=100, textvariable=self.boxplot_title_wrap_width, width=5).grid(row=2, column=1, padx=2)
        ttk.Label(font_frame, text="chars").grid(row=2, column=2, padx=2, sticky='w')

        # Y-axis label customization
        ylabel_frame = ttk.LabelFrame(panel, text="Y-axis Label", padding=5)
        ylabel_frame.pack(fill='x', pady=4)
        self.boxplot_ylabel_text = tk.StringVar(value='Relative Abundance (%)')
        ttk.Label(ylabel_frame, text="Label text:").pack(side='left', padx=(0, 5))
        ttk.Entry(ylabel_frame, textvariable=self.boxplot_ylabel_text, width=30).pack(side='left', expand=True, fill='x')

        # X-axis tick rotation controls
        rotation_frame = ttk.LabelFrame(panel, text="X-axis Tick Rotation", padding=5)
        rotation_frame.pack(fill='x', pady=4)
        self.boxplot_rotate_xticks = tk.BooleanVar(value=True)
        self.boxplot_xtick_rotation = tk.IntVar(value=45)
        ttk.Checkbutton(rotation_frame, text="Rotate X-tick labels", variable=self.boxplot_rotate_xticks).grid(row=0, column=0, sticky='w', padx=2)
        ttk.Label(rotation_frame, text="Rotation angle:").grid(row=0, column=1, padx=8, sticky='w')
        ttk.Spinbox(rotation_frame, from_=0, to=90, textvariable=self.boxplot_xtick_rotation, width=6).grid(row=0, column=2, padx=2)

        # Excel export control (MOVED TO END)
        save_frame = ttk.LabelFrame(panel, text="💾 Save Options", padding=5)
        save_frame.pack(fill='x', pady=4)
        self.boxplot_save_excel = tk.BooleanVar(value=True)
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.boxplot_save_excel).pack(anchor='w')

    def create_bargraph_panel(self, parent):
        """Create Bar graph parameter panel with independent controls."""
        panel = ttk.LabelFrame(parent, text="Bar Graphs", padding=10)
        panel.pack(fill='both', expand=True, padx=5, pady=2)

        ttk.Checkbutton(panel, text="Generate bar graphs",
                        variable=self.viz_selected['bargraph']).pack(anchor='w')

        # Keep layout mode at the top so it is always visible.
        mode_frame = ttk.LabelFrame(panel, text="Layout Mode", padding=5)
        mode_frame.pack(fill='x', pady=4)
        self.bargraph_display_mode = tk.StringVar(value='separate')
        ttk.Radiobutton(mode_frame, text='Separate (one metabolite per plot)', variable=self.bargraph_display_mode, value='separate').pack(anchor='w')
        ttk.Radiobutton(mode_frame, text='Grouped (all metabolites in one grouped plot)', variable=self.bargraph_display_mode, value='grouped').pack(anchor='w')

        grouped_title_frame = ttk.LabelFrame(panel, text="Grouped Plot Title (Grouped mode only)", padding=5)
        grouped_title_frame.pack(fill='x', pady=4)
        self.bargraph_grouped_title = tk.StringVar(value='')
        ttk.Label(grouped_title_frame, text="Title text:").pack(side='left', padx=(0, 5))
        grouped_title_entry = ttk.Entry(grouped_title_frame, textvariable=self.bargraph_grouped_title, width=45)
        grouped_title_entry.pack(side='left', fill='x', expand=True)

        def _toggle_grouped_title_state(*_args):
            state = 'normal' if self.bargraph_display_mode.get().strip().lower() == 'grouped' else 'disabled'
            try:
                grouped_title_entry.configure(state=state)
            except Exception:
                pass

        self.bargraph_display_mode.trace_add('write', _toggle_grouped_title_state)
        _toggle_grouped_title_state()

        boost_frame = ttk.LabelFrame(panel, text="Grouped Low-Value Boost (Grouped mode only)", padding=5)
        boost_frame.pack(fill='x', pady=4)
        self.bargraph_low_boost_enabled = tk.BooleanVar(value=False)
        ttk.Checkbutton(boost_frame, text="Enable low-value boost", variable=self.bargraph_low_boost_enabled).grid(row=0, column=0, columnspan=4, sticky='w', padx=2)
        ttk.Label(boost_frame, text="Boost values below:").grid(row=1, column=0, sticky='w', padx=2)
        self.bargraph_low_boost_threshold = tk.DoubleVar(value=0.25)
        ttk.Spinbox(boost_frame, from_=0.0, to=10.0, increment=0.05, textvariable=self.bargraph_low_boost_threshold, width=8).grid(row=1, column=1, padx=2, sticky='w')
        ttk.Label(boost_frame, text="Scale by:").grid(row=1, column=2, sticky='w', padx=8)
        self.bargraph_low_boost_factor = tk.DoubleVar(value=2.0)
        ttk.Spinbox(boost_frame, from_=1.0, to=20.0, increment=0.5, textvariable=self.bargraph_low_boost_factor, width=8).grid(row=1, column=3, padx=2, sticky='w')
        ttk.Label(boost_frame, text="Applies equal scaling to all bars in each low-value metabolite cluster.",
              font=('Arial', 8), foreground='#666').grid(row=2, column=0, columnspan=4, sticky='w', padx=2, pady=(3, 0))

        comp_metab_frame = ttk.LabelFrame(panel, text="🔍 Group & Annotation Configuration", padding=5)
        comp_metab_frame.pack(fill='x', pady=4)

        left_column = ttk.Frame(comp_metab_frame)
        left_column.pack(side='left', fill='both', expand=True, padx=(0, 5))

        right_column = ttk.Frame(comp_metab_frame)
        right_column.pack(side='left', fill='both', expand=True, padx=(5, 0))

        ttk.Label(left_column, text="Groups to Display", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.bargraph_selected_groups = None
        self.bargraph_groups_status_label = ttk.Label(left_column, text="All groups will be displayed",
                                font=('Arial', 8), foreground='#666')
        self.bargraph_groups_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(left_column, text="Set Groups...",
                   command=self._configure_bargraph_groups).pack(anchor='w')

        ttk.Label(right_column, text="Annotation Comparisons", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.bargraph_annotate_comparisons = None
        self.bargraph_annotate_status_label = ttk.Label(right_column, text="All will be annotated",
                                font=('Arial', 8), foreground='#666')
        self.bargraph_annotate_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(right_column, text="Configure Annotation...",
                   command=self._configure_bargraph_annotations).pack(anchor='w')

        comp_select_frame = ttk.LabelFrame(panel, text="Pairwise Comparison Selection", padding=5)
        comp_select_frame.pack(fill='x', pady=4)
        self.bargraph_selected_comparisons = None
        self.bargraph_comp_status_label = ttk.Label(comp_select_frame, text="All comparisons will be plotted",
                                                    font=('Arial', 8), foreground='#666')
        self.bargraph_comp_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(comp_select_frame, text="Configure Comparisons...",
                   command=self._configure_bargraph_comparisons).pack(anchor='w')

        custom_list_frame = ttk.LabelFrame(panel, text="Custom Metabolite List", padding=5)
        custom_list_frame.pack(fill='x', pady=4)
        self.bargraph_custom_list = tk.StringVar()
        list_row = ttk.Frame(custom_list_frame)
        list_row.pack(fill='x', pady=2)
        ttk.Entry(list_row, textvariable=self.bargraph_custom_list, width=40).pack(side='left', expand=True, fill='x', padx=(0,5))
        ttk.Button(list_row, text="Browse", command=lambda: self.browse_metabolite_list('bargraph')).pack(side='right')
        self.bargraph_use_custom_only = tk.BooleanVar(value=False)
        ttk.Checkbutton(custom_list_frame, text='Use custom list ONLY (ignore p-value/FC filters)',
                        variable=self.bargraph_use_custom_only).pack(anchor='w', pady=(2,0))

        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        ttk.Label(params_frame, text="Top N metabolites:").grid(row=0, column=0, sticky='w', padx=2)
        self.bargraph_top_n = tk.IntVar(value=10)
        ttk.Spinbox(params_frame, from_=1, to=100, textvariable=self.bargraph_top_n, width=8).grid(row=0, column=1, padx=2)
        self.bargraph_no_limit = tk.BooleanVar(value=False)
        ttk.Checkbutton(params_frame, text="No limit (plot all)", variable=self.bargraph_no_limit).grid(row=0, column=2, padx=10)
        self.bargraph_annotate = tk.BooleanVar(value=True)
        ttk.Checkbutton(params_frame, text="Show significance stars", variable=self.bargraph_annotate).grid(row=1, column=0, columnspan=3, sticky='w', pady=5)

        sig_frame = ttk.LabelFrame(panel, text="Significance Filters")
        sig_frame.pack(fill='x', pady=4)
        ttk.Label(sig_frame, text="P-value <").grid(row=0, column=0, sticky='w', padx=2)
        self.bargraph_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(sig_frame, textvariable=self.bargraph_p_thresh, width=8).grid(row=0, column=1, padx=2)
        ttk.Label(sig_frame, text="|Fold Change| ≥").grid(row=0, column=2, sticky='w', padx=8)
        self.bargraph_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(sig_frame, textvariable=self.bargraph_fc_thresh, width=8).grid(row=0, column=3, padx=2)
        self.bargraph_skip_fc = tk.BooleanVar(value=True)
        def _on_bargraph_skip_fc_change():
            try:
                state = 'disabled' if self.bargraph_skip_fc.get() else 'normal'
                for w in sig_frame.grid_slaves(row=0, column=3):
                    try:
                        w.configure(state=state)
                    except Exception:
                        pass
            except Exception:
                pass
        ttk.Checkbutton(sig_frame, text="Skip FC cutoff", variable=self.bargraph_skip_fc, command=_on_bargraph_skip_fc_change).grid(row=0, column=4, padx=4, sticky='w')
        _on_bargraph_skip_fc_change()
        ttk.Label(sig_frame, text="Mode:").grid(row=1, column=0, sticky='w', padx=2)
        self.bargraph_filter_mode = tk.StringVar(value='any')
        ttk.Radiobutton(sig_frame, text='Any', value='any', variable=self.bargraph_filter_mode).grid(row=1, column=1, sticky='w')
        ttk.Radiobutton(sig_frame, text='All', value='all', variable=self.bargraph_filter_mode).grid(row=1, column=2, sticky='w')
        ttk.Radiobutton(sig_frame, text='Specific', value='specific', variable=self.bargraph_filter_mode).grid(row=1, column=3, sticky='w', padx=0)
        ttk.Label(sig_frame, text="Comparison:", font=('Arial', 8)).grid(row=1, column=4, sticky='w', padx=0)
        self.bargraph_specific_comparison = tk.StringVar(value="")
        specific_entry = ttk.Entry(sig_frame, textvariable=self.bargraph_specific_comparison, width=15)
        specific_entry.grid(row=1, column=5, sticky='ew', padx=0)
        self._create_tooltip(specific_entry, "Enter comparison (e.g., 'GroupA|GroupB' or 'AD|Control')\nUsed when Mode='Specific'")

        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        self.bargraph_fig_width = tk.DoubleVar(value=3.0)
        self.bargraph_fig_height = tk.DoubleVar(value=3.0)
        self.bargraph_fig_dpi = tk.IntVar(value=240)
        ttk.Label(size_frame, text="Width:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.bargraph_fig_width, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.bargraph_fig_height, width=6).grid(row=0, column=3, padx=2)
        ttk.Label(size_frame, text="DPI:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.bargraph_fig_dpi, width=6).grid(row=0, column=5, padx=2)

        font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        font_frame.pack(fill='x', pady=4)
        self.bargraph_xlabel_fontsize = tk.IntVar(value=12)
        self.bargraph_ylabel_fontsize = tk.IntVar(value=10)
        self.bargraph_title_fontsize = tk.IntVar(value=14)
        self.bargraph_tick_fontsize = tk.IntVar(value=10)
        self.bargraph_legend_fontsize = tk.IntVar(value=10)

        ttk.Label(font_frame, text="X-label:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.bargraph_xlabel_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(font_frame, text="Y-label:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.bargraph_ylabel_fontsize, width=5).grid(row=0, column=3, padx=2)
        ttk.Label(font_frame, text="Title:").grid(row=1, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.bargraph_title_fontsize, width=5).grid(row=1, column=1, padx=2)
        ttk.Label(font_frame, text="Tick:").grid(row=1, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.bargraph_tick_fontsize, width=5).grid(row=1, column=3, padx=2)
        ttk.Label(font_frame, text="Legend:").grid(row=2, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.bargraph_legend_fontsize, width=5).grid(row=2, column=3, padx=2)

        ttk.Label(font_frame, text="Title wrap:").grid(row=3, column=0, padx=2, sticky='w')
        self.bargraph_title_wrap_width = tk.IntVar(value=25)
        ttk.Spinbox(font_frame, from_=20, to=100, textvariable=self.bargraph_title_wrap_width, width=5).grid(row=3, column=1, padx=2)
        ttk.Label(font_frame, text="chars").grid(row=3, column=2, padx=2, sticky='w')

        ylabel_frame = ttk.LabelFrame(panel, text="Y-axis Label", padding=5)
        ylabel_frame.pack(fill='x', pady=4)
        self.bargraph_ylabel_text = tk.StringVar(value='Relative Abundance (%)')
        ttk.Label(ylabel_frame, text="Label text:").pack(side='left', padx=(0, 5))
        ttk.Entry(ylabel_frame, textvariable=self.bargraph_ylabel_text, width=30).pack(side='left', expand=True, fill='x')

        rotation_frame = ttk.LabelFrame(panel, text="X-axis Tick Rotation", padding=5)
        rotation_frame.pack(fill='x', pady=4)
        self.bargraph_rotate_xticks = tk.BooleanVar(value=True)
        self.bargraph_xtick_rotation = tk.IntVar(value=45)
        ttk.Checkbutton(rotation_frame, text="Rotate X-tick labels", variable=self.bargraph_rotate_xticks).grid(row=0, column=0, sticky='w', padx=2)
        ttk.Label(rotation_frame, text="Rotation angle:").grid(row=0, column=1, padx=8, sticky='w')
        ttk.Spinbox(rotation_frame, from_=0, to=90, textvariable=self.bargraph_xtick_rotation, width=6).grid(row=0, column=2, padx=2)

        save_frame = ttk.LabelFrame(panel, text="💾 Save Options", padding=5)
        save_frame.pack(fill='x', pady=4)
        self.bargraph_save_excel = tk.BooleanVar(value=True)
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.bargraph_save_excel).pack(anchor='w')

        info = (
            "Bar graphs have their own independent controls in this tab.\n"
            "Bars show mean with SD error bars and significance brackets/stars."
        )
        ttk.Label(panel, text=info, foreground='#555', wraplength=650, justify='left').pack(anchor='w', pady=4)
        ttk.Label(panel, text="Tip: Use Separate for one metabolite per figure, or Grouped for all selected metabolites in one figure.",
                  font=('Arial', 8, 'italic'), foreground='#666').pack(anchor='w', pady=(2, 0))

    def create_roc_panel(self, parent):
        """Create ROC parameter panel - REORGANIZED VERSION."""
        panel = ttk.LabelFrame(parent, text="ROC Curves", padding=10)
        panel.pack(fill='both', expand=True, padx=5, pady=2)

        ttk.Checkbutton(panel, text="Generate ROC curves",
                        variable=self.viz_selected['roc']).pack(anchor='w')

        # Comparison and Metabolite List controls side-by-side (MOVED UP for visibility)
        comp_metab_frame = ttk.LabelFrame(panel, text="🔍 Comparison & Metabolite Configuration", padding=5)
        comp_metab_frame.pack(fill='x', pady=4)
        
        # Create two columns
        left_col = ttk.Frame(comp_metab_frame)
        left_col.pack(side='left', fill='both', expand=True, padx=(0,5))
        right_col = ttk.Frame(comp_metab_frame)
        right_col.pack(side='right', fill='both', expand=True, padx=(5,0))
        
        # Left: Comparison Selection
        ttk.Label(left_col, text="Comparison Selection", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.roc_selected_comparisons = None
        self.roc_comp_status_label = ttk.Label(left_col, text="All comparisons will be plotted", 
                                            font=('Arial', 8), foreground='#666')
        self.roc_comp_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(left_col, text="Configure Comparisons...", 
                command=self._configure_roc_comparisons).pack(anchor='w')
        
        # Right: Per-Comparison Metabolite Lists
        ttk.Label(right_col, text="Per-Comparison Lists", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.roc_metabolite_lists = {}
        self.roc_skip_unlisted = tk.BooleanVar(value=False)
        info_label = ttk.Label(right_col, text="Upload different lists for comparisons", 
                            font=('Arial', 8, 'italic'), foreground='#666')
        info_label.pack(anchor='w', pady=(0, 5))
        ttk.Checkbutton(right_col, text="Skip comparisons without lists", 
                    variable=self.roc_skip_unlisted).pack(anchor='w', pady=(0, 5))
        ttk.Button(right_col, text="Manage Metabolite Lists...", 
                command=self._configure_roc_metabolite_lists).pack(anchor='w')
        
        # Custom Metabolite List - Top Right
        custom_list_frame = ttk.LabelFrame(panel, text="Custom Metabolite List", padding=5)
        custom_list_frame.pack(fill='x', pady=4)
        self.roc_custom_list = tk.StringVar()
        list_row = ttk.Frame(custom_list_frame)
        list_row.pack(fill='x', pady=2)
        ttk.Entry(list_row, textvariable=self.roc_custom_list, width=40).pack(side='left', expand=True, fill='x', padx=(0,5))
        ttk.Button(list_row, text="Browse", command=lambda: self.browse_metabolite_list('roc')).pack(side='right')
        self.roc_use_custom_only = tk.BooleanVar(value=False)
        ttk.Checkbutton(custom_list_frame, text='Use custom list ONLY (ignore p-value/FC filters)', 
                    variable=self.roc_use_custom_only).pack(anchor='w', pady=(2,0))

        # Group Selection (MOVED UP for visibility)
        group_frame = ttk.LabelFrame(panel, text="🔍 Group Selection", padding=5)
        group_frame.pack(fill='x', pady=4)
        ttk.Label(group_frame, text="Groups to Display", font=('Arial', 9, 'bold')).pack(anchor='w', pady=(0,5))
        self.roc_selected_groups = None  # Will store list of group names to display
        self.roc_groups_status_label = ttk.Label(group_frame, text="All groups will be displayed", 
                                font=('Arial', 8), foreground='#666')
        self.roc_groups_status_label.pack(anchor='w', pady=(0, 5))
        ttk.Button(group_frame, text="Set Groups...", 
                   command=self._configure_roc_groups).pack(anchor='w')

        params_frame = ttk.Frame(panel)
        params_frame.pack(fill='x', pady=5)
        self.roc_all_pairs = tk.BooleanVar(value=True)
        ttk.Checkbutton(params_frame, text="All pairwise comparisons", variable=self.roc_all_pairs).grid(row=0, column=0, columnspan=2, sticky='w')
        self.roc_include_combined = tk.BooleanVar(value=False)
        ttk.Checkbutton(params_frame, text="Include combined ROC curve", variable=self.roc_include_combined).grid(row=0, column=2, columnspan=2, sticky='w', padx=(15,0))
        ttk.Label(params_frame, text="Max metabolites per ROC:").grid(row=1, column=0, sticky='w', padx=2)
        self.roc_max_metabolites = tk.IntVar(value=100)
        ttk.Spinbox(params_frame, from_=1, to=500, textvariable=self.roc_max_metabolites, width=8).grid(row=1, column=1, padx=2)
        ttk.Label(params_frame, text="Min AUC:").grid(row=2, column=0, sticky='w', padx=2)
        self.roc_min_auc = tk.DoubleVar(value=0.5)
        ttk.Spinbox(params_frame, from_=0.0, to=1.0, increment=0.05, textvariable=self.roc_min_auc, width=8).grid(row=2, column=1, padx=2)

        # Figure size
        size_frame = ttk.LabelFrame(panel, text="Figure Size (inches & DPI)")
        size_frame.pack(fill='x', pady=4)
        self.roc_fig_width = tk.DoubleVar(value=8.0)
        self.roc_fig_height = tk.DoubleVar(value=6.8)
        self.roc_fig_dpi = tk.IntVar(value=260)
        ttk.Label(size_frame, text="Width:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.roc_fig_width, width=6).grid(row=0, column=1, padx=2)
        ttk.Label(size_frame, text="Height:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=2.0, to=20.0, increment=0.5, textvariable=self.roc_fig_height, width=6).grid(row=0, column=3, padx=2)
        ttk.Label(size_frame, text="DPI:").grid(row=0, column=4, padx=8, sticky='w')
        ttk.Spinbox(size_frame, from_=72, to=600, textvariable=self.roc_fig_dpi, width=6).grid(row=0, column=5, padx=2)

        # Significance filters
        sig_frame = ttk.LabelFrame(panel, text="Significance Filters")
        sig_frame.pack(fill='x', pady=4)
        ttk.Label(sig_frame, text="P-value <").grid(row=0, column=0, sticky='w', padx=2)
        self.roc_p_thresh = tk.DoubleVar(value=0.05)
        ttk.Entry(sig_frame, textvariable=self.roc_p_thresh, width=8).grid(row=0, column=1, padx=2)
        ttk.Label(sig_frame, text="|Fold Change| ≥").grid(row=0, column=2, sticky='w', padx=8)
        self.roc_fc_thresh = tk.DoubleVar(value=2.0)
        ttk.Entry(sig_frame, textvariable=self.roc_fc_thresh, width=8).grid(row=0, column=3, padx=2)
        self.roc_skip_fc = tk.BooleanVar(value=True)
        def _on_roc_skip_fc_change():
            try:
                state = 'disabled' if self.roc_skip_fc.get() else 'normal'
                for w in sig_frame.grid_slaves(row=0, column=3):
                    try:
                        w.configure(state=state)
                    except Exception:
                        pass
            except Exception:
                pass
        ttk.Checkbutton(sig_frame, text="Skip FC cutoff", variable=self.roc_skip_fc, command=_on_roc_skip_fc_change).grid(row=0, column=4, padx=4, sticky='w')
        _on_roc_skip_fc_change()
        ttk.Label(sig_frame, text="Mode:").grid(row=1, column=0, sticky='w', padx=2)
        self.roc_filter_mode = tk.StringVar(value='any')
        ttk.Radiobutton(sig_frame, text='Any', value='any', variable=self.roc_filter_mode).grid(row=1, column=1, sticky='w')
        ttk.Radiobutton(sig_frame, text='All', value='all', variable=self.roc_filter_mode).grid(row=1, column=2, sticky='w')
        ttk.Radiobutton(sig_frame, text='Specific', value='specific', variable=self.roc_filter_mode).grid(row=1, column=3, sticky='w', padx=0)
        ttk.Label(sig_frame, text="Comparison:", font=('Arial', 8)).grid(row=1, column=4, sticky='w', padx=0)
        self.roc_specific_comparison = tk.StringVar(value="")
        specific_entry = ttk.Entry(sig_frame, textvariable=self.roc_specific_comparison, width=15)
        specific_entry.grid(row=1, column=5, sticky='ew', padx=0)
        self._create_tooltip(specific_entry, "Enter comparison (e.g., 'GroupA|GroupB' or 'AD|Control')\nUsed when Mode='Specific'")

        # Font size controls
        font_frame = ttk.LabelFrame(panel, text="Font Sizes", padding=5)
        font_frame.pack(fill='x', pady=4)
        self.roc_xlabel_fontsize = tk.IntVar(value=12)
        self.roc_ylabel_fontsize = tk.IntVar(value=12)
        self.roc_title_fontsize = tk.IntVar(value=14)
        self.roc_tick_fontsize = tk.IntVar(value=10)
        self.roc_legend_fontsize = tk.IntVar(value=10)
        
        ttk.Label(font_frame, text="X-label:").grid(row=0, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.roc_xlabel_fontsize, width=5).grid(row=0, column=1, padx=2)
        ttk.Label(font_frame, text="Y-label:").grid(row=0, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.roc_ylabel_fontsize, width=5).grid(row=0, column=3, padx=2)
        ttk.Label(font_frame, text="Title:").grid(row=1, column=0, padx=2, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.roc_title_fontsize, width=5).grid(row=1, column=1, padx=2)
        ttk.Label(font_frame, text="Tick:").grid(row=1, column=2, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.roc_tick_fontsize, width=5).grid(row=1, column=3, padx=2)
        ttk.Label(font_frame, text="Legend:").grid(row=1, column=4, padx=8, sticky='w')
        ttk.Spinbox(font_frame, from_=4, to=24, textvariable=self.roc_legend_fontsize, width=5).grid(row=1, column=5, padx=2)

        # Excel export control (MOVED TO END)
        save_frame = ttk.LabelFrame(panel, text="💾 Save Options", padding=5)
        save_frame.pack(fill='x', pady=4)
        self.roc_save_excel = tk.BooleanVar(value=True)
        ttk.Checkbutton(save_frame, text="Save Excel Files (CSV)", variable=self.roc_save_excel).pack(anchor='w')
        self.roc_excel_only = tk.BooleanVar(value=False)
        ttk.Checkbutton(save_frame, text="Excel only (AUC values, no plots) - for testing", 
                        variable=self.roc_excel_only).pack(anchor='w')

    def choose_viz_output_dir(self):
        """Choose output directory for visualizations."""
        directory = filedialog.askdirectory(title="Select Visualization Output Directory")
        if directory:
            self.viz_output_dir.set(directory)

    def open_viz_output_folder(self):
        """Open the visualization output folder in file explorer."""
        output_dir = self.viz_output_dir.get()
        if os.path.exists(output_dir):
            os.startfile(output_dir)
        else:
            messagebox.showwarning("Directory Not Found", f"Output directory does not exist:\n{output_dir}")

    def _confirm_venn_filter_alignment(self) -> bool:
        """Require explicit user confirmation that Venn filters match statistics choices."""
        pval = self.venn_p_thresh.get() if hasattr(self, 'venn_p_thresh') else 0.05
        skip_fc = bool(hasattr(self, 'venn_skip_fc') and self.venn_skip_fc.get())
        fc = self.venn_fc_thresh.get() if hasattr(self, 'venn_fc_thresh') else 2.0

        allmol_enabled = bool(hasattr(self, 'venn_generate_all_molecules') and self.venn_generate_all_molecules.get())
        presence_type = self.venn_allmol_min_presence_type.get() if hasattr(self, 'venn_allmol_min_presence_type') else 'count'
        presence_count = self.venn_allmol_min_presence_count.get() if hasattr(self, 'venn_allmol_min_presence_count') else 3
        presence_percent = self.venn_allmol_min_presence_percent.get() if hasattr(self, 'venn_allmol_min_presence_percent') else 50.0

        filtered_line = f"Filtered Venn: p < {pval}, FC cutoff skipped" if skip_fc else f"Filtered Venn: p < {pval}, FC >= {fc}"

        if allmol_enabled:
            if str(presence_type).lower() == 'percentage':
                allmol_line = f"All Molecules presence: Percentage >= {presence_percent}% non-zero/non-NaN per group"
            else:
                allmol_line = f"All Molecules presence: Count >= {presence_count} non-zero/non-NaN per group"
        else:
            allmol_line = "All Molecules Venn: disabled"

        msg = (
            "Venn Filter Confirmation Required\n\n"
            "Count/Percentage thresholds are critical and must match your statistics filtering choices.\n"
            "Please confirm these settings are EXACTLY what you intend to use before proceeding.\n\n"
            f"{filtered_line}\n"
            f"{allmol_line}\n\n"
            "Proceed with Venn generation?"
        )
        return bool(messagebox.askyesno("Confirm Venn Filtering", msg, icon='warning'))

    def generate_selected_viz_plots(self):
        """Generate only the selected visualization plots."""
        selected_plots = [plot_type for plot_type, var in self.viz_selected.items() if var.get()]
        
        if not selected_plots:
            messagebox.showwarning("No Plots Selected", "Please select at least one plot type to generate.")
            return

        if 'venn' in selected_plots and not self._confirm_venn_filter_alignment():
            self.log_viz_message("⚠️ Venn generation cancelled: filter confirmation not accepted.")
            return
        
        self.log_viz_message(f"Starting generation of selected plots: {', '.join(selected_plots)}")
        # TODO: Implement plot generation logic
        threading.Thread(target=self._generate_viz_plots, args=(selected_plots,), daemon=True).start()

    def generate_all_viz_plots(self):
        """Generate all available visualization plots."""
        all_plots = list(self.viz_selected.keys())

        if 'venn' in all_plots and not self._confirm_venn_filter_alignment():
            self.log_viz_message("⚠️ Venn generation cancelled: filter confirmation not accepted.")
            return

        self.log_viz_message(f"Starting generation of all plots: {', '.join(all_plots)}")
        # TODO: Implement plot generation logic
        threading.Thread(target=self._generate_viz_plots, args=(all_plots,), daemon=True).start()

    def stop_viz_plots(self):
        """Stop the current visualization plot generation."""
        self.viz_cancel_flag.set()
        self.log_viz_message("Plot generation stopped by user")
        self.viz_progress_label.config(text="Stopped")

    def _normalize_verified_assignments(self, assignments_dict):
        """Convert verified assignments from UI format to analysis format.
        
        Converts from:
            {'Log2 Fold Change': 'abd', 'P-Value': 'cdf', '-Log10 P-Value': 'gft', ...}
        
        To:
            {'log2FC': 'abd', 'pvalue': 'cdf', 'neglog10pvalue': 'gft', ...}
        """
        if not assignments_dict:
            logger.info(f"🔍 _normalize_verified_assignments: assignments_dict is None/empty")
            return None
        
        # Mapping from UI column type names to normalized keys used by analysis code
        key_mapping = {
            'Log2 Fold Change': 'log2FC',
            'P-Value': 'pvalue',
            '-Log10 P-Value': 'neglog10pvalue',
            'Fold Change': 'FC',
            'Feature ID': 'feature_id',
            'HMDB ID': 'hmdb_id',
            'KEGG ID': 'kegg_id',
            'ChEBI ID': 'chebi_id',
            'InChIKey': 'inchikey',
            'Super Class': 'super_class',
            'Class': 'class',
        }
        
        normalized = {}
        for ui_key, col_name in assignments_dict.items():
            if col_name:  # Skip None values
                normalized_key = key_mapping.get(ui_key, ui_key.lower().replace(' ', '_'))
                normalized[normalized_key] = col_name
        
        if normalized:
            logger.info(f"🔍 _normalize_verified_assignments: Converted {len(assignments_dict)} assignments to {len(normalized)} normalized: {normalized}")
        
        return normalized if normalized else None

    def _generate_viz_plots(self, plot_types):
        """Internal method to generate visualization plots."""
        try:
            self.viz_cancel_flag.clear()

            # **IMPORTANT**: Check if column verification is needed
            # For pairwise data (with _vs_ columns), automated detection works
            # For global data (without _vs_ columns), verification is required
            complete_df_check = self.get_current_complete_df()
            if complete_df_check is not None:
                pairwise_cols = [c for c in complete_df_check.columns if '_vs_' in c]
                is_pairwise = len(pairwise_cols) > 0
                
                # Column verification removed - automated detection handles all cases
                if is_pairwise:
                    self.log_viz_message(f"✅ Using automated pairwise column detection (detected {len(pairwise_cols)} _vs_ columns)")
            else:
                self.log_viz_message("⚠️ Warning: Could not check data type for column verification")

            # Verify group configuration
            if not getattr(self, 'viz_group_mapping', None):
                self.log_viz_message("❌ Error: No group configuration found. Please configure groups first.")
                return

            # Acquire data
            complete_df = self.get_current_complete_df()
            logger.info(f"📊 Retrieved complete_df: {len(complete_df) if complete_df is not None else 0} rows")
            if complete_df is None or complete_df.empty:
                self.log_viz_message("❌ Error: No data available. Please load or compute statistics first.")
                return

            # Determine groups / samples
            groups = self.get_viz_groups()
            sample_to_group = self.get_viz_sample_to_group_mapping()
            sample_cols = list(sample_to_group.keys())

            min_required = self._get_min_group_size()
            filtered_map, viz_counts, viz_excluded = self._filter_groups_by_min_samples(sample_to_group, min_required)
            if viz_excluded:
                excluded_info = ', '.join(f"{grp} (n={cnt})" for grp, cnt in viz_excluded.items())
                self.log_viz_message(f"⚠️ Excluding groups below minimum sample threshold ({min_required}): {excluded_info}")
            sample_to_group = filtered_map
            if not sample_to_group:
                self.log_viz_message(f"❌ No groups with at least {min_required} samples available for visualization. Aborting plots.")
                return
            sample_cols = list(sample_to_group.keys())
            groups = self.ordered_groups(list(dict.fromkeys(sample_to_group.values())))
            included_counts = {grp: viz_counts.get(grp, 0) for grp in groups}
            if len(groups) < 2:
                self.log_viz_message(f"⚠️ Only one group meets the minimum threshold ({min_required}). Pairwise-dependent plots may be skipped or empty.")

            # Output directory
            output_dir = self.viz_output_dir.get()
            os.makedirs(output_dir, exist_ok=True)

            # Normalize verified assignments from UI format to analysis format
            normalized_verified = self._normalize_verified_assignments(self.verified_assignments)
            if self.verified_assignments and normalized_verified:
                self.log_viz_message(f"✅ Using verified column assignments: {normalized_verified}")
                logger.info(f"🔍 DEBUG: normalized_verified = {normalized_verified}")
            else:
                logger.info(f"🔍 DEBUG: self.verified_assignments = {self.verified_assignments}, normalized_verified = {normalized_verified}")

            # Ensure parameter objects exist & refresh from GUI controls if newer system present
            if hasattr(self, '_update_viz_params') and hasattr(self, 'viz_params'):
                try:
                    self._update_viz_params()
                except Exception as e:
                    self.log_viz_message(f"Parameter update warning: {e}")

            self.log_viz_message(f"Starting generation of {len(plot_types)} plot types")
            self.log_viz_message(f"Output directory: {output_dir}")
            group_report = ', '.join(f"{grp}:{included_counts.get(grp, 0)}" for grp in groups)
            self.log_viz_message(f"Groups: {group_report} ({len(sample_cols)} samples)")

            total = max(1, len(plot_types))
            for idx, ptype in enumerate(plot_types):
                if self.viz_cancel_flag.is_set():
                    break

                # Progress update
                pct = int((idx / total) * 100)
                self.viz_progress['value'] = pct
                self.viz_progress_label.config(text=f"Generating {ptype}...")
                self.log_viz_message(f"Generating {ptype} plots...")

                try:
                    if ptype == 'volcano':
                        self._generate_volcano_plots(complete_df, groups, sample_cols, sample_to_group, output_dir, normalized_verified)
                    elif ptype == 'pca':
                        self._generate_pca_plots(complete_df, groups, sample_cols, sample_to_group, output_dir)
                        self.log_viz_message(f"✅ {ptype} completed")
                    elif ptype == 'venn':
                        self._generate_venn_plots(complete_df, groups, sample_cols, sample_to_group, output_dir)
                        self.log_viz_message(f"✅ {ptype} completed")
                    elif ptype == 'boxplot':
                        self._generate_boxplot_plots(complete_df, groups, sample_cols, sample_to_group, output_dir, normalized_verified)
                        self.log_viz_message(f"✅ {ptype} completed")
                    elif ptype == 'bargraph':
                        self._generate_bargraph_plots(complete_df, groups, sample_cols, sample_to_group, output_dir, normalized_verified)
                        self.log_viz_message(f"✅ {ptype} completed")
                    elif ptype == 'heatmap':
                        self._generate_heatmap_plots(complete_df, groups, sample_cols, sample_to_group, output_dir, normalized_verified)
                        self.log_viz_message(f"✅ {ptype} completed")
                    elif ptype == 'roc':
                        self._generate_roc_plots(complete_df, groups, sample_cols, sample_to_group, output_dir)
                        self.log_viz_message(f"✅ {ptype} completed")
                    else:
                        self.log_viz_message(f"⚠️ Unknown plot type: {ptype}")
                        continue
                except Exception as e:
                    import traceback
                    self.log_viz_message(f"❌ {ptype} failed: {e}")
                    self.log_viz_message(traceback.format_exc())

            if not self.viz_cancel_flag.is_set():
                self.viz_progress['value'] = 100
                self.viz_progress_label.config(text="All plots completed!")
                self.log_viz_message("🎉 All selected plots generated successfully!")
                self.log_viz_message(f"📁 Check output folder: {output_dir}")
        except Exception as e:
            import traceback
            self.log_viz_message(f"❌ Visualization generation crashed: {e}")
            self.log_viz_message(traceback.format_exc())

    # --- Updated helper generators using unified service API --- #
    def _build_common_context(self, complete_df, groups, sample_cols, sample_to_group, outdir, verified_assignments=None):
        from main_script.metabolites_visualization import CommonVizContext
        
        # DEBUG: Log incoming DataFrame size
        logger.info(f"🔍 _build_common_context ENTRY: complete_df has {len(complete_df)} rows, outdir={outdir}")
        
        # Pre-filter complete_df for custom metabolite lists to ensure combine mode respects custom lists
        plot_type = None
        if 'heatmaps' in outdir:
            plot_type = 'heatmap'
        elif 'boxplots' in outdir:
            plot_type = 'boxplot'
        elif 'bargraphs' in outdir:
            plot_type = 'bargraph'
        elif 'roc' in outdir:
            plot_type = 'roc'
        elif 'venn' in outdir.lower():
            plot_type = 'venn'
        
        logger.info(f"🔍 _build_common_context: Detected plot_type = '{plot_type}'")
        
        if plot_type and plot_type != 'venn' and hasattr(self, 'viz_params') and plot_type in self.viz_params:
            params = self.viz_params[plot_type]
            if (hasattr(params, 'use_custom_only') and params.use_custom_only and 
                hasattr(params, 'include_metabolites') and params.include_metabolites):
                # Filter complete_df to only include custom metabolites - use multi-column matching
                # Import the matching function
                from main_script.metabolites_visualization import match_metabolites_multi_column

                # Improved ID column selection
                id_col = None
                for candidate in ['Gene', 'Protein', 'Name']:
                    if candidate in complete_df.columns:
                        id_col = candidate
                        break
                if id_col is None:
                    # Find first non-numeric column (feature column)
                    for col in complete_df.columns:
                        if not pd.api.types.is_numeric_dtype(complete_df[col]):
                            id_col = col
                            break
                if id_col is None and 'metabolite_id' in complete_df.columns:
                    id_col = 'metabolite_id'
                if id_col is None:
                    id_col = complete_df.columns[0]

                logger.info(f"🔍 Pre-filtering complete_df with custom list (use_custom_only mode)...")
                logger.info(f"   Before filter: {len(complete_df)} metabolites")

                # Use multi-column matching instead of simple .isin()
                match_mask = match_metabolites_multi_column(complete_df, params.include_metabolites, id_col)
                complete_df = complete_df[match_mask].copy()

                logger.info(f"   After filter: {len(complete_df)} metabolites")

        # Ensure visualization respects minimum group sample threshold even for direct calls
        min_required = self._get_min_group_size()
        filtered_map, _, _ = self._filter_groups_by_min_samples(sample_to_group, min_required)
        if filtered_map:
            sample_to_group = filtered_map
            sample_cols = [col for col in sample_cols if col in sample_to_group]
            filtered_groups = list(dict.fromkeys(sample_to_group.values()))
            if filtered_groups:
                groups = self.ordered_groups(filtered_groups)
        
        # If groups is still empty, use configured group definitions (required for accuracy)
        if not groups:
            # Use configured group definitions - this is always configured in Statistics tab
            if hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
                groups = list(self.viz_group_definitions.values())
                logger.info(f"✅ Using configured group definitions: {groups}")
        
        # Debug: Log all available pairwise comparisons and their stat columns
        pairwise_cols = [col for col in complete_df.columns if '_vs_' in col]
        print(f"DEBUG: Total pairwise columns: {len(pairwise_cols)}")
        if pairwise_cols:
            print(f"DEBUG: Sample pairwise columns: {pairwise_cols[:10]}")
            logger.info(f"🔍 Total pairwise columns in data: {len(pairwise_cols)}")
            # Group columns by comparison prefix
            comp_stats = {}
            for col in pairwise_cols:
                # Find the comparison prefix (before stat suffix)
                if '_log2FC' in col:
                    prefix = col.replace('_log2FC', '')
                elif '_p_adj' in col:
                    prefix = col.replace('_p_adj', '')
                elif '_adj_p' in col:
                    prefix = col.replace('_adj_p', '')
                elif '_FC' in col and not '_log2FC' in col:
                    prefix = col.replace('_FC', '')
                elif '_neg_log10_' in col:
                    prefix = col.replace('_neg_log10_p_adj', '').replace('_neg_log10_adj_p', '')
                else:
                    continue  # Skip unknown stat columns
                if prefix not in comp_stats:
                    comp_stats[prefix] = []
                comp_stats[prefix].append(col)
            # Log each comparison and its available stat columns
            for comp in sorted(comp_stats.keys()):
                print(f"DEBUG: {comp}: {sorted(comp_stats[comp])}")
                logger.info(f"  🔍 {comp}: {sorted(comp_stats[comp])}")
        else:
            print("DEBUG: No pairwise columns found")
            logger.info("🔍 No pairwise columns found in imported data")
        
        # Build color map with fallback
        try:
            color_map = self.build_color_map()
        except Exception:
            # Fallback simple palette
            palette = {}
            try:
                import seaborn as sns
                base = sns.color_palette('tab10')
                palette = {g: base[i % len(base)] for i, g in enumerate(groups)}
            except Exception:
                palette = {g: '#%02x%02x%02x' % (int(255/(i+1)), 50, 150) for i, g in enumerate(groups)}
            color_map = palette
        # If user prefers raw p-values (unchecked 'Use adjusted p-values'),
        # try to enrich the imported Complete Results with any raw p-value
        # columns that may exist on other sheets in the same Excel file.
        try:
            prefer_adj = getattr(self, 'use_adj_p_var', tk.BooleanVar(value=True)).get()
        except Exception:
            prefer_adj = True
        
        # AUTO-DETECT: Check what p-value columns actually exist in the data
        # This helps when the user changed the checkbox but the data was created with different settings
        has_adj_p_cols = any(col.lower().endswith('_adj_p') or col.lower().endswith('_adj_p_value') for col in complete_df.columns)
        has_raw_p_cols = any(col.lower().endswith('_pvalue') or col.lower().endswith('_p_value') or col.lower().endswith('_p') 
                             for col in complete_df.columns if not (col.lower().endswith('_adj_p') or col.lower().endswith('_adj_p_value')))
        
        # If only one type exists, prefer that type regardless of checkbox setting
        if has_raw_p_cols and not has_adj_p_cols:
            prefer_adj = False
            logger.info(f"🔍 Auto-detected raw p-values in data (no adjusted p-values found); preferring raw p-values")
        elif has_adj_p_cols and not has_raw_p_cols:
            prefer_adj = True
            logger.info(f"🔍 Auto-detected adjusted p-values in data (no raw p-values found); preferring adjusted p-values")

        # Best-effort: if the user asked to prefer raw p-values and we do NOT
    # already have raw p-value columns in Complete Results, scan other sheets
    # and copy matching pairwise raw p-value columns into the complete_df
    # keyed by metabolite identifier ('metabolite_id' or 'Name').
    # Note: avoid importing derived metrics like *_neg_log10_p and require
    # a comparison prefix pattern like <G1>_vs_<G2>_pvalue.
        if not prefer_adj and not has_raw_p_cols:
            try:
                import re
                # Determine Excel path: prefer the explicit import file control,
                # else try the previously imported file attribute.
                excel_path = None
                if hasattr(self, 'viz_import_file') and isinstance(self.viz_import_file, tk.StringVar):
                    excel_path = self.viz_import_file.get()
                if not excel_path and hasattr(self, 'imported_complete_df_path'):
                    excel_path = getattr(self, 'imported_complete_df_path')
                if not excel_path and hasattr(self, 'viz_import_file') and isinstance(self.viz_import_file, str):
                    excel_path = self.viz_import_file
                if excel_path and os.path.exists(excel_path):
                    try:
                        xl = pd.ExcelFile(excel_path)
                        id_col = 'metabolite_id' if 'metabolite_id' in complete_df.columns else ('Name' if 'Name' in complete_df.columns else None)
                        if id_col is not None:
                            # Keep track of new columns added
                            added = []
                            for sheet in xl.sheet_names:
                                # Skip the Complete Results sheet to avoid redundant work
                                if sheet.lower().strip() == 'complete results':
                                    continue
                                try:
                                    df_sheet = xl.parse(sheet_name=sheet)
                                except Exception:
                                    continue
                                if id_col not in df_sheet.columns:
                                    continue
                                # Find candidate raw p-value columns present on this sheet
                                for col in df_sheet.columns:
                                    col_lower = str(col).lower()
                                    if col in complete_df.columns:
                                        continue
                                    # Only accept pairwise p-value columns with a comparison prefix
                                    # Examples: PC3_2D_vs_DU145_3D_pvalue, PC3_2D_vs_DU145_3D_p_value
                                    # Exclude derived metrics like *_neg_log10_p
                                    if 'neg_log10' in col_lower:
                                        continue
                                    if re.search(r'.+_vs_.+_(p_value|pvalue|p)$', col_lower):
                                        # Map by id_col
                                        try:
                                            mapping = df_sheet.set_index(id_col)[col]
                                            # Align into complete_df
                                            complete_df = complete_df.set_index(id_col)
                                            complete_df[col] = complete_df.index.map(mapping)
                                            complete_df = complete_df.reset_index()
                                            added.append(col)
                                        except Exception:
                                            # If mapping fails, skip this column
                                            complete_df = complete_df.reset_index() if id_col in complete_df.index.names else complete_df
                                            continue
                            if added:
                                self.log_viz_message(f"Imported raw p-value columns from other sheets: {', '.join(added)}")
                    except Exception as e:
                        self.log_viz_message(f"Warning: failed to scan Excel for raw p-value columns: {e}")
            except Exception:
                pass

        # Lipid mode detection - check user's mode selection
        is_lipid = False
        try:
            # User override from viz_mode radio buttons
            user_mode = None
            if hasattr(self, 'viz_mode'):
                user_mode = self.viz_mode.get()
            if user_mode == 'lipid':
                is_lipid = True
            elif user_mode in ['metabolite', 'custom']:
                is_lipid = False
        except Exception:
            pass

        # DEBUG: Log outgoing DataFrame size
        logger.info(f"🔍 _build_common_context EXIT: complete_df has {len(complete_df)} rows (outdir={outdir})")

        # Determine ID column from multiple sources
        id_col = None
        # Priority 1: stat_column_assignments (from saved config)
        if getattr(self, 'stat_column_assignments', None):
            id_col = self.stat_column_assignments.get('id_column')
        # Priority 2: memory_store (from Statistics tab)
        if not id_col and hasattr(self, 'memory_store') and self.memory_store:
            # Check for lipid class ID column first (more specific)
            if 'id_column_class' in self.memory_store:
                id_col = self.memory_store['id_column_class']
            elif 'id_column' in self.memory_store:
                id_col = self.memory_store['id_column']
        # Priority 3: Auto-detect from dataframe columns (for lipid class data)
        if not id_col and 'Class' in complete_df.columns:
            # If 'Class' column exists and we're likely using lipid class data, use it as ID
            # Check if this looks like lipid class data (has Class column but not LipidID)
            if 'LipidID' not in complete_df.columns or complete_df['LipidID'].isna().all():
                id_col = 'Class'
        
        return CommonVizContext(
            complete_df=complete_df,
            groups=groups,
            sample_cols=sample_cols,
            sample_to_group=sample_to_group,
            output_dir=outdir,
            color_map=color_map,
            # Propagate whether to prefer adjusted p-values (default True)
            use_adj_p=prefer_adj,
            preferred_group_order=[g.strip() for g in self.viz_preferred_group_order.get().split(',') if g.strip()] if hasattr(self, 'viz_preferred_group_order') and self.viz_preferred_group_order.get().strip() else None,
            is_lipid_mode=is_lipid,
            verified_assignments=verified_assignments,
            stat_column_assignments=getattr(self, 'stat_column_assignments', None),
            id_column=id_col
        )

    def _generate_pca_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir):
        from main_script.metabolites_visualization import run_pca_analysis
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'pca'))
        os.makedirs(ctx.output_dir, exist_ok=True)
        params = None
        if hasattr(self, 'viz_params') and 'pca' in self.viz_params:
            params = self.viz_params['pca']
        else:
            from main_script.metabolites_visualization import PCAParams
            params = PCAParams()
        result = run_pca_analysis(ctx, params)
        if result.files_created:
            self.log_viz_message(f"PCA: {len(result.files_created)} files saved")
            # If lipid feature metadata columns were filtered out, log them
            if hasattr(ctx, '_lipid_feature_columns_removed') and ctx._lipid_feature_columns_removed:
                removed_list = ', '.join(ctx._lipid_feature_columns_removed[:10])
                extra = '' if len(ctx._lipid_feature_columns_removed) <= 10 else f" (+{len(ctx._lipid_feature_columns_removed)-10} more)"
                self.log_viz_message(f"PCA (Lipid mode): excluded metadata columns: {removed_list}{extra}")
            # Additional clarity: report lipid class outputs directory if present
            try:
                if getattr(params, 'include_lipid_class', False):
                    class_subdir = getattr(params, 'class_subdir_name', 'pca_class')
                    class_dir_path = os.path.join(ctx.output_dir, class_subdir)
                    if os.path.isdir(class_dir_path):
                        class_files = [f for f in os.listdir(class_dir_path) if f.lower().startswith('pca_lipid_class') or f.lower().endswith('.png') or f.lower().endswith('.csv')]
                        if class_files:
                            self.log_viz_message(f"PCA (Class): {len(class_files)} class-level files in {class_dir_path}")
            except Exception:
                pass
        for err in result.errors:
            self.log_viz_message(f"PCA error: {err}")

    def _generate_volcano_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir, verified_assignments=None):
        from main_script.metabolites_visualization import run_volcano_analysis
        
        # DEBUG: Log the DataFrame size BEFORE building context
        self.log_viz_message(f"🔍 DEBUG Volcano: complete_df has {len(complete_df)} rows BEFORE _build_common_context")
        
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'volcano'), verified_assignments)
        
        # DEBUG: Log the DataFrame size AFTER building context
        self.log_viz_message(f"🔍 DEBUG Volcano: ctx.complete_df has {len(ctx.complete_df)} rows AFTER _build_common_context")
        
        os.makedirs(ctx.output_dir, exist_ok=True)
        if hasattr(self, 'viz_params') and 'volcano' in self.viz_params:
            params = self.viz_params['volcano']
        else:
            from main_script.metabolites_visualization import VolcanoParams
            # Support legacy variable names
            pval = getattr(self, 'volcano_p_thresh', getattr(self, 'volcano_p_threshold', tk.DoubleVar(value=0.05))).get() if hasattr(self, 'volcano_p_thresh') or hasattr(self, 'volcano_p_threshold') else 0.05
            fc = getattr(self, 'volcano_fc_thresh', getattr(self, 'volcano_fc_threshold', tk.DoubleVar(value=2.0))).get() if hasattr(self, 'volcano_fc_thresh') or hasattr(self, 'volcano_fc_threshold') else 2.0
            params = VolcanoParams(p_threshold=pval, fc_threshold=fc)
        # Debug: report which p-value & fold-change thresholds (and figure sizing) will be used
        try:
            self.log_viz_message(f"DEBUG: Volcano params → p_threshold={getattr(params, 'p_threshold', None)}, fc_threshold={getattr(params, 'fc_threshold', None)}, fig_width={getattr(params, 'fig_width', getattr(self, 'volcano_fig_width', tk.DoubleVar(value=6.0)).get())}, fig_height={getattr(params, 'fig_height', getattr(self, 'volcano_fig_height', tk.DoubleVar(value=5.0)).get())}, dpi={getattr(params, 'fig_dpi', getattr(self, 'volcano_fig_dpi', tk.IntVar(value=250)).get())}")
            # Also print to stdout for immediate console traceability
            print(f"DEBUG: Volcano params -> p_threshold={getattr(params, 'p_threshold', None)}, fc_threshold={getattr(params, 'fc_threshold', None)}")
        except Exception:
            pass
        result = run_volcano_analysis(ctx, params)
        if getattr(result, 'summary', None):
            self.log_viz_message(result.summary)
        if result.files_created:
            self.log_viz_message(f"Volcano: {len(result.files_created)} files saved")
        for err in result.errors:
            self.log_viz_message(f"Volcano error: {err}")
        if result.files_created:
            self.log_viz_message("✅ volcano completed")
        elif result.errors:
            self.log_viz_message("❌ volcano failed")
        else:
            self.log_viz_message("⚠️ volcano completed (no plots generated)")

    def _generate_boxplot_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir, verified_assignments=None):
        from main_script.metabolites_visualization import run_boxplot_analysis
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'boxplots'), verified_assignments)
        os.makedirs(ctx.output_dir, exist_ok=True)
        if hasattr(self, 'viz_params') and 'boxplot' in self.viz_params:
            params = self.viz_params['boxplot']
        else:
            from main_script.metabolites_visualization import BoxplotParams
            params = BoxplotParams()
        # Propagate x-tick rotation preferences to the dataframe for plotting
        try:
            rotate_on = getattr(params, 'rotate_xticks', True)
            rotation = getattr(params, 'xtick_rotation', 45)
            setattr(ctx.complete_df, '_viz_rotate_xticks', bool(rotate_on))
            setattr(ctx.complete_df, '_viz_xtick_rotation', int(rotation))
        except Exception:
            pass
        result = run_boxplot_analysis(ctx, params)
        if result.files_created:
            self.log_viz_message(f"Boxplots: {len(result.files_created)} files saved")
        for err in result.errors:
            self.log_viz_message(f"Boxplot error: {err}")
        if result.files_created and not result.errors:
            self.log_viz_message("✅ boxplot completed")
        elif result.files_created and result.errors:
            self.log_viz_message("⚠️ boxplot completed with errors")
        elif result.errors:
            self.log_viz_message("❌ boxplot failed")
        else:
            self.log_viz_message("⚠️ boxplot completed (no plots generated)")

    def _generate_bargraph_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir, verified_assignments=None):
        from main_script.metabolites_visualization import run_bargraph_analysis
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'bargraphs'), verified_assignments)
        os.makedirs(ctx.output_dir, exist_ok=True)
        if hasattr(self, 'viz_params') and 'bargraph' in self.viz_params:
            params = self.viz_params['bargraph']
        else:
            from main_script.metabolites_visualization import BargraphParams
            params = BargraphParams()
        try:
            rotate_on = getattr(params, 'rotate_xticks', True)
            rotation = getattr(params, 'xtick_rotation', 45)
            setattr(ctx.complete_df, '_viz_rotate_xticks', bool(rotate_on))
            setattr(ctx.complete_df, '_viz_xtick_rotation', int(rotation))
        except Exception:
            pass
        result = run_bargraph_analysis(ctx, params)
        if result.files_created:
            self.log_viz_message(f"Bar Graphs: {len(result.files_created)} files saved")
        for err in result.errors:
            self.log_viz_message(f"Bar graph error: {err}")
        if result.files_created and not result.errors:
            self.log_viz_message("✅ bargraph completed")
        elif result.files_created and result.errors:
            self.log_viz_message("⚠️ bargraph completed with errors")
        elif result.errors:
            self.log_viz_message("❌ bargraph failed")
        else:
            self.log_viz_message("⚠️ bargraph completed (no plots generated)")

    def _generate_heatmap_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir, verified_assignments=None):
        from main_script.metabolites_visualization import run_heatmap_analysis
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'heatmaps'), verified_assignments)
        os.makedirs(ctx.output_dir, exist_ok=True)
        if hasattr(self, 'viz_params') and 'heatmap' in self.viz_params:
            params = self.viz_params['heatmap']
        else:
            from main_script.metabolites_visualization import HeatmapParams
            params = HeatmapParams()
        
        # Check if auto_size is enabled - only override dimensions if True
        auto_size_enabled = getattr(params, 'auto_size', True)
        
        # Log for debugging
        self.log_viz_message(f"Heatmap auto_size: {auto_size_enabled}, Width: {getattr(params, 'fig_width', 'N/A')}, Height: {getattr(params, 'fig_height', 'N/A')}")
        
        # Calculate dynamic figure size ONLY if auto_size is enabled
        if auto_size_enabled:
            # Calculate number of samples
            num_samples = len(sample_cols)
            
            # Estimate number of metabolites to plot
            # This is approximate since the actual filtering happens in the visualization module
            max_metabolites = getattr(params, 'max_metabolites', 50)
            if hasattr(params, 'include_metabolites') and params.include_metabolites:
                num_metabolites = len(params.include_metabolites)
            else:
                # Estimate based on available data and filters
                num_metabolites = min(max_metabolites if max_metabolites > 0 else len(complete_df), len(complete_df))
            
            # Find longest metabolite name for width calculation
            metabolite_col = 'Name' if 'Name' in complete_df.columns else complete_df.columns[0]
            if hasattr(params, 'include_metabolites') and params.include_metabolites:
                # Use only the metabolites that will be plotted
                filtered_df = complete_df[complete_df[metabolite_col].isin(params.include_metabolites)]
            else:
                filtered_df = complete_df
            
            if len(filtered_df) > 0:
                longest_name = max(len(str(name)) for name in filtered_df[metabolite_col].head(num_metabolites))
            else:
                longest_name = 20  # fallback
            
            # Calculate dimensions
            # Width: base width + space for samples + space for longest name
            base_width = 3  # inches for colorbar and margins
            sample_width = num_samples * 0.15  # 0.15 inches per sample
            name_width = longest_name * 0.08  # 0.08 inches per character
            calculated_width = max(6, base_width + sample_width + name_width)
            
            # Height: base height + space for metabolites
            base_height = 1.5  # inches for title and margins
            metabolite_height = num_metabolites * 0.1  # 0.1 inches per metabolite
            calculated_height = max(3, base_height + metabolite_height)
            
            # Apply calculated sizes
            params.fig_width = min(calculated_width, 20)  # Cap at reasonable maximum
            params.fig_height = min(calculated_height, 15)  # Cap at reasonable maximum
            self.log_viz_message(f"Auto-calculated dimensions: {params.fig_width:.1f} x {params.fig_height:.1f} inches")
        else:
            # User has disabled auto-size, keep their manual settings
            self.log_viz_message(f"Using manual dimensions: {params.fig_width:.1f} x {params.fig_height:.1f} inches")
        
        result = run_heatmap_analysis(ctx, params)
        if result.files_created:
            self.log_viz_message(f"Heatmaps: {len(result.files_created)} files saved")
        for err in result.errors:
            self.log_viz_message(f"Heatmap error: {err}")
        if result.files_created and not result.errors:
            self.log_viz_message("✅ heatmap completed")
        elif result.files_created and result.errors:
            self.log_viz_message("⚠️ heatmap completed with errors")
        elif result.errors:
            self.log_viz_message("❌ heatmap failed")
        else:
            self.log_viz_message("⚠️ heatmap completed (no plots generated)")

    def _generate_roc_plots(self, complete_df, groups, sample_cols, sample_to_group, output_dir):
        from main_script.metabolites_visualization import run_roc_analysis
        ctx = self._build_common_context(complete_df, groups, sample_cols, sample_to_group, os.path.join(output_dir, 'roc'))
        os.makedirs(ctx.output_dir, exist_ok=True)
        if hasattr(self, 'viz_params') and 'roc' in self.viz_params:
            params = self.viz_params['roc']
        else:
            from main_script.metabolites_visualization import ROCParams
            params = ROCParams()
        result = run_roc_analysis(ctx, params)
        if result.files_created:
            self.log_viz_message(f"ROC: {len(result.files_created)} files saved")
        for err in result.errors:
            self.log_viz_message(f"ROC error: {err}")

    def create_metabolite_selection_panel(self, parent):
        """Create metabolite selection panel for custom lists."""
        panel = ttk.LabelFrame(parent, text="Metabolite Selection", padding=10)
        panel.pack(fill='x', padx=5, pady=2)
        
        # Custom metabolite list for specific plots
        ttk.Label(panel, text="Custom metabolite list:").pack(anchor='w')
        
        # File selection
        file_row = ttk.Frame(panel)
        file_row.pack(fill='x', pady=2)
        
        self.viz_metabolite_list_file = tk.StringVar()
        ttk.Entry(file_row, textvariable=self.viz_metabolite_list_file, width=25).pack(side='left', expand=True, fill='x')
        ttk.Button(file_row, text="Browse", command=self.browse_metabolite_list).pack(side='right', padx=(2, 0))
        
        # Text area for manual entry
        ttk.Label(panel, text="Or enter metabolite names (one per line):").pack(anchor='w', pady=(5, 0))
        self.viz_metabolite_text = scrolledtext.ScrolledText(panel, height=4, width=30)
        self.viz_metabolite_text.pack(fill='x', pady=2)
        
        # Apply to which plots
        ttk.Label(panel, text="Apply custom list to:").pack(anchor='w', pady=(5, 0))
        self.viz_custom_list_boxplot = tk.BooleanVar(value=False)
        self.viz_custom_list_bargraph = tk.BooleanVar(value=False)
        self.viz_custom_list_heatmap = tk.BooleanVar(value=False)
        self.viz_custom_list_roc = tk.BooleanVar(value=False)
        
        ttk.Checkbutton(panel, text="Boxplots", variable=self.viz_custom_list_boxplot).pack(anchor='w')
        ttk.Checkbutton(panel, text="Bar Graphs", variable=self.viz_custom_list_bargraph).pack(anchor='w')
        ttk.Checkbutton(panel, text="Heatmaps", variable=self.viz_custom_list_heatmap).pack(anchor='w')
        ttk.Checkbutton(panel, text="ROC curves", variable=self.viz_custom_list_roc).pack(anchor='w')

    def browse_metabolite_list(self, plot_type: Optional[str] = None):
        """Browse for a metabolite list file.

        Unified handler used by: general custom list (no plot_type) and per-plot list selectors
        (plot_type in {'boxplot','bargraph','heatmap','roc'}). Prevents duplicate method name errors.
        """
        filename = filedialog.askopenfilename(
            title=f"Select Metabolite List{' for ' + plot_type.title() if plot_type else ''}",
            filetypes=[("Excel files", "*.xlsx"), ("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not filename:
            return
        # Per-plot assignments
        if plot_type == 'boxplot':
            self.boxplot_custom_list.set(filename)
        elif plot_type == 'bargraph':
            self.bargraph_custom_list.set(filename)
        elif plot_type == 'heatmap':
            self.heatmap_custom_list.set(filename)
        elif plot_type == 'roc':
            self.roc_custom_list.set(filename)
        else:
            # Generic metabolite list (metabolite selection panel)
            if hasattr(self, 'viz_metabolite_list_file'):
                self.viz_metabolite_list_file.set(filename)
        # Attempt to load contents into text area if generic panel present
        if hasattr(self, 'viz_metabolite_text') and plot_type is None:
            metabolites = self.load_metabolite_list(filename)
            if metabolites:
                content = '\n'.join(metabolites)
                self.viz_metabolite_text.delete(1.0, tk.END)
                self.viz_metabolite_text.insert(1.0, content)
            else:
                messagebox.showerror("Error", "Failed to load metabolite list from file.")
                return
        self.log_viz_message(f"Selected metabolite list{(' for ' + plot_type) if plot_type else ''}: {filename}")

    def configure_viz_groups(self):
        """Open group configuration dialog with full group management."""
        # Create group configuration window
        group_window = tk.Toplevel(self.root)
        group_window.title("Configure Visualization Groups")
        group_window.geometry("900x500")
        group_window.minsize(700, 450)  # Much reduced minimum height
        group_window.transient(self.root)
        group_window.grab_set()
        
        # Create main canvas and scrollbar for entire dialog
        main_canvas = tk.Canvas(group_window, highlightthickness=0)
        main_scrollbar = ttk.Scrollbar(group_window, orient="vertical", command=main_canvas.yview)
        main_frame = ttk.Frame(main_canvas)
        
        main_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        )
        
        main_canvas.create_window((0, 0), window=main_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=main_scrollbar.set)
        
        # Enable mouse wheel scrolling for entire dialog
        def _on_main_mousewheel(event):
            try:
                if main_canvas.winfo_exists():
                    main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except tk.TclError:
                pass
        
        main_canvas.bind("<MouseWheel>", _on_main_mousewheel)
        
        main_canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        main_scrollbar.pack(side="right", fill="y", pady=10)
        
        # Instructions with buttons on the right
        inst_frame = ttk.Frame(main_frame)
        inst_frame.pack(fill='x', pady=(0, 10))
        
        # Left side - instructions
        inst_left = ttk.Frame(inst_frame)
        inst_left.pack(side='left', fill='x', expand=True)
        ttk.Label(inst_left, text="Configure Groups and Sample Assignments", 
                 font=('Arial', 12, 'bold')).pack(anchor='w')
        ttk.Label(inst_left, text="Step 1: Define your groups, Step 2: Assign samples to groups", 
                 font=('Arial', 9)).pack(anchor='w')
        
        # Right side - Cancel and Apply buttons
        inst_right = ttk.Frame(inst_frame)
        inst_right.pack(side='right', padx=(10, 0))
        
        def apply_groups():
            # Keep definitions in sync with any in-dialog label edits.
            self._sync_viz_group_definitions_from_entries()
            valid_labels = set(self.viz_group_definitions.values())

            # Build group mapping using display names
            group_mapping = {}
            for col, var in self.viz_sample_group_vars.items():
                selected = var.get()
                if selected in valid_labels and selected != 'Ignore':
                    group_mapping[col] = selected
            
            # Track which columns are currently being used as samples
            current_samples = set(group_mapping.keys())
            
            # Update visualization group mapping
            self.viz_group_mapping = group_mapping
            # Keep preferred order synchronized with the current configured labels order.
            try:
                if hasattr(self, 'viz_preferred_group_order'):
                    ordered_labels = list(self.viz_group_definitions.values())
                    self.viz_preferred_group_order.set(','.join(ordered_labels))
            except Exception:
                pass
            # Preserve configured order when listing groups
            groups = self.ordered_groups(list(dict.fromkeys(group_mapping.values())))
            
            # Update status
            status_groups = ', '.join(groups) if groups else 'None (all Ignore)'
            self.viz_group_status.set(f"Groups configured: {status_groups} ({len(group_mapping)} samples)")
            self.update_group_color_controls(groups)
            self.log_viz_message(f"Groups configured: {len(groups)} groups, {len(group_mapping)} samples")
            try:
                self._save_viz_config()
            except Exception:
                pass
            group_window.destroy()
        
        ttk.Button(inst_right, text="Apply Groups", command=apply_groups).pack(side='right', padx=2)
        ttk.Button(inst_right, text="Cancel", command=group_window.destroy).pack(side='right', padx=2)
        
        # Get current data
        complete_df = self.get_current_complete_df()
        if complete_df is None:
            messagebox.showwarning("No Data", "Please load statistical results first.")
            group_window.destroy()
            return
        
        # PRIORITY 1: Use verified sample columns from Visualization's own Verify Columns step
        # All feature detection should be done in Verify Columns step - trust that assignment
        sample_cols = []
        if hasattr(self, 'memory_store') and self.memory_store:
            # Check for columns verified in Visualization tab's Verify Columns
            if 'viz_verified_sample_cols' in self.memory_store and self.memory_store['viz_verified_sample_cols']:
                sample_cols = [col for col in self.memory_store['viz_verified_sample_cols'] 
                              if col in complete_df.columns]
                if sample_cols:
                    self.log_viz_message(f"📊 Using {len(sample_cols)} verified sample columns from Verify Columns")
            # PRIORITY 2: Check for sample columns from Statistics tab (memory mode)
            elif 'sample_to_group' in self.memory_store and self.memory_store['sample_to_group']:
                sample_cols = [col for col in self.memory_store['sample_to_group'].keys() 
                              if col in complete_df.columns]
                if sample_cols:
                    self.log_viz_message(f"📊 Using {len(sample_cols)} verified sample columns from Statistics tab")
            # PRIORITY 3: Check for preprocessed sample columns
            elif 'preprocessed_sample_cols' in self.memory_store:
                sample_cols = [col for col in self.memory_store['preprocessed_sample_cols'] 
                              if col in complete_df.columns]
                if sample_cols:
                    self.log_viz_message(f"📊 Using {len(sample_cols)} preprocessed sample columns")
        
        # PRIORITY 4: Use detected_sample_cols if available
        if not sample_cols and hasattr(self, 'detected_sample_cols') and self.detected_sample_cols:
            sample_cols = [col for col in self.detected_sample_cols if col in complete_df.columns]
            if sample_cols:
                self.log_viz_message(f"📊 Using {len(sample_cols)} detected sample columns")
        
        # PRIORITY 5: Fallback to auto-detection only if no verified columns available
        if not sample_cols:
            self.log_viz_message("⚠️ No verified sample columns found. Please use '🔍 Verify Columns' first.")
            self.log_viz_message("   Auto-detection disabled to prevent incorrect column assignments.")
            # Return early - user must verify columns first for file imports
            if self.viz_data_source.get() == 'file':
                messagebox.showwarning(
                    "Columns Not Verified",
                    "Please click '🔍 Verify Columns' first to define which columns are sample data.\n\n"
                    "This ensures correct group assignment without auto-detection errors."
                )
                group_window.destroy()
                return
        
        # If sample columns are not found, do NOT block the user.
        # Proceed to allow manual configuration and/or pairwise-based group detection.
        groups_from_pairwise = []
        if not sample_cols:
            # Attempt to extract groups from pairwise statistical columns
            pairwise_groups = set()
            for col in complete_df.columns:
                if '_vs_' in col:
                    # Accept any pairwise stat suffix, not just log2FC
                    base = col
                    for suf in ['_log2FC', '_FC', '_pvalue', '_p_value', '_adj_p', '_p_adj', '_neg_log10_p', '_neg_log10_adj_p']:
                        base = base.replace(suf, '')
                    parts = base.split('_vs_')
                    if len(parts) == 2:
                        pairwise_groups.add(parts[0])
                        pairwise_groups.add(parts[1])
            if pairwise_groups:
                groups_from_pairwise = self.ordered_groups(list(pairwise_groups))
                self.log_viz_message(f"ℹ️ Proceeding without sample columns. Detected groups from pairwise columns: {groups_from_pairwise}")
                # Initialize definitions if not set
                if not hasattr(self, 'viz_group_definitions') or not self.viz_group_definitions:
                    self.viz_group_definitions = {}
                    for i, group_name in enumerate(groups_from_pairwise, 1):
                        self.viz_group_definitions[f"Group{i}"] = group_name
                # Keep mapping empty (no raw samples), user can still Apply to confirm
                self.viz_group_mapping = {}
                # Update visuals
                self.viz_group_status.set(f"Groups (no samples): {', '.join(groups_from_pairwise)}")
                self.update_group_color_controls(groups_from_pairwise)
            else:
                self.log_viz_message("ℹ️ Proceeding without sample columns and no pairwise-based group detection.")
        # Track all available sample candidates for exclusion tracking
        self.viz_all_sample_candidates = set(sample_cols)
        
        # Initialize group definitions with loaded groups or defaults
        if not hasattr(self, 'viz_group_definitions') or not self.viz_group_definitions:
            # PRIORITY 1: Try to sync from Statistics tab if available
            if hasattr(self, 'group_definitions') and self.group_definitions:
                # Use Statistics group definitions
                self.viz_group_definitions = self.group_definitions.copy()
                self.log_viz_message(f"📊 Loaded group definitions from Statistics tab: {list(self.viz_group_definitions.values())}")
            # PRIORITY 2: Try to get groups from current loaded data
            elif (groups := self.get_viz_groups()) and len(groups) > 0:
                # Use actual loaded groups and preserve configured order
                self.viz_group_definitions = {}
                for i, group_name in enumerate(self.ordered_groups(groups), 1):
                    group_id = f"Group{i}"
                    self.viz_group_definitions[group_id] = group_name
            else:
                # PRIORITY 3: Try to detect common group patterns from column names
                detected_groups = {}
                group_patterns = {
                    'Control': ['ctrl', 'control', 'con', 'c_'],
                    'Ortho': ['ortho', 'orthopedic', 'orth', 'or_'],
                    'TBI': ['tbi', 'brain', 'injury', 'tb_'],
                    'Treatment': ['treat', 'trt', 'tx', 'drug']
                }
                
                for col in sample_cols:
                    col_lower = col.lower()
                    for group_name, patterns in group_patterns.items():
                        if any(pattern in col_lower for pattern in patterns):
                            detected_groups[group_name] = group_name
                            break
                
                if detected_groups:
                    # Create group IDs from detected groups
                    self.viz_group_definitions = {}
                    for i, group_name in enumerate(self.ordered_groups(list(detected_groups.keys())), 1):
                        group_id = f"Group{i}"
                        self.viz_group_definitions[group_id] = group_name
                    self.log_viz_message(f"Auto-detected groups from column names: {', '.join(detected_groups.keys())}")
                else:
                    # PRIORITY 4: Fall back to defaults
                    self.viz_group_definitions = {
                        'Group1': 'Control',
                        'Group2': 'Disease', 
                        'Group3': 'Treatment'
                    }
        
        # Group definition section - USING STATISTICS TAB STYLE
        group_def_frame = ttk.LabelFrame(main_frame, text="Group IDs & Labels", padding=10)
        group_def_frame.pack(fill='x', pady=(0, 10))
        
        # Initialize group count if not present
        if not hasattr(self, 'group_count'):
            self.group_count = len(self.viz_group_definitions)
        
        # Add/Remove group buttons at the top
        buttons_frame = tk.Frame(group_def_frame, bg='#f0f0f0')
        buttons_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Button(buttons_frame, text='+ Add Group', command=self.add_group,
                  bg='#27ae60', fg='white', font=('Arial', 9, 'bold'), pady=4, padx=15).pack(side='left', padx=5)
        tk.Button(buttons_frame, text='- Remove', command=self.remove_group,
                  bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'), pady=4, padx=15).pack(side='left', padx=5)
        
        # Canvas for group IDs (scrollable list)
        groups_canvas = tk.Canvas(group_def_frame, bg='#f0f0f0', highlightthickness=0, height=120)
        groups_scrollbar = ttk.Scrollbar(group_def_frame, orient='vertical', command=groups_canvas.yview)
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
        
        # Refresh group UI (will create the Group 1:, Group 2:, etc entries)
        self.group_id_vars = {}
        self.refresh_group_ui()
        
        # Auto-assign button at bottom of Step 1
        auto_assign_frame = ttk.Frame(main_frame)
        auto_assign_frame.pack(fill='x', pady=(10, 0))
        ttk.Button(auto_assign_frame, text="🎯 Auto-Assign by Pattern", 
                  command=lambda: self.auto_assign_viz_groups()).pack(side='left')
        
        # Sample assignment section
        assign_frame = ttk.LabelFrame(main_frame, text="Step 2: Assign Samples to Groups", padding=10)
        assign_frame.pack(fill='both', expand=True)
        
        # Filter out verified statistics columns from sample display
        verified_stats_cols = set()
        if hasattr(self, 'verified_assignments') and self.verified_assignments:
            verified_stats_cols = set(self.verified_assignments.values())
        
        # Filter sample_cols to exclude verified statistics columns
        filtered_sample_cols = [col for col in sample_cols if col not in verified_stats_cols]
        
        if len(filtered_sample_cols) < len(sample_cols):
            excluded = len(sample_cols) - len(filtered_sample_cols)
            self.log_viz_message(f"✓ Excluding {excluded} verified statistics column(s) from sample assignment")
        
        sample_cols = filtered_sample_cols
        
        # Direct frame for sample assignments (no internal scrollbar needed)
        scrollable_frame = ttk.Frame(assign_frame)
        scrollable_frame.pack(fill='both', expand=True)
        
        # Header
        header_frame = ttk.Frame(scrollable_frame)
        header_frame.pack(fill='x', pady=(0, 5))
        ttk.Label(header_frame, text="Sample Column", font=('Arial', 10, 'bold')).pack(side='left', padx=20)
        ttk.Label(header_frame, text="Assigned Group", font=('Arial', 10, 'bold')).pack(side='left', padx=(200, 0))
        ttk.Label(header_frame, text="Action", font=('Arial', 10, 'bold')).pack(side='right', padx=20)
        
        # Sample assignment variables
        self.viz_sample_group_vars = {}
        # Keep explicit references to combobox widgets so we can update their values reliably
        self.viz_sample_comboboxes = {}
        
        # Sample assignments
        for col in sample_cols:
            row_frame = ttk.Frame(scrollable_frame)
            row_frame.pack(fill='x', padx=10, pady=2)
            
            # Sample column name (left side, wider)
            ttk.Label(row_frame, text=col, width=50).pack(side='left', padx=(10, 50))
            
            # Get current assignment - try multiple sources in priority order
            current_assignment = None
            
            # PRIORITY 1: Check if Statistics has this sample assigned
            if hasattr(self, 'sample_group_vars') and col in self.sample_group_vars:
                stats_assignment = self.sample_group_vars[col].get()
                # Check if this assignment exists in viz_group_definitions (by label or ID)
                if stats_assignment in self.viz_group_definitions.values():
                    current_assignment = stats_assignment
                elif stats_assignment in self.viz_group_definitions.keys():
                    # Convert group ID to group label
                    current_assignment = self.viz_group_definitions[stats_assignment]
            
            # PRIORITY 2: Check existing viz_group_mapping
            if not current_assignment and hasattr(self, 'viz_group_mapping') and col in self.viz_group_mapping:
                current_assignment = self.viz_group_mapping[col]
            
            # PRIORITY 3: Default to 'Ignore' (user should explicitly assign groups)
            if current_assignment and current_assignment in self.viz_group_definitions.values():
                initial_value = current_assignment
            elif current_assignment == 'Ignore':
                initial_value = 'Ignore'
            else:
                initial_value = 'Ignore'  # Default all samples to Ignore
            
            var = tk.StringVar(value=initial_value)
            # Auto-save when sample group assignment changes
            var.trace_add('write', lambda *args: self._viz_config_changed())
            self.viz_sample_group_vars[col] = var
            combo = ttk.Combobox(row_frame, textvariable=var, width=20, state='readonly')
            # Add 'Ignore' as first option in dropdown
            combo['values'] = ['Ignore'] + list(self.viz_group_definitions.values())
            combo.pack(side='left', padx=(0, 50))
            # store combobox reference
            self.viz_sample_comboboxes[col] = combo
            
            # Delete sample button (right side)
            ttk.Button(row_frame, text="🗑️", width=3,
                      command=lambda c=col, rf=row_frame: self.delete_viz_sample_column(c, rf)).pack(side='right', padx=10)
        
        # Initialize the display
        self.refresh_viz_group_display()

        # NOTE: Do NOT load config from disk here - use current in-memory state
        # Config is loaded once at tab initialization. Loading here would overwrite
        # user changes made during the session (like deleted groups).
    
    def add_group(self):
        """Add a new group to the group definitions"""
        self.group_count += 1

        # Keep key format consistent with existing definitions (Group1 vs Group 1).
        has_spaced = any(str(k).startswith('Group ') for k in self.viz_group_definitions.keys())
        new_group_id = f'Group {self.group_count}' if has_spaced else f'Group{self.group_count}'
        new_group_label = f'Group{self.group_count}'

        self.viz_group_definitions[new_group_id] = new_group_label
        self.group_id_vars[new_group_id] = tk.StringVar(value=new_group_label)

        # Immediate live update in the open dialog.
        self.refresh_viz_group_display()
        self.log_viz_message(f"✅ Added {new_group_id}")
        try:
            self._save_viz_config()
        except Exception:
            pass
    
    def remove_group(self):
        """Remove the last group from group definitions"""
        if self.group_count <= 2:  # Don't allow less than 2 groups
            messagebox.showwarning("Minimum Groups", "At least 2 groups are required for visualization.")
            return

        # Support both key styles used across the codebase.
        candidates = [f'Group{self.group_count}', f'Group {self.group_count}']
        last_group_id = next((gid for gid in candidates if gid in self.viz_group_definitions), None)
        if not last_group_id:
            # Fallback: remove last inserted key to avoid stale count drift.
            keys = list(self.viz_group_definitions.keys())
            if not keys:
                return
            last_group_id = keys[-1]

        removed_label = self.viz_group_definitions.get(last_group_id, last_group_id)
        del self.viz_group_definitions[last_group_id]
        if last_group_id in self.group_id_vars:
            del self.group_id_vars[last_group_id]

        # If sample rows were assigned to removed label, set them to Ignore immediately.
        if hasattr(self, 'viz_sample_group_vars') and self.viz_sample_group_vars:
            for _, var in self.viz_sample_group_vars.items():
                try:
                    if var.get() == removed_label:
                        var.set('Ignore')
                except Exception:
                    pass

        self.group_count = max(2, len(self.viz_group_definitions))
        # Immediate live update in the open dialog.
        self.refresh_viz_group_display()
        self.log_viz_message(f"❌ Removed {last_group_id}")
        try:
            self._save_viz_config()
        except Exception:
            pass
    
    def refresh_group_ui(self):
        """Refresh the group UI with current group definitions"""
        # Clear existing group entries
        if not hasattr(self, 'groups_scrollable_frame'):
            return
        
        for widget in self.groups_scrollable_frame.winfo_children():
            widget.destroy()
        
        # Recreate group entries
        for i, (group_id, default_label) in enumerate(self.viz_group_definitions.items()):
            id_frame = tk.Frame(self.groups_scrollable_frame, bg='#f0f0f0')
            id_frame.pack(fill='x', padx=3, pady=2)
            tk.Label(id_frame, text=f'{group_id}:', bg='#f0f0f0', width=10).pack(side='left')
            
            if group_id not in self.group_id_vars:
                self.group_id_vars[group_id] = tk.StringVar(value=default_label)
            
            entry_var = self.group_id_vars[group_id]
            tk.Entry(id_frame, textvariable=entry_var, font=('Arial', 9), width=20).pack(side='left', padx=(5,5))
            
            # Trace for auto-updating labels
            if not hasattr(self, '_group_label_trace_ids'):
                self._group_label_trace_ids = {}
            if group_id not in self._group_label_trace_ids:
                def _on_label_change(*args, gid=group_id):
                    if getattr(self, '_suppress_group_label_trace', False) or getattr(self, '_shutting_down', False):
                        return
                    try:
                        new_label = self.group_id_vars[gid].get().strip()
                        if new_label:
                            old_label = self.viz_group_definitions.get(gid, gid)
                            self.viz_group_definitions[gid] = new_label
                            
                            # Update all sample assignment dropdowns
                            display_labels = [self.viz_group_definitions[g] for g in self.viz_group_definitions.keys()]
                            if hasattr(self, 'viz_sample_comboboxes'):
                                for sample_col, combo in self.viz_sample_comboboxes.items():
                                    if combo and combo.winfo_exists():
                                        combo['values'] = ['Ignore'] + display_labels
                                        current_selection = self.viz_sample_group_vars.get(sample_col)
                                        if current_selection and current_selection.get() == old_label:
                                            current_selection.set(new_label)
                    except Exception:
                        pass
                    
                    try:
                        if hasattr(self, '_group_label_change_after'):
                            self.root.after_cancel(self._group_label_change_after)
                        self._group_label_change_after = self.root.after(200, lambda: self._save_viz_config() if hasattr(self, '_save_viz_config') else None)
                    except Exception:
                        pass
                
                trace_id = entry_var.trace_add('write', _on_label_change)
                self._group_label_trace_ids[group_id] = trace_id
        
        # Update scroll region
        if hasattr(self, 'groups_canvas'):
            try:
                self.groups_canvas.update_idletasks()
                self.groups_scrollable_frame.update_idletasks()
                width = self.groups_scrollable_frame.winfo_reqwidth()
                height = self.groups_scrollable_frame.winfo_reqheight()
                self.groups_canvas.configure(scrollregion=(0, 0, width, height))
            except Exception:
                pass

    def verify_viz_columns(self):
        """Verify and assign columns for visualization using unified dialog.
        
        This allows users to explicitly define sample vs feature columns,
        which will be used downstream in Configure Groups without any filtering.
        """
        import threading
        from gui.shared.column_assignment import show_column_assignment_dialog, ColumnDetector
        
        def _verify_thread():
            try:
                # Get current data
                complete_df = self.get_current_complete_df()
                if complete_df is None:
                    self.root.after(0, lambda: messagebox.showwarning(
                        "No Data", 
                        "Please load data first using 'Load & Analyze' button."
                    ))
                    return
                
                mode = self.viz_mode.get()
                
                # Check if data was autoloaded from Statistics tab and already verified
                # If so, skipverification and use existing verified columns
                if hasattr(self, 'memory_store') and self.memory_store:
                    # Check for verified columns from Statistics tab
                    if mode == 'lipid':
                        pos_verified = self.memory_store.get('verified_pos_lipid_sample_cols', [])
                        neg_verified = self.memory_store.get('verified_neg_lipid_sample_cols', [])
                    else:
                        pos_verified = self.memory_store.get('verified_pos_sample_cols', [])
                        neg_verified = self.memory_store.get('verified_neg_sample_cols', [])
                    
                    # If we have verified columns from Statistics, use them directly
                    if pos_verified or neg_verified:
                        verified_sample_cols = list(set(pos_verified + neg_verified))
                        self.memory_store['viz_verified_sample_cols'] = verified_sample_cols
                        self.detected_sample_cols = verified_sample_cols
                        self._columns_verified = True
                        
                        self.log_viz_message("\n✅ Using pre-verified columns from Statistics tab!")
                        self.log_viz_message(f"   • Sample columns: {len(verified_sample_cols)}")
                        self.log_viz_message(f"\nNext: Click '⚙️ Configure Groups' to assign samples to groups.\n")
                        
                        # Enable Configure Groups button
                        self.root.after(0, lambda: self.configure_groups_btn.config(state='normal'))
                        return
                
                # Determine tab type for dialog
                if mode == 'lipid':
                    tab_type = 'statistics_lipid'
                else:
                    tab_type = 'statistics_metabolite'
                
                self.log_viz_message("\n🔍 Verifying columns...")
                
                # Use the same shared detector as Statistics tab for consistent behavior.
                auto_sample_cols = ColumnDetector.detect_sample_columns(complete_df)
                
                # Show column assignment dialog
                result = show_column_assignment_dialog(
                    parent=self.root,
                    df=complete_df,
                    tab_type=tab_type,
                    auto_calculate=False,
                    dialog_title=f"Visualization - Verify Columns ({mode.title()} Mode)",
                    detected_sample_cols=auto_sample_cols,
                )
                
                if not result:
                    self.log_viz_message("❌ Verification cancelled\n")
                    return
                
                # Store verified columns
                verified_sample_cols = result.get('sample_cols', [])
                verified_assignments = result.get('assignments', {})
                feature_id_col = result.get('feature_id_col')
                
                # Store in memory for use by Configure Groups
                if not hasattr(self, 'memory_store') or self.memory_store is None:
                    self.memory_store = {}
                
                self.memory_store['viz_verified_sample_cols'] = verified_sample_cols
                self.memory_store['viz_verified_assignments'] = verified_assignments
                if feature_id_col:
                    self.memory_store['viz_feature_id_col'] = feature_id_col
                    self._normalize_feature_id_column(complete_df, feature_id_col)
                    if self.viz_data_source.get() == 'file' and hasattr(self, 'imported_complete_df'):
                        self.imported_complete_df = complete_df
                    elif hasattr(self, 'memory_store') and self.memory_store:
                        stat_results = self.memory_store.get('statistical_test_results')
                        if isinstance(stat_results, dict) and 'enhanced_metabolites' in stat_results:
                            stat_results['enhanced_metabolites'] = complete_df
                    if hasattr(self, 'statistical_test_results') and isinstance(getattr(self, 'statistical_test_results', None), dict) and 'enhanced_metabolites' in self.statistical_test_results:
                        self.statistical_test_results['enhanced_metabolites'] = complete_df
                
                # Also store as detected_sample_cols for backward compatibility
                self.detected_sample_cols = verified_sample_cols
                
                # Log success
                self.log_viz_message(f"\n✅ Column verification complete!")
                self.log_viz_message(f"   • Verified sample columns: {len(verified_sample_cols)}")
                if feature_id_col:
                    self.log_viz_message(f"   • Feature ID column: {feature_id_col}")
                self.log_viz_message(f"\nNext: Click '⚙️ Configure Groups' to assign samples to groups.\n")
                
                # Enable Configure Groups button
                self.root.after(0, lambda: self.configure_groups_btn.config(state='normal'))
                
                # Update any status indicators
                self._columns_verified = True
                
            except Exception as e:
                import traceback
                self.log_viz_message(f"\n❌ Verification error: {str(e)}")
                self.log_viz_message(traceback.format_exc())
                self.root.after(0, lambda: messagebox.showerror(
                    "Verification Error", 
                    f"Failed to verify columns: {str(e)}"
                ))
        
        # Run in background thread
        thread = threading.Thread(target=_verify_thread, daemon=True)
        thread.start()

    def _configure_groups_wrapper(self):
        """Wrapper for configure_viz_groups that updates state after configuration."""
        self.configure_viz_groups()
        # After the dialog closes, check if groups are configured
        self._update_configuration_state()

    def _update_configuration_state(self):
        """Update the configuration state and enable/disable generate buttons."""
        # Check if groups are configured
        self.groups_configured = bool(getattr(self, 'viz_group_mapping', None)) and len(self.viz_group_mapping) > 0
        
        # Check if stat columns are configured
        # Require at least one comparison with "ready" or "partial" status (has pvalue and/or FC)
        comparisons = self.stat_column_assignments.get('comparisons', {})
        if comparisons:
            # Check if at least one comparison has pvalue or FC assigned
            has_valid_comparison = any(
                comp.get('status') in ['ready', 'partial'] or comp.get('pvalue') or comp.get('log2fc') or comp.get('fc')
                for comp in comparisons.values()
            )
            self.stat_cols_configured = has_valid_comparison
        else:
            self.stat_cols_configured = False
        
        # Update status indicators
        if hasattr(self, 'groups_status_label'):
            self.groups_status_label.config(text="✅" if self.groups_configured else "⚠️")
        if hasattr(self, 'stat_cols_status_label'):
            self.stat_cols_status_label.config(text="✅" if self.stat_cols_configured else "⚠️")
        
        # Update config status message
        if self.groups_configured and self.stat_cols_configured:
            self.config_status.set("✅ Ready to generate plots")
        elif not self.groups_configured:
            self.config_status.set("⚠️ Step 1: Configure groups first")
        else:
            self.config_status.set("⚠️ Step 2: Configure stat columns")
        
        # Enable/disable generate buttons
        can_generate = self.groups_configured and self.stat_cols_configured
        state = 'normal' if can_generate else 'disabled'
        if hasattr(self, 'generate_selected_btn'):
            self.generate_selected_btn.config(state=state)
        if hasattr(self, 'generate_all_btn'):
            self.generate_all_btn.config(state=state)

    def _get_all_possible_comparisons(self):
        """Get all possible pairwise comparisons from configured groups."""
        if not hasattr(self, 'viz_group_mapping') or not self.viz_group_mapping:
            return []
        
        groups = list(dict.fromkeys(self.viz_group_mapping.values()))
        groups = self.ordered_groups(groups)
        
        comparisons = []
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                comparisons.append((groups[i], groups[j]))
        
        return comparisons

    def _auto_detect_stat_columns(self, df, g1, g2):
        """Auto-detect statistical columns for a comparison.
        
        Returns dict with keys: pvalue, pvalue_adj, log2fc, fc, neglogp
        """
        prefix_variants = [
            f"{g1}_vs_{g2}",
            f"{g2}_vs_{g1}",
            f"{g1}vs{g2}",
            f"{g2}vs{g1}",
        ]
        
        result = {
            'pvalue': None,
            'pvalue_adj': None,
            'log2fc': None,
            'fc': None,
            'neglogp': None,
            'status': 'missing'
        }
        
        columns_lower = {c.lower(): c for c in df.columns}
        
        for prefix in prefix_variants:
            prefix_lower = prefix.lower()
            
            # Look for p-value columns
            pvalue_suffixes = ['_pvalue', '_p_value', '_p', '_rawp']
            for suffix in pvalue_suffixes:
                test_col = f"{prefix_lower}{suffix}"
                if test_col in columns_lower and not result['pvalue']:
                    result['pvalue'] = columns_lower[test_col]
            
            # Look for adjusted p-value columns
            padj_suffixes = ['_adj_p', '_p_adj', '_padj', '_fdr', '_adj_pvalue', '_pvalue_adj', '_adjpvalue']
            for suffix in padj_suffixes:
                test_col = f"{prefix_lower}{suffix}"
                if test_col in columns_lower and not result['pvalue_adj']:
                    result['pvalue_adj'] = columns_lower[test_col]
            
            # Look for log2FC columns
            log2fc_suffixes = ['_log2fc', '_log2_fc', '_log2foldchange']
            for suffix in log2fc_suffixes:
                test_col = f"{prefix_lower}{suffix}"
                if test_col in columns_lower and not result['log2fc']:
                    result['log2fc'] = columns_lower[test_col]
            
            # Look for FC columns (non-log)
            fc_suffixes = ['_fc', '_foldchange', '_fold_change']
            for suffix in fc_suffixes:
                test_col = f"{prefix_lower}{suffix}"
                if test_col in columns_lower and not result['fc']:
                    # Make sure it's not the log2fc column
                    actual_col = columns_lower[test_col]
                    if 'log2' not in actual_col.lower():
                        result['fc'] = actual_col
            
            # Look for neg log p columns
            neglog_suffixes = ['_neg_log10_p', '_neglogp', '_neg_log10_pvalue', '_neglog10p', '_neg_log10_adj_p']
            for suffix in neglog_suffixes:
                test_col = f"{prefix_lower}{suffix}"
                if test_col in columns_lower and not result['neglogp']:
                    result['neglogp'] = columns_lower[test_col]
        
        # Determine status
        has_pvalue = result['pvalue'] or result['pvalue_adj']
        has_fc = result['log2fc'] or result['fc']
        
        if has_pvalue and has_fc:
            result['status'] = 'ready'
        elif has_pvalue or has_fc:
            result['status'] = 'partial'
        else:
            result['status'] = 'missing'
        
        return result

    def _auto_configure_stat_columns(self, complete_df, groups):
        """Automatically configure stat columns based on auto-detection.
        
        This runs silently when data is loaded to set up stat_column_assignments.
        """
        try:
            from itertools import combinations
            
            # Find ID column - only set if a clear candidate exists
            # Don't default to first non-numeric column to avoid incorrect assumptions
            id_col = None
            for candidate in ['metabolite_id', 'Name', 'Metabolite', 'LipidMolec', 'Compound', 'ID']:
                if candidate in complete_df.columns:
                    id_col = candidate
                    self.log_viz_message(f"ℹ️ Auto-detected ID column: {id_col} (can be changed in Configure Stat Columns)", 'INFO')
                    break
            
            if id_col is None:
                # Try fallback: find first non-numeric, non-statistical column
                for col in complete_df.columns:
                    # Skip statistical columns
                    if any(marker in str(col).lower() for marker in ['_vs_', '_log2fc', '_adj_p', '_p_value', '_pvalue', '_fc', 'fold']):
                        continue
                    # Skip purely numeric columns
                    try:
                        pd.to_numeric(complete_df[col], errors='raise')
                        continue  # Skip numeric columns
                    except (ValueError, TypeError):
                        # This is a non-numeric column - good candidate for ID
                        id_col = col
                        self.log_viz_message(f"ℹ️ Using fallback ID column: {id_col} (can be changed in Configure Stat Columns)", 'INFO')
                        break
            
            if id_col is None:
                self.log_viz_message("⚠️ Could not auto-detect ID column. Please use 'Configure Stat Columns' dialog to specify.", 'WARNING')
                # Still continue with comparisons, but leave ID as None
                # This allows volcano plots and other analyses that don't strictly require IDs
            
            # Generate comparisons
            ordered_groups = self.ordered_groups(groups)
            comparisons = []
            for i in range(len(ordered_groups)):
                for j in range(i + 1, len(ordered_groups)):
                    comparisons.append((ordered_groups[i], ordered_groups[j]))
            
            if not comparisons:
                return
            
            # Auto-detect columns for each comparison
            self.stat_column_assignments = {
                'id_column': id_col,
                'comparisons': {}
            }
            
            ready_count = 0
            partial_count = 0
            missing_count = 0
            
            for g1, g2 in comparisons:
                detected = self._auto_detect_stat_columns(complete_df, g1, g2)
                
                # Prioritize adj_pvalue over pvalue
                pval = detected['pvalue_adj'] or detected['pvalue']
                log2fc = detected['log2fc']
                fc = detected['fc']
                
                has_pval = bool(pval)
                has_fc = bool(log2fc) or bool(fc)
                
                if has_pval and has_fc:
                    status = 'ready'
                    ready_count += 1
                elif has_pval or has_fc:
                    status = 'partial'
                    partial_count += 1
                else:
                    status = 'missing'
                    missing_count += 1
                
                self.stat_column_assignments['comparisons'][(g1, g2)] = {
                    'pvalue': pval,
                    'log2fc': log2fc,
                    'fc': fc,
                    'status': status
                }
            
            # Mark as configured if at least one comparison is ready
            if ready_count > 0:
                self.stat_cols_configured = True
                self.log_viz_message(f"✅ Auto-configured stat columns: {ready_count} ready, {partial_count} partial, {missing_count} missing", 'INFO')
            else:
                self.stat_cols_configured = False
                self.log_viz_message(f"⚠️ Stat columns need manual configuration (0 ready, {partial_count} partial, {missing_count} missing)", 'WARNING')
        
        except Exception as e:
            self.log_viz_message(f"⚠️ Auto-configure stat columns failed: {str(e)}", 'WARNING')
            self.stat_cols_configured = False

    def _is_statistical_feature_column(self, col_name):
        """Check if a column is a statistical feature that should be excluded from sample columns.
        
        Returns True for columns like r_squared, adj_r_squared, Control_Mean, PD_Mean, 
        n_Control, f_statistic, f_pvalue, status, etc.
        
        Does NOT exclude comparison-based columns like Control_vs_PD_pvalue, Control_vs_PD_log2FC.
        """
        col_lower = col_name.lower()
        
        # If it's a comparison column (contains _vs_), it's a valid statistical result column
        if '_vs_' in col_lower:
            return False
        
        # Exact match patterns for standalone statistical metadata
        exact_patterns = [
            'r_squared', 'adj_r_squared', 'r2', 'adj_r2',
            'f_statistic', 'f_pvalue', 'f_stat', 'f_p',
            'status', 'stat_status', 'result_status'
        ]
        
        # Prefix/suffix patterns (case-insensitive)
        prefix_patterns = ['n_', 'count_', 'num_']
        suffix_patterns = ['_mean', '_median', '_std', '_sd', '_sem', '_n', '_count']
        
        # Check exact matches
        if col_lower in exact_patterns:
            return True
        
        # Check prefix patterns (e.g., n_Control, n_PD, count_Treatment)
        if any(col_lower.startswith(prefix) for prefix in prefix_patterns):
            return True
        
        # Check suffix patterns (e.g., Control_Mean, PD_Median, Treatment_std)
        if any(col_lower.endswith(suffix) for suffix in suffix_patterns):
            return True
        
        # Check for comparison-based statistical columns (e.g., Control_vs_PD_pvalue)
        if '_vs_' in col_lower:
            # Extract suffix after the comparison
            parts = col_lower.split('_vs_')
            if len(parts) >= 2:
                # Check if it ends with a statistical suffix
                remainder = parts[-1]
                stat_suffixes = ['pvalue', 'p_value', 'p', 'adj_p', 'padj', 'fdr', 
                                'fc', 'log2fc', 'neglogp', 'significant', 'sig']
                if any(remainder.endswith(f'_{suff}') or remainder.split('_')[-1] == suff 
                      for suff in stat_suffixes):
                    return True
        
        return False

    def configure_stat_columns(self):
        """Open dialog to configure statistical columns for each comparison."""
        if not self.groups_configured:
            messagebox.showwarning("Groups Required", 
                "Please configure groups first (Step 1) before configuring stat columns.")
            return
        
        complete_df = self.get_current_complete_df()
        if complete_df is None:
            messagebox.showwarning("No Data", "Please load data first.")
            return
        
        comparisons = self._get_all_possible_comparisons()
        if not comparisons:
            messagebox.showwarning("No Comparisons", 
                "No comparisons available. Please configure at least 2 groups.")
            return
        
        # Create dialog with wider layout for better readability
        dialog = tk.Toplevel(self.root)
        dialog.title("Configure Statistical Columns")
        dialog.geometry("1100x600")
        dialog.minsize(900, 500)
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Get feature columns (non-sample columns, excluding statistical features)
        sample_cols = set(self.viz_group_mapping.keys()) if hasattr(self, 'viz_group_mapping') else set()
        all_feature_cols = [c for c in complete_df.columns if c not in sample_cols]
        
        # Filter out statistical feature columns
        feature_cols = [c for c in all_feature_cols if not self._is_statistical_feature_column(c)]
        
        # ID column candidates: Show ALL non-sample columns (no filtering)
        id_candidates = ['--Select--'] + feature_cols
        
        # Main frame with scrollbar
        main_canvas = tk.Canvas(dialog, highlightthickness=0)
        scrollbar = ttk.Scrollbar(dialog, orient="vertical", command=main_canvas.yview)
        main_frame = ttk.Frame(main_canvas)
        
        # Configure canvas to expand the content window
        def on_frame_configure(event):
            main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        
        def on_canvas_configure(event):
            # Update the inner frame width to match canvas width
            canvas_width = event.width
            main_canvas.itemconfig(canvas_window, width=canvas_width)
        
        main_frame.bind("<Configure>", on_frame_configure)
        canvas_window = main_canvas.create_window((0, 0), window=main_frame, anchor="nw")
        main_canvas.bind("<Configure>", on_canvas_configure)
        main_canvas.configure(yscrollcommand=scrollbar.set)
        
        def _on_mousewheel(event):
            main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        main_canvas.bind("<MouseWheel>", _on_mousewheel)
        
        main_canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        scrollbar.pack(side="right", fill="y", pady=10)
        
        # Header
        header_frame = ttk.Frame(main_frame)
        header_frame.pack(fill='x', pady=(0, 10))
        ttk.Label(header_frame, text="Configure Statistical Columns", 
                 font=('Arial', 14, 'bold')).pack(side='left')
        
        # ID Column selection
        id_frame = ttk.LabelFrame(main_frame, text="Step 1: Identify Metabolite/Lipid ID Column", padding=10)
        id_frame.pack(fill='x', pady=5)
        
        id_var = tk.StringVar(value=self.stat_column_assignments.get('id_column') or '--Select--')
        ttk.Label(id_frame, text="ID Column:").pack(side='left', padx=5)
        id_combo = ttk.Combobox(id_frame, textvariable=id_var, values=id_candidates, 
                               state='readonly', width=40)
        id_combo.pack(side='left', padx=5)
        
        # Comparisons section
        comp_frame = ttk.LabelFrame(main_frame, text=f"Step 2: Assign Columns for {len(comparisons)} Comparisons", padding=10)
        comp_frame.pack(fill='both', expand=True, pady=5)
        
        # Column headers (increased widths for better readability)
        header_row = ttk.Frame(comp_frame)
        header_row.pack(fill='x', pady=(0, 5))
        ttk.Label(header_row, text="Comparison", width=22, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        ttk.Label(header_row, text="P-Value", width=26, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        ttk.Label(header_row, text="NegLogP (optional)", width=26, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        ttk.Label(header_row, text="Log2FC", width=26, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        ttk.Label(header_row, text="FC", width=22, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        ttk.Label(header_row, text="Status", width=12, font=('Arial', 9, 'bold')).pack(side='left', padx=3)
        
        # Dropdown options for stat columns
        stat_col_options = ['--None--'] + feature_cols
        
        # Store variable references
        comparison_vars = {}
        status_labels = {}
        
        def update_status(comp_key):
            """Update status indicator for a comparison (FC required)."""
            vars_dict = comparison_vars[comp_key]
            pval = vars_dict['pvalue'].get()
            fc = vars_dict['fc'].get()
            # log2FC optional if FC present; we can compute from FC later if missing
            has_pval = pval and pval != '--None--'
            has_fc = fc and fc != '--None--'
            if has_pval and has_fc:
                status_labels[comp_key].config(text="✅ Ready", foreground='green')
            elif has_pval and not has_fc:
                status_labels[comp_key].config(text="⚠️ FC required", foreground='orange')
            elif has_fc and not has_pval:
                status_labels[comp_key].config(text="⚠️ P-value required", foreground='orange')
            else:
                status_labels[comp_key].config(text="❌ Missing", foreground='red')
        
        for g1, g2 in comparisons:
            comp_key = (g1, g2)
            row = ttk.Frame(comp_frame)
            row.pack(fill='x', pady=2)
            
            # Comparison label (increased width for better visibility)
            ttk.Label(row, text=f"{g1} vs {g2}", width=22).pack(side='left', padx=3)
            
            # Auto-detect for this comparison
            detected = self._auto_detect_stat_columns(complete_df, g1, g2)
            
            # DEBUG: Log detection results
            print(f"[VIZ] Auto-detect for {g1} vs {g2}:")
            print(f"      pvalue: {detected['pvalue']}")
            print(f"      pvalue_adj: {detected['pvalue_adj']}")
            print(f"      log2fc: {detected['log2fc']}")
            print(f"      fc: {detected['fc']}")
            print(f"      neglogp: {detected['neglogp']}")
            print(f"      status: {detected['status']}")
            
            # P-value dropdown (prioritize adj_pvalue over pvalue)
            pval_detected = detected.get('pvalue_adj') or detected.get('pvalue')
            pval_var = tk.StringVar(value=pval_detected if pval_detected else '--None--')
            pval_combo = ttk.Combobox(row, textvariable=pval_var, values=stat_col_options, 
                                     state='readonly', width=24)
            pval_combo.pack(side='left', padx=3)

            # NegLogP dropdown (optional, increased width)
            neglogp_detected = detected.get('neglogp')
            neglogp_var = tk.StringVar(value=neglogp_detected if neglogp_detected else '--None--')
            neglogp_combo = ttk.Combobox(row, textvariable=neglogp_var, values=stat_col_options,
                                        state='readonly', width=24)
            neglogp_combo.pack(side='left', padx=3)

            # Log2FC dropdown (increased width)
            log2fc_detected = detected.get('log2fc')
            log2fc_var = tk.StringVar(value=log2fc_detected if log2fc_detected else '--None--')
            log2fc_combo = ttk.Combobox(row, textvariable=log2fc_var, values=stat_col_options, 
                                       state='readonly', width=24)
            log2fc_combo.pack(side='left', padx=3)
            
            # FC dropdown (mandatory, increased width)
            fc_detected = detected.get('fc')
            fc_var = tk.StringVar(value=fc_detected if fc_detected else '--None--')
            fc_combo = ttk.Combobox(row, textvariable=fc_var, values=stat_col_options, 
                                   state='readonly', width=20)
            fc_combo.pack(side='left', padx=3)
            
            # Status indicator (increased width)
            status_label = ttk.Label(row, text="", width=12)
            status_label.pack(side='left', padx=3)
            status_labels[comp_key] = status_label
            
            # Store vars
            comparison_vars[comp_key] = {
                'pvalue': pval_var,
                'neglogp': neglogp_var,
                'log2fc': log2fc_var,
                'fc': fc_var
            }
            
            # Bind updates
            pval_var.trace_add('write', lambda *args, k=comp_key: update_status(k))
            neglogp_var.trace_add('write', lambda *args, k=comp_key: update_status(k))
            log2fc_var.trace_add('write', lambda *args, k=comp_key: update_status(k))
            fc_var.trace_add('write', lambda *args, k=comp_key: update_status(k))
            
            # Initial status
            update_status(comp_key)
        
        # Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill='x', pady=10)
        
        def auto_detect_all():
            """Re-run auto-detection for all comparisons."""
            for (g1, g2), vars_dict in comparison_vars.items():
                detected = self._auto_detect_stat_columns(complete_df, g1, g2)
                # Prioritize adj_pvalue over pvalue
                vars_dict['pvalue'].set(detected['pvalue_adj'] or detected['pvalue'] or '--None--')
                vars_dict['neglogp'].set(detected['neglogp'] or '--None--')
                vars_dict['log2fc'].set(detected['log2fc'] or '--None--')
                vars_dict['fc'].set(detected['fc'] or '--None--')
            self.log_viz_message("Auto-detection refreshed for all comparisons")
        
        def apply_config():
            """Apply the stat column configuration."""
            # Validate ID column
            id_col = id_var.get()
            if id_col == '--Select--':
                messagebox.showwarning("ID Column Required", 
                    "Please select the metabolite/lipid ID column.")
                return
            
            # Build assignments
            self.stat_column_assignments = {
                'id_column': id_col,
                'comparisons': {}
            }
            
            ready_count = 0
            partial_count = 0
            missing_count = 0
            
            for (g1, g2), vars_dict in comparison_vars.items():
                pval = vars_dict['pvalue'].get()
                neglogp = vars_dict['neglogp'].get()
                log2fc = vars_dict['log2fc'].get()
                fc = vars_dict['fc'].get()
                
                pval = pval if pval != '--None--' else None
                neglogp = neglogp if neglogp != '--None--' else None
                log2fc = log2fc if log2fc != '--None--' else None
                fc = fc if fc != '--None--' else None

                # If FC missing but log2FC present, synthesize FC
                if not fc and log2fc and log2fc in complete_df.columns:
                    synth_fc_col = f"{g1}_vs_{g2}_FC_calc"
                    try:
                        if synth_fc_col not in complete_df.columns:
                            import numpy as _np
                            complete_df[synth_fc_col] = _np.power(2.0, complete_df[log2fc].astype(float))
                        fc = synth_fc_col
                        self.log_viz_message(f"   Synthesized FC from log2FC for {g1}_vs_{g2}: {synth_fc_col}")
                    except Exception:
                        pass

                has_pval = bool(pval)
                has_fc = bool(fc)

                if has_pval and has_fc:
                    status = 'ready'
                    ready_count += 1
                elif has_pval and not has_fc:
                    status = 'partial'
                    partial_count += 1
                elif has_fc and not has_pval:
                    status = 'partial'
                    partial_count += 1
                else:
                    status = 'missing'
                    missing_count += 1
                    self.log_viz_message(f"   ... missing pvalue/FC for {g1}_vs_{g2}")

                self.stat_column_assignments['comparisons'][(g1, g2)] = {
                    'pvalue': pval,
                    'neglogp': neglogp,
                    'log2fc': log2fc,
                    'fc': fc,
                    'status': status
                }
            
            self.log_viz_message(f"Stat columns configured: {ready_count} ready, {partial_count} partial, {missing_count} missing")
            
            # Mark as configured if at least one comparison is ready
            self.stat_cols_configured = ready_count > 0
            self._update_configuration_state()
            
            try:
                self._save_viz_config()
            except Exception:
                pass
            
            dialog.destroy()
        
        ttk.Button(btn_frame, text="🔄 Auto-Detect All", command=auto_detect_all).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Apply", command=apply_config).pack(side='right', padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side='right', padx=5)

    def auto_assign_viz_groups(self):
        """Auto-assign visualization groups based on common naming patterns."""
        if not hasattr(self, 'viz_sample_group_vars') or not self.viz_sample_group_vars:
            messagebox.showwarning('No Columns', 'No sample columns detected yet.')
            return

        # Keep auto-assign in sync with current labels typed in Group IDs & Labels.
        self._sync_viz_group_definitions_from_entries()
        
        pattern_window = tk.Toplevel(self.root)
        pattern_window.title('Auto-Assign Visualization Groups by Pattern')
        pattern_window.geometry('600x500')
        pattern_window.transient(self.root)
        pattern_window.grab_set()
        
        # Top frame for title and buttons
        top_frame = ttk.Frame(pattern_window)
        top_frame.pack(fill='x', padx=10, pady=10)
        
        ttk.Label(top_frame, text='Auto-Assignment Patterns', 
                 font=('Arial', 14, 'bold')).pack(side='left')
        
        # Buttons in top right
        btn_frame = ttk.Frame(top_frame)
        btn_frame.pack(side='right')
        
        # Create canvas with scrollbar for pattern definitions
        canvas = tk.Canvas(pattern_window, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(pattern_window, orient='vertical', command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            '<Configure>',
            lambda e: canvas.configure(scrollregion=canvas.bbox('all'))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pattern definition frame inside scrollable area
        pattern_frame = ttk.LabelFrame(scrollable_frame, text='Define Patterns', padding=10)
        pattern_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        ttk.Label(pattern_frame, text='Enter keywords/patterns for each group (one per line):').pack(anchor='w', pady=5)
        
        pattern_vars = {}
        for group_id, group_name in self.viz_group_definitions.items():
            group_frame = ttk.LabelFrame(pattern_frame, text=f'{group_id}: {group_name}', padding=5)
            group_frame.pack(fill='x', pady=5)
            
            pattern_text = tk.Text(group_frame, height=3, font=('Arial', 9))
            pattern_text.pack(fill='x', pady=5)
            pattern_vars[group_id] = pattern_text
            
            # Load saved patterns if available, otherwise use defaults
            saved_patterns = getattr(self, 'auto_assign_patterns', {}).get(group_id, '')
            if saved_patterns:
                # Use saved patterns
                pattern_text.insert(tk.END, saved_patterns)
            else:
                # Add default patterns based on common naming (only if no saved patterns)
                name_lower = group_name.lower()
                if any(k in name_lower for k in ['control', 'baseline', 'base', 'sham', 'healthy']):
                    pattern_text.insert(tk.END, 'CTRL\nControl\nCtrl')
                elif any(k in name_lower for k in ['disease', 'ortho', 'orthopedic', 'headache', 'patient']):
                    pattern_text.insert(tk.END, 'Ortho\nORTHO\nDisease')
                elif any(k in name_lower for k in ['tbi', 'injury', 'brain', 'trauma']):
                    pattern_text.insert(tk.END, 'TBI\nBrain\nInjury')
                elif 'treatment' in name_lower:
                    pattern_text.insert(tk.END, 'TRT\nTx\nTreatment')

        # Auto-discover prefixes from sample column names
        try:
            sample_cols = list(self.viz_sample_group_vars.keys())
            prefix_counts = {}
            for col_name in sample_cols:
                # Token = up to first underscore OR first digit block
                token = col_name.split('_')[0]
                if token and len(token) <= 30:
                    prefix_counts[token] = prefix_counts.get(token, 0) + 1
            
            # Sort by frequency (desc)
            sorted_tokens = sorted(prefix_counts.items(), key=lambda x: (-x[1], x[0]))

            ordered_group_items = list(self.viz_group_definitions.items())

            def _find_group_id_by_label_keywords(keywords):
                for gid, gname in ordered_group_items:
                    gl = str(gname).strip().lower()
                    if any(k in gl for k in keywords):
                        return gid
                return None

            first_gid = ordered_group_items[0][0] if ordered_group_items else None
            second_gid = ordered_group_items[1][0] if len(ordered_group_items) > 1 else first_gid
            third_gid = ordered_group_items[2][0] if len(ordered_group_items) > 2 else (ordered_group_items[-1][0] if ordered_group_items else None)

            control_gid = _find_group_id_by_label_keywords(['control', 'baseline', 'base', 'sham', 'healthy']) or first_gid
            disease_gid = _find_group_id_by_label_keywords(['disease', 'ortho', 'headache', 'patient']) or second_gid
            tbi_gid = _find_group_id_by_label_keywords(['tbi', 'injury', 'brain', 'trauma']) or second_gid
            treatment_gid = _find_group_id_by_label_keywords(['treatment', 'trt', 'drug', 'tx']) or third_gid
            
            # Map tokens heuristically to groups
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
                if token_upper.startswith(('CTR', 'CTL', 'CONT')) or token_upper == 'C':
                    target_group_id = control_gid
                elif token_upper.startswith(('BASE', 'BASLIN', 'BASAL')):
                    target_group_id = control_gid
                elif token_upper.startswith(('ORTH', 'ORT', 'OR', 'DIS')):
                    target_group_id = disease_gid
                elif token_upper.startswith(('HEAD', 'HA')):
                    target_group_id = disease_gid
                elif token_upper.startswith(('TBI', 'TB', 'BRAIN', 'INJ')):
                    target_group_id = tbi_gid
                elif token_upper in {'TRT','TX','TREAT','DRUG'}:
                    target_group_id = treatment_gid
                
                # If we found a target group, append token to its pattern box
                if target_group_id and target_group_id in pattern_vars:
                    tv = pattern_vars[target_group_id]
                    existing_text = tv.get('1.0', tk.END).strip()
                    if existing_text:
                        tv.insert(tk.END, f'\n{token}')
                    else:
                        tv.insert(tk.END, token)
            
            # Log what was auto-detected
            detected_tokens = ', '.join(f"{t}({c})" for t,c in sorted_tokens[:10])
            self.log_viz_message(f"Auto-detected sample name prefixes: {detected_tokens}")
            
        except Exception as e:
            self.log_viz_message(f"Warning: Prefix auto-detect failed: {e}")

        # Pack canvas and scrollbar
        canvas.pack(side='left', fill='both', expand=True, padx=(10, 0), pady=5)
        scrollbar.pack(side='right', fill='y', pady=5)
        
        def apply_patterns():
            self._sync_viz_group_definitions_from_entries()

            # Save patterns to config before applying
            if not hasattr(self, 'auto_assign_patterns'):
                self.auto_assign_patterns = {}
            for group_id, pattern_text_widget in pattern_vars.items():
                patterns = pattern_text_widget.get('1.0', tk.END).strip()
                self.auto_assign_patterns[group_id] = patterns
            
            assignments_made = 0
            
            # First, set all columns to 'Ignore' as default
            for col_name, group_var in self.viz_sample_group_vars.items():
                group_var.set('Ignore')
            
            # Then assign only those matching patterns
            for col_name, group_var in self.viz_sample_group_vars.items():
                matched = False
                
                # Collect all pattern candidates from all groups
                all_pattern_candidates = []
                for group_id, pattern_text in pattern_vars.items():
                    patterns = [p.strip() for p in pattern_text.get('1.0', tk.END).splitlines() if p.strip()]
                    for pattern in patterns:
                        if pattern.lower():
                            all_pattern_candidates.append((group_id, pattern, len(pattern)))
                
                # Sort by length descending (longest patterns first for more specific matching)
                # This ensures "HFD_TBI" matches before "HFD", preventing incorrect assignments
                all_pattern_candidates.sort(key=lambda x: -x[2])
                
                col_lower = col_name.lower()
                for group_id, pattern, _ in all_pattern_candidates:
                    if pattern.lower() in col_lower:
                        # Assign using group name
                        group_name = self.viz_group_definitions.get(group_id, group_id)
                        group_var.set(group_name)
                        assignments_made += 1
                        matched = True
                        break  # Stop after first match
                
                # If nothing matched, it stays as 'Ignore' (already set above)

            # Post-assignment diagnostics
            group_counts = {}
            for gv in self.viz_sample_group_vars.values():
                grp = gv.get()
                group_counts[grp] = group_counts.get(grp, 0) + 1
            
            summary = ', '.join(f"{gid}:{cnt}" for gid, cnt in group_counts.items())
            self.log_viz_message(f"Auto-assignment completed: {assignments_made} columns assigned")
            self.log_viz_message(f"Group distribution: {summary}")
            
            pattern_window.destroy()
        
        ttk.Button(btn_frame, text='Apply Patterns', command=apply_patterns).pack(side='right', padx=5)
        ttk.Button(btn_frame, text='Cancel', command=pattern_window.destroy).pack(side='right', padx=5)

    def delete_viz_sample_column(self, col_name, row_frame):
        """Delete a sample column from visualization group assignments."""
        if messagebox.askyesno("Confirm Delete", f"Remove sample column '{col_name}' from group assignments?"):
            try:
                # Remove from sample group variables
                if hasattr(self, 'viz_sample_group_vars') and col_name in self.viz_sample_group_vars:
                    del self.viz_sample_group_vars[col_name]
                
                # Remove from visualization group mapping
                if hasattr(self, 'viz_group_mapping') and col_name in self.viz_group_mapping:
                    del self.viz_group_mapping[col_name]
                
                # Remove the row from the UI
                row_frame.destroy()
                # Remove stored combobox reference if present
                try:
                    if hasattr(self, 'viz_sample_comboboxes') and col_name in self.viz_sample_comboboxes:
                        del self.viz_sample_comboboxes[col_name]
                except Exception:
                    pass

                self.log_viz_message(f"Removed sample column: {col_name}")
                try:
                    self._save_viz_config()
                except Exception:
                    pass
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to remove column: {str(e)}")
                self.log_viz_message(f"Error removing column {col_name}: {str(e)}")

    def on_viz_group_select(self, event):
        """Handle group selection - NO LONGER USED WITH NEW UI."""
        pass

    # --- Visualization persistence helpers ---
    def _viz_config_file(self):
        try:
            base_dir = os.path.dirname(os.path.abspath(__file__))
        except Exception:
            base_dir = os.getcwd()
        return resolve_runtime_config_path('visualization_config.json', {base_dir})

    def _manual_viz_config_file(self):
        """Path for manual save/load separate from auto-save."""
        try:
            base_dir = os.path.dirname(os.path.abspath(__file__))
        except Exception:
            base_dir = os.getcwd()
        return resolve_runtime_config_path('visualization_config_manual.json', {base_dir})

    def _gather_viz_config(self):
        import matplotlib.colors as mcolors
        
        def serialize_for_json(obj):
            """Recursively convert tuples to lists/strings for JSON serialization."""
            if isinstance(obj, dict):
                result = {}
                for key, value in obj.items():
                    # Convert tuple keys to strings
                    if isinstance(key, tuple):
                        key_str = "|".join(str(k) for k in key)
                        result[key_str] = serialize_for_json(value)
                    else:
                        result[str(key)] = serialize_for_json(value)
                return result
            elif isinstance(obj, (list, tuple)):
                return [serialize_for_json(item) for item in obj]
            else:
                return obj
        
        # Gather current sample-to-group assignments from UI
        current_mapping = {}
        if hasattr(self, 'viz_sample_group_vars'):
            for col, var in self.viz_sample_group_vars.items():
                current_mapping[col] = var.get()
        stat_assign = getattr(self, 'stat_column_assignments', None)
        
        # Serialize stat_column_assignments properly
        if stat_assign:
            stat_assign = serialize_for_json(stat_assign)
        
        id_col = stat_assign.get('id_column') if stat_assign and stat_assign.get('id_column') else None
        return {
            'viz_group_definitions': getattr(self, 'viz_group_definitions', {}),
            'viz_group_mapping': current_mapping,
            'viz_preferred_group_order': self.viz_preferred_group_order.get() if hasattr(self, 'viz_preferred_group_order') else '',
            'viz_color_map': {g: (mcolors.to_hex(c) if not isinstance(c, str) else c)
                              for g, c in getattr(self, 'viz_color_map', {}).items()},
            'stat_column_assignments': stat_assign,
            'id_column': id_col,
            'auto_assign_patterns': getattr(self, 'auto_assign_patterns', {}),
            'venn_pairs_list': getattr(self, 'venn_pairs_list', []),
            'venn_specs': getattr(self, 'venn_specs', []),
            'venn_allmol_pairs_list': getattr(self, 'venn_allmol_pairs_list', []),
            'venn_allmol_specs': getattr(self, 'venn_allmol_specs', []),
            'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
        }

    def reset_visualization_defaults(self):
        """Reset visualization-specific settings to sensible defaults."""
        try:
            # Reset viz params to defaults
            from main_script.metabolites_visualization import VolcanoParams, BoxplotParams, BargraphParams, PCAParams, HeatmapParams, ROCParams
            self.viz_params['volcano'] = VolcanoParams()
            self.viz_params['boxplot'] = BoxplotParams()
            self.viz_params['bargraph'] = BargraphParams()
            self.viz_params['pca'] = PCAParams()
            self.viz_params['heatmap'] = HeatmapParams()
            self.viz_params['roc'] = ROCParams()
            # Clear custom groups & mapping
            self.viz_group_definitions = {}
            self.viz_group_mapping = {}
            self.viz_color_map = {}
            # Clear Venn configurations (they depend on groups)
            self.venn_pairs_list = []
            self.venn_specs = []
            self.venn_allmol_pairs_list = []
            self.venn_allmol_specs = []
            # Clear display variables if they exist
            if hasattr(self, 'venn_pairs_var'):
                self.venn_pairs_var.set('')
            if hasattr(self, 'venn_specs_var'):
                self.venn_specs_var.set('')
            # Update UI widgets
            try:
                if hasattr(self, 'refresh_viz_group_display'):
                    self.refresh_viz_group_display()
            except Exception:
                pass
            self.log_viz_message('Visualization settings reset to defaults.')
        except Exception as e:
            self.log_viz_message(f'Failed to reset visualization defaults: {e}')

    def _apply_viz_config(self, data: dict):
        """Apply visualization config with detailed logging."""
        loaded_items = []
        gd = data.get('viz_group_definitions') or {}
        gm = data.get('viz_group_mapping') or {}
        cm = data.get('viz_color_map') or {}
        preferred_order = data.get('viz_preferred_group_order') or ''
        stat_assign = data.get('stat_column_assignments') or None
        
        # Convert string keys back to tuples for stat_column_assignments
        if stat_assign and isinstance(stat_assign, dict):
            stat_assign_restored = {}
            for key, value in stat_assign.items():
                if isinstance(key, str) and '|' in key and key != 'id_column':
                    # Convert string "g1|g2" back to tuple (g1, g2)
                    parts = key.split('|', 1)
                    if len(parts) == 2:
                        stat_assign_restored[(parts[0], parts[1])] = value
                    else:
                        stat_assign_restored[key] = value
                else:
                    stat_assign_restored[key] = value
            stat_assign = stat_assign_restored
        
        id_column = data.get('id_column') or (stat_assign.get('id_column') if stat_assign else None)

        if gd:
            self.viz_group_definitions = gd
            loaded_items.append(f"Group definitions: {len(gd)} groups ({', '.join(gd.values())})")
        if gm:
            self.viz_group_mapping = gm
            loaded_items.append(f"Sample-to-group mapping: {len(gm)} samples")
        if cm:
            self.viz_color_map = {g: cm[g] for g in cm}
            loaded_items.append(f"Color map: {', '.join([f'{g}={c}' for g, c in cm.items()])}")
        if preferred_order and hasattr(self, 'viz_preferred_group_order'):
            self.viz_preferred_group_order.set(preferred_order)
            loaded_items.append(f"Preferred order: {preferred_order}")
        if stat_assign:
            self.stat_column_assignments = stat_assign
            loaded_items.append(f"Stat assignments: {len(stat_assign.get('comparisons', {}))} comparisons")
        if id_column and self.stat_column_assignments:
            self.stat_column_assignments['id_column'] = id_column
            loaded_items.append(f"ID column: {id_column}")
        
        # Load auto-assign patterns
        aap = data.get('auto_assign_patterns') or {}
        if aap:
            self.auto_assign_patterns = aap
            loaded_items.append(f"Auto-assign patterns: {len(aap)} groups")
        
        # Get current valid groups for Venn validation
        current_groups = set()
        if gd:
            current_groups = set(gd.values())
        elif hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
            current_groups = set(self.viz_group_definitions.values())
        elif hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            current_groups = set(self.viz_group_mapping.values())
        
        def _validate_venn_pair(pair, valid_groups):
            """Check if both groups in a pair exist in current dataset."""
            if isinstance(pair, (list, tuple)) and len(pair) >= 2:
                return pair[0] in valid_groups and pair[1] in valid_groups
            return False
        
        def _validate_venn_spec(spec, valid_groups):
            """Check if all groups in a spec exist in current dataset."""
            if isinstance(spec, (list, tuple)):
                return all(g in valid_groups for g in spec)
            return False
        
        # Load Venn configuration - validate against current groups
        venn_pairs = data.get('venn_pairs_list') or []
        if venn_pairs:
            raw_pairs = [tuple(p) if isinstance(p, list) else p for p in venn_pairs]
            if current_groups:
                # Filter to only pairs that match current groups
                valid_pairs = [p for p in raw_pairs if _validate_venn_pair(p, current_groups)]
                invalid_count = len(raw_pairs) - len(valid_pairs)
                if invalid_count > 0:
                    logger.info(f"Skipped {invalid_count} Venn pairs not matching current groups: {current_groups}")
                self.venn_pairs_list = valid_pairs
                if valid_pairs:
                    loaded_items.append(f"Filtered Venn pairs: {len(valid_pairs)}")
            else:
                # No groups defined yet, clear saved pairs
                self.venn_pairs_list = []
                logger.info("Cleared saved Venn pairs (no groups defined yet)")
        
        venn_specs = data.get('venn_specs') or []
        if venn_specs:
            if current_groups:
                valid_specs = [s for s in venn_specs if _validate_venn_spec(s, current_groups)]
                invalid_count = len(venn_specs) - len(valid_specs)
                if invalid_count > 0:
                    logger.info(f"Skipped {invalid_count} Venn specs not matching current groups")
                self.venn_specs = valid_specs
                if valid_specs:
                    loaded_items.append(f"Filtered Venn specs: {len(valid_specs)}")
            else:
                self.venn_specs = []
        
        # Load All Molecules Venn configuration - same validation
        allmol_pairs = data.get('venn_allmol_pairs_list') or []
        if allmol_pairs:
            raw_pairs = [tuple(p) if isinstance(p, list) else p for p in allmol_pairs]
            if current_groups:
                valid_pairs = [p for p in raw_pairs if _validate_venn_pair(p, current_groups)]
                self.venn_allmol_pairs_list = valid_pairs
                if valid_pairs:
                    loaded_items.append(f"All Molecules Venn pairs: {len(valid_pairs)}")
            else:
                self.venn_allmol_pairs_list = []
        
        allmol_specs = data.get('venn_allmol_specs') or []
        if allmol_specs:
            if current_groups:
                valid_specs = [s for s in allmol_specs if _validate_venn_spec(s, current_groups)]
                self.venn_allmol_specs = valid_specs
                if valid_specs:
                    loaded_items.append(f"All Molecules Venn specs: {len(valid_specs)}")
            else:
                self.venn_allmol_specs = []

        if loaded_items:
            logger.info("Loaded visualization config:")
            for item in loaded_items:
                logger.info(f"  ✓ {item}")

        try:
            if hasattr(self, 'viz_group_listbox'):
                self.refresh_viz_group_display()
            if hasattr(self, 'viz_sample_comboboxes') and self.viz_sample_comboboxes:
                names = list(self.viz_group_definitions.values())
                for col, combo in self.viz_sample_comboboxes.items():
                    try:
                        combo['values'] = names
                        if col in self.viz_group_mapping and self.viz_group_mapping[col] in names:
                            self.viz_sample_group_vars[col].set(self.viz_group_mapping[col])
                    except Exception:
                        pass
            if hasattr(self, 'viz_color_frame') and self.viz_color_frame:
                groups = list(self.viz_group_definitions.values())
                if groups:
                    self.update_group_color_controls(groups)
            # Restore ID column selection widget if present
            if id_column and hasattr(self, 'viz_id_column_var'):
                try:
                    self.viz_id_column_var.set(id_column)
                except Exception:
                    pass
        except Exception as e:
            if hasattr(self, 'log_viz_message'):
                self.log_viz_message(f"Warning during config apply: {e}")

        # Recompute workflow state
        try:
            self._update_configuration_state()
        except Exception:
            pass

    def save_viz_config_manual(self):
        """Manually save current visualization configuration to distinct file."""
        path = self._manual_viz_config_file()
        try:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(self._gather_viz_config(), f, indent=2)
            self._remove_readonly_attribute(path)
            self.log_viz_message(f"💾 Manual config saved: {os.path.basename(path)}")
        except Exception as e:
            self.log_viz_message(f"Manual save failed: {e}")

    def load_viz_config_manual(self):
        """Load visualization configuration from manual file."""
        path = self._manual_viz_config_file()
        if not os.path.exists(path):
            messagebox.showinfo("Load Config", "No manual visualization config found yet.")
            return
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self._apply_viz_config(data)
            self.log_viz_message(f"📥 Manual config loaded: {os.path.basename(path)}")
        except Exception as e:
            self.log_viz_message(f"Manual load failed: {e}")

    def _save_viz_config(self):
        path = self._viz_config_file()
        try:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(self._gather_viz_config(), f, indent=2)
            # Remove read-only attribute
            self._remove_readonly_attribute(path)
        except Exception as e:
            self.log_viz_message(f"Failed to save visualization config: {e}")

    def _load_viz_config(self):
        path = self._viz_config_file()
        if os.path.exists(path):
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self._apply_viz_config(data)
                self.log_viz_message('Loaded previous visualization configuration.')
            except Exception as e:
                self.log_viz_message(f'Failed to load visualization config: {e}')

    def refresh_viz_group_display(self):
        """Refresh the group display and combobox options."""
        self._sync_viz_group_definitions_from_entries()

        # Ensure viz_group_definitions exists
        if not hasattr(self, 'viz_group_definitions') or not self.viz_group_definitions:
            if not hasattr(self, 'viz_group_definitions'):
                self.viz_group_definitions = {}
            return
        
        # Refresh the group UI (scrollable list with "Group 1:", "Group 2:", etc.)
        if hasattr(self, 'groups_scrollable_frame'):
            self.refresh_group_ui()
        
        # Update combobox values for sample assignments
        group_names = list(self.viz_group_definitions.values())
        if hasattr(self, 'viz_sample_group_vars'):
            try:
                if hasattr(self, 'viz_sample_comboboxes') and self.viz_sample_comboboxes:
                    for col, combo in self.viz_sample_comboboxes.items():
                        try:
                            combo['values'] = ['Ignore'] + group_names
                            # Keep current value if it's still valid
                            cur = self.viz_sample_group_vars.get(col).get() if col in self.viz_sample_group_vars else None
                            if cur not in (['Ignore'] + group_names):
                                if col in self.viz_sample_group_vars:
                                    self.viz_sample_group_vars[col].set('Ignore')
                        except Exception:
                            try:
                                if col in self.viz_sample_group_vars:
                                    self.viz_sample_group_vars[col].set('Ignore')
                            except Exception:
                                pass
            except Exception as e:
                self.log_viz_message(f"Warning: Could not update all comboboxes: {e}")

    def _sync_viz_group_definitions_from_entries(self):
        """Sync currently typed group labels from dialog entries into definitions."""
        try:
            if not hasattr(self, 'viz_group_definitions') or not self.viz_group_definitions:
                return
            if not hasattr(self, 'group_id_vars') or not self.group_id_vars:
                return

            updated = False
            for gid in list(self.viz_group_definitions.keys()):
                if gid not in self.group_id_vars:
                    continue
                raw_val = self.group_id_vars[gid].get()
                new_label = str(raw_val).strip() if raw_val is not None else ''
                if not new_label:
                    continue
                if self.viz_group_definitions.get(gid) != new_label:
                    self.viz_group_definitions[gid] = new_label
                    updated = True

            if updated and hasattr(self, 'viz_sample_comboboxes') and hasattr(self, 'viz_sample_group_vars'):
                labels = list(self.viz_group_definitions.values())
                valid = ['Ignore'] + labels
                for col, combo in self.viz_sample_comboboxes.items():
                    try:
                        combo['values'] = valid
                        if col in self.viz_sample_group_vars:
                            cur = self.viz_sample_group_vars[col].get()
                            if cur not in valid:
                                self.viz_sample_group_vars[col].set('Ignore')
                    except Exception:
                        pass
        except Exception:
            pass

    def _populate_pca_custom_group_checkboxes(self):
        """No-op method for backward compatibility - PCA group checkboxes removed."""
        pass

    def get_viz_groups(self):
        """Get the list of visualization groups.
        
        Priority:
        1. Groups from viz_group_mapping (if samples are assigned to groups)
        2. Groups extracted from actual data column names (e.g., "Control_vs_PD")
        3. Groups from viz_group_definitions (may be from saved config)
        """
        # Priority 1: Get from sample-to-group mapping (most reliable)
        if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            return self.ordered_groups(list(dict.fromkeys(self.viz_group_mapping.values())))
        
        # Priority 2: Extract groups from actual data column names 
        # This handles the case where old config has wrong groups
        groups_from_data = self._get_groups_from_data_columns()
        if groups_from_data:
            return self.ordered_groups(groups_from_data)
        
        # Priority 3: Fallback to viz_group_definitions (from config or UI)
        if hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
            return list(self.viz_group_definitions.values())
        return []
    
    def _get_groups_from_data_columns(self):
        """Extract group names from data column names (e.g., 'Control_vs_PD' -> ['Control', 'PD'])."""
        try:
            df = self._build_complete_df()
            if df is None or df.empty:
                return []
            
            groups = set()
            for col in df.columns:
                # Look for comparison columns (e.g., Control_vs_PD_log2FC)
                if '_vs_' in col:
                    # Extract the comparison part (before any suffix like _log2FC, _pvalue)
                    parts = col.split('_vs_')
                    if len(parts) >= 2:
                        group1 = parts[0].strip()
                        # Group2 might have suffixes, extract just the group name
                        group2_full = parts[1]
                        # Common suffixes to strip
                        for suffix in ['_log2FC', '_pvalue', '_p_adj', '_FC', '_foldchange', '_mean', '_std']:
                            if group2_full.endswith(suffix):
                                group2_full = group2_full[:-len(suffix)]
                                break
                        # Also handle cases like "PD_log2FC" -> "PD"
                        group2 = group2_full.split('_')[0] if '_' in group2_full else group2_full
                        
                        if group1:
                            groups.add(group1)
                        if group2:
                            groups.add(group2)
            
            return list(groups) if groups else []
        except Exception as e:
            logger.debug(f"Could not extract groups from data columns: {e}")
            return []

    def ordered_groups(self, groups: list[str]) -> list[str]:
        """Return groups in a stable, user-respecting order.

        Preference order:
        1. If preferred group order is set, honor it first.
        2. Else if viz_group_definitions exists, follow its value order.
        3. Else if viz_group_mapping exists, follow order of first appearance in mapping values.
        4. Else fall back to the input list order.
        """
        if not groups:
            return []
        # If explicit preferred order exists, use it first.
        try:
            if hasattr(self, 'viz_preferred_group_order'):
                order_str = self.viz_preferred_group_order.get().strip()
                if order_str:
                    preferred = [g.strip() for g in order_str.split(',') if g.strip()]
                    ordered = [g for g in preferred if g in groups]
                    ordered += [g for g in groups if g not in ordered]
                    return ordered
        except Exception:
            pass
        # If user-defined group definitions present, follow that ordering
        if hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
            defs = list(self.viz_group_definitions.values())
            ordered = [g for g in defs if g in groups]
            # Append any remaining groups not in definitions in original order
            ordered += [g for g in groups if g not in ordered]
            return ordered
        # If mapping exists, preserve the first-seen order from mapping
        if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            seen = []
            for sample, grp in self.viz_group_mapping.items():
                if grp in groups and grp not in seen:
                    seen.append(grp)
            seen += [g for g in groups if g not in seen]
            return seen
        # Default: return as provided
        return list(groups)

    def get_current_viz_groups(self):
        """Get current visualization groups with consistent ordering."""
        groups = self.get_viz_groups()
        if not groups and hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
            groups = list(self.viz_group_definitions.values())
        return self.ordered_groups(groups) if groups else []

    def get_viz_sample_to_group_mapping(self):
        """Get the sample to group mapping for visualization.
        
        Filters out samples assigned to 'Ignore' as they should be excluded from analysis.
        """
        if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            # Filter out samples assigned to 'Ignore'
            filtered_mapping = {
                col: group 
                for col, group in self.viz_group_mapping.items() 
                if group != 'Ignore'
            }
            return filtered_mapping
        return {}

    def import_stats_groups(self):
        """Import group assignments from the statistics tab."""
        if hasattr(self, 'sample_group_vars') and self.sample_group_vars:
            # Parse existing group assignments from statistics tab
            self.viz_group_mapping = self._parse_group_assignments()
            groups = self.ordered_groups(list(dict.fromkeys(self.viz_group_mapping.values())))
            
            # Update status
            self.viz_group_status.set(f"Imported from Statistics: {', '.join(groups)} ({len(self.viz_group_mapping)} samples)")
            self.update_group_color_controls(groups)
            self.log_viz_message(f"Imported group assignments from Statistics tab: {len(groups)} groups")
        else:
            messagebox.showwarning("No Statistics Groups", 
                                 "No group assignments found in Statistics tab. Please configure groups there first.")

    def update_group_color_controls(self, groups):
        """Update the group color selection controls."""
        # Initialize color map if it doesn't exist
        if not hasattr(self, 'viz_color_map'):
            self.viz_color_map = {}
        
        # Preserve existing color mappings and clean up obsolete ones
        old_color_map = self.viz_color_map.copy()
        self.viz_color_map = {}
        
        # Default colors for consistent assignment
        default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
        # Specific requested defaults (case-insensitive matching)
        requested_defaults = {
            'du145_2d': '#1f77b4',
            'du145_3d': '#ff7f0e',
            'pc3_2d': '#2ca02c',
            'pc3_3d': '#d62728'
        }
        
        # Assign colors to groups, preserving existing ones
        for i, group in enumerate(self.ordered_groups(groups)):
            if group in old_color_map:
                # Preserve existing color
                self.viz_color_map[group] = old_color_map[group]
            else:
                # If the group name matches a requested default, prefer that
                group_key = group.lower() if isinstance(group, str) else group
                if group_key in requested_defaults:
                    self.viz_color_map[group] = requested_defaults[group_key]
                else:
                    # Assign new color from palette
                    self.viz_color_map[group] = default_colors[i % len(default_colors)]
        
        # Clear existing controls
        for widget in self.viz_color_frame.winfo_children():
            widget.destroy()
        
        # Create horizontal layout for color controls
        color_row = ttk.Frame(self.viz_color_frame)
        color_row.pack(fill='x', pady=2)
        
        # Create color controls for each group horizontally
        for group in self.ordered_groups(groups):
            group_frame = ttk.Frame(color_row)
            group_frame.pack(side='left', padx=2)
            
            ttk.Label(group_frame, text=f"{group}:", font=('Arial', 8)).pack()
            
            color_button = tk.Button(group_frame, text="●", font=('Arial', 12), 
                                   bg=self.viz_color_map[group], width=2, height=1,
                                   command=lambda g=group: self.choose_group_color(g))
            color_button.pack(pady=1)
        
        # Log the color assignments for debugging
        color_assignments = ', '.join([f"{g}:{self.viz_color_map[g]}" for g in self.ordered_groups(groups)])
        self.log_viz_message(f"Color controls updated: {color_assignments}")

    def choose_group_color(self, group_name):
        """Open color chooser with both visual picker and hex input option."""
        import tkinter.colorchooser as colorchooser
        import matplotlib.colors as mcolors
        from tkinter import simpledialog
        
        current_color = self.viz_color_map.get(group_name, '#000000')
        
        # Create dialog
        method_dialog = tk.Toplevel(self.root)
        method_dialog.title(f"Choose Color for {group_name}")
        method_dialog.geometry("400x350")
        method_dialog.resizable(False, False)
        method_dialog.transient(self.root)
        method_dialog.grab_set()
        
        result = {'color': current_color}
        
        main_frame = ttk.Frame(method_dialog, padding=20)
        main_frame.pack(fill='both', expand=True)
        
        ttk.Label(main_frame, text=f"Select color for: {group_name}", font=('Arial', 11, 'bold')).pack(pady=(0, 15))
        
        # Current color preview
        preview_frame = ttk.Frame(main_frame)
        preview_frame.pack(pady=(0, 15))
        ttk.Label(preview_frame, text="Current:", font=('Arial', 9)).pack(side='left', padx=(0, 5))
        preview_canvas = tk.Canvas(preview_frame, width=80, height=35, bg=current_color, 
                                   highlightthickness=1, highlightbackground='black')
        preview_canvas.pack(side='left')
        
        # Visual Color Picker Button
        picker_frame = ttk.Frame(main_frame)
        picker_frame.pack(pady=(0, 15), fill='x')
        
        def use_visual_picker():
            color = colorchooser.askcolor(
                title=f"Pick Color for {group_name}",
                initialcolor=result['color'],
                parent=method_dialog
            )
            if color and color[1]:
                result['color'] = color[1]
                preview_canvas.configure(bg=color[1])
                hex_var.set(color[1])
        
        ttk.Button(picker_frame, text="🎨 Visual Color Picker", 
                  command=use_visual_picker, width=25).pack()
        
        # Hex input section with textbox inline
        hex_frame = ttk.LabelFrame(main_frame, text="Or Enter Hex Code", padding=10)
        hex_frame.pack(fill='x', pady=(0, 15))
        
        hex_input_row = ttk.Frame(hex_frame)
        hex_input_row.pack(fill='x')
        
        ttk.Label(hex_input_row, text="Hex:", font=('Arial', 10)).pack(side='left', padx=(0, 5))
        hex_var = tk.StringVar(value=current_color)
        hex_entry = ttk.Entry(hex_input_row, textvariable=hex_var, width=15, font=('Arial', 10))
        hex_entry.pack(side='left', padx=(0, 5))
        
        def preview_hex():
            try:
                hex_val = hex_var.get().strip()
                if not hex_val.startswith('#'):
                    hex_val = '#' + hex_val
                # Validate hex color
                mcolors.to_hex(hex_val)
                result['color'] = hex_val
                preview_canvas.configure(bg=hex_val)
            except (ValueError, AttributeError):
                messagebox.showwarning("Invalid Color", 
                                     f"'{hex_var.get()}' is not a valid hex color code.\n"
                                     "Please use format: #RRGGBB (e.g., #00BA38)",
                                     parent=method_dialog)
        
        ttk.Button(hex_input_row, text="Preview", command=preview_hex).pack(side='left')
        
        ttk.Label(hex_frame, text="Example: #00BA38, #FF5733, #1F77B4", 
                 font=('Arial', 8), foreground='gray').pack(anchor='w', pady=(5, 0))
        
        # OK and Cancel buttons side by side
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill='x')
        
        def apply_color():
            preview_hex()  # Validate current hex input
            method_dialog.destroy()
        
        def cancel():
            result['color'] = None
            method_dialog.destroy()
        
        ttk.Button(btn_frame, text="OK", command=apply_color, width=15).pack(side='left', padx=(0, 5))
        ttk.Button(btn_frame, text="Cancel", command=cancel, width=15).pack(side='left')
        
        # Bind Enter key to apply
        hex_entry.bind('<Return>', lambda e: apply_color())
        
        # Wait for dialog to close
        self.root.wait_window(method_dialog)
        
        # Apply the color if selected
        if result['color']:
            self.viz_color_map[group_name] = result['color']
            self.update_single_color_button(group_name, result['color'])
            self.log_viz_message(f"Updated color for {group_name} to {result['color']}")
            try:
                self._save_viz_config()
            except Exception as e:
                self.log_viz_message(f"Warning: Could not auto-save color settings: {e}")

    def update_single_color_button(self, group_name, color):
        """Update a single color button without recreating all controls."""
        try:
            # Find the color button for this group and update its color
            for widget in self.viz_color_frame.winfo_children():
                if isinstance(widget, ttk.Frame):
                    for child in widget.winfo_children():
                        if isinstance(child, ttk.Frame):
                            label_widget = None
                            button_widget = None
                            for grandchild in child.winfo_children():
                                if isinstance(grandchild, ttk.Label):
                                    label_widget = grandchild
                                elif isinstance(grandchild, tk.Button):
                                    button_widget = grandchild
                            
                            # Check if this is the right group
                            if label_widget and button_widget and label_widget.cget('text') == f"{group_name}:":
                                button_widget.config(bg=color)
                                return
        except Exception as e:
            # Fallback to full recreation if individual update fails
            self.log_viz_message(f"Warning: Could not update individual button, recreating all: {e}")
            current_groups = []
            if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
                current_groups = self.ordered_groups(list(dict.fromkeys(self.viz_group_mapping.values())))
            elif hasattr(self, 'viz_group_definitions') and self.viz_group_definitions:
                current_groups = list(self.viz_group_definitions.values())
            
            if current_groups:
                self.update_group_color_controls(current_groups)

    def log_viz_message(self, message, tag='INFO'):
        """Add a message to visualization log with timestamp."""
        timestamp = time.strftime('%H:%M:%S')
        tagged_message = f"[{timestamp}] {tag}: {message}\n"
        if hasattr(self, 'viz_log'):
            self.viz_log.insert(tk.END, tagged_message)
            self.viz_log.see(tk.END)
            self.root.update_idletasks()
        else:
            print(tagged_message.strip())
    
    @staticmethod
    def _remove_readonly_attribute(filepath):
        """Remove read-only attribute from a file (Windows-specific fix)"""
        try:
            import os
            import stat
            os.chmod(filepath, stat.S_IWRITE)
        except Exception as e:
            logger.debug(f"Could not remove read-only attribute from {filepath}: {e}")
    
    def _normalize_col(self, col):
        """Normalize a column name for comparison (lowercase, no spaces/underscores)."""
        if col is None:
            return ''
        return str(col).lower().replace('_', '').replace(' ', '')

    @staticmethod
    def _normalize_feature_id_value(value):
        """Normalize feature identifiers so integer-like numeric values stay text."""
        if value is None:
            return value
        try:
            if pd.isna(value):
                return value
        except Exception:
            pass

        if isinstance(value, (int, np.integer)):
            return str(int(value))

        if isinstance(value, (float, np.floating)):
            if float(value).is_integer():
                return str(int(value))
            return str(value)

        text = str(value).strip()
        if not text:
            return text

        if re.fullmatch(r'[-+]?\d+\.0+', text):
            try:
                return str(int(float(text)))
            except Exception:
                return text

        return text

    def _normalize_feature_id_column(self, df, feature_id_col):
        """Normalize a feature ID column in-place when possible."""
        if df is None or not feature_id_col or feature_id_col not in df.columns:
            return False

        try:
            df[feature_id_col] = df[feature_id_col].map(self._normalize_feature_id_value)
            return True
        except Exception as exc:
            logger.debug(f"Could not normalize feature ID column '{feature_id_col}': {exc}")
            return False
    
    def _is_lipid_feature_col(self, col):
        """Check if a column is a canonical lipid feature column."""
        if col is None:
            return False
        
        # List of canonical lipid feature columns (normalized)
        # These are metadata/annotation columns, NOT sample intensity columns
        lipid_features = [
            'lipidid', 'lipid_id', 'lipidname', 'lipid_name',
            'class', 'subclass', 'lipid', 'lipidgroup', 'lipidclass',
            'charge', 'calcmz', 'basert', 'obsrt', 'ppmdiff',
            'mz', 'retention', 'rt', 'mass'
        ]
        
        normalized = self._normalize_col(col)
        return any(feature in normalized for feature in lipid_features)
    
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
    
    def _get_min_group_size_percent(self) -> float:
        """Return minimum percentage threshold."""
        try:
            value = float(self.min_samples_percent_var.get()) if hasattr(self, 'min_samples_percent_var') else 50.0
        except Exception:
            value = 50.0
        return max(0.0, min(100.0, value))

    def _filter_groups_by_min_samples(self, sample_to_group: dict[str, str], min_required: int):
        """Filter groups based on a minimum samples-per-group threshold.

        Returns a tuple: (filtered_map, included_counts, excluded_counts)
        - filtered_map: dict of sample -> group including only groups that meet threshold
        - included_counts: dict of group -> count for included groups
        - excluded_counts: dict of group -> count for excluded groups
        """
        try:
            from collections import Counter
            import math

            if not sample_to_group:
                return {}, {}, {}

            # Count current samples per group
            counts = Counter(sample_to_group.values())

            # Determine thresholding mode
            mode = self._get_min_group_size_type() if hasattr(self, '_get_min_group_size_type') else 'absolute'
            threshold = int(min_required)

            if mode == 'percentage':
                pct = self._get_min_group_size_percent() if hasattr(self, '_get_min_group_size_percent') else 50.0
                # Use the largest group size as reference for percentage threshold
                max_n = max(counts.values()) if counts else 0
                threshold = max(1, math.ceil((pct / 100.0) * max_n))

            # Partition groups by threshold
            included = {g for g, n in counts.items() if n >= threshold}
            excluded = {g for g, n in counts.items() if n < threshold}

            filtered_map = {s: g for s, g in sample_to_group.items() if g in included}
            included_counts = {g: counts[g] for g in included}
            excluded_counts = {g: counts[g] for g in excluded}

            return filtered_map, included_counts, excluded_counts
        except Exception:
            # On any error, do not filter
            return dict(sample_to_group), {}, {}

    def browse_viz_import_file(self):
        """Browse for statistical results Excel file."""
        # Only allow browsing when file mode is selected
        if self.viz_data_source.get() != 'file':
            messagebox.showwarning("File Mode Required", 
                                 "Please select 'Import statistical results Excel file' first.")
            return
            
        filename = filedialog.askopenfilename(
            title="Select Statistical Results Excel File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            self.viz_import_file.set(filename)

    def update_viz_import_state(self):
        """Enable/disable file import controls based on data source selection."""
        def _safe_config(widget_name, **kwargs):
            """Configure widget only if it exists and is still valid."""
            widget = getattr(self, widget_name, None)
            if widget is None:
                return
            try:
                if widget.winfo_exists():
                    widget.config(**kwargs)
            except Exception:
                pass

        if self.viz_data_source.get() == 'file':
            # Enable file import controls
            _safe_config('viz_file_entry', state='normal')
            _safe_config('viz_browse_btn', state='normal')
            # File mode: user must verify columns before Configure Groups
            _safe_config('verify_cols_btn', state='disabled')  # Enable after Load & Analyze
            _safe_config('configure_groups_btn', state='disabled')  # Enable after Verify Columns
        else:
            # Disable file import controls
            _safe_config('viz_file_entry', state='disabled')
            _safe_config('viz_browse_btn', state='disabled')
            # Clear the file path when switching to memory mode
            if hasattr(self, 'viz_import_file'):
                self.viz_import_file.set("")
            # Memory mode: enable verify columns to allow re-verification if needed
            _safe_config('verify_cols_btn', state='normal')  # Allow re-verification
            _safe_config('configure_groups_btn', state='normal')  # Allow direct config

    def load_viz_data(self):
        """Load and analyze data from Excel file or current session."""
        try:
            if self.viz_data_source.get() == 'file':
                # File mode - load Excel file
                if not self.viz_import_file.get():
                    messagebox.showwarning("No File", "Please select an Excel file first.")
                    return
                
                # Clear old session data when loading from file to prevent contamination
                if hasattr(self, 'statistical_test_results'):
                    self.statistical_test_results = None
                if hasattr(self, 'statistical_test_results_class'):
                    self.statistical_test_results_class = None
                
                # Load Excel file
                self.log_viz_message("Loading Excel file...")
                excel_path = self.viz_import_file.get()
                
                # Auto-load first sheet from Excel file (no hardcoded sheet name)
                try:
                    xls = pd.ExcelFile(excel_path)
                    if not xls.sheet_names:
                        messagebox.showerror("Load Error", "Excel file contains no sheets")
                        return
                    
                    first_sheet = xls.sheet_names[0]
                    df = pd.read_excel(excel_path, sheet_name=first_sheet)
                    self.imported_complete_df = df
                    self.log_viz_message(f"Loaded sheet '{first_sheet}': {len(df)} rows, {len(df.columns)} columns")
                    
                    # Enable Verify Columns button for file imports
                    self.verify_cols_btn.config(state='normal')
                    self.log_viz_message("\n📋 Next: Click '🔍 Verify Columns' to define sample columns.")
                    
                except Exception as e:
                    messagebox.showerror("Load Error", f"Could not load Excel file:\n{str(e)}")
                    return
                
                # If in lipid mode, optionally load class Complete Results
                try:
                    mode = self.viz_mode.get() if hasattr(self, 'viz_mode') else 'metabolite'
                except Exception:
                    mode = 'metabolite'
                
                # Check if lipid mode
                is_lipid = (mode == 'lipid')
                
                if is_lipid:
                    # Check if additional class-level sheet exists (flexible naming)
                    try:
                        # Already have xls from above
                        if len(xls.sheet_names) > 1:
                            class_sheet_candidates = [s for s in xls.sheet_names if 'class' in s.lower() and s != first_sheet]
                            if class_sheet_candidates:
                                self.log_viz_message(f"ℹ️  File contains class-level sheet: '{class_sheet_candidates[0]}'")
                                self.log_viz_message("   Use 'Import Lipid Class Data' button to load class-level data for visualizations")
                    except Exception as e:
                        # Silently ignore if sheet doesn't exist or can't be loaded
                        pass
                        
            else:
                # Memory mode - use session data from Statistics tab
                # Clear imported file data when switching back to memory mode
                if hasattr(self, 'imported_complete_df'):
                    self.imported_complete_df = None
                
                self.log_viz_message("Using data from current session...")
                
                # In memory mode, columns are already verified from Statistics tab
                # Enable Configure Groups directly and keep Verify Columns available
                # so users can re-verify if needed.
                self.configure_groups_btn.config(state='normal')
                self.verify_cols_btn.config(state='normal')
            
            # Update status and controls for both modes
            self.update_viz_data_status()
            self.update_group_color_controls(self.get_current_groups())
            
        except Exception as e:
            messagebox.showerror("Load Error", f"Failed to load data:\n{str(e)}")
            self.log_viz_message(f"Load failed: {str(e)}", 'ERROR')

    def update_viz_data_status(self):
        """Update visualization status by logging to the right panel log only."""
        try:
            # Get current data
            complete_df = self.get_current_complete_df()
            groups = self.get_current_groups()
            sample_cols, sample_to_group = self.get_sample_mapping()
            
            if complete_df is None:
                self.log_viz_message("❌ No statistical data available. Load data from Statistics tab or import Excel file", 'INFO')
                return
            
            # DEBUG: Detailed column classification to diagnose missing pairwise / padj columns
            all_cols = list(complete_df.columns)
            pairwise_cols = [c for c in all_cols if '_vs_' in c]
            padj_cols = [c for c in all_cols if c.endswith('_p_adj') or c.endswith('_adj')]
            overall_cols = [c for c in all_cols if c.startswith('overall_')]
            sample_like = [c for c in all_cols if c in sample_cols]
            # Log high-level summary first
            self.log_viz_message(f"✅ Data loaded: {len(complete_df)} metabolites, {len(all_cols)} columns", 'INFO')
            self.log_viz_message(f"🔎 DEBUG Columns: samples={len(sample_like)} pairwise={len(pairwise_cols)} padj={len(padj_cols)} overall={len(overall_cols)}", 'INFO')
            # If expected pairwise columns appear absent, hint probable cause
            if len(pairwise_cols) == 0:
                self.log_viz_message("🛠 DEBUG: No '_vs_' columns detected. Possible causes: (1) Two-way ANOVA not run, (2) group assignments missing during stats, (3) export used stripped frame.", 'WARNING')
            if len(padj_cols) == 0:
                self.log_viz_message("🛠 DEBUG: No adjusted p-value columns present. If FDR was enabled, verify 'use_adj_p_var' and sufficient non-NaN p-values.", 'WARNING')
            # Log a trimmed preview of first few pairwise/padj columns for confirmation
            if pairwise_cols:
                self.log_viz_message(f"🧪 DEBUG Pairwise head: {', '.join(pairwise_cols[:6])}", 'INFO')
            if padj_cols:
                self.log_viz_message(f"🧪 DEBUG Padj head: {', '.join(padj_cols[:6])}", 'INFO')
            self.log_viz_message(f"📊 Groups ({len(groups)}): {', '.join(groups)}", 'INFO')
            
            # Group sample counts
            if sample_to_group:
                group_counts = {}
                for sample, group in sample_to_group.items():
                    group_counts[group] = group_counts.get(group, 0) + 1
                
                for group, count in group_counts.items():
                    self.log_viz_message(f"   • {group}: {count} samples", 'INFO')
                
                # AUTO-POPULATE VISUALIZATION GROUPS FROM LOADED DATA
                if groups:
                    # Create group mapping for visualization
                    self.viz_group_mapping = sample_to_group.copy()
                    
                    # Update group definitions from loaded groups
                    self.viz_group_definitions = {}
                    for i, group_name in enumerate(self.ordered_groups(groups), 1):
                        group_id = f"Group{i}"
                        self.viz_group_definitions[group_id] = group_name
                    
                    # Update visualization status and colors
                    self.viz_group_status.set(f"Auto-loaded: {', '.join(groups)} ({len(sample_to_group)} samples)")
                    self.update_group_color_controls(groups)
                    self.log_viz_message(f"Auto-populated {len(groups)} groups from loaded data: {', '.join(groups)}", 'INFO')
                    
                    # Update PCA custom group checkboxes
                    self._populate_pca_custom_group_checkboxes()
                    
                    # Mark groups as configured and reset stat columns for new data
                    self.groups_configured = True
                    self.stat_cols_configured = False
                    self.stat_column_assignments = {}
                    
                    # Auto-configure stat columns with detection
                    self._auto_configure_stat_columns(complete_df, groups)
                    self._update_configuration_state()
                    
            else:
                self.log_viz_message("⚠️ No sample-to-group mapping available", 'WARNING')
        
        except Exception as e:
            self.log_viz_message(f"❌ Error analyzing data: {str(e)}", 'ERROR')

    def get_current_complete_df(self):
        """Get the current complete DataFrame from session or imported file."""
        # If user chose file, prefer imported data
        if self.viz_data_source.get() == 'file' and hasattr(self, 'imported_complete_df'):
            return self.imported_complete_df

        # Prefer lipid class results when in lipid mode and available
        try:
            mode = self.viz_mode.get() if hasattr(self, 'viz_mode') else 'metabolite'
        except Exception:
            mode = 'metabolite'
        if mode == 'lipid':
            # If class statistical results exist, use their Complete Results
            if hasattr(self, 'statistical_test_results_class') and self.statistical_test_results_class and \
               'enhanced_metabolites' in self.statistical_test_results_class and \
               isinstance(self.statistical_test_results_class['enhanced_metabolites'], pd.DataFrame) and \
               not self.statistical_test_results_class['enhanced_metabolites'].empty:
                return self.statistical_test_results_class['enhanced_metabolites']
            # No automatic prompt - user can manually import via the "Import Lipid Class Data" button
        # Fallback to standard build from metabolite results
        return self.build_complete_results_df()

    def load_metabolite_list(self, file_path):
        """Load metabolite list from file with multi-column support for better matching.
        
        Returns:
            dict with 'names' list and optional ID columns, or None on error
            Format: {
                'names': [list of metabolite names],
                'pubchem_ids': [list of PubChem IDs] (if column exists),
                'hmdb_ids': [list of HMDB IDs] (if column exists),
                'cas_ids': [list of CAS IDs] (if column exists),
                'gene_symbols': [list of gene symbols] (if column exists),
                'accessions': [list of protein accessions e.g., UniProt] (if column exists)
            }
        """
        if not file_path or not os.path.exists(file_path):
            return None
        
        try:
            # Load file into DataFrame
            if file_path.lower().endswith('.xlsx'):
                df = pd.read_excel(file_path)
            elif file_path.lower().endswith('.csv'):
                df = pd.read_csv(file_path)
            else:
                # Assume text file with one metabolite per line (simple name list)
                with open(file_path, 'r', encoding='utf-8') as f:
                    names = [line.strip() for line in f if line.strip()]
                return {'names': names}
            
            if df.empty:
                return None
            
            # Normalize column names for case-insensitive matching
            df.columns = df.columns.str.strip()
            col_lower_map = {col.lower(): col for col in df.columns}
            
            # Result dictionary
            result = {}
            
            # Find Name column (first column by default, or specific name)
            name_col_candidates = ['name', 'metabolite', 'compound', 'molecule', 'metabolite name']
            name_col = None
            for candidate in name_col_candidates:
                if candidate in col_lower_map:
                    name_col = col_lower_map[candidate]
                    break
            
            if name_col is None:
                # Use first column as names
                name_col = df.columns[0]
            
            result['names'] = df[name_col].dropna().astype(str).str.strip().tolist()
            
            # Find PubChem ID column
            pubchem_candidates = ['pubchem id', 'pubchem', 'pubchem_id', 'pubchem cid', 'pubchem_cid', 'cid']
            pubchem_col = None
            for candidate in pubchem_candidates:
                if candidate in col_lower_map:
                    pubchem_col = col_lower_map[candidate]
                    break
            if pubchem_col:
                result['pubchem_ids'] = df[pubchem_col].dropna().astype(str).str.strip().tolist()
            
            # Find HMDB ID column
            hmdb_candidates = ['hmdb id', 'hmdb', 'hmdb_id', 'hmdbid']
            hmdb_col = None
            for candidate in hmdb_candidates:
                if candidate in col_lower_map:
                    hmdb_col = col_lower_map[candidate]
                    break
            if hmdb_col:
                result['hmdb_ids'] = df[hmdb_col].dropna().astype(str).str.strip().tolist()
            
            # Find CAS ID column
            cas_candidates = ['cas id', 'cas', 'cas_id', 'casid', 'cas number', 'cas_number']
            # Find Gene symbol column
            gene_candidates = ['gene', 'symbol', 'gene_symbol', 'genesymbol', 'gene symbol']
            gene_col = None
            for candidate in gene_candidates:
                if candidate in col_lower_map:
                    gene_col = col_lower_map[candidate]
                    break
            if gene_col:
                result['gene_symbols'] = df[gene_col].dropna().astype(str).str.strip().str.upper().tolist()

            # Find protein accession column (UniProt/Swiss-Prot/GenPept)
            acc_candidates = [
                'accession', 'uniprot', 'uniprot_id', 'uniprot accession', 'uniprot accession id',
                'swiss-prot', 'swissprot', 'genpept', 'protein accession', 'protein_accession', 'protein id', 'protein_id',
                'genpept/uniprot/swiss-prot accession'
            ]
            acc_col = None
            for candidate in acc_candidates:
                if candidate in col_lower_map:
                    acc_col = col_lower_map[candidate]
                    break
            if acc_col:
                result['accessions'] = df[acc_col].dropna().astype(str).str.strip().str.upper().tolist()
            cas_col = None
            for candidate in cas_candidates:
                if candidate in col_lower_map:
                    cas_col = col_lower_map[candidate]
                    break
            if cas_col:
                result['cas_ids'] = df[cas_col].dropna().astype(str).str.strip().tolist()
            
            # Debug logging
            self.log_viz_message(f"📋 Loaded custom list: {len(result['names'])} metabolites", 'INFO')
            if 'pubchem_ids' in result:
                self.log_viz_message(f"   ✅ PubChem IDs: {len(result['pubchem_ids'])} entries", 'INFO')
            if 'hmdb_ids' in result:
                self.log_viz_message(f"   ✅ HMDB IDs: {len(result['hmdb_ids'])} entries", 'INFO')
            if 'cas_ids' in result:
                self.log_viz_message(f"   ✅ CAS IDs: {len(result['cas_ids'])} entries", 'INFO')
            if 'gene_symbols' in result:
                self.log_viz_message(f"   ✅ Gene symbols: {len(result['gene_symbols'])} entries", 'INFO')
            if 'accessions' in result:
                self.log_viz_message(f"   ✅ Accessions: {len(result['accessions'])} entries", 'INFO')
            
            return result
            
        except Exception as e:
            self.log_viz_message(f"Failed to load metabolite list: {str(e)}", 'ERROR')
            import traceback
            traceback.print_exc()
            return None
    
    def match_metabolites_multi_column(self, df, metabolite_list_dict, id_col='metabolite_id'):
        """Match metabolites using multiple columns for better accuracy.
        
        Args:
            df: DataFrame with metabolite data
            metabolite_list_dict: Dictionary from load_metabolite_list with names and optional ID columns
            id_col: Primary ID column in df (usually 'metabolite_id' or 'Name')
        
        Returns:
            Boolean mask indicating which rows match the metabolite list
        """
        if metabolite_list_dict is None or not isinstance(metabolite_list_dict, dict):
            return pd.Series([False] * len(df), index=df.index)
        
        # Start with all False
        match_mask = pd.Series([False] * len(df), index=df.index)
        matched_count = {'name': 0, 'pubchem': 0, 'hmdb': 0, 'cas': 0, 'gene': 0, 'accession': 0}
        
        # Match by Name (case-insensitive)
        if 'names' in metabolite_list_dict:
            name_set = {str(n).lower().strip() for n in metabolite_list_dict['names']}
            name_matches = df[id_col].astype(str).str.lower().str.strip().isin(name_set)
            match_mask |= name_matches
            matched_count['name'] = name_matches.sum()
        
        # Match by PubChem ID if available in both list and dataframe
        if 'pubchem_ids' in metabolite_list_dict:
            pubchem_col = None
            for col in ['PubChem_CID', 'PubChem', 'pubchem_cid', 'pubchem']:
                if col in df.columns:
                    pubchem_col = col
                    break
            
            if pubchem_col:
                pubchem_set = {str(pid).strip() for pid in metabolite_list_dict['pubchem_ids']}
                pubchem_matches = df[pubchem_col].astype(str).str.strip().isin(pubchem_set)
                new_matches = pubchem_matches & ~match_mask
                match_mask |= pubchem_matches
                matched_count['pubchem'] = new_matches.sum()
        
        # Match by HMDB ID if available
        if 'hmdb_ids' in metabolite_list_dict:
            hmdb_col = None
            for col in ['HMDB_ID', 'HMDB', 'hmdb_id', 'hmdb']:
                if col in df.columns:
                    hmdb_col = col
                    break
            
            if hmdb_col:
                hmdb_set = {str(hid).strip() for hid in metabolite_list_dict['hmdb_ids']}
                hmdb_matches = df[hmdb_col].astype(str).str.strip().isin(hmdb_set)
                new_matches = hmdb_matches & ~match_mask
                match_mask |= hmdb_matches
                matched_count['hmdb'] = new_matches.sum()
        
        # Match by CAS ID if available
        if 'cas_ids' in metabolite_list_dict:
            cas_col = None
            for col in ['CAS', 'cas', 'CAS_ID', 'cas_id']:
                if col in df.columns:
                    cas_col = col
                    break
            
            if cas_col:
                cas_set = {str(cid).strip() for cid in metabolite_list_dict['cas_ids']}
                cas_matches = df[cas_col].astype(str).str.strip().isin(cas_set)
                new_matches = cas_matches & ~match_mask
                match_mask |= cas_matches
                matched_count['cas'] = new_matches.sum()

        # Match by Gene symbol if available
        if 'gene_symbols' in metabolite_list_dict:
            gene_col = None
            for col in ['Gene', 'gene', 'Symbol', 'symbol', 'Gene_Symbol', 'GeneSymbol', 'gene_symbol', 'genesymbol']:
                if col in df.columns:
                    gene_col = col
                    break
            if gene_col:
                gene_set = {str(gs).strip().upper() for gs in metabolite_list_dict['gene_symbols']}
                gene_matches = df[gene_col].astype(str).str.strip().str.upper().isin(gene_set)
                new_matches = gene_matches & ~match_mask
                match_mask |= gene_matches
                matched_count['gene'] = new_matches.sum()

        # Match by protein accession if available
        if 'accessions' in metabolite_list_dict:
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
                acc_set = {str(a).strip().upper() for a in metabolite_list_dict['accessions']}
                acc_matches = df[acc_col].astype(str).str.strip().str.upper().isin(acc_set)
                new_matches = acc_matches & ~match_mask
                match_mask |= acc_matches
                matched_count['accession'] = new_matches.sum()
        
        # Debug logging
        total_matched = match_mask.sum()
        self.log_viz_message(f"🔍 Multi-column matching results:", 'INFO')
        self.log_viz_message(f"   Total matched: {total_matched}/{len(df)} metabolites", 'INFO')
        self.log_viz_message(f"   Matched by Name: {matched_count['name']}", 'INFO')
        if matched_count['pubchem'] > 0:
            self.log_viz_message(f"   Additional matched by PubChem ID: {matched_count['pubchem']}", 'INFO')
        if matched_count['hmdb'] > 0:
            self.log_viz_message(f"   Additional matched by HMDB ID: {matched_count['hmdb']}", 'INFO')
        if matched_count['cas'] > 0:
            self.log_viz_message(f"   Additional matched by CAS ID: {matched_count['cas']}", 'INFO')
        if matched_count['gene'] > 0:
            self.log_viz_message(f"   Additional matched by Gene symbol: {matched_count['gene']}", 'INFO')
        if matched_count['accession'] > 0:
            self.log_viz_message(f"   Additional matched by Accession: {matched_count['accession']}", 'INFO')
        
        return match_mask

    def run_selected_visualizations(self):
        """Run only the selected visualization types."""
        selected = [name for name, var in self.viz_selected.items() if var.get()]
        if not selected:
            messagebox.showwarning("No Selection", "Please select at least one visualization type.")
            return
        self._run_visualizations(selected)

    def run_all_visualizations(self):
        """Run all visualization types."""
        all_types = list(self.viz_selected.keys())
        self._run_visualizations(all_types)

    def stop_visualizations(self):
        """Stop running visualizations."""
        self.viz_cancel_flag.set()
        self.log_viz_message("Cancellation requested...", 'WARN')

    def _run_visualizations(self, viz_types):
        """Execute visualizations in a background thread."""
        # Validate data is ready
        try:
            complete_df = self.build_complete_results_df()
            if complete_df is None:
                messagebox.showerror("Data Not Ready", "No statistical results available. Please run statistical analysis first.")
                return
        except Exception as e:
            messagebox.showerror("Data Error", f"Failed to prepare visualization data:\n{str(e)}")
            return
        
        # Start visualization thread
        self.viz_cancel_flag.clear()
        thread = threading.Thread(target=self._visualization_worker, args=(viz_types,))
        thread.daemon = True
        thread.start()

    def _visualization_worker(self, viz_types):
        """Background worker for visualization generation."""
        try:
            self.log_viz_message("Starting visualization generation...")
            
            # Prepare context
            complete_df = self.build_complete_results_df()
            groups = self.get_current_groups()
            sample_cols, sample_to_group = self.get_sample_mapping()
            
            from main_script.metabolites_visualization import CommonVizContext, validate_visualization_ready
            
            # Validate
            is_ready, message = validate_visualization_ready(complete_df, groups, sample_cols)
            if not is_ready:
                self.log_viz_message(f"Validation failed: {message}", 'ERROR')
                return
            
            # Determine ID column from multiple sources
            id_col = None
            # Priority 1: stat_column_assignments (from saved config)
            if getattr(self, 'stat_column_assignments', None):
                id_col = self.stat_column_assignments.get('id_column')
            # Priority 2: memory_store (from Statistics tab)
            if not id_col and hasattr(self, 'memory_store') and self.memory_store:
                # Check for lipid class ID column first (more specific)
                if 'id_column_class' in self.memory_store:
                    id_col = self.memory_store['id_column_class']
                elif 'id_column' in self.memory_store:
                    id_col = self.memory_store['id_column']
            # Priority 3: Auto-detect from dataframe columns (for lipid class data)
            if not id_col and 'Class' in complete_df.columns:
                # If 'Class' column exists and we're likely using lipid class data, use it as ID
                # Check if this looks like lipid class data (has Class column but not LipidID)
                if 'LipidID' not in complete_df.columns or complete_df['LipidID'].isna().all():
                    id_col = 'Class'
            
            # Create context
            context = CommonVizContext(
                complete_df=complete_df,
                groups=groups,
                sample_cols=sample_cols,
                sample_to_group=sample_to_group,
                output_dir=self.viz_output_dir.get(),
                color_map=self.build_color_map(),
                preferred_group_order=[g.strip() for g in self.viz_preferred_group_order.get().split(',') if g.strip()] if hasattr(self, 'viz_preferred_group_order') and self.viz_preferred_group_order.get().strip() else None,
                id_column=id_col
            )
            
            # Ensure output directory exists
            os.makedirs(context.output_dir, exist_ok=True)
            
            # Run each selected visualization type
            total_types = len(viz_types)
            for i, viz_type in enumerate(viz_types):
                if self.viz_cancel_flag.is_set():
                    self.log_viz_message("Visualization cancelled by user", 'WARN')
                    break
                
                # Update progress bar and status
                progress_value = (i / total_types) * 100
                self.root.after(0, lambda v=progress_value: self.viz_progress.configure(value=v))
                self.root.after(0, lambda t=viz_type, c=i+1, total=total_types: self.viz_progress_label.configure(text=f"Generating {t} ({c}/{total})..."))
                
                self.log_viz_message(f"Generating {viz_type} ({i+1}/{total_types})...")
                
                try:
                    result = self._run_single_visualization(context, viz_type)
                    self.log_viz_message(result.summary, 'SUCCESS' if not result.errors else 'WARN')
                    
                    for error in result.errors:
                        self.log_viz_message(f"Error: {error}", 'ERROR')
                        
                except Exception as e:
                    self.log_viz_message(f"Failed to generate {viz_type}: {str(e)}", 'ERROR')
            
            # Update final progress
            self.root.after(0, lambda: self.viz_progress.configure(value=100))
            self.root.after(0, lambda: self.viz_progress_label.configure(text="Complete!"))
            self.log_viz_message("Visualization generation finished!")
            
        except Exception as e:
            self.log_viz_message(f"Visualization worker failed: {str(e)}", 'ERROR')
            import traceback
            self.log_viz_message(traceback.format_exc(), 'ERROR')

    def _run_single_visualization(self, context, viz_type):
        """Run a single visualization type."""
        from main_script.metabolites_visualization import (
            run_pca_analysis, run_volcano_analysis, run_boxplot_analysis, run_bargraph_analysis,
            run_heatmap_analysis, run_roc_analysis
        )
        
        # Update parameters from GUI
        self._update_viz_params()
        
        if viz_type == 'pca':
            return run_pca_analysis(context, self.viz_params['pca'])
        elif viz_type == 'volcano':
            return run_volcano_analysis(context, self.viz_params['volcano'])
        elif viz_type == 'boxplot':
            return run_boxplot_analysis(context, self.viz_params['boxplot'])
        elif viz_type == 'bargraph':
            return run_bargraph_analysis(context, self.viz_params['bargraph'])
        elif viz_type == 'heatmap':
            return run_heatmap_analysis(context, self.viz_params['heatmap'])
        elif viz_type == 'roc':
            return run_roc_analysis(context, self.viz_params['roc'])
        else:
            raise ValueError(f"Unknown visualization type: {viz_type}")

    def _update_viz_params(self):
        # ROC filtering controls
        if hasattr(self, 'roc_filter_mode'):
            self.viz_params['roc'].filter_mode = self.roc_filter_mode.get()
        # Parse specific comparison textbox
        if hasattr(self, 'roc_specific_comparison') and self.roc_specific_comparison.get().strip():
            comp_txt = self.roc_specific_comparison.get().strip()
            if '|' in comp_txt:
                parts = comp_txt.split('|', 1)
                if len(parts) == 2:
                    g1, g2 = parts[0].strip(), parts[1].strip()
                    if g1 and g2:
                        self.viz_params['roc'].filter_pairs = [(g1, g2)]
        if hasattr(self, 'roc_filter_pairs') and self.roc_filter_pairs.get().strip():
            pairs_txt = self.roc_filter_pairs.get().strip()
            parsed = []
            for token in pairs_txt.split(';'):
                if '|' in token:
                    a,b = token.split('|',1)
                    a = a.strip(); b = b.strip()
                    if a and b:
                        parsed.append((a,b))
            self.viz_params['roc'].filter_pairs = parsed if parsed else None
        if hasattr(self, 'roc_use_custom_only'):
            self.viz_params['roc'].use_custom_only = self.roc_use_custom_only.get()
        if hasattr(self, 'roc_excel_only'):
            self.viz_params['roc'].excel_only = self.roc_excel_only.get()
        """Update parameter objects from GUI values."""
        # PCA parameters
        self.viz_params['pca'].components = self.pca_components.get()
        self.viz_params['pca'].plot_3d = self.pca_3d.get()
        self.viz_params['pca'].interactive_3d = self.pca_interactive.get()
        self.viz_params['pca'].scree = self.pca_scree.get()
        self.viz_params['pca'].loadings = self.pca_loadings.get()
        self.viz_params['pca'].loadings_top_k = self.pca_loadings_k.get()
        if hasattr(self, 'pca_fig_width'):
            self.viz_params['pca'].fig_width = self.pca_fig_width.get()
            self.viz_params['pca'].fig_height = self.pca_fig_height.get()
            self.viz_params['pca'].fig_dpi = self.pca_fig_dpi.get()
        # Point sizes
        if hasattr(self, 'pca_point_size_2d'):
            self.viz_params['pca'].point_size_2d = self.pca_point_size_2d.get()
        if hasattr(self, 'pca_point_size_3d'):
            self.viz_params['pca'].point_size_3d = self.pca_point_size_3d.get()
        # 3D viewing angles
        if hasattr(self, 'pca_view_azim'):
            self.viz_params['pca'].view_azim = self.pca_view_azim.get()
        if hasattr(self, 'pca_view_elev'):
            self.viz_params['pca'].view_elev = self.pca_view_elev.get()
        # Save options
        if hasattr(self, 'pca_save_2d'):
            self.viz_params['pca'].save_2d = self.pca_save_2d.get()
        if hasattr(self, 'pca_save_3d'):
            self.viz_params['pca'].save_3d = self.pca_save_3d.get()
        if hasattr(self, 'pca_save_excel'):
            self.viz_params['pca'].save_excel = self.pca_save_excel.get()
        # Font sizes
        if hasattr(self, 'pca_xlabel_fontsize'):
            self.viz_params['pca'].xlabel_fontsize = self.pca_xlabel_fontsize.get()
            self.viz_params['pca'].ylabel_fontsize = self.pca_ylabel_fontsize.get()
            self.viz_params['pca'].title_fontsize = self.pca_title_fontsize.get()
            self.viz_params['pca'].tick_fontsize = self.pca_tick_fontsize.get()
        if hasattr(self, 'pca_legend_fontsize'):
            self.viz_params['pca'].legend_fontsize = self.pca_legend_fontsize.get()
        if hasattr(self, 'pca_show_legend'):
            self.viz_params['pca'].show_legend = bool(self.pca_show_legend.get())
        if hasattr(self, 'pca_output_format'):
            self.viz_params['pca'].output_format = self.pca_output_format.get()
        # Custom group selection
        if hasattr(self, 'pca_custom_group_vars'):
            selected = [g for g, var in self.pca_custom_group_vars.items() if var.get()]
            self.viz_params['pca'].custom_groups = selected if selected else None
        # Comparison selection
        if hasattr(self, 'pca_selected_comparisons'):
            self.viz_params['pca'].selected_comparisons = self.pca_selected_comparisons
        
        # Volcano parameters
        self.viz_params['volcano'].p_threshold = self.volcano_p_thresh.get()
        # If user requested to skip FC cutoff, set threshold to 0 so metabolites_visualization treats it as 'no FC cutoff'
        if hasattr(self, 'volcano_skip_fc') and self.volcano_skip_fc.get():
            self.viz_params['volcano'].fc_threshold = 0
        else:
            self.viz_params['volcano'].fc_threshold = self.volcano_fc_thresh.get()
        self.viz_params['volcano'].annotate = self.volcano_annotate.get()
        self.viz_params['volcano'].annotate_top_n = self.volcano_top_n.get()
        self.viz_params['volcano'].annot_fontsize = self.volcano_annot_fontsize.get()
        self.viz_params['volcano'].output_format = self.volcano_output_format.get()
        
        # Debug: Log volcano parameters
        logger.info(f"📊 Volcano params from GUI:")
        logger.info(f"   ├─ p_threshold = {self.viz_params['volcano'].p_threshold}")
        logger.info(f"   ├─ fc_threshold = {self.viz_params['volcano'].fc_threshold} {'(SKIPPED)' if self.viz_params['volcano'].fc_threshold == 0 else ''}")
        logger.info(f"   ├─ annotate = {self.viz_params['volcano'].annotate}")
        logger.info(f"   └─ annotate_top_n = {self.viz_params['volcano'].annotate_top_n}")
        
        # Lipid class option removed for volcano
        # if hasattr(self, 'viz_include_lipid_class') and 'volcano' in self.viz_include_lipid_class:
        #     self.viz_params['volcano'].include_lipid_class = self.viz_include_lipid_class['volcano'].get()
        self.viz_params['volcano'].include_lipid_class = False  # Force disable for volcano
        if hasattr(self, 'volcano_fig_width'):
            self.viz_params['volcano'].fig_width = self.volcano_fig_width.get()
            self.viz_params['volcano'].fig_height = self.volcano_fig_height.get()
            self.viz_params['volcano'].fig_dpi = self.volcano_fig_dpi.get()
        # Point sizes
        if hasattr(self, 'volcano_point_size_sig'):
            self.viz_params['volcano'].point_size_sig = self.volcano_point_size_sig.get()
        if hasattr(self, 'volcano_point_size_nonsig'):
            self.viz_params['volcano'].point_size_nonsig = self.volcano_point_size_nonsig.get()
        # Font sizes
        if hasattr(self, 'volcano_xlabel_fontsize'):
            self.viz_params['volcano'].xlabel_fontsize = self.volcano_xlabel_fontsize.get()
            self.viz_params['volcano'].ylabel_fontsize = self.volcano_ylabel_fontsize.get()
            self.viz_params['volcano'].title_fontsize = self.volcano_title_fontsize.get()
            self.viz_params['volcano'].tick_fontsize = self.volcano_tick_fontsize.get()
            self.viz_params['volcano'].count_fontsize = self.volcano_count_fontsize.get()
            self.viz_params['volcano'].total_fontsize = self.volcano_total_fontsize.get()
        if hasattr(self, 'volcano_legend_fontsize'):
            self.viz_params['volcano'].legend_fontsize = self.volcano_legend_fontsize.get()
        # Count background style
        if hasattr(self, 'volcano_count_background'):
            self.viz_params['volcano'].count_background = self.volcano_count_background.get()
        # Excel export control
        if hasattr(self, 'volcano_save_excel'):
            self.viz_params['volcano'].save_excel = self.volcano_save_excel.get()
        # Comparison selection
        if hasattr(self, 'volcano_selected_comparisons'):
            self.viz_params['volcano'].selected_comparisons = self.volcano_selected_comparisons
        
        # Boxplot parameters
        self.viz_params['boxplot'].top_n = self.boxplot_top_n.get()
        self.viz_params['boxplot'].no_limit = self.boxplot_no_limit.get()
        self.viz_params['boxplot'].annotate = self.boxplot_annotate.get()
        self.viz_params['boxplot'].include_metabolites = None
        # Lipid class option removed for boxplot
        # if hasattr(self, 'viz_include_lipid_class') and 'boxplot' in self.viz_include_lipid_class:
        #     self.viz_params['boxplot'].include_lipid_class = self.viz_include_lipid_class['boxplot'].get()
        self.viz_params['boxplot'].include_lipid_class = False  # Force disable for boxplot
        if hasattr(self, 'boxplot_p_thresh'):
            self.viz_params['boxplot'].p_threshold = self.boxplot_p_thresh.get()
        # Honor skip-FC option for boxplots
        if hasattr(self, 'boxplot_skip_fc') and self.boxplot_skip_fc.get():
            self.viz_params['boxplot'].fc_threshold = 0
        elif hasattr(self, 'boxplot_fc_thresh'):
            self.viz_params['boxplot'].fc_threshold = self.boxplot_fc_thresh.get()
        if hasattr(self, 'boxplot_filter_mode'):
            self.viz_params['boxplot'].filter_mode = self.boxplot_filter_mode.get()
        # Parse specific comparison textbox
        if hasattr(self, 'boxplot_specific_comparison') and self.boxplot_specific_comparison.get().strip():
            comp_txt = self.boxplot_specific_comparison.get().strip()
            if '|' in comp_txt:
                parts = comp_txt.split('|', 1)
                if len(parts) == 2:
                    g1, g2 = parts[0].strip(), parts[1].strip()
                    if g1 and g2:
                        self.viz_params['boxplot'].filter_pairs = [(g1, g2)]
        if hasattr(self, 'boxplot_filter_pairs') and self.boxplot_filter_pairs.get().strip():
            pairs_txt = self.boxplot_filter_pairs.get().strip()
            parsed = []
            for token in pairs_txt.split(';'):
                if '|' in token:
                    a,b = token.split('|',1)
                    a = a.strip(); b = b.strip()
                    if a and b:
                        parsed.append((a,b))
            self.viz_params['boxplot'].filter_pairs = parsed if parsed else None
        if hasattr(self, 'boxplot_use_custom_only'):
            self.viz_params['boxplot'].use_custom_only = self.boxplot_use_custom_only.get()
        if hasattr(self, 'boxplot_fig_width'):
            self.viz_params['boxplot'].fig_width = self.boxplot_fig_width.get()
            self.viz_params['boxplot'].fig_height = self.boxplot_fig_height.get()
            self.viz_params['boxplot'].fig_dpi = self.boxplot_fig_dpi.get()
        # Font sizes
        if hasattr(self, 'boxplot_xlabel_fontsize'):
            self.viz_params['boxplot'].xlabel_fontsize = self.boxplot_xlabel_fontsize.get()
            self.viz_params['boxplot'].ylabel_fontsize = self.boxplot_ylabel_fontsize.get()
            self.viz_params['boxplot'].title_fontsize = self.boxplot_title_fontsize.get()
            self.viz_params['boxplot'].tick_fontsize = self.boxplot_tick_fontsize.get()
        # Title wrap width
        if hasattr(self, 'boxplot_title_wrap_width'):
            self.viz_params['boxplot'].title_wrap_width = self.boxplot_title_wrap_width.get()
        # Y-axis label text
        if hasattr(self, 'boxplot_ylabel_text'):
            ylabel_text = self.boxplot_ylabel_text.get().strip()
            if ylabel_text:  # Only set if non-empty
                self.viz_params['boxplot'].ylabel_text = ylabel_text
        # X-tick rotation controls
        if hasattr(self, 'boxplot_rotate_xticks'):
            self.viz_params['boxplot'].rotate_xticks = self.boxplot_rotate_xticks.get()
        if hasattr(self, 'boxplot_xtick_rotation'):
            try:
                self.viz_params['boxplot'].xtick_rotation = int(self.boxplot_xtick_rotation.get())
            except Exception:
                # Keep whatever default in params if parsing fails
                pass
        # Excel export control
        if hasattr(self, 'boxplot_save_excel'):
            self.viz_params['boxplot'].save_excel = self.boxplot_save_excel.get()
        # Comparison selection
        if hasattr(self, 'boxplot_selected_comparisons'):
            self.viz_params['boxplot'].selected_comparisons = self.boxplot_selected_comparisons
        # Annotation comparison selection
        if hasattr(self, 'boxplot_annotate_comparisons'):
            self.viz_params['boxplot'].annotate_comparisons = self.boxplot_annotate_comparisons
        # Group selection - which groups to display
        if hasattr(self, 'boxplot_selected_groups'):
            self.viz_params['boxplot'].selected_groups = self.boxplot_selected_groups
        # Filter comparison for All/Specific modes
        if hasattr(self, 'boxplot_filter_comparison'):
            self.viz_params['boxplot'].filter_comparison = self.boxplot_filter_comparison
        # Load custom metabolite list if provided
        if self.boxplot_custom_list.get():
            metabolites = self.load_metabolite_list(self.boxplot_custom_list.get())
            self.viz_params['boxplot'].include_metabolites = metabolites
            if hasattr(self, 'boxplot_use_custom_only') and self.boxplot_use_custom_only.get():
                logger.info("   🎯 Boxplot: 'Use custom list ONLY' mode ENABLED - ignoring p-value/FC thresholds")
            else:
                logger.info("   🔀 Boxplot: Custom list will be intersected with p-value/FC filtered metabolites")

        # Bargraph parameters (independent from boxplot)
        if 'bargraph' in self.viz_params:
            self.viz_params['bargraph'].include_metabolites = None
            if hasattr(self, 'bargraph_top_n'):
                self.viz_params['bargraph'].top_n = self.bargraph_top_n.get()
            if hasattr(self, 'bargraph_no_limit'):
                self.viz_params['bargraph'].no_limit = self.bargraph_no_limit.get()
            if hasattr(self, 'bargraph_annotate'):
                self.viz_params['bargraph'].annotate = self.bargraph_annotate.get()
            if hasattr(self, 'bargraph_p_thresh'):
                self.viz_params['bargraph'].p_threshold = self.bargraph_p_thresh.get()
            if hasattr(self, 'bargraph_skip_fc') and self.bargraph_skip_fc.get():
                self.viz_params['bargraph'].fc_threshold = 0
            elif hasattr(self, 'bargraph_fc_thresh'):
                self.viz_params['bargraph'].fc_threshold = self.bargraph_fc_thresh.get()
            if hasattr(self, 'bargraph_filter_mode'):
                self.viz_params['bargraph'].filter_mode = self.bargraph_filter_mode.get()
            if hasattr(self, 'bargraph_specific_comparison') and self.bargraph_specific_comparison.get().strip():
                comp_txt = self.bargraph_specific_comparison.get().strip()
                if '|' in comp_txt:
                    parts = comp_txt.split('|', 1)
                    if len(parts) == 2:
                        g1, g2 = parts[0].strip(), parts[1].strip()
                        if g1 and g2:
                            self.viz_params['bargraph'].filter_pairs = [(g1, g2)]
            if hasattr(self, 'bargraph_use_custom_only'):
                self.viz_params['bargraph'].use_custom_only = self.bargraph_use_custom_only.get()
            if hasattr(self, 'bargraph_fig_width'):
                self.viz_params['bargraph'].fig_width = self.bargraph_fig_width.get()
                self.viz_params['bargraph'].fig_height = self.bargraph_fig_height.get()
                self.viz_params['bargraph'].fig_dpi = self.bargraph_fig_dpi.get()
            if hasattr(self, 'bargraph_xlabel_fontsize'):
                self.viz_params['bargraph'].xlabel_fontsize = self.bargraph_xlabel_fontsize.get()
                self.viz_params['bargraph'].ylabel_fontsize = self.bargraph_ylabel_fontsize.get()
                self.viz_params['bargraph'].title_fontsize = self.bargraph_title_fontsize.get()
                self.viz_params['bargraph'].tick_fontsize = self.bargraph_tick_fontsize.get()
            if hasattr(self, 'bargraph_legend_fontsize'):
                self.viz_params['bargraph'].legend_fontsize = self.bargraph_legend_fontsize.get()
            if hasattr(self, 'bargraph_title_wrap_width'):
                self.viz_params['bargraph'].title_wrap_width = self.bargraph_title_wrap_width.get()
            if hasattr(self, 'bargraph_ylabel_text'):
                ylabel_text = self.bargraph_ylabel_text.get().strip()
                if ylabel_text:
                    self.viz_params['bargraph'].ylabel_text = ylabel_text
            if hasattr(self, 'bargraph_rotate_xticks'):
                self.viz_params['bargraph'].rotate_xticks = self.bargraph_rotate_xticks.get()
            if hasattr(self, 'bargraph_xtick_rotation'):
                try:
                    self.viz_params['bargraph'].xtick_rotation = int(self.bargraph_xtick_rotation.get())
                except Exception:
                    pass
            if hasattr(self, 'bargraph_save_excel'):
                self.viz_params['bargraph'].save_excel = self.bargraph_save_excel.get()
            if hasattr(self, 'bargraph_selected_comparisons'):
                self.viz_params['bargraph'].selected_comparisons = self.bargraph_selected_comparisons
            if hasattr(self, 'bargraph_annotate_comparisons'):
                self.viz_params['bargraph'].annotate_comparisons = self.bargraph_annotate_comparisons
            if hasattr(self, 'bargraph_selected_groups'):
                self.viz_params['bargraph'].selected_groups = self.bargraph_selected_groups
            if hasattr(self, 'bargraph_display_mode'):
                self.viz_params['bargraph'].display_mode = self.bargraph_display_mode.get().strip().lower()
            if hasattr(self, 'bargraph_grouped_title'):
                self.viz_params['bargraph'].grouped_title = self.bargraph_grouped_title.get().strip()
            if hasattr(self, 'bargraph_low_boost_enabled'):
                self.viz_params['bargraph'].low_value_boost_enabled = self.bargraph_low_boost_enabled.get()
            if hasattr(self, 'bargraph_low_boost_threshold'):
                try:
                    self.viz_params['bargraph'].low_value_boost_threshold = float(self.bargraph_low_boost_threshold.get())
                except Exception:
                    pass
            if hasattr(self, 'bargraph_low_boost_factor'):
                try:
                    self.viz_params['bargraph'].low_value_boost_factor = max(1.0, float(self.bargraph_low_boost_factor.get()))
                except Exception:
                    pass
            if hasattr(self, 'bargraph_custom_list') and self.bargraph_custom_list.get():
                metabolites = self.load_metabolite_list(self.bargraph_custom_list.get())
                self.viz_params['bargraph'].include_metabolites = metabolites
                if hasattr(self, 'bargraph_use_custom_only') and self.bargraph_use_custom_only.get():
                    logger.info("   🎯 Bargraph: 'Use custom list ONLY' mode ENABLED - ignoring p-value/FC thresholds")
                else:
                    logger.info("   🔀 Bargraph: Custom list will be intersected with p-value/FC filtered metabolites")
        
        # Heatmap parameters
        self.viz_params['heatmap'].max_metabolites = self.heatmap_max.get()
        if hasattr(self, 'heatmap_show_fc_divider'):
            self.viz_params['heatmap'].show_fc_divider = self.heatmap_show_fc_divider.get()
        if hasattr(self, 'heatmap_combined'):
            self.viz_params['heatmap'].combined = self.heatmap_combined.get()
        if hasattr(self, 'heatmap_combined_mode'):
            self.viz_params['heatmap'].combined_mode = self.heatmap_combined_mode.get()
        if hasattr(self, 'heatmap_p_thresh'):
            self.viz_params['heatmap'].p_threshold = self.heatmap_p_thresh.get()
        if hasattr(self, 'heatmap_p_thresh'):
            self.viz_params['heatmap'].p_threshold = self.heatmap_p_thresh.get()
        # Honor skip-FC option for heatmaps
        if hasattr(self, 'heatmap_skip_fc') and self.heatmap_skip_fc.get():
            self.viz_params['heatmap'].fc_threshold = 0
        elif hasattr(self, 'heatmap_fc_thresh'):
            self.viz_params['heatmap'].fc_threshold = self.heatmap_fc_thresh.get()
        if hasattr(self, 'heatmap_filter_mode'):
            self.viz_params['heatmap'].filter_mode = self.heatmap_filter_mode.get()
        # Parse specific comparison textbox
        if hasattr(self, 'heatmap_specific_comparison') and self.heatmap_specific_comparison.get().strip():
            comp_txt = self.heatmap_specific_comparison.get().strip()
            if '|' in comp_txt:
                parts = comp_txt.split('|', 1)
                if len(parts) == 2:
                    g1, g2 = parts[0].strip(), parts[1].strip()
                    if g1 and g2:
                        self.viz_params['heatmap'].filter_pairs = [(g1, g2)]
        if hasattr(self, 'heatmap_filter_pairs') and self.heatmap_filter_pairs.get().strip():
            pairs_txt = self.heatmap_filter_pairs.get().strip()
            parsed = []
            for token in pairs_txt.split(';'):
                if '|' in token:
                    a,b = token.split('|',1)
                    a = a.strip(); b = b.strip()
                    if a and b:
                        parsed.append((a,b))
            self.viz_params['heatmap'].filter_pairs = parsed if parsed else None
        if hasattr(self, 'heatmap_use_custom_only'):
            self.viz_params['heatmap'].use_custom_only = self.heatmap_use_custom_only.get()
        # Divider line options
        if hasattr(self, 'heatmap_show_fc_divider'):
            self.viz_params['heatmap'].show_fc_divider = self.heatmap_show_fc_divider.get()
        if hasattr(self, 'heatmap_show_sample_divider'):
            # Invert because parameter is "no_col_split" but GUI shows "show_sample_divider"
            self.viz_params['heatmap'].no_col_split = not self.heatmap_show_sample_divider.get()
        # Color scale options
        if hasattr(self, 'heatmap_use_fixed_scale'):
            self.viz_params['heatmap'].use_fixed_scale = self.heatmap_use_fixed_scale.get()
        if hasattr(self, 'heatmap_auto_scale'):
            self.viz_params['heatmap'].auto_scale = self.heatmap_auto_scale.get()
        if hasattr(self, 'heatmap_vmin'):
            self.viz_params['heatmap'].vmin = self.heatmap_vmin.get()
        if hasattr(self, 'heatmap_vmax'):
            self.viz_params['heatmap'].vmax = self.heatmap_vmax.get()
        if hasattr(self, 'heatmap_fig_width'):
            self.viz_params['heatmap'].fig_width = self.heatmap_fig_width.get()
            self.viz_params['heatmap'].fig_height = self.heatmap_fig_height.get()
            self.viz_params['heatmap'].fig_dpi = self.heatmap_fig_dpi.get()
        if hasattr(self, 'heatmap_auto_size'):
            self.viz_params['heatmap'].auto_size = self.heatmap_auto_size.get()
        # Layout controls
        if hasattr(self, 'heatmap_cluster'):
            self.viz_params['heatmap'].cluster = self.heatmap_cluster.get()
        if hasattr(self, 'heatmap_show_colorbar'):
            self.viz_params['heatmap'].show_colorbar = self.heatmap_show_colorbar.get()
        if hasattr(self, 'heatmap_dendro_width_pct'):
            try:
                self.viz_params['heatmap'].dendrogram_width_ratio = float(self.heatmap_dendro_width_pct.get()) / 100.0
            except Exception:
                self.viz_params['heatmap'].dendrogram_width_ratio = 0.18
        if hasattr(self, 'heatmap_cbar_height_inches'):
            try:
                self.viz_params['heatmap'].colorbar_height_inches = float(self.heatmap_cbar_height_inches.get())
            except Exception:
                self.viz_params['heatmap'].colorbar_height_inches = 0.6
        # Font sizes - simplified
        if hasattr(self, 'heatmap_feature_fontsize'):
            self.viz_params['heatmap'].feature_fontsize = self.heatmap_feature_fontsize.get()
        if hasattr(self, 'heatmap_sample_fontsize'):
            self.viz_params['heatmap'].sample_fontsize = self.heatmap_sample_fontsize.get()
        if hasattr(self, 'heatmap_title_fontsize'):
            self.viz_params['heatmap'].title_fontsize = self.heatmap_title_fontsize.get()
        # Excel export control
        if hasattr(self, 'heatmap_save_excel'):
            self.viz_params['heatmap'].save_excel = self.heatmap_save_excel.get()
        if hasattr(self, 'heatmap_output_format'):
            self.viz_params['heatmap'].output_format = self.heatmap_output_format.get()
        # Comparison selection
        if hasattr(self, 'heatmap_selected_comparisons'):
            self.viz_params['heatmap'].selected_comparisons = self.heatmap_selected_comparisons
        # Per-comparison metabolite lists
        if hasattr(self, 'heatmap_metabolite_lists'):
            # Convert file paths to metabolite lists
            metabolite_lists_dict = {}
            for pair, filepath in self.heatmap_metabolite_lists.items():
                if filepath and os.path.exists(filepath):
                    metabolites = self.load_metabolite_list(filepath)
                    if metabolites:
                        metabolite_lists_dict[pair] = metabolites
            self.viz_params['heatmap'].metabolite_lists = metabolite_lists_dict
        # Skip unlisted comparisons
        if hasattr(self, 'heatmap_skip_unlisted'):
            self.viz_params['heatmap'].skip_unlisted_comparisons = self.heatmap_skip_unlisted.get()
        # Load custom metabolite list if provided
        if hasattr(self, 'heatmap_custom_list') and self.heatmap_custom_list.get():
            metabolites = self.load_metabolite_list(self.heatmap_custom_list.get())
            self.viz_params['heatmap'].include_metabolites = metabolites
            if hasattr(self, 'heatmap_use_custom_only') and self.heatmap_use_custom_only.get():
                logger.info("   🎯 Heatmap: 'Use custom list ONLY' mode ENABLED - ignoring p-value/FC thresholds")
            else:
                logger.info("   🔀 Heatmap: Custom list will be intersected with p-value/FC filtered metabolites")
        
        # ROC parameters
        self.viz_params['roc'].all_pairs = self.roc_all_pairs.get()
        self.viz_params['roc'].min_auc = self.roc_min_auc.get()
        self.viz_params['roc'].max_metabolites = self.roc_max_metabolites.get()
        if hasattr(self, 'roc_include_combined'):
            self.viz_params['roc'].include_combined = self.roc_include_combined.get()
        # Lipid class option removed for ROC
        # if hasattr(self, 'viz_include_lipid_class') and 'roc' in self.viz_include_lipid_class:
        #     self.viz_params['roc'].include_lipid_class = self.viz_include_lipid_class['roc'].get()
        self.viz_params['roc'].include_lipid_class = False  # Force disable for ROC
        if hasattr(self, 'roc_p_thresh'):
            self.viz_params['roc'].p_threshold = self.roc_p_thresh.get()
        if hasattr(self, 'roc_p_thresh'):
            self.viz_params['roc'].p_threshold = self.roc_p_thresh.get()
        # Honor skip-FC option for ROC
        if hasattr(self, 'roc_skip_fc') and self.roc_skip_fc.get():
            self.viz_params['roc'].fc_threshold = 0
        elif hasattr(self, 'roc_fc_thresh'):
            self.viz_params['roc'].fc_threshold = self.roc_fc_thresh.get()
        if hasattr(self, 'roc_fig_width'):
            self.viz_params['roc'].fig_width = self.roc_fig_width.get()
            self.viz_params['roc'].fig_height = self.roc_fig_height.get()
            self.viz_params['roc'].fig_dpi = self.roc_fig_dpi.get()
        # Font sizes
        if hasattr(self, 'roc_xlabel_fontsize'):
            self.viz_params['roc'].xlabel_fontsize = self.roc_xlabel_fontsize.get()
            self.viz_params['roc'].ylabel_fontsize = self.roc_ylabel_fontsize.get()
            self.viz_params['roc'].title_fontsize = self.roc_title_fontsize.get()
            self.viz_params['roc'].tick_fontsize = self.roc_tick_fontsize.get()
        if hasattr(self, 'roc_legend_fontsize'):
            self.viz_params['roc'].legend_fontsize = self.roc_legend_fontsize.get()
        # Excel export control
        if hasattr(self, 'roc_save_excel'):
            self.viz_params['roc'].save_excel = self.roc_save_excel.get()
        # Comparison selection
        if hasattr(self, 'roc_selected_comparisons'):
            self.viz_params['roc'].selected_comparisons = self.roc_selected_comparisons
        # Per-comparison metabolite lists
        if hasattr(self, 'roc_metabolite_lists'):
            # Convert file paths to metabolite lists
            metabolite_lists_dict = {}
            for pair, filepath in self.roc_metabolite_lists.items():
                if filepath and os.path.exists(filepath):
                    metabolites = self.load_metabolite_list(filepath)
                    if metabolites:
                        metabolite_lists_dict[pair] = metabolites
            self.viz_params['roc'].metabolite_lists = metabolite_lists_dict
        # Skip unlisted comparisons
        if hasattr(self, 'roc_skip_unlisted'):
            self.viz_params['roc'].skip_unlisted_comparisons = self.roc_skip_unlisted.get()
        
        # Custom metabolite lists
        custom_metabolites = None
        if self.viz_metabolite_text.get(1.0, tk.END).strip():
            custom_metabolites = [line.strip() for line in self.viz_metabolite_text.get(1.0, tk.END).strip().split('\n') if line.strip()]
        
        if custom_metabolites:
            if self.viz_custom_list_boxplot.get():
                self.viz_params['boxplot'].include_metabolites = custom_metabolites
            if hasattr(self, 'viz_custom_list_bargraph') and self.viz_custom_list_bargraph.get():
                self.viz_params['bargraph'].include_metabolites = custom_metabolites
            if self.viz_custom_list_heatmap.get():
                self.viz_params['heatmap'].include_metabolites = custom_metabolites
            if self.viz_custom_list_roc.get():
                self.viz_params['roc'].metabolites = custom_metabolites
        self.viz_params['roc'].max_metabolites = self.roc_max_metabolites.get()
        # Load custom metabolite list if provided
        if self.roc_custom_list.get():
            metabolites = self.load_metabolite_list(self.roc_custom_list.get())
            self.viz_params['roc'].metabolites = metabolites
            if hasattr(self, 'roc_use_custom_only') and self.roc_use_custom_only.get():
                logger.info("   🎯 ROC: 'Use custom list ONLY' mode ENABLED - ignoring p-value/FC thresholds")
            else:
                logger.info("   🔀 ROC: Custom list will be intersected with p-value/FC filtered metabolites")

    def build_complete_results_df(self):
        """Build a Complete Results DataFrame for visualization."""
        # Check if using imported file
        if hasattr(self, 'imported_complete_df') and self.viz_data_source.get() == 'file':
            logger.info(f"📂 Using imported file: {len(self.imported_complete_df)} rows")
            return self.imported_complete_df
        
        # Use session data: first try memory_store (shared from Statistics tab)
        logger.info(f"📊 Building from session data...")
        
        # Try memory_store first (populated by Statistics tab after completion)
        if hasattr(self, 'memory_store') and self.memory_store:
            if 'statistical_test_results' in self.memory_store and self.memory_store['statistical_test_results']:
                stat_results = self.memory_store['statistical_test_results']
                if isinstance(stat_results, dict) and 'enhanced_metabolites' in stat_results:
                    complete_df = stat_results['enhanced_metabolites']
                    if isinstance(complete_df, pd.DataFrame) and not complete_df.empty:
                        logger.info(f"✅ Found enhanced_metabolites in memory_store: {len(complete_df)} rows")
                        return complete_df.copy()
        
        # Fallback to tab attribute (backward compatibility)
        if not hasattr(self, 'statistical_test_results') or not self.statistical_test_results:
            logger.warning("⚠️ No statistical_test_results found in memory_store or tab!")
            return None
        
        # Use enhanced_metabolites which contains complete results with stats
        if 'enhanced_metabolites' in self.statistical_test_results:
            complete_df = self.statistical_test_results['enhanced_metabolites'].copy()
            logger.info(f"✅ Found enhanced_metabolites (tab attribute): {len(complete_df)} rows")
            # DEBUG: log column diagnostic snapshot before any further cleaning
            try:
                col_list = list(complete_df.columns)
                pairwise_cols = [c for c in col_list if '_vs_' in c]
                padj_cols = [c for c in col_list if c.endswith('_p_adj') or c.endswith('_adj')]
                logger.info(f"🔎 DEBUG (enhanced_metabolites) total_cols={len(col_list)} pairwise={len(pairwise_cols)} padj={len(padj_cols)}")
                if pairwise_cols:
                    logger.info(f"🧪 DEBUG sample pairwise cols: {pairwise_cols[:8]}")
                if padj_cols:
                    logger.info(f"🧪 DEBUG sample padj cols: {padj_cols[:8]}")
            except Exception:
                pass
            if not complete_df.empty:
                return complete_df
        else:
            logger.warning("⚠️ 'enhanced_metabolites' key not found in statistical_test_results!")
        
        return None

    def get_current_groups(self):
        """Get current group definitions."""
        # First check visualization-specific groups
        if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            return self.ordered_groups(list(dict.fromkeys(self.viz_group_mapping.values())))
        
        # Fallback to statistics tab groups
        if hasattr(self, 'group_definitions') and self.group_definitions:
            return list(self.group_definitions.values())
        return []

    def get_sample_mapping(self):
        """Get sample column to group mapping."""
        from main_script.metabolites_visualization import identify_sample_columns
        
        # First check if we have visualization-specific group mapping
        if hasattr(self, 'viz_group_mapping') and self.viz_group_mapping:
            complete_df = self.build_complete_results_df()
            if complete_df is None:
                return [], {}
            
            # Use visualization group mapping directly
            sample_cols = list(self.viz_group_mapping.keys())
            sample_to_group = self.viz_group_mapping.copy()
            
            # Filter to columns that actually exist in data
            existing_cols = [col for col in sample_cols if col in complete_df.columns]
            sample_to_group = {col: group for col, group in sample_to_group.items() if col in existing_cols}
            
            return existing_cols, sample_to_group
        
        # Try to get sample mapping from memory_store (populated by Statistics tab)
        if hasattr(self, 'memory_store') and self.memory_store and 'sample_to_group' in self.memory_store:
            sample_to_group = self.memory_store['sample_to_group']
            if sample_to_group:
                complete_df = self.build_complete_results_df()
                if complete_df is not None:
                    sample_cols = [col for col in sample_to_group.keys() if col in complete_df.columns]
                    return sample_cols, {col: sample_to_group[col] for col in sample_cols}
        
        # Fallback to statistics tab group assignments
        if not hasattr(self, 'sample_group_vars') or not self.sample_group_vars:
            return [], {}
        
        # Get current complete results
        complete_df = self.build_complete_results_df()
        if complete_df is None:
            return [], {}
        
        # Parse current group assignments from statistics tab
        group_mapping = self._parse_group_assignments()
        
        # Build pattern map from current assignments
        groups = self.get_current_groups()
        pattern_map = {}
        
        # Create reverse mapping: group -> sample patterns
        for sample, group in group_mapping.items():
            if group not in pattern_map:
                pattern_map[group] = []
            # Use the sample name itself as the pattern
            pattern_map[group].append(sample)
        
        # Use the identification function to get proper sample columns
        sample_cols, sample_to_group = identify_sample_columns(complete_df, groups, pattern_map)
        
        return sample_cols, sample_to_group

    def build_color_map(self):
        """Build color mapping for groups using custom colors if available."""
        groups = self.get_current_groups()
        if not groups:
            return {}
        
        # Use visualization custom colors if available
        color_map = {}
        if hasattr(self, 'viz_color_map') and self.viz_color_map:
            for group in groups:
                if group in self.viz_color_map:
                    color_map[group] = self.viz_color_map[group]
                else:
                    # Fallback to default
                    palette = sns.color_palette('tab10')
                    idx = list(groups).index(group)
                    color_map[group] = palette[idx % len(palette)]
        elif hasattr(self, 'viz_custom_colors'):
            # Legacy support for old viz_custom_colors
            for group in groups:
                if group in self.viz_custom_colors:
                    # Convert hex to RGB tuple for matplotlib
                    hex_color = self.viz_custom_colors[group]
                    if hex_color.startswith('#'):
                        hex_color = hex_color[1:]
                    rgb = tuple(int(hex_color[i:i+2], 16)/255.0 for i in (0, 2, 4))
                    color_map[group] = rgb
                else:
                    palette = sns.color_palette('tab10')
                    idx = list(groups).index(group)
                    color_map[group] = palette[idx % len(palette)]
        else:
            # Default palette
            palette = sns.color_palette('tab10')
            color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups)}
        
        return color_map
    # Normalization methods removed - this tab only handles visualization of already-normalized data from Statistics tab

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
        self.stats_log_lines = getattr(self, 'stats_log_lines', 23)

        # Ensure stats config change handler exists before widgets bind to it
        if not hasattr(self, '_stats_config_changed'):
            def _stats_config_changed(log: str | None = None):
                # Minimal placeholder; full version defined later when pattern window used
                if log and hasattr(self, 'stats_log'):
                    ts = time.strftime('%H:%M:%S')
                    self.stats_log.insert(tk.END, f"[{ts}] {log}\n")
                    self.stats_log.see(tk.END)
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

        # Create Statistics tab
        self.stats_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.stats_tab, text='📈 Statistics')

        # Main container with scrollbar
        main_container = tk.Frame(self.stats_tab, bg='#f0f0f0')
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
        self.stats_canvas.bind('<MouseWheel>', _on_mousewheel)
        scrollable_frame.bind('<MouseWheel>', _on_mousewheel)

        container = tk.Frame(scrollable_frame, bg='#f0f0f0')
        container.pack(fill='both', expand=True, padx=10, pady=10)
        tk.Label(container, text='📊 Statistics & Normalization', font=('Arial', 16, 'bold'), bg='#f0f0f0').pack(pady=(0, 10))

        body = tk.Frame(container, bg='#f0f0f0')
        body.pack(fill='both', expand=True)
        # keep a reference for runtime layout updates
        self.stats_body = body

        # Column layout: allow columns to expand freely without minimum size constraints
        # This allows the content to fill the screen properly and avoids blank space
        body.grid_columnconfigure(0, weight=1)  # Left column (Configuration)
        body.grid_columnconfigure(1, weight=3)  # Middle column (Group Assignments)
        body.grid_columnconfigure(2, weight=2)  # Right column (Log)
        # Set minimum height to ensure configuration content is visible
        body.grid_rowconfigure(0, weight=1, minsize=800)

        # LEFT COLUMN - Completely rebuilt with proper scrolling
        cfg = tk.LabelFrame(body, text='⚙️ Configuration', bg='#f0f0f0')
        cfg.grid(row=0, column=0, sticky='nsew', padx=(0, 5))
        
        # Canvas with scrollbar for left column
        cfg_canvas = tk.Canvas(cfg, bg='#f0f0f0', highlightthickness=0)
        cfg_scrollbar = ttk.Scrollbar(cfg, orient="vertical", command=cfg_canvas.yview)
        cfg_scrollable_frame = tk.Frame(cfg_canvas, bg='#f0f0f0')
        
        cfg_scrollable_frame.bind(
            "<Configure>",
            lambda e: cfg_canvas.configure(scrollregion=cfg_canvas.bbox("all"))
        )
        
        cfg_canvas_window = cfg_canvas.create_window((0, 0), window=cfg_scrollable_frame, anchor="nw")
        cfg_canvas.configure(yscrollcommand=cfg_scrollbar.set)
        
        def configure_cfg_scroll(event):
            cfg_canvas.configure(scrollregion=cfg_canvas.bbox("all"))
            cfg_canvas.itemconfig(cfg_canvas_window, width=event.width)
        
        cfg_canvas.bind('<Configure>', configure_cfg_scroll)
        
        # Pack scrollbar first (right side), then canvas (fills remaining space)
        cfg_scrollbar.pack(side="right", fill="y")
        cfg_canvas.pack(side="left", fill="both", expand=True)
        
        def _on_cfg_mousewheel(event):
            cfg_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        cfg_canvas.bind("<MouseWheel>", _on_cfg_mousewheel)
        cfg_scrollable_frame.bind("<MouseWheel>", _on_cfg_mousewheel)
        
        # ========== IMPORT SECTION (First) ==========
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        
        import_frame = tk.LabelFrame(cfg_scrollable_frame, text='📂 Data Import', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        import_frame.pack(fill='x', padx=5, pady=(5, 10))
        
        # Data Mode Selection (Metabolite vs Lipid)
        mode_frame = tk.Frame(import_frame, bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=5, pady=(5, 5))
        tk.Label(mode_frame, text='Data Mode:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(side='left', padx=(0, 10))
        self.statistics_data_mode = tk.StringVar(value='metabolite')
        self.statistics_data_mode.trace_add('write', lambda *a: self._stats_config_changed(log=f"Data mode: {self.statistics_data_mode.get()}"))
        tk.Radiobutton(mode_frame, text='Metabolite', variable=self.statistics_data_mode,
                       value='metabolite', bg='#f0f0f0', command=self.on_statistics_mode_change).pack(side='left', padx=5)
        tk.Radiobutton(mode_frame, text='Lipid', variable=self.statistics_data_mode,
                       value='lipid', bg='#f0f0f0', command=self.on_statistics_mode_change).pack(side='left', padx=5)
        
        tk.Button(import_frame, text='📂 Import Excel Data (Pos/Neg)', command=self.import_statistics_excel,
                  bg='#16a085', fg='white', **btn_style).pack(fill='x', padx=5, pady=5)
        
        # ========== STATISTICS SECTION ==========
        # Note: Normalization is done in the Statistics tab. This tab imports normalized data for visualization.
        tests = tk.LabelFrame(cfg_scrollable_frame, text='📊 Statistical Tests', bg='#f0f0f0', font=('Arial', 10, 'bold'))
        tests.pack(fill='x', padx=5, pady=(0, 10))
        
        # Test type selection
        self.stat_test_type = tk.StringVar(value='overall')
        
        overall_frame = tk.Frame(tests, bg='#f0f0f0')
        overall_frame.pack(fill='x', padx=5, pady=(5, 2))
        tk.Radiobutton(overall_frame, text='Overall (>2 groups):', variable=self.stat_test_type,
                       value='overall', bg='#f0f0f0', command=self.on_test_type_change).pack(side='left')
        
        self.stat_overall_test = tk.StringVar(value='anova')
        self.stat_overall_test.trace_add('write', lambda *a: self._stats_config_changed(log=f"Overall test: {self.stat_overall_test.get()}"))
        self.overall_combo = ttk.Combobox(tests, values=['anova', 'kruskal'], textvariable=self.stat_overall_test, state='readonly')
        self.overall_combo.pack(fill='x', padx=15, pady=(0, 5))
        
        pairwise_frame = tk.Frame(tests, bg='#f0f0f0')
        pairwise_frame.pack(fill='x', padx=5, pady=2)
        tk.Radiobutton(pairwise_frame, text='Pairwise:', variable=self.stat_test_type,
                       value='pairwise', bg='#f0f0f0', command=self.on_test_type_change).pack(side='left')
        
        self.stat_pairwise_test = tk.StringVar(value='welch')
        self.stat_pairwise_test.trace_add('write', lambda *a: self._stats_config_changed(log=f"Pairwise test: {self.stat_pairwise_test.get()}"))
        self.pairwise_combo = ttk.Combobox(tests, values=['welch', 'mannwhitney', 'rots', 'limma'], textvariable=self.stat_pairwise_test, state='readonly')
        self.pairwise_combo.pack(fill='x', padx=15, pady=(0, 5))
        self.pairwise_combo.config(state='disabled')
        
        # Base Group
        base_frame = tk.Frame(tests, bg='#f0f0f0')
        base_frame.pack(fill='x', padx=5, pady=(8, 2))
        tk.Label(base_frame, text='Base Group (optional):', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w')
        tk.Label(base_frame, text='Compare all groups vs this base only', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w')
        self.stat_base_group = tk.StringVar(value='')
        self.base_group_combo = ttk.Combobox(base_frame, values=[''], textvariable=self.stat_base_group, state='readonly')
        self.base_group_combo.pack(fill='x', pady=(2, 5))
        self.base_group_combo.bind('<<ComboboxSelected>>', lambda e: self._stats_config_changed(log=f"Base group: {self.stat_base_group.get() or '[None]'}"))
        
        # Custom Comparisons
        custom_comp_frame = tk.Frame(tests, bg='#f0f0f0')
        custom_comp_frame.pack(fill='x', padx=5, pady=(8, 2))
        tk.Label(custom_comp_frame, text='Custom Comparisons (optional):', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w')
        tk.Label(custom_comp_frame, text='E.g., "Group1-Group2,Group3-Group4"', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w')
        self.custom_comparisons_var = tk.StringVar(value='')
        self.custom_comparisons_var.trace_add('write', lambda *a: self._stats_config_changed())
        tk.Entry(custom_comp_frame, textvariable=self.custom_comparisons_var, font=('Arial', 9)).pack(fill='x', pady=(2, 2))
        tk.Label(custom_comp_frame, text='Leave empty for all pairwise comparisons', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', pady=(0, 5))
        
        # FDR Scope
        fdr_frame = tk.Frame(tests, bg='#f0f0f0')
        fdr_frame.pack(fill='x', padx=5, pady=(6, 2))
        tk.Label(fdr_frame, text='FDR Scope:', bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w', pady=(0, 2))
        self.fdr_scope_var = tk.StringVar(value='per-comparison')
        self.fdr_scope_var.trace_add('write', lambda *a: self._on_fdr_scope_changed())
        fdr_radio_frame = tk.Frame(fdr_frame, bg='#f0f0f0')
        fdr_radio_frame.pack(fill='x', pady=(0, 5))
        for txt, val in [('Per-Comparison', 'per-comparison'), ('Per-Metabolite', 'per-metabolite')]:
            tk.Radiobutton(fdr_radio_frame, text=txt, variable=self.fdr_scope_var, value=val, bg='#f0f0f0').pack(side='left', padx=2)
        
        # Add help text for FDR scope
        fdr_help = tk.Label(fdr_frame, text='⚠️ Only use when comparing 3+ groups \n (returns identical p-values with 2 groups)', 
                           bg='#fff3cd', fg='#856404', font=('Arial', 8), wraplength=450, justify='left', padx=5, pady=3)
        # Don't pack initially - will be shown only when per-metabolite is selected
        self.fdr_scope_warning = fdr_help  # Store reference for dynamic updates
        
        # Group Order (optional)
        group_order_frame = tk.Frame(tests, bg='#f0f0f0')
        group_order_frame.pack(fill='x', padx=5, pady=(4, 2))
        tk.Label(group_order_frame, text='Group Order (optional):', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w')
        tk.Label(group_order_frame, text='Comma separated labels e.g. PC3_2D,PC3_3D,DU145_2D', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w')
        self.statistics_group_order_var = tk.StringVar(value='')
        self.statistics_group_order_var.trace_add('write', lambda *a: self._stats_config_changed(log=f"Group order updated"))
        tk.Entry(group_order_frame, textvariable=self.statistics_group_order_var, font=('Arial', 9)).pack(fill='x', pady=(2, 2))
        tk.Label(group_order_frame, text='If empty, order follows group definitions listing.', bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', pady=(0, 4))

        # Alpha and adjusted p-values
        alpha_frame = tk.Frame(tests, bg='#f0f0f0')
        alpha_frame.pack(fill='x', padx=5, pady=(6, 5))
        tk.Label(alpha_frame, text='α:', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
        self.alpha_var = tk.StringVar(value='0.05')
        self.alpha_var.trace_add('write', lambda *a: self._stats_config_changed())
        tk.Entry(alpha_frame, textvariable=self.alpha_var, width=8, font=('Arial', 9)).pack(side='left', padx=5)
        
        self.use_adj_p_var = tk.BooleanVar(value=True)
        self.use_adj_p_var.trace_add('write', lambda *a: self._stats_config_changed())
        tk.Checkbutton(alpha_frame, text='adjust p-values', variable=self.use_adj_p_var,
                       bg='#f0f0f0', font=('Arial', 9)).pack(side='left', padx=(10, 0))

        # MIDDLE COLUMN - Completely rebuilt with proper height distribution
        grp_mgmt = tk.LabelFrame(body, text='👥 Group Management', bg='#f0f0f0')
        grp_mgmt.grid(row=0, column=1, sticky='nsew', padx=5)
        self.stats_grp_mgmt = grp_mgmt
        
        # Configure grid: Group IDs gets 25% height, Assignment gets 75% height
        # Remove minsize constraints to allow proper expansion
        grp_mgmt.grid_rowconfigure(0, weight=1)  # Group IDs - smaller, compact
        grp_mgmt.grid_rowconfigure(1, weight=3)  # Assignment - larger, takes most space
        grp_mgmt.grid_columnconfigure(0, weight=1)
        
        # ===== TOP SECTION: Group IDs & Labels (25% of space) =====
        group_ids_frame = tk.LabelFrame(grp_mgmt, text='Group IDs & Labels', bg='#f0f0f0')
        group_ids_frame.grid(row=0, column=0, sticky='nsew', padx=5, pady=(5, 2))
        
        # Initialize group data - use defaults only if no saved config exists
        if not hasattr(self, 'group_definitions') or not self.group_definitions:
            self.group_definitions = {'Group1': 'Control', 'Group2': 'Disease', 'Group3': 'Treatment', 'Group4': 'Other'}
            self.group_count = 4
        if hasattr(self, 'base_group_combo'):
            self.base_group_combo['values'] = [''] + list(self.group_definitions.keys())
        
        # Canvas with scrollbar for group IDs
        self.groups_canvas = tk.Canvas(group_ids_frame, bg='#f0f0f0', highlightthickness=0, height=80)
        groups_scrollbar = ttk.Scrollbar(group_ids_frame, orient='vertical', command=self.groups_canvas.yview)
        self.groups_scrollable_frame = tk.Frame(self.groups_canvas, bg='#f0f0f0')
        
        self.groups_scrollable_frame.bind(
            "<Configure>",
            lambda e: self.groups_canvas.configure(scrollregion=self.groups_canvas.bbox("all"))
        )
        
        groups_canvas_window = self.groups_canvas.create_window((0, 0), window=self.groups_scrollable_frame, anchor='nw')
        self.groups_canvas.configure(yscrollcommand=groups_scrollbar.set)
        
        def configure_groups_scroll(event):
            self.groups_canvas.configure(scrollregion=self.groups_canvas.bbox('all'))
            self.groups_canvas.itemconfig(groups_canvas_window, width=event.width)
        
        self.groups_canvas.bind('<Configure>', configure_groups_scroll)
        
        self.groups_canvas.pack(side='left', fill='both', expand=True, padx=(5, 0), pady=5)
        groups_scrollbar.pack(side='right', fill='y')
        
        def _on_groups_mousewheel(event):
            self.groups_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        self.groups_canvas.bind('<MouseWheel>', _on_groups_mousewheel)
        self.groups_scrollable_frame.bind('<MouseWheel>', _on_groups_mousewheel)
        
        # Refresh group UI
        self.group_id_vars = {}
        self.refresh_group_ui()
        
        # Buttons for adding/removing groups
        group_buttons_frame = tk.Frame(group_ids_frame, bg='#f0f0f0')
        group_buttons_frame.pack(fill='x', padx=5, pady=(0, 5))
        button_opts = {'font': ('Arial', 9, 'bold'), 'height': 1, 'padx': 10, 'pady': 4}
        tk.Button(group_buttons_frame, text='+ Add Group', command=self.add_group,
                  bg='#27ae60', fg='white', width=60, **button_opts).pack(fill='x', padx=4, pady=2)
        tk.Button(group_buttons_frame, text='- Remove Group', command=self.remove_group,
                  bg='#e74c3c', fg='white', width=60, **button_opts).pack(fill='x', padx=4, pady=2)
        
        # ===== BOTTOM SECTION: Sample Assignment (75% of space) =====
        assign_frame = tk.LabelFrame(grp_mgmt, text='Sample → Group Assignment', bg='#f0f0f0')
        assign_frame.grid(row=1, column=0, sticky='nsew', padx=5, pady=(2, 5))
        
        # Header row
        header_frame = tk.Frame(assign_frame, bg='#f0f0f0')
        header_frame.pack(fill='x', padx=5, pady=(5, 2))
        tk.Label(header_frame, text='Sample Column', bg='#f0f0f0', font=('Arial', 9, 'bold'), width=20).pack(side='left')
        tk.Label(header_frame, text='→', bg='#f0f0f0', font=('Arial', 9, 'bold'), width=3).pack(side='left')
        tk.Label(header_frame, text='Group', bg='#f0f0f0', font=('Arial', 9, 'bold'), width=12).pack(side='left')
        
        # Canvas with scrollbar for assignments
        assign_canvas = tk.Canvas(assign_frame, bg='#f0f0f0', highlightthickness=0)
        assign_scrollbar = ttk.Scrollbar(assign_frame, orient='vertical', command=assign_canvas.yview)
        self.assign_scrollable_frame = tk.Frame(assign_canvas, bg='#f0f0f0')
        
        self.assign_scrollable_frame.bind(
            "<Configure>",
            lambda e: assign_canvas.configure(scrollregion=assign_canvas.bbox("all"))
        )
        
        canvas_window = assign_canvas.create_window((0, 0), window=self.assign_scrollable_frame, anchor='nw')
        assign_canvas.configure(yscrollcommand=assign_scrollbar.set)
        
        def configure_assign_scroll(event):
            assign_canvas.configure(scrollregion=assign_canvas.bbox('all'))
            assign_canvas.itemconfig(canvas_window, width=event.width)
        
        assign_canvas.bind('<Configure>', configure_assign_scroll)
        
        assign_canvas.pack(side='left', fill='both', expand=True, padx=(5, 0), pady=5)
        assign_scrollbar.pack(side='right', fill='y')
        
        def _on_assign_mousewheel(event):
            assign_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        assign_canvas.bind('<MouseWheel>', _on_assign_mousewheel)
        self.assign_scrollable_frame.bind('<MouseWheel>', _on_assign_mousewheel)
        
        self.assign_canvas = assign_canvas
        self.assign_canvas_window = canvas_window
        self.sample_group_vars = {}
        
        # Initial placeholder message
        initial_msg = tk.Label(self.assign_scrollable_frame,
                               text='Complete data normalization first to populate sample columns',
                               bg='#f0f0f0', fg='#7f8c8d', font=('Arial', 10, 'italic'))
        initial_msg.pack(pady=20)
        
        # Control buttons at bottom
        ctrl_frame = tk.Frame(assign_frame, bg='#f0f0f0')
        ctrl_frame.pack(fill='x', padx=5, pady=(0, 5))
        CTRL_BTN_WIDTH = 80
        tk.Button(ctrl_frame, text='🔄 Refresh Groups', command=self.refresh_group_assignments,
                  bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                  width=CTRL_BTN_WIDTH, padx=12, pady=5).pack(fill='x', padx=2, pady=2)
        tk.Button(ctrl_frame, text='🎯 Auto-Assign by Pattern', command=self.auto_assign_groups,
                  bg='#e67e22', fg='white', font=('Arial', 9, 'bold'),
                  width=CTRL_BTN_WIDTH, padx=12, pady=5).pack(fill='x', padx=2, pady=2)
        
        # Data Cleanup section (renamed from Pre-Normalization Filtering)
        min_group_section = tk.LabelFrame(assign_frame, text='🧹 Data Cleanup', 
                                         bg='#f0f0f0', font=('Arial', 9, 'bold'))
        min_group_section.pack(fill='x', padx=5, pady=(10, 5))
        
        # Filter timing option (when to apply minimum sample filtering)
        if not hasattr(self, 'filter_timing_var'):
            self.filter_timing_var = tk.StringVar(value='before')  # Default: filter before normalization
        
        timing_frame = tk.Frame(min_group_section, bg='#f0f0f0')
        timing_frame.pack(fill='x', padx=5, pady=(5, 8))

        # Place the label at the top, radios stacked vertically underneath (works on narrow screens)
        tk.Label(timing_frame, text='Apply filtering:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w', padx=(0, 10))

        radio_stack = tk.Frame(timing_frame, bg='#f0f0f0')
        radio_stack.pack(fill='x', padx=10, pady=(4, 0))

        tk.Radiobutton(radio_stack, text='Before Normalization',
                       variable=self.filter_timing_var, value='before',
                       bg='#f0f0f0', font=('Arial', 9),
                       command=self._on_filter_timing_change).pack(anchor='w', pady=2)
        tk.Radiobutton(radio_stack, text='After Normalization',
                       variable=self.filter_timing_var, value='after',
                       bg='#f0f0f0', font=('Arial', 9),
                       command=self._on_filter_timing_change).pack(anchor='w', pady=2)
        
        # Explanation text that updates based on selection
        self.filter_timing_explanation = tk.Label(min_group_section, 
                text="", bg='#f0f0f0', font=('Arial', 8), fg='#555', 
                wraplength=500, justify='left')
        self.filter_timing_explanation.pack(anchor='w', padx=5, pady=(0, 8))
        
        # Update explanation text initially
        self._update_filter_timing_explanation()
        
        # Initialize variables first (defaults)
        if not hasattr(self, 'min_samples_type_var'):
            # Default to Percentage (user preference)
            self.min_samples_type_var = tk.StringVar(value='percentage')   # 'absolute' or 'percentage'
        if not hasattr(self, 'min_samples_per_group_var'):
            # sensible absolute default (still available if user switches)
            self.min_samples_per_group_var = tk.IntVar(value=2)
        else:
            try:
                current_min = int(self.min_samples_per_group_var.get())
                self.min_samples_per_group_var.set(max(1, current_min))
            except Exception:
                self.min_samples_per_group_var.set(2)
        if not hasattr(self, 'min_samples_percent_var'):
            # sensible percent default
            self.min_samples_percent_var = tk.DoubleVar(value=70.0)
        else:
            try:
                pct = float(self.min_samples_percent_var.get())
                self.min_samples_percent_var.set(max(0.0, min(100.0, pct)))
            except Exception:
                self.min_samples_percent_var.set(50.0)


        def _log_min_group_threshold(*_args):
            try:
                if self.min_samples_type_var.get() == 'absolute':
                    threshold = int(self.min_samples_per_group_var.get())
                    threshold = max(1, threshold)
                    if self.min_samples_per_group_var.get() != threshold:
                        self.min_samples_per_group_var.set(threshold)
                    self._stats_config_changed(log=f"Minimum per-group samples: {threshold} (absolute)")
                else:
                    percent = float(self.min_samples_percent_var.get())
                    percent = max(0, min(100, percent))
                    if self.min_samples_percent_var.get() != percent:
                        self.min_samples_percent_var.set(percent)
                    self._stats_config_changed(log=f"Minimum per-group samples: {percent}% (percentage)")
            except Exception:
                pass

        self.min_samples_per_group_var.trace_add('write', _log_min_group_threshold)
        self.min_samples_type_var.trace_add('write', _log_min_group_threshold)
        self.min_samples_percent_var.trace_add('write', _log_min_group_threshold)
        
        # Threshold Type header
        tk.Label(min_group_section, text='Threshold Type:', bg='#f0f0f0', font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(0, 2))
        
        # Radio buttons stacked vertically
        radio_frame = tk.Frame(min_group_section, bg='#f0f0f0')
        radio_frame.pack(fill='x', padx=10, pady=(0, 5))
        
        # Absolute Count radio + spinbox
        abs_frame = tk.Frame(radio_frame, bg='#f0f0f0')
        abs_frame.pack(fill='x', pady=2)
        tk.Radiobutton(abs_frame, text='Absolute Count:', variable=self.min_samples_type_var, value='absolute',
                      bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w')
        count_input_frame = tk.Frame(abs_frame, bg='#f0f0f0')
        count_input_frame.pack(fill='x', padx=(20, 0), pady=(2, 0))
        min_spin = ttk.Spinbox(count_input_frame, from_=1, to=50, textvariable=self.min_samples_per_group_var, width=8)
        min_spin.pack(side='left', padx=(0, 5))
        tk.Label(count_input_frame, text='samples', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
        self._create_tooltip(min_spin, "Minimum number of replicates per group (applied per-metabolite before normalization)")
        
        # Percentage radio + spinbox
        pct_frame = tk.Frame(radio_frame, bg='#f0f0f0')
        pct_frame.pack(fill='x', pady=2)
        tk.Radiobutton(pct_frame, text='Percentage:', variable=self.min_samples_type_var, value='percentage',
                      bg='#f0f0f0', font=('Arial', 9)).pack(anchor='w')
        percent_input_frame = tk.Frame(pct_frame, bg='#f0f0f0')
        percent_input_frame.pack(fill='x', padx=(20, 0), pady=(2, 0))
        percent_spin = ttk.Spinbox(percent_input_frame, from_=0, to=100, increment=5, 
                                   textvariable=self.min_samples_percent_var, width=8)
        percent_spin.pack(side='left', padx=(0, 5))
        tk.Label(percent_input_frame, text='% group n', bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
        self._create_tooltip(percent_spin, "Minimum percentage of group samples required (applied per-metabolite before normalization)")
        
        # Update scroll regions
        try:
            self._update_assignment_scroll_region()
            self._update_main_scroll_region()
        except Exception:
            pass
        # RIGHT COLUMN
        results = tk.LabelFrame(body, text='📊 Statistics Log', bg='#f0f0f0')
        results.grid(row=0, column=2, sticky='nsew', padx=(5, 0))
        
        # Configure grid: Fixed rows for buttons/progress (weight=0), expandable row for log (weight=1)
        results.grid_rowconfigure(0, weight=0)  # Check Data button - fixed
        results.grid_rowconfigure(1, weight=0)  # Action buttons - fixed
        results.grid_rowconfigure(2, weight=0)  # Progress label - fixed
        results.grid_rowconfigure(3, weight=0)  # Progress bar - fixed
        results.grid_rowconfigure(4, weight=1)  # Log - expandable
        results.grid_columnconfigure(0, weight=1)
        
        # Add Check Data Availability button at top
        btn_style = {'font': ('Arial', 9, 'bold'), 'relief': 'raised', 'bd': 2, 'pady': 3}
        tk.Button(results, text='🔍 Check Available Data', command=self.check_memory_store_status,
                  bg='#34495e', fg='white', **btn_style).grid(row=0, column=0, sticky='ew', padx=5, pady=(5, 5))
        
        # ========== STATISTICAL ANALYSIS ACTIONS (at the top, side by side) ==========
        action_buttons_frame = tk.Frame(results, bg='#f0f0f0')
        action_buttons_frame.grid(row=1, column=0, sticky='ew', padx=5, pady=(0, 5))
        action_buttons_frame.grid_columnconfigure(0, weight=1)
        action_buttons_frame.grid_columnconfigure(1, weight=1)
        
        tk.Button(action_buttons_frame, text='🧪 Run Statistical Tests', command=self.run_statistical_tests,
                  bg='#e74c3c', fg='white', **btn_style).grid(row=0, column=0, sticky='ew', padx=(0, 3))
        tk.Button(action_buttons_frame, text='📈 Export Statistical Results', command=self.export_statistical_results,
                  bg='#f39c12', fg='white', **btn_style).grid(row=0, column=1, sticky='ew', padx=(3, 0))
        
        # Add progress bar and label above the log
        self.stats_progress_label = tk.Label(results, text="", bg='#f0f0f0', font=('Arial', 9))
        self.stats_progress_label.grid(row=2, column=0, sticky='ew', padx=5, pady=(0, 0))
        self.stats_progress_label.grid_remove()  # Hide initially
        
        self.stats_progress = ttk.Progressbar(results, mode='indeterminate', length=400)
        self.stats_progress.grid(row=3, column=0, sticky='ew', padx=5, pady=(2, 2))
        self.stats_progress.grid_remove()  # Hide initially
        
        # Statistics log - expands to fill remaining space
        self.stats_log = scrolledtext.ScrolledText(results, font=('Courier', 9), wrap=tk.WORD)
        self.stats_log.grid(row=4, column=0, sticky='nsew', padx=5, pady=(2, 5))
        self.stats_log.insert(tk.END, '📊 Statistics Log Ready\n')
        self.stats_log.insert(tk.END, 'Ready for statistical analysis operations...\n\n')
        
        try:
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
        
        # Initialize flag to track if statistics config has been loaded
        self._stats_config_loaded = False
        
        # Note: Config will be auto-loaded when the tab is first opened (see on_tab_changed)
        # This prevents loading the config during initialization when widgets aren't ready
        self.stats_log.insert(tk.END, 'Statistics configuration will be auto-loaded when tab is opened.\n')

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
            self.stats_progress_label.grid()
            self.stats_progress.grid()
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
            self.stats_progress.grid_remove()
            self.stats_progress_label.grid_remove()
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
        """Handle data mode change (Metabolite <-> Lipid) in Statistics tab."""
        mode = self.statistics_data_mode.get()
        if hasattr(self, 'stats_log'):
            self.stats_log.insert(tk.END, f"\n🔄 Data mode changed to: {mode.upper()}\n")
            self.stats_log.see(tk.END)

    def import_statistics_excel(self):
        """Load Excel file with Pos/Neg sheets for Statistics tab based on selected data mode."""
        try:
            mode = self.statistics_data_mode.get() if hasattr(self, 'statistics_data_mode') else 'metabolite'
            
            title = 'Select Lipid Data Excel' if mode == 'lipid' else 'Select ID Annotated Excel'
            path = filedialog.askopenfilename(title=title, filetypes=[('Excel Files','*.xlsx *.xls')])
            if not path:
                return
            if not hasattr(self, 'stats_log'):
                messagebox.showerror('Error', 'Statistics log not initialized.')
                return
            
            self.stats_log.insert(tk.END, f"\n===== Importing {mode.upper()} Excel =====\n{os.path.basename(path)}\n")
            self.stats_log.see(tk.END)
            
            # Show progress bar
            self.show_stats_progress(f"Importing {mode} Excel file...")
            
            xl = pd.ExcelFile(path)
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
                for candidate in ['Positive_Lipid_Class', 'Pos_Lipid_Class', 'Positive_Class']:
                    if candidate in xl.sheet_names:
                        pos_class_sheet = candidate
                        break
                for candidate in ['Negative_Lipid_Class', 'Neg_Lipid_Class', 'Negative_Class']:
                    if candidate in xl.sheet_names:
                        neg_class_sheet = candidate
                        break
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
                for label, df in [('Positive Lipids', pos_df), ('Negative Lipids', neg_df)]:
                    if df is None:
                        continue
                    try:
                        # Use robust feature detection
                        sample_cols = []
                        for col in df.columns:
                            # treat as feature if normalized name matches canonical lipid features
                            if self._is_lipid_feature_col(col):
                                continue
                            # numeric columns not identified as features are treated as sample intensity columns
                            if pd.api.types.is_numeric_dtype(df[col]):
                                sample_cols.append(col)
                        self.stats_log.insert(tk.END, f'{label}: detected {len(sample_cols)} sample cols.\n')
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
                
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
                            if pd.api.types.is_numeric_dtype(df[col]):
                                sample_cols.append(col)
                        self.stats_log.insert(tk.END, f'{label}: detected {len(sample_cols)} sample cols.\n')
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
            else:
                # Use existing metabolite detection
                from main_script.metabolite_statistics_analysis import detect_feature_and_sample_columns
                for label, df in [('Positive', pos_df), ('Negative', neg_df)]:
                    if df is None:
                        continue
                    try:
                        feature_cols, sample_cols = detect_feature_and_sample_columns(df)
                        self.stats_log.insert(tk.END, f'{label}: detected {len(feature_cols)} feature cols, {len(sample_cols)} sample cols.\n')
                        for c in sample_cols:
                            if c not in seen:
                                seen.add(c)
                                union_sample_cols.append(c)
                    except Exception as e:
                        self.stats_log.insert(tk.END, f'{label}: detection error {e}\n')
            
            if union_sample_cols:
                # Ensure listbox exists
                if not hasattr(self, 'sample_cols_list'):
                    self.sample_cols_list = tk.Listbox(self.assign_scrollable_frame)
                self.sample_cols_list.delete(0, tk.END)
                for c in union_sample_cols:
                    self.sample_cols_list.insert(tk.END, c)
                self.populate_sample_assignments(union_sample_cols)
                self.stats_log.insert(tk.END, f'Union sample columns loaded: {len(union_sample_cols)} columns.\n')
            else:
                self.stats_log.insert(tk.END, 'No sample columns auto-detected from imported file.\n')
            
            self.stats_log.insert(tk.END, 'Ready for normalization. Configure groups and click "Run Normalization & Merge".\n')
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
            self.stats_log.insert(tk.END, f'3. Click "Run Normalization & Merge" (no merging will occur)\n')
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
                                'subclass', 'lipid', 'lipidgroup']
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
        """Update the scroll region for the sample assignment panel."""
        # Support both legacy and current attribute names ('assign_canvas' used
        # in newer code). Prefer current attributes: self.assign_canvas and
        # self.assign_scrollable_frame.
        canvas_attr = getattr(self, 'assign_canvas', None) or getattr(self, 'assignment_canvas', None)
        frame_attr = getattr(self, 'assign_scrollable_frame', None) or getattr(self, 'assignment_scrollable_frame', None)
        if canvas_attr is None or frame_attr is None:
            return
        try:
            canvas_attr.update_idletasks()
            frame_attr.update_idletasks()
            width = frame_attr.winfo_reqwidth()
            height = frame_attr.winfo_reqheight()
            canvas_attr.configure(scrollregion=(0, 0, width, height))
        except Exception:
            # Best-effort: ignore failures to avoid crashing the UI
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
                                 'Use "Run Normalization & Merge" button to process your data.')
            return
        
        # Use the normalized data for sample column detection
        self.populate_sample_assignments_from_normalized_data()
    
    # REMOVED: run_statistics_pipeline() - normalization is done in Statistics tab, not Visualization tab


