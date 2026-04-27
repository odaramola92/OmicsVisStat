"""
Help Tab - User guide and documentation for OmicsVisStat.
"""
import tkinter as tk
from tkinter import ttk, scrolledtext
import logging

from gui.shared.base_tab import BaseTab

logger = logging.getLogger(__name__)


class HelpTab(BaseTab):
    """Tab for help and documentation."""

    def __init__(self, parent, data_manager):
        """Initialize the help tab."""
        super().__init__(parent, data_manager)
        self.setup_ui()
        print("[OK] Help Tab initialized")

    def setup_ui(self):
        """Set up the comprehensive Help tab with sub-tabs for each main tab."""
        title_frame = tk.Frame(self.frame, bg="#27ae60", height=60)
        title_frame.pack(fill="x", pady=(0, 5))
        title_frame.pack_propagate(False)

        tk.Label(
            title_frame,
            text="User Guide & Help",
            font=("Arial", 16, "bold"),
            bg="#27ae60",
            fg="white",
        ).pack(pady=15)

        content_frame = tk.Frame(self.frame)
        content_frame.pack(fill="both", expand=True, padx=5, pady=5)

        help_notebook = ttk.Notebook(content_frame)
        help_notebook.pack(fill="both", expand=True)

        help_contents = self._get_default_help_content()
        for tab_name, content in help_contents.items():
            self._create_help_subtab(help_notebook, tab_name, content)

        resources_frame = ttk.LabelFrame(self.frame, text="Additional Resources", padding=15)
        resources_frame.pack(fill="x", padx=5, pady=(0, 5))

        tk.Label(
            resources_frame,
            text="Documentation files available in the docs/ folder",
            font=("Arial", 9),
            bg="#ecf0f1",
            justify="left",
        ).pack(anchor="w")

    def _create_help_subtab(self, notebook, tab_name, content):
        """Create a help sub-tab with formatted content."""
        tab_frame = ttk.Frame(notebook)
        notebook.add(tab_frame, text=" " + tab_name + " ")

        text_area = scrolledtext.ScrolledText(
            tab_frame,
            wrap=tk.WORD,
            font=("Arial", 10),
            bg="#ffffff",
            fg="#2c3e50",
            padx=15,
            pady=15,
        )
        text_area.pack(fill="both", expand=True, padx=5, pady=5)

        text_area.tag_config("title", font=("Arial", 14, "bold"), foreground="#2c3e50")
        text_area.tag_config("heading", font=("Arial", 12, "bold"), foreground="#3498db")
        text_area.tag_config("subheading", font=("Arial", 10, "bold"), foreground="#27ae60")
        text_area.tag_config("bullet", font=("Arial", 10), foreground="#555")
        text_area.tag_config("warning", font=("Arial", 10), foreground="#e74c3c")
        text_area.tag_config("tip", font=("Arial", 10), foreground="#f39c12")

        text_area.insert("1.0", content)
        text_area.config(state="disabled")

    def _get_default_help_content(self):
        """Provide help content for OmicsVisStat."""
        return {
            "Statistics": """STATISTICS TAB

Run normalization and statistical testing for metabolite/lipid/custom datasets.

OVERVIEW:
The Statistics tab follows a guided 5-step workflow with optional covariate-adjusted analysis.

WORKFLOW:
1. Step 1: Select Data Mode
   - Metabolite, Lipid, or Custom
2. Step 2: Import and Verify Columns
   - Import Excel Pos/Neg data
   - Verify feature and sample columns
3. Step 3: Configure Groups
   - Assign samples to groups (manual or pattern-based)
   - Add/remove groups and rename labels
   - Set replicate filtering (absolute count or percentage)
   - Optional base group, custom comparisons, and group order
4. Step 4: Normalization
   - Select normalization methods as an ordered chain
   - Run normalization + normality checks
   - Optional: variability filtering, imputation, PCA-based outlier removal
5. Step 5: Statistical Tests
   - Overall and pairwise testing workflows
   - Pairwise p-value adjustment and FDR scope controls
   - Optional Two-Way ANOVA setup and non-parametric two-way settings
   - Optional ROTS parameter configuration

OPTIONAL STEP 5B:
- Covariate adjustment (for Age, Sex, BMI, etc.)

EXPORTS AND OUTPUT:
- Export normalized data
- Export statistical results
- Open output folder from the right panel

PRACTICAL TIPS:
- Verify sample names between Pos and Neg sheets before running
- Configure groups before normalization/testing
- Use consistent group labels for cleaner downstream plots
- Keep logs for reproducibility
""",

            "Machine Learning": """MACHINE LEARNING TAB

Train and evaluate classification models with configurable preprocessing, validation, and figure generation.

OVERVIEW:
The ML tab is organized into Steps 1-4 and uses the same data/group mindset as the Statistics tab.

WORKFLOW:
1. Step 1: Import and Verify Data
   - Import Excel input and verify columns
   - Set working folder for outputs
2. Step 2: Configure Groups
   - Configure group IDs/labels
   - Add/remove groups
   - Use pattern-based auto-assignment when needed
3. Step 3: ML Configuration
   - Select one or multiple models
   - Configure test size, CV folds, scaling, and class weighting
   - Configure linear-model regularization options
   - Configure robustness testing (repeats, seed, feature stability)
   - Optional advanced validation/tuning:
     hyperparameter tuning, repeated CV, nested CV, SVM calibration,
     permutation runs, imputation, and feature selection
   - Optional feature filters:
     replicate filter, endogenous-only, HMDB-present, p-value threshold
   - Optional pairwise p-value column verification for pairwise ML workflows
4. Step 4: Run and Review
   - Run analysis
   - Test model combinations
   - Review logs and generated outputs

FIGURE SUPPORT:
- Auto-generate ROC/model comparison/top-metabolite figures
- Manual figure controls and top-N settings available

PRACTICAL TIPS:
- Complete column verification before enabling p-value-based filters
- Use repeated CV/nested CV for more stable performance estimates
- Keep a fixed seed for reproducible runs
""",

            "Visualization": """VISUALIZATION TAB

Generate publication-ready plots from current-session statistics or imported statistical result files.

OVERVIEW:
This tab enforces required setup (groups + stat columns) before plot generation.

WORKFLOW:
1. Data Source
   - Use statistics from current session, or import a statistics Excel file
   - Mode selection: Metabolite, Lipid, or Custom
   - Load and analyze, then verify columns when required
2. Required Configuration
   - Configure Groups
   - Configure Stat Columns
3. Plot Selection and Generation
   - Generate selected plots or all plots
   - Stop running jobs if needed
   - Open output folder directly

PLOT PANELS AVAILABLE:
- PCA
- Volcano
- Venn
- Boxplots
- Bar Graphs
- Heatmaps
- ROC
- Custom List

COMMON CONTROLS ACROSS PANELS:
- Figure width/height/DPI
- Font sizes
- Legend/title options
- Output format options (for example PNG/SVG where supported)
- Panel-specific settings (thresholds, annotations, 2D/3D PCA options, etc.)

PRACTICAL TIPS:
- Finish both required configuration steps before plotting
- Keep output directory set before batch plot generation
- Import group definitions from Statistics config for consistent labeling
""",

            "Utility": """UTILITY TAB

Use specialized tools grouped into multiple utility sub-tabs.

SUB-TABS AVAILABLE:
- Chord Diagram
- Venn Diagram
- Pie Chart
- Heatmap
- Effect Plot
- Regression Utility
- Linear Regression Plot
- Glycan Classification

COMMON PATTERN:
Most utility tools follow this sequence:
1. Upload data
2. Map or verify required columns
3. Configure tool-specific settings
4. Generate output and export

TOOL SUMMARY:
- Chord Diagram: show relationships between compounds and classifications.
- Venn Diagram: compare overlaps between 2 to 6 sets.
- Pie Chart: compare class composition across groups.
- Heatmap: generate grouped heatmaps with normalization and color options.
- Effect Plot: plot effect sizes with confidence intervals.
- Regression Utility: run per-category outcome analysis with optional covariates.
- Linear Regression Plot: build regression plots from uploaded data.
- Glycan Classification: classify glycan entries and export processed outputs.

PRACTICAL TIPS:
- Complete mapping/verification before running generation buttons
- Use consistent naming for groups and samples across files
- If a missing-library error appears, install the required package and rerun
"""
        }
