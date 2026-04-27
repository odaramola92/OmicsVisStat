# User Guide - OmicsVisStat

## Overview

This guide reflects the current application layout:

- Statistics tab
- Machine Learning tab
- Visualization tab
- Utility tab
- Help tab

Use this guide for the core end-user workflow from data import to exports.

---

## Quick Start

1. Launch the app using run_gui.py or your packaged executable.
2. Go to Statistics tab and complete data setup.
3. Run normalization and statistical tests.
4. Export results.
5. Use Visualization tab for plot generation.
6. Use Machine Learning and Utility tabs as needed.

---

## Tab Summary

### 1. Statistics Tab

Primary analysis pipeline with guided steps.

Step 1: Select Data Mode
- Metabolite
- Lipid
- Custom

Step 2: Import and Verify Columns
- Import Excel Pos/Neg data
- Verify sample and feature columns

Step 3: Configure Groups
- Configure groups and sample assignment
- Add/remove groups
- Optional pattern-based auto-assignment
- Optional base group, custom comparisons, group order
- Replicate filtering settings (absolute or percentage)

Step 4: Normalization
- Select normalization methods in order
- Run normalization and normality testing
- Optional post-normalization filters and imputation
- Optional PCA-based sample outlier removal

Step 5: Statistical Tests
- Overall and pairwise testing modes
- Pairwise p-value adjustment options
- FDR scope controls
- Optional Two-Way ANOVA setup
- Optional non-parametric two-way settings
- Optional ROTS parameter configuration

Optional Step 5b: Covariate Adjustment
- Run statistics with covariate adjustment where needed

Outputs
- Export normalized data
- Export statistical results
- Open output folder

### 2. Machine Learning Tab

ML workflow with data verification, group setup, model configuration, and export.

High-level flow
1. Import and verify data
2. Configure groups
3. Configure ML settings
4. Run analysis and export

See MACHINE_LEARNING_GUIDE.md for full details.

### 3. Visualization Tab

Generates publication-ready plots using current-session statistics data or imported statistical results.

Required setup
1. Choose data source and mode
2. Load/analyze data
3. Verify columns
4. Configure groups
5. Configure stat columns

Plot panels available
- PCA
- Volcano
- Venn
- Boxplots
- Bar Graphs
- Heatmaps
- ROC
- Custom List

Generation controls
- Generate selected plots
- Generate all plots
- Stop current plotting job
- Open output folder

### 4. Utility Tab

Collection of dedicated tools:
- Chord Diagram
- Venn Diagram
- Pie Chart
- Heatmap
- Effect Plot
- Regression Utility
- Linear Regression Plot
- Glycan Classification

See UTILITY_GUIDE.md for per-tool steps.

### 5. Help Tab

In-app quick-reference content aligned to current tab behavior.

---

## Data Format Notes

General expectations
- Features in rows
- Samples in columns
- At least one feature identifier column
- Numeric sample columns for analysis

Recommended
- Consistent sample naming across sheets
- Consistent group naming
- Clean missing values before analysis

---

## Typical End-to-End Workflow

1. Import data and verify columns in Statistics tab.
2. Configure groups and filtering settings.
3. Run normalization and review logs.
4. Run statistical tests and export results.
5. In Visualization tab, load results/current session data and configure required mappings.
6. Generate selected/all plots and export files.
7. Optionally run Machine Learning analyses and Utility tools.

---

## Reproducibility Checklist

- Save exported normalized/statistical outputs.
- Keep group labels and comparison settings documented.
- Keep a copy of config files used during runs.
- Record software version/date in study notes.

---

## Troubleshooting

If something fails:
1. Check TROUBLESHOOTING.md
2. Verify column mapping and group assignments
3. Confirm required files and dependencies are present
4. Re-run with a smaller subset to isolate issues
