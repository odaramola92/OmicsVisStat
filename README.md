# OmicsVisStat

OmicsVisStat: Statistical and Visualization Platform for Omics Data Analysis.

## Overview

OmicsVisStat is a GUI application for statistical analysis, machine learning, and figure generation across metabolomics, lipidomics, and other omics datasets.

The current interface includes:
- Statistics
- Machine Learning
- Visualization
- Utility
- Help

## Key Features

### Statistics
- One-way ANOVA with post-hoc tests
- Two-way ANOVA with interaction effects
- Non-parametric tests
- Pairwise comparisons with multiple testing correction
- Covariate adjustment
- Normalization, filtering, and imputation workflows
- ROTS and limma-style analysis support where configured

### Machine Learning
- Classification models such as Random Forest, SVM, Gradient Boosting, and Logistic Regression
- Stratified cross-validation
- Model comparison and performance metrics
- Feature importance for biomarker discovery
- PCA and LDA workflows

### Visualization
- PCA plots
- Volcano plots
- Heatmaps
- Venn diagrams
- Boxplots and bar graphs
- ROC plots
- Custom list-based plots
- Export-ready figure controls

### Utility
- Chord diagrams
- Venn diagrams
- Pie charts
- Heatmaps
- Effect plots
- Regression utility
- Linear regression plots
- Glycan classification

## Installation

### Option 1: Standalone Executable (Windows)

1. Download the release ZIP package.
2. Extract it to your desired location.
3. Run `OmicsVisStat.exe`.

### Option 2: From Source

#### Prerequisites
- Python 3.8 or higher
- pip package manager

#### Installation Steps

1. Clone or download this repository.
2. Install dependencies:
```bash
pip install -r requirements.txt
```
3. Run the application:
```bash
python run_gui.py
```

### Building from Source

To build the executable:
```bash
pyinstaller statistics.spec
```

The build produces:
- `dist/OmicsVisStat/OmicsVisStat.exe`

## Quick Start

1. Launch the application.
2. Use the Statistics tab to import and verify data.
3. Configure groups and analysis settings.
4. Run normalization and statistical tests.
5. Use Visualization, Machine Learning, and Utility tabs as needed.
6. Export results and figures.

## Data Format

Use a wide table format with:
- One row per feature
- One column per sample
- Optional metadata columns on the left
- Numeric measurement columns for analysis

Example:

| Feature_ID | Sample_1 | Sample_2 | Sample_3 |
| --- | --- | --- | --- |
| Metab_001 | 123.45 | 145.67 | 167.89 |
| Metab_002 | 234.56 | 256.78 | 278.90 |

## Documentation

The full documentation is in the docs folder:

- [USER_GUIDE.md](docs/USER_GUIDE.md)
- [MACHINE_LEARNING_GUIDE.md](docs/MACHINE_LEARNING_GUIDE.md)
- [UTILITY_GUIDE.md](docs/UTILITY_GUIDE.md)
- [METABOLOMICS_GUIDE.md](docs/METABOLOMICS_GUIDE.md)
- [CUSTOM_OMICS_GUIDE.md](docs/CUSTOM_OMICS_GUIDE.md)
- [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)

## System Requirements

- Windows 10/11, macOS 10.14+, or Linux
- 4 GB RAM minimum
- 8 GB RAM recommended for large datasets
- 1280x720 display minimum

## Citation

If you use this tool in your research, please cite the project and release metadata in your manuscript or repository record.

## License

This project is licensed under the MIT License.

## Version History

### Version 1.0.0
- Initial public release
- Statistics workflows
- Machine learning workflows
- Visualization workflows
- Utility sub-tools

---

Always validate results with domain knowledge and your study design.
