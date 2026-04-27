# Utility Tab Guide

## Overview

The Utility tab contains multiple independent sub-tools for targeted analysis and figure generation.

Current sub-tabs:
- Chord Diagram
- Venn Diagram
- Pie Chart
- Heatmap
- Effect Plot
- Regression Utility
- Linear Regression Plot
- Glycan Classification

---

## Shared Best Practice

Most utility tools follow this sequence:
1. Upload data
2. Map or verify required columns
3. Configure tool-specific settings
4. Generate output and export

If generation fails, re-check:
- Input file type and sheet content
- Column mappings
- Group/pattern configuration where required
- Required dependency availability

---

## Chord Diagram

Purpose:
- Show relationships between compounds and classification categories.

Typical steps:
1. Upload Excel file
2. Map compound and classification columns
3. Optionally map log2FC and weight/value columns
4. Configure label and figure settings
5. Generate diagram

---

## Venn Diagram

Purpose:
- Compare overlaps between 2 to 6 sets.

Typical steps:
1. Select number of sets
2. For each set, define name/color and data source
3. Enter data by paste or upload
4. Configure figure size options
5. Generate diagram

---

## Pie Chart

Purpose:
- Compare class composition across configured sample groups.

Typical steps:
1. Upload data
2. Map class and sample columns
3. Configure groups and pattern-based sample assignment
4. Configure labels/font/size/save options
5. Generate charts

---

## Heatmap

Purpose:
- Generate grouped heatmaps with normalization and color options.

Typical steps:
1. Upload data
2. Map identifier and sample columns
3. Configure groups and patterns
4. Configure normalization, color scheme, and save options
5. Generate heatmap

---

## Effect Plot

Purpose:
- Plot effect sizes with confidence intervals and significance coloring.

Typical steps:
1. Upload Excel/CSV
2. Map required columns:
   - grouping variable
   - effect/estimate
   - significance column
   - CI lower/upper columns
3. Configure sorting/label/font/selection options
4. Generate and save plot

---

## Regression Utility

Purpose:
- Per-category outcome analysis with optional covariate adjustment.

Typical steps:
1. Upload data
2. Map sample ID, category columns, covariates, and metabolite columns
3. Configure method and p-adjust settings
4. Run analysis and export outputs

---

## Linear Regression Plot

Purpose:
- Build continuous regression plots from uploaded data or Regression Utility data.

Typical steps:
1. Select data source
2. Map X/Y and optional ID fields
3. Configure plot settings
4. Set output folder and file prefix
5. Generate and save plot

---

## Glycan Classification

Purpose:
- Classify glycan entries and export processed outputs.

Typical steps:
1. Upload Excel
2. Select glycan feature ID and sample columns
3. Review classification summary preview
4. Set output folder
5. Run classification and export

---

## Troubleshooting Tips

- Confirm numeric columns are truly numeric.
- Use consistent naming in sample/group labels.
- Avoid mixed-type columns in quantitative fields.
- Re-run column verification after changing input files.

For broader issues, see TROUBLESHOOTING.md.

---

Updated: April 27, 2026
