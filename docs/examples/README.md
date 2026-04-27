# Example Datasets

This folder is reserved for optional example datasets.

## Current Status

At the moment, no example dataset files are included in this folder.

If you add examples later, place them here and keep this README in sync.

## Data Format Template

All examples follow this general structure:

```excel
Feature_ID   | [Feature_Metadata...] | Sample_001 | Sample_002 | Sample_003 | ... | Sample_N
-------------|------------------------|------------|------------|------------|-----|---------
Metab_001    | ...                    | 123.45     | 145.67     | 167.89     | ... | 189.01
Metab_002    | ...                    | 234.56     | 256.78     | 278.90     | ... | 290.12
...          | ...                    | ...        | ...        | ...        | ... | ...
```

### Column Types

**Required:**
- A feature identifier column (e.g., `Feature_ID`, `metabolite`, `name`)
- Sample columns: Numeric data values (one column per sample)

**Optional:**
- Feature metadata columns (m/z, RT, adduct, lipid class, etc.)

### Creating Your Own Data

To format your data for the tool:

1. **Organize in Excel/CSV**:
   - Features in rows
   - Samples in columns
   - First row = column headers

2. **Required columns**:
   - Sample IDs (unique)
   - Group labels (consistent)
   - Numeric data

3. **Missing values**:
   - Leave blank or use "NA"
   - Be consistent

4. **Save format**:
   - Excel: `.xlsx` (recommended)
   - CSV: `.csv` with UTF-8 encoding

## Using Your Own Test Data

Until example files are added, use your own test data and compare format to the template below.

## Example Workflows

### Workflow 1: Metabolomics Discovery
```
1. Load: metabolomics_example.xlsx
2. Filter: Remove features with >50% missing
3. Impute: Half-minimum method
4. Normalize: Log2 transform
5. Test: One-way ANOVA (Disease vs Control)
6. Visualize: Volcano plot, PCA, heatmap
```

### Workflow 2: Time-Series Analysis
```
1. Load: transcriptomics_example.xlsx
2. Configure: Two-way ANOVA (Treatment × Time)
3. Test: Main effects and interaction
4. Visualize: Line plots over time, heatmap by cluster
```

### Workflow 3: Biomarker Selection
```
1. Load: proteomics_example.xlsx
2. Impute: KNN imputation
3. Adjust: Control for age and sex
4. Test: Pairwise comparison (Case vs Control)
5. Select: Features with FDR < 0.05 and |FC| > 1.5
6. Visualize: ROC curves, box plots
```

## Creating Custom Examples

Want to contribute an example dataset?

**Requirements:**
- Anonymized or simulated data
- Clear biological context
- Proper format (matches template)
- Documented expected results
- Small size (< 5 MB)

**Include:**
- Brief description
- Study design
- Sample size and groups
- Features measured
- Suggested analyses

## Notes

If you contribute example data:
- Use simulated or properly anonymized data only.
- Keep files small and easy to open.
- Add a short description of intended workflow coverage.

📝 **Best Practices**:
- Test with examples before using your data
- Compare your data format to examples
- Use examples to learn analysis workflows
- Adapt workflows to your specific study

## Troubleshooting Without Example Files

If your data will not load:
1. Compare your file to the template structure in this README.
2. Check USER_GUIDE.md for mapping expectations.
3. Check TROUBLESHOOTING.md for specific error cases.

Common issues:
- Extra header rows
- Inconsistent group names
- Text in numeric columns
- Special characters in headers

---

## Additional Resources

- [User Guide](../USER_GUIDE.md): Detailed instructions
- [Metabolomics Guide](../METABOLOMICS_GUIDE.md): Metabolite-specific workflows
- [Custom Omics Guide](../CUSTOM_OMICS_GUIDE.md): Other omics types
- [Troubleshooting](../TROUBLESHOOTING.md): Common problems

---

Note: This folder currently contains documentation only.
