# Troubleshooting Guide

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Data Loading Problems](#data-loading-problems)
3. [Statistical Analysis Errors](#statistical-analysis-errors)
4. [Visualization Issues](#visualization-issues)
5. [Performance Problems](#performance-problems)
6. [Common Error Messages](#common-error-messages)
7. [FAQ](#faq)

---

## Installation Issues

### Problem: Application won't start (Windows .exe)

**Symptoms:**
- Double-clicking .exe does nothing
- Application starts then immediately closes
- "Windows protected your PC" message

**Solutions:**

1. **Windows SmartScreen warning:**
   ```
   - Click "More info"
   - Click "Run anyway"
   - This is normal for unsigned executables
   ```

2. **Missing dependencies:**
   ```
   - Download and install Visual C++ Redistributable
   - Available from Microsoft website
   - Choose appropriate version (x64 for 64-bit Windows)
   ```

3. **Antivirus blocking:**
   ```
   - Check antivirus quarantine
   - Add application to exceptions
   - Temporarily disable to test
   ```

4. **Permission issues:**
   ```
   - Right-click .exe → Properties
   - Check "Unblock" if present
   - Run as Administrator
   ```

### Problem: Python script won't run

**Symptoms:**
- `ModuleNotFoundError`
- `ImportError`
- Script crashes on launch

**Solutions:**

1. **Verify Python version:**
   ```bash
   python --version
   # Should be 3.8 or higher
   ```

2. **Install all dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Check specific missing modules:**
   ```bash
   pip install pandas numpy scipy matplotlib seaborn scikit-learn statsmodels openpyxl
   ```

4. **Use virtual environment:**
   ```bash
   python -m venv venv
   venv\Scripts\activate  # Windows
   source venv/bin/activate  # Linux/Mac
   pip install -r requirements.txt
   ```

5. **Update pip:**
   ```bash
   python -m pip install --upgrade pip
   ```

---

## Data Loading Problems

### Problem: File won't load

**Error:** "Failed to load data file"

**Solutions:**

1. **Check file format:**
   ```
   Supported: .xlsx, .xls, .csv, .tsv
   Not supported: .xlsm (macro-enabled), .txt (unless tab-delimited)
   ```

2. **Verify file is not corrupted:**
   ```
   - Open file in Excel or text editor
   - Check if data displays correctly
   - Save as new file
   ```

3. **Check file path:**
   ```
   - Avoid special characters in filename
   - Avoid extremely long paths
   - No network drives (copy locally first)
   ```

4. **File permissions:**
   ```
   - Ensure file is not open in another program
   - Check read permissions
   - Close Excel if file is open there
   ```

### Problem: Data loaded but columns not recognized

**Symptoms:**
- All columns treated as data
- Group column not detected
- Sample IDs missing

**Solutions:**

1. **Verify data structure:**
   ```
   Required format:
   - First row: Column headers
   - First column: Sample IDs
   - Include a group/condition column
   - Remaining columns: numeric data
   ```

2. **Check for hidden characters:**
   ```
   - Extra spaces in headers
   - Non-breaking spaces
   - Invisible characters
   - Clean headers in Excel: TRIM() function
   ```

3. **Manually assign columns:**
   ```
   - Use column assignment interface
   - Specify which column is Sample ID
   - Specify which column is Group
   - Select data columns
   ```

### Problem: Missing values not handled correctly

**Symptoms:**
- "NaN" appearing in results
- Analysis fails with missing data error
- Unexpected number of missing values

**Solutions:**

1. **Standardize missing value encoding:**
   ```
   Recognized formats:
   - Empty cells (blank)
   - "NA"
   - "N/A"
   - "nan"
   - "NaN"
   
   Not recognized:
   - "." (period alone)
   - "0" (if zero has biological meaning)
   - "-" (hyphen)
   - Text like "below detection"
   ```

2. **Pre-process in Excel:**
   ```
   - Find and replace all missing indicators with blank
   - Or use "NA"
   - Ensure consistency
   ```

3. **Check data types:**
   ```
   - Numeric columns should only contain numbers
   - Remove any text in data columns
   - Check for mixed content (e.g., "< 0.001")
   ```

### Problem: Groups not detected correctly

**Symptoms:**
- Wrong number of groups
- Groups have strange names
- Some samples not assigned to groups

**Solutions:**

1. **Standardize group names:**
   ```
   Good examples:
   - Control, Treatment
   - Healthy, Disease
   - Time_0, Time_6, Time_24
   
   Bad examples:
   - control, Control, CONTROL (inconsistent case)
   - Group 1, Group_1, Group1 (mixed formats)
   - Extra spaces: "Control " vs "Control"
   ```

2. **Check for hidden characters:**
   ```
   In Excel:
   - Select group column
   - Data → Remove Duplicates (to see unique values)
   - Look for subtle differences
   ```

3. **Verify all samples have group labels:**
   ```
   - No blank cells in group column
   - Each sample assigned to exactly one group
   ```

---

## Statistical Analysis Errors

### Problem: "Insufficient samples" error

**Error:** "Not enough samples for statistical analysis"

**Minimum requirements:**
- One-way ANOVA: n ≥ 3 per group, 3+ groups
- Two-way ANOVA: n ≥ 2 per combination
- t-test: n ≥ 3 per group (recommended)
- Pairwise comparisons: n ≥ 3 per group

**Solutions:**

1. **Verify sample counts:**
   ```
   - Check how many samples per group
   - Look for samples with missing group labels
   - Confirm no samples were accidentally filtered
   ```

2. **Adjust analysis:**
   ```
   - If only 2 samples per group: Cannot perform test
   - Consider non-parametric alternatives
   - Combine timepoints or conditions if appropriate
   ```

### Problem: "All features fail statistical test"

**Symptoms:**
- No significant results
- All p-values = 1.0 or NaN
- Analysis completes but no significant features

**Solutions:**

1. **Check data quality:**
   ```
   - View data in table: Are values varying?
   - Check if normalization removed all variance
   - Verify data was loaded correctly
   ```

2. **Verify groups are different:**
   ```
   - Generate PCA plot
   - If no separation, groups may be truly similar
   - Check if group labels are correct
   ```

3. **Review thresholds:**
   ```
   - Try less stringent p-value (p < 0.1)
   - Remove fold change filter temporarily
   - Check if raw p-values are significant (before FDR)
   ```

4. **Check for zero variance features:**
   ```
   - Filter features with no variation
   - Remove constant features
   - Check if log-transform created infinities
   ```

### Problem: "Singular matrix" or "Linear dependency" error

**Error:** Common in two-way ANOVA or regression

**Causes:**
- Unbalanced design (some factor combinations missing)
- Collinear covariates
- Too few samples for number of factors

**Solutions:**

1. **Check design balance:**
   ```
   For two-way ANOVA:
   - Ensure all factor combinations have samples
   - Example: If testing Treatment × Time
     Must have: Control_T0, Control_T1, Drug_T0, Drug_T1
   ```

2. **Remove redundant covariates:**
   ```
   - Don't include perfectly correlated variables
   - E.g., Age and AgeGroup derived from Age
   ```

3. **Simplify model:**
   ```
   - Use one-way ANOVA instead
   - Analyze factors separately
   - Reduce number of covariates
   ```

### Problem: All features significant (suspiciously many)

**Symptoms:**
- > 50% features significant
- All p-values near zero
- Results seem too good to be true

**Likely causes:**

1. **Batch effect:**
   ```
   Check:
   - Are groups confounded with batch?
   - Generate PCA colored by batch
   - If batch = group, results are invalid
   
   Solution:
   - Reorganize data to break confounding
   - Apply batch correction
   - Re-collect data with balanced design
   ```

2. **Group labels swapped or incorrect:**
   ```
   - Verify each sample's group assignment
   - Check original study records
   - Look for systematic errors
   ```

3. **Data leakage:**
   ```
   - Normalization was done including group information
   - Should normalize within samples, not across groups
   ```

---

## Visualization Issues

### Problem: Plot not displaying

**Symptoms:**
- Blank plot window
- "Figure not generated" error
- Window opens but is empty

**Solutions:**

1. **Check data selection:**
   ```
   - Ensure enough data points
   - Verify groups are selected
   - Check if features are filtered out
   ```

2. **Restart application:**
   ```
   - Close and reopen
   - Reload data
   - Try plot again
   ```

3. **Update graphics backend:**
   ```python
   # If running from Python
   import matplotlib
   matplotlib.use('TkAgg')  # Or 'Qt5Agg'
   ```

### Problem: Heatmap is blank or all one color

**Causes:**
- No variance in data
- All values identical
- Incorrect scaling

**Solutions:**

1. **Check data variance:**
   ```
   - Verify features have different values
   - Look at data preview
   ```

2. **Adjust scaling:**
   ```
   - Try different normalization (Z-score, min-max)
   - Check color scale limits
   - Use row normalization (per feature)
   ```

3. **Filter low-variance features:**
   ```
   - Remove constant features
   - Keep only top variable features
   ```

### Problem: PCA plot shows no separation

**Causes:**
- Groups are truly similar
- Need more components
- Batch effects dominating

**Solutions:**

1. **Try different components:**
   ```
   - Plot PC2 vs PC3 instead of PC1 vs PC2
   - Check variance explained by each PC
   - May need higher PCs for separation
   ```

2. **Supervised method:**
   ```
   - Use Machine Learning tab classification (supervised)
   - Compare multiple models and CV settings
   ```

3. **Feature selection:**
   ```
   - Use only significant features
   - Or top 500-1000 most variable
   - Reduces noise
   ```

4. **Check for batch effects:**
   ```
   - Color by batch instead of group
   - If batch dominates, apply correction
   ```

### Problem: Labels overlap or are unreadable

**Solutions:**

1. **Adjust plot size:**
   ```
   - Increase figure dimensions
   - Use larger plot window
   ```

2. **Reduce number of labels:**
   ```
   - Label only top N features
   - Remove labels (show in legend instead)
   ```

3. **Adjust font size:**
   ```
   - Decrease font size in settings
   - Rotate labels if needed
   ```

4. **Export and edit:**
   ```
   - Export as PDF or SVG
   - Edit in vector graphics editor (Inkscape, Illustrator)
   ```

---

## Performance Problems

### Problem: Application is very slow

**Symptoms:**
- Analysis takes very long
- Application freezes
- High CPU or memory usage

**Solutions:**

1. **Reduce data size:**
   ```
   - Filter features first (e.g., remove low variance)
   - Subsample features for exploration
   - Remove QC samples after QC checks
   ```

2. **Close other programs:**
   ```
   - Free up RAM
   - Close unnecessary applications
   ```

3. **Optimize settings:**
   ```
   For PCA:
   - Use subset of samples for preview
   - Reduce number of components
   
   For clustering:
   - Use fewer features
   - Simplify distance metric
   ```

4. **Check data type:**
   ```
   - Ensure data columns are numeric
   - Remove text columns from data section
   ```

### Problem: Application crashes or freezes

**When:**
- During statistical analysis
- Generating complex plots
- With large datasets

**Solutions:**

1. **Increase memory (if Python):**
   ```bash
   # Run with more memory
   python -Xmx4g run_gui.py
   ```

2. **Break up analysis:**
   ```
   - Analyze subsets of features
   - Process in batches
   - Combine results afterward
   ```

3. **Simplify visualization:**
   ```
   - Reduce number of features in heatmap
   - Use 2D instead of 3D plots
   - Lower resolution for export
   ```

4. **Check for infinite loops:**
   ```
   - Look for warning messages
   - Check log files
   - Report to support
   ```

---

## Common Error Messages

### "ValueError: could not convert string to float"

**Cause:** Non-numeric data in data columns

**Solution:**
```
1. Check data columns for text
2. Look for:
   - Currency symbols ($)
   - Percentages with %
   - Commas in numbers (1,000)
   - Text annotations
3. Clean data in Excel before loading
```

### "KeyError: 'Group'" or "Column not found"

**Cause:** Expected column missing

**Solution:**
```
1. Verify column name exactly matches
2. Check for:
   - Extra spaces in column name
   - Different case (Group vs group)
   - Renamed column
3. Manually specify column in settings
```

### "MemoryError" or "Cannot allocate memory"

**Cause:** Dataset too large for available RAM

**Solutions:**
```
1. Close other applications
2. Reduce dataset size:
   - Fewer samples
   - Fewer features
   - Filter before analysis
3. Use computer with more RAM
4. Process in chunks
```

### "LinAlgError: Singular matrix"

**Cause:** Mathematical impossibility (no unique solution)

**Common reasons:**
```
1. Perfectly correlated variables
2. More variables than samples (p >> n)
3. Missing factor combinations in design
```

**Solutions:**
```
1. Remove redundant features
2. Use regularization
3. Simplify statistical model
4. Add more samples
```

### "Warning: p-value adjusted to 1.0 for all tests"

**Cause:** Very stringent multiple testing correction

**Interpretation:**
```
- No features pass FDR threshold
- May indicate:
  * True null (no differences)
  * Underpowered study
  * High noise
```

**Options:**
```
1. Report raw p-values with caution
2. Use less stringent FDR threshold
3. Focus on effect sizes
4. Increase sample size for future studies
```

---

## FAQ

### Q: How do I know if my p-values are reliable?

**A:** Check these:
```
✓ Assumptions met (normality for parametric tests)
✓ Sufficient sample size (n ≥ 3 per group minimum, n ≥ 10 recommended)
✓ No outliers driving significance
✓ Multiple testing corrected (FDR)
✓ Effect sizes are meaningful (not just statistically significant)
✓ Results make biological sense
```

### Q: Should I remove outliers?

**A:** Depends on cause:
```
Technical outliers (instrument error):
- Yes, remove and document

Biological outliers (real biological variation):
- Generally keep
- Report sensitivity analysis (with and without)
- Use robust methods (non-parametric tests)

Unknown cause:
- Investigate first
- Check sample quality metrics
- Consider technical replicates
- Document decision either way
```

### Q: What p-value threshold should I use?

**A:** Depends on context:
```
Exploratory (hypothesis generation):
- p < 0.05 (uncorrected) or FDR < 0.1

Standard hypothesis testing:
- FDR < 0.05 (recommended)

Validation/Confirmation:
- p < 0.05 (Bonferroni corrected)

Clinical applications:
- Very stringent (p < 0.001 or lower)

Note: Always report:
- Exact p-values (not just < 0.05)
- Correction method used
- Effect sizes
```

### Q: How many samples do I need?

**A:** Rule of thumb:
```
Minimum viable:
- n = 3-5 per group (pilot only)

Standard study:
- n = 10-20 per group

Well-powered study:
- n = 30-50 per group

Clinical/Validation:
- n = 100+ per group

Calculate properly:
- Use power analysis
- Consider effect size
- Account for dropouts
- Plan for validation cohort
```

### Q: Should I use parametric or non-parametric tests?

**A:** Decision tree:
```
1. Check sample size:
   - Small (n < 30): Consider non-parametric
   - Large (n ≥ 30): Parametric often robust

2. Check distribution:
   - Normal: Parametric
   - Skewed: Non-parametric or transform
   - Heavy tails: Non-parametric

3. Check variance:
   - Equal variance: Parametric
   - Unequal variance: Welch's t-test or non-parametric

4. When in doubt:
   - Run both
   - If similar conclusions: Report parametric
   - If different: Report both, explain discrepancy
```

### Q: What's the difference between p-value and FDR?

**A:**
```
P-value:
- Probability of seeing this result if null hypothesis true
- For single test
- Type I error rate per test

FDR (False Discovery Rate):
- Proportion of false positives among significant results
- For multiple tests
- Controls overall false discovery rate

Example:
- Test 10,000 features
- 100 significant at p < 0.05
- Without correction: Expect ~500 false positives
- With FDR < 0.05: Expect ~5 false positives among 100

Always use FDR when testing multiple features!
```

### Q: How do I choose the right normalization method?

**A:** By data type:
```
RNA-seq:
- DESeq2, TMM, or edgeR (external)
- If already normalized: log2 transform

Microarray:
- Quantile normalization
- RMA (for Affymetrix)

Proteomics:
- Median normalization
- Log2 transform
- Quantile for batch correction

Metabolomics:
- TIC or PQN
- Internal standard if available
- Log2 transform

General:
- Z-score for visualization
- Median centering for batch correction
- Check distributions before/after
```

### Q: What if my experiment has no significant results?

**A:** Consider:
```
1. Is it a true null result?
   - Groups may actually be similar
   - Effect size too small to detect

2. Statistical power:
   - Perform post-hoc power analysis
   - May need more samples

3. Data quality:
   - High technical noise
   - Check QC metrics
   - Consider technical replicates

4. Analysis approach:
   - Try different normalization
   - Check for batch effects
   - Use multivariate methods

5. Biological interpretation:
   - No change can be important finding
   - Report negative results (avoid publication bias)
```

---

## Getting Additional Help

### Before contacting support:

1. **Check documentation:**
   - User Guide
   - This troubleshooting guide
   - Omics-specific guides

2. **Review error message:**
   - Take screenshot
   - Copy exact error text
   - Note when it occurs

3. **Test with example data:**
   - Does error occur with examples?
   - Helps isolate data vs. software issue

### When reporting an issue:

Include:
```
1. Operating system and version
2. Application version
3. Exact steps to reproduce
4. Error message (screenshot)
5. Data format description (not actual data)
6. What you expected to happen
7. What actually happened
```

### Example data request:

If asked to share data:
```
- Anonymize sample IDs
- Remove sensitive information
- Include only first few samples and features
- Keep same format and structure
- Preserve the error-causing issue
```

---

## Reporting Bugs

### Bug report template:

```
**Description:**
Clear, concise description of the bug

**Steps to reproduce:**
1. Launch application
2. Load data file
3. Click on Statistics tab
4. ...

**Expected behavior:**
What you expected to happen

**Actual behavior:**
What actually happened

**Screenshots:**
If applicable

**Environment:**
- OS: Windows 10
- Version: 1.0.0
- Python version (if applicable): 3.9

**Additional context:**
Any other relevant information
```

---

**Still need help?**
- Check the [User Guide](USER_GUIDE.md)
- Review omics-specific guides
- Contact support with detailed information
- Join user community (if available)

Remember: Most issues are due to data formatting. Double-check your input data structure before reporting bugs!
