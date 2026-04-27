# Custom Omics Data Analysis Guide

## Table of Contents

1. [Overview](#overview)
2. [Supported Data Types](#supported-data-types)
3. [Data Preparation](#data-preparation)
4. [Generic Workflows](#generic-workflows)
5. [Omics-Specific Considerations](#omics-specific-considerations)
6. [Multi-Omics Integration](#multi-omics-integration)
7. [Case Studies](#case-studies)
8. [Best Practices](#best-practices)

---

## Overview

This guide demonstrates how to use OmicsVisStat for any omics data type beyond metabolomics/lipidomics, including:

- **Transcriptomics** (RNA-seq, microarray)
- **Proteomics** (mass spectrometry, antibody arrays)
- **Genomics** (SNPs, variants, CNVs)
- **Epigenomics** (DNA methylation, histone modifications)
- **Microbiomics** (16S rRNA, metagenomics)
- **Glycomics** (glycan profiles)
- **And any custom data type**

The tool provides flexible, platform-agnostic statistical analysis and visualization suitable for any high-dimensional biological dataset.

---

## Supported Data Types

### 1. Transcriptomics

**Data types:**
- RNA-seq count data (after normalization)
- Microarray intensity values
- Single-cell RNA-seq (aggregated)
- NanoString gene expression

**Typical format:**
```
Sample_ID | Group | Treatment | Gene_1 | Gene_2 | ... | Gene_N
```

**Considerations:**
- Data should be normalized (TPM, FPKM, log2(counts))
- Check for batch effects
- Consider log transformation if not already applied
- Filter low-expressed genes

### 2. Proteomics

**Data types:**
- Label-free quantification (LFQ)
- TMT/iTRAQ intensities
- SWATH/DIA quantification
- Antibody array data
- Western blot quantification

**Typical format:**
```
Sample_ID | Group | Protein_1 | Protein_2 | ... | Protein_N
```

**Considerations:**
- High proportion of missing values (impute appropriately)
- Log2 transform recommended
- Normalize by total protein or loading control
- Account for batch effects (TMT batches)

### 3. Genomics

**Data types:**
- SNP genotypes (coded as 0, 1, 2)
- Copy number variations (CNV ratios)
- Structural variants (presence/absence)
- Allele frequencies

**Typical format:**
```
Sample_ID | Group | Population | SNP_1 | SNP_2 | ... | SNP_N
```

**Considerations:**
- Hardy-Weinberg equilibrium testing
- Population stratification correction
- Linkage disequilibrium (LD) filtering
- Multiple testing burden (many variants)

### 4. Epigenomics

**Data types:**
- DNA methylation (β-values, M-values)
- Histone modifications (ChIP-seq peaks)
- ATAC-seq accessibility scores
- Bisulfite sequencing

**Typical format:**
```
Sample_ID | Group | CellType | CpG_1 | CpG_2 | ... | CpG_N
```

**Considerations:**
- Use M-values for statistics (β-values for visualization)
- Cell type heterogeneity correction
- Probe filtering (cross-reactive, SNPs)
- Regional analysis (DMR identification)

### 5. Microbiomics

**Data types:**
- 16S rRNA OTU/ASV counts
- Metagenomic relative abundances
- Fungal ITS
- Virome data

**Typical format:**
```
Sample_ID | Group | Location | OTU_1 | OTU_2 | ... | OTU_N
```

**Considerations:**
- Compositional data (relative abundances)
- Rarefaction or normalization (CSS, TSS)
- Zero-inflation (many zeros)
- Alpha/beta diversity measures
- Taxonomic level aggregation

### 6. Glycomics

**Data types:**
- HPLC glycan peaks
- Mass spectrometry glycan profiles
- Glycan array binding data

**Typical format:**
```
Sample_ID | Group | Glycan_1 | Glycan_2 | ... | Glycan_N
```

**Considerations:**
- Similar to metabolomics preprocessing
- Normalize to total glycan content
- Structural isomers may be unresolved
- Derived traits (branching, sialylation)

---

## Data Preparation

### Universal Data Format

The tool accepts any omics data in this general structure:

```
[Metadata columns] | [Data columns]
```

**Minimum requirements:**
- Sample identifiers (unique per row)
- Group/condition labels
- Numeric data values

**Example:**
```excel
Sample_ID  | Group    | Batch | Age | Sex | Feature_1 | Feature_2 | ... | Feature_N
Sample_01  | Control  | 1     | 25  | M   | 5.23      | 8.91      | ... | 3.45
Sample_02  | Control  | 1     | 30  | F   | 5.12      | 8.67      | ... | 3.23
Sample_03  | Disease  | 2     | 28  | M   | 7.89      | 10.45     | ... | 5.67
Sample_04  | Disease  | 2     | 33  | F   | 8.12      | 11.23     | ... | 6.01
```

### Column Types

**Metadata columns** (not analyzed):
- Sample_ID, Subject_ID
- Group, Condition, Treatment
- Time, Timepoint
- Batch, Run
- Age, Sex, BMI
- Any other clinical/experimental variables

**Data columns** (analyzed):
- Genes, proteins, metabolites, OTUs, etc.
- All numeric values
- Column names should be unique and descriptive

### Data Formatting Tips

1. **Headers**
   - First row contains column names
   - No special characters except underscore (_)
   - Avoid spaces (use underscores)

2. **Sample IDs**
   - Unique for each sample
   - No duplicates
   - Consistent format

3. **Group labels**
   - Consistent naming (e.g., "Control" not "Ctrl" or "control")
   - No spelling variations
   - Short and descriptive

4. **Numeric data**
   - Decimal format (use . not ,)
   - No text in data columns
   - Missing values: blank or "NA"

5. **File format**
   - Excel (.xlsx) recommended
   - CSV/TSV also supported
   - UTF-8 encoding for special characters

---

## Generic Workflows

### Workflow 1: Basic Group Comparison

**Objective**: Compare two or more groups across all features

**Steps:**

#### 1. Load Data
```
Import data using the current tab workflow (Statistics or Machine Learning tab)
- Select your data file
- Verify columns are detected correctly
```

#### 2. Assign Columns
```
Settings:
- Sample ID: First column
- Group column: Column with group labels
- Data columns: All feature columns
```

#### 3. Preprocessing
```
Missing values:
- View summary
- Choose imputation method (mean, median, KNN)
- Apply

Normalization:
- Log2 transform (if not already done)
- Z-score (for visualization)
- Or custom method
```

#### 4. Statistical Testing
```
For 2 groups:
- Statistics → Pairwise Comparison
- Select both groups
- Choose test: t-test or Mann-Whitney
- Apply FDR correction

For 3+ groups:
- Statistics → One-Way ANOVA
- Select all groups
- Choose post-hoc test: Tukey HSD
- Set significance: p < 0.05, FDR < 0.05
```

#### 5. Results Review
```
- Check number of significant features
- Review effect sizes
- Note fold changes or group differences
```

#### 6. Visualization
```
- Volcano plot (2 groups)
- PCA plot (all samples)
- Heatmap (significant features)
- Box plots (top features)
```

#### 7. Export
```
- Export statistical results (Excel)
- Export figures (PDF, PNG)
```

### Workflow 2: Multi-Factorial Design

**Objective**: Analyze effects of multiple factors (e.g., Treatment × Time)

**Data structure:**
```
Sample_ID | Factor1 | Factor2 | Features...
S001      | Ctrl    | T0      | ...
S002      | Ctrl    | T1      | ...
S003      | Drug    | T0      | ...
S004      | Drug    | T1      | ...
```

**Alternative structure:**
```
Sample_ID | Group      | Features...
S001      | Ctrl_T0    | ...
S002      | Ctrl_T1    | ...
S003      | Drug_T0    | ...
S004      | Drug_T1    | ...
```

**Steps:**

#### 1. Configure Factor Mapping
```
- Define Factor 1 (e.g., Treatment)
- Define Factor 2 (e.g., Time)
- Map groups to factor combinations
- Verify mapping is correct
```

#### 2. Run Two-Way ANOVA
```
Statistics → Two-Way ANOVA
- Select all groups
- Configure options:
  ☑ Test main effects
  ☑ Test interaction
  ☑ Post-hoc comparisons
  ☑ Effect sizes
```

#### 3. Interpret Results
```
Main effect of Factor 1:
- Overall effect of treatment

Main effect of Factor 2:
- Overall effect of time

Interaction (Factor 1 × Factor 2):
- Does treatment effect differ over time?

Post-hoc:
- Which specific combinations differ?
```

#### 4. Visualization
```
- PCA colored by Factor 1, shaped by Factor 2
- Heatmap with Factor 1 and Factor 2 annotations
- Line plots (Time vs. Feature, grouped by Treatment)
- Interaction plots
```

### Workflow 3: Time-Series Analysis

**Objective**: Identify features changing over time

**Data structure:**
```
Sample_ID | Subject_ID | Time | Group | Features...
```

**Steps:**

#### 1. Organize Data
```
- Ensure Subject_ID links repeated measures
- Time as numeric or ordered factor
- Group indicates treatment/condition
```

#### 2. Statistical Testing
```
Options:

A. Simple approach (independent timepoints):
   - One-way ANOVA with Time as factor
   - Post-hoc for specific timepoint comparisons

B. Advanced approach (repeated measures):
   - Two-way ANOVA: Group × Time
   - Account for subject correlation
   - Test for Group × Time interaction

C. Trend analysis:
   - Linear, quadratic, or cubic trends
   - Significant time trends
   - Different trends between groups
```

#### 3. Clustering
```
- Cluster features by temporal profile
- Identify co-regulated patterns:
  * Early responders
  * Late responders
  * Transient changes
  * Sustained changes
```

#### 4. Visualization
```
- Line plots: Time vs. Feature (mean ± SEM)
- Separate lines per group
- Heatmap: Features × Timepoints
- Profile plots for clusters
```

### Workflow 4: Batch Effect Correction

**Objective**: Remove technical variation while preserving biological signal

**Steps:**

#### 1. Assess Batch Effects
```
Visualization → PCA
- Color by Batch
- Check for batch separation
- If visible: proceed with correction
```

#### 2. Load Covariate Data
```
If not in main data:
- Prepare file: Sample_ID | Batch | OtherCovariates
- Load as covariate file
```

#### 3. Apply Correction
```
Preprocessing → Covariate Adjustment
- Select Batch as covariate
- Choose method:
  * Linear regression (simple)
  * ComBat (recommended for batch)
  * Residualization

- Apply correction
```

#### 4. Verify Correction
```
Visualization → PCA
- Color by Batch again
- Batch separation should be reduced
- Group separation should remain
```

#### 5. Proceed with Analysis
```
- Use corrected data for statistical testing
- Report correction method in methods
```

---

## Omics-Specific Considerations

### Transcriptomics (RNA-seq)

**Preprocessing:**
```
1. Start with normalized counts (TPM, FPKM, or DESeq2 normalized)
2. Filter low-expressed genes:
   - Remove genes with < 1 count in > 80% samples
   - Or keep genes with CPM > 1 in ≥ smallest group size
3. Log2 transform: log2(x + 1)
4. No additional normalization needed if already TMM/DESeq2 normalized
```

**Statistical considerations:**
```
- RNA-seq counts: Use external tool (DESeq2, edgeR) first, import results
- Microarray: Can analyze directly in this tool
- Batch effects common: Correct before testing
- Multiple testing: FDR essential (many genes tested)
```

**Visualization:**
```
- PCA: Check for outliers and batch effects
- Heatmap: Top variable genes (e.g., top 1000 by variance)
- Volcano plot: Log2FC vs. -log10(p-value)
- Box plots: Validation genes
```

**Interpretation:**
```
- Consider biological pathways (use pathway tools)
- Validate with qPCR
- Check literature for known markers
```

### Proteomics

**Preprocessing:**
```
1. Log2 transform (intensity values)
2. Missing value imputation:
   - If MNAR: Use minimum value imputation
   - If MAR: Use KNN (k=5)
   - Consider hybrid approach
3. Normalization:
   - Median normalization (per sample)
   - Or quantile normalization
4. Batch correction if multiple TMT batches
```

**Statistical considerations:**
```
- High missingness: May need to filter (e.g., present in 70% of samples)
- Imputation choice affects results: Test sensitivity
- Effect sizes important (p-values can be inflated with imputation)
- Consider protein families (shared peptides)
```

**Visualization:**
```
- PCA: Identify outliers
- Heatmap: Complete cases only (for cleaner visualization)
- Volcano plot: Mark proteins with high missingness
- Correlation with RNA-seq (if available)
```

### Genomics (SNPs/CNVs)

**Preprocessing:**
```
SNPs (coded 0, 1, 2 for AA, AB, BB):
1. Hardy-Weinberg equilibrium test (remove if p < 1e-6)
2. Filter by MAF (minor allele frequency > 0.05)
3. LD pruning if needed (r² < 0.5)
4. No transformation needed

CNVs (copy number ratios):
1. Log2 transform if ratios
2. Segment values if continuous
3. Normalize to diploid (2 copies)
```

**Statistical considerations:**
```
- Many variants: Very stringent FDR (q < 0.01 or lower)
- Population structure: Include PCs as covariates
- Family structure: Use appropriate test
- Gene-based or region-based testing may be more powerful
```

**Visualization:**
```
- PCA for population structure
- Manhattan plot for GWAS-like results
- Heatmap of top variants
```

### Epigenomics (DNA Methylation)

**Preprocessing:**
```
1. Convert β-values to M-values for statistics:
   M = log2(β / (1 - β))
2. Filter probes:
   - Cross-reactive probes
   - SNP-containing probes
   - Sex chromosomes (if mixed sex)
3. Normalize (if not already):
   - Illumina: Noob, SWAN, or Funnorm
   - Bisulfite-seq: normalize coverage
4. Cell type deconvolution (if heterogeneous tissue)
```

**Statistical considerations:**
```
- Use M-values for testing (more homoscedastic)
- Convert back to β-values for visualization
- Multiple testing: FDR < 0.05
- Consider regional analysis (DMR) not just single CpGs
- Cell type proportions as covariates
```

**Visualization:**
```
- Heatmap: Use β-values (0-1 scale)
- PCA: Can use either M or β-values
- Volcano plot: Use M-value differences
- Genomic tracks: Visualize DMRs in genome browser
```

### Microbiomics

**Preprocessing:**
```
1. Compositional data requires special handling:
   - Total sum scaling (TSS): Relative abundances
   - Or CSS (cumulative sum scaling)
   - CLR (centered log-ratio) transformation
   - DO NOT use log without proper transformation

2. Filter rare taxa:
   - Remove OTUs present in < 10% samples
   - Or with < 0.01% relative abundance

3. Rarefaction (optional, controversial):
   - Downsample to minimum read depth
   - Or use normalization instead
```

**Statistical considerations:**
```
- Zero-inflation: Many features with zeros
- Non-parametric tests often more appropriate
- DESeq2 can be used (treats as counts)
- Consider diversity metrics (Shannon, Simpson)
- Permutational tests (PERMANOVA) for community composition
```

**Visualization:**
```
- PCA or PCoA (principal coordinates analysis)
- Taxonomic bar charts (stacked)
- Heatmap (top abundant taxa)
- Alpha diversity box plots
- Beta diversity ordination
```

---

## Multi-Omics Integration

### Approach 1: Sequential Analysis

**Workflow:**
```
1. Analyze each omics dataset separately in the tool
2. Export significant features from each
3. Integrate results manually:
   - Venn diagrams (overlap)
   - Concordance analysis
   - Joint pathway enrichment
```

**Use case:** Simple integration, different sample sets

### Approach 2: Correlation-Based Integration

**Workflow:**
```
1. Ensure matched samples across omics
2. Identify significant features in each omics
3. Calculate cross-omics correlations:
   - Gene-protein correlations
   - Gene-metabolite correlations
   - Protein-metabolite correlations
4. Visualize correlation networks
5. Identify multi-omics hubs
```

**Use case:** Same samples, hypothesis generation

### Approach 3: Concatenated Analysis

**Workflow:**
```
1. Normalize each omics separately
2. Combine into single data matrix:
   [Samples × (Genes + Proteins + Metabolites)]
3. Analyze as single dataset:
   - PCA (shows overall structure)
   - Clustering (identifies multi-omics patterns)
4. Interpret which omics drives separation
```

**Use case:** Exploratory analysis, pattern identification

### Integration Best Practices

✅ **DO:**
- Normalize within each omics type separately
- Use same sample IDs across datasets
- Analyze each omics separately first
- Consider biological relationships (gene→protein→metabolite)
- Validate correlations
- Use pathway/network context

❌ **DON'T:**
- Directly combine raw data from different platforms
- Treat all omics equally (different noise levels)
- Over-interpret correlations without validation
- Ignore causality (correlation ≠ causation)
- Forget to account for multiple testing

---

## Case Studies

### Case Study 1: Transcriptomics - Drug Response

**Scenario:**
- Cells treated with drug vs. control
- RNA-seq gene expression
- Goal: Identify drug mechanism

**Data:**
```
Sample_ID | Treatment | Replicate | Gene_1 | Gene_2 | ... | Gene_20000
```

**Analysis:**

1. **Load and preprocess**
   - Load normalized counts (TPM)
   - Log2 transform
   - Filter low expression (< 1 TPM in > 80% samples)

2. **PCA**
   - Check for outliers
   - Verify treatment separation

3. **Differential expression**
   - Pairwise: Drug vs. Control
   - FDR < 0.05, |log2FC| > 1

4. **Visualization**
   - Volcano plot
   - Heatmap of top 100 DE genes
   - Box plots for known pathway members

5. **Export**
   - DE gene list → Pathway enrichment tool
   - Identify affected biological processes

### Case Study 2: Proteomics - Biomarker Panel

**Scenario:**
- Serum samples from disease vs. healthy
- Targeted proteomics (50 proteins)
- Goal: Build diagnostic classifier

**Data:**
```
Sample_ID | Group | Age | Sex | Protein_1 | ... | Protein_50
```

**Analysis:**

1. **Preprocessing**
   - Log2 transform
   - Impute missing (median per group)
   - Z-score normalization

2. **Covariate adjustment**
   - Adjust for Age and Sex
   - Use adjusted values

3. **Feature selection**
   - t-test for each protein
   - FDR < 0.05
   - Effect size |d| > 0.8
   - → 15 significant proteins

4. **Multivariate analysis**
   - Supervised ML tab classification with 15 proteins
   - 10-fold cross-validation
   - Accuracy: 85%, Q² = 0.68

5. **Top biomarkers**
   - VIP scores > 1.5 → 8 proteins
   - Box plots for visualization
   - ROC curves for each

6. **Export**
   - Biomarker panel for validation

### Case Study 3: Microbiome - Diet Intervention

**Scenario:**
- Gut microbiome before/after diet change
- 16S rRNA sequencing
- Goal: Identify responsive taxa

**Data:**
```
Sample_ID | Subject_ID | Timepoint | Diet | OTU_1 | ... | OTU_500
```

**Analysis:**

1. **Preprocessing**
   - Filter rare OTUs (< 0.1% abundance, < 10% prevalence)
   - TSS normalization (relative abundance)
   - CLR transformation

2. **Alpha diversity**
   - Calculate Shannon index per sample
   - Compare Before vs. After (paired t-test)

3. **Beta diversity**
   - PCA of CLR-transformed data
   - PERMANOVA (not in tool, use external)

4. **Differential abundance**
   - Paired comparison: Before vs. After
   - Wilcoxon signed-rank test
   - FDR < 0.05

5. **Visualization**
   - Taxonomic bar charts (external tool)
   - PCA colored by Time, shaped by Diet
   - Box plots of significant taxa

6. **Functional inference**
   - Map taxa to functional pathways
   - Test pathway enrichment

---

## Best Practices

### General Guidelines

✅ **DO:**
- Understand your data type and its characteristics
- Apply appropriate preprocessing for your omics
- Check assumptions of statistical tests
- Use FDR correction for multiple testing
- Report effect sizes with p-values
- Visualize data at each step
- Validate findings when possible
- Document all analysis steps

❌ **DON'T:**
- Apply methods blindly without understanding
- Skip quality control steps
- Ignore batch effects
- Use p < 0.05 as the only criterion
- Cherry-pick results
- Forget about biological context
- Over-interpret exploratory findings
- Skip normalization

### Quality Control Checklist

Before statistical testing:
```
☑ Data loaded correctly (rows = samples, columns = features)
☑ Group labels are correct and consistent
☑ No duplicate sample IDs
☑ Missing values handled appropriately
☑ Outliers investigated (biological vs. technical)
☑ Batch effects assessed and corrected if needed
☑ Data normalized/transformed as appropriate
☑ Distributions checked (histograms, Q-Q plots)
```

### Reporting Standards

In your methods, report:
```
1. Data preprocessing:
   - Software and version
   - Normalization method
   - Transformation applied
   - Missing value handling
   - Filtering criteria
   - Batch correction method

2. Statistical analysis:
   - Test used and justification
   - Multiple testing correction method
   - Significance thresholds
   - Effect size measures
   - Software/tool used

3. Visualization:
   - Plot types
   - Scaling methods
   - Color schemes
```

### Sample Size Guidelines

Minimum recommendations:
```
Exploratory/Pilot: n ≥ 5 per group
Hypothesis-testing: n ≥ 10 per group
Validation: n ≥ 30 per group
Clinical: n ≥ 50-100 per group

Consider:
- Effect size (smaller effect needs larger n)
- Variability in your system
- Multiple testing burden
- Dropout rate
- Independent validation cohort
```

---

## Advanced Topics

### Handling High Dimensionality

**When**: p (features) >> n (samples)

**Strategies:**
```
1. Feature selection:
   - Filter by variance (keep top X%)
   - Filter by significance (univariate)
   - Use regularization (LASSO, Elastic Net)

2. Dimension reduction:
   - PCA (unsupervised)
   - PLS (supervised)
   - Keep top components

3. Multiple testing:
   - Very stringent FDR
   - Consider Bonferroni for critical findings
```

### Dealing with Confounders

**Common confounders:**
- Age, sex, BMI
- Batch, run date
- Tissue heterogeneity
- Sample collection variations

**Approaches:**
```
1. Study design:
   - Match confounders across groups
   - Randomize sample order
   - Balance batches

2. Statistical adjustment:
   - Include as covariates in model
   - Stratified analysis
   - Propensity score matching

3. Computational:
   - Batch correction algorithms
   - Cell type deconvolution
   - Surrogate variable analysis
```

### Power Analysis

**Before study:**
- Determine required sample size
- Based on:
  * Expected effect size
  * Desired power (typically 0.8)
  * Significance level (0.05)
  * Variability (from pilot data)

**After study:**
- Post-hoc power for non-significant findings
- Helps interpret negative results

---

## Troubleshooting

### Common Issues

**Issue**: Groups not separating in PCA
```
Possible causes:
- Weak biological effect
- High technical noise
- Batch effects dominating
- Need more samples

Solutions:
- Use supervised method (Machine Learning tab classification)
- Increase stringency (top variable features only)
- Correct batch effects
- Filter noise
```

**Issue**: Too many/too few significant features
```
Too many (> 50%):
- Check if groups are mislabeled
- Verify preprocessing
- May indicate batch effect
- Consider stricter FDR

Too few (< 1%):
- May be true (no difference)
- Check if enough samples
- Verify data quality
- Consider less stringent threshold
```

**Issue**: Results don't make biological sense
```
Steps:
1. Verify group labels are correct
2. Check for sample mix-ups
3. Review preprocessing steps
4. Look for confounders
5. Consult domain expert
6. Validate top features independently
```

---

## Resources

### Statistical Methods
- PCA: Principal Component Analysis
- Supervised classification: Machine Learning tab models with cross-validation
- ANOVA: Analysis of Variance
- FDR: False Discovery Rate (Benjamini-Hochberg)

### External Tools for Omics
- **Transcriptomics**: DESeq2, edgeR, limma
- **Proteomics**: MaxQuant, Proteome Discoverer
- **Genomics**: PLINK, GCTA, SNPTEST
- **Epigenomics**: minfi, RnBeads, methylKit
- **Microbiomics**: QIIME2, Mothur, DADA2, phyloseq

### Further Reading
- Statistical methods for omics data
- Experimental design for high-throughput data
- Multiple testing correction methods
- Data normalization techniques

---

**Next Steps:**
- [User Guide](USER_GUIDE.md) for general features
- [Metabolomics Guide](METABOLOMICS_GUIDE.md) for metabolite-specific features
- [Troubleshooting](TROUBLESHOOTING.md) for common issues

---

**Need Help?**
- Check example datasets in `docs/examples/`
- Review error messages carefully
- Ensure data format matches requirements
- Contact support with specific questions
