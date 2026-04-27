# Metabolomics & Lipidomics User Guide

## Table of Contents

1. [Overview](#overview)
2. [Metabolomics-Specific Features](#metabolomics-specific-features)
3. [Data Preparation](#data-preparation)
4. [Metabolite Identification](#metabolite-identification)
5. [Statistical Workflows](#statistical-workflows)
6. [Visualization for Metabolomics](#visualization-for-metabolomics)
7. [Lipidomics Analysis](#lipidomics-analysis)
8. [Pathway Integration](#pathway-integration)
9. [Common Metabolomics Scenarios](#common-metabolomics-scenarios)
10. [Best Practices](#best-practices)

---

## Overview

This guide covers specialized workflows for metabolomics and lipidomics data analysis. OmicsVisStat provides enhanced features for metabolite data, including:

- Metabolite-specific data formatting
- Lipid class annotation and analysis
- Integration-ready outputs for pathway analysis
- Specialized normalization for metabolomics
- Lipid nomenclature handling

---

## Metabolomics-Specific Features

### Enhanced Data Import

The tool recognizes metabolomics-specific columns:

**Standard metabolomics format:**
```
Sample_ID | Group | Compound_Name | KEGG_ID | HMDB_ID | Metabolite_1 | Metabolite_2 | ...
```

**Recognized identifier columns:**
- Compound names (common names, IUPAC)
- KEGG IDs (C00001, C00002, etc.)
- HMDB IDs (HMDB0000001, etc.)
- PubChem IDs
- ChEBI IDs
- InChI Keys

### Metabolite Annotation

The tool can handle various metabolite naming conventions:

- **Untargeted metabolomics**: m/z values, retention times
- **Targeted metabolomics**: Compound names
- **Semi-targeted**: Putative identifications with confidence levels

---

## Data Preparation

### Format Your Metabolomics Data

#### Untargeted Metabolomics (LC-MS/GC-MS)

**Example structure:**
```excel
Sample_ID    | Group       | Age | Sex | mz_123.456_1.23 | mz_234.567_2.34 | ...
Sample_001   | Control     | 25  | M   | 1234567.89      | 2345678.90      | ...
Sample_002   | Control     | 30  | F   | 1345678.90      | 2456789.01      | ...
Sample_003   | Disease     | 28  | M   | 1456789.01      | 2567890.12      | ...
```

**Column naming conventions:**
- **m/z_RetentionTime format**: mz_123.456_1.23
- **Feature ID**: Feature_001, Feature_002, etc.
- Include metadata columns (Age, Sex, BMI) for covariate adjustment

#### Targeted Metabolomics

**Example structure:**
```excel
Sample_ID    | Group    | Glucose | Lactate | Pyruvate | Alanine | ...
Sample_001   | Control  | 5.2     | 1.8     | 0.15     | 0.32    | ...
Sample_002   | Control  | 4.9     | 1.5     | 0.12     | 0.28    | ...
Sample_003   | Treated  | 6.8     | 2.4     | 0.22     | 0.45    | ...
```

**Requirements:**
- Metabolite names should be standardized
- Units should be consistent (μM, mM, etc.)
- Include QC samples if available

### Data Preprocessing

#### Step 1: Quality Control

**Filter features:**
```
Criteria:
- RSD in QC samples < 30%
- Detection rate > 80% in at least one group
- Blank ratio > 3 (signal in samples vs. blanks)
```

**How to apply:**
1. Load data with QC samples labeled
2. Click "QC Filtering" in preprocessing section
3. Set thresholds
4. Review filtered features

#### Step 2: Missing Value Imputation

**Recommended methods for metabolomics:**

1. **Half minimum** (for MNAR - Missing Not At Random)
   - Use when values are below detection limit
   - Replaces with half of the minimum detected value
   - Common in untargeted metabolomics

2. **KNN imputation** (for MAR - Missing At Random)
   - Use for random technical missingness
   - k = 5 or 10 neighbors
   - Better preserves data structure

3. **Group minimum**
   - Replace with minimum value within the same group
   - Useful for targeted metabolomics

**Application:**
```
1. Go to "Missing Values" section
2. View missing data summary
3. Select method based on missingness pattern:
   - > 50% missing: Consider removing feature
   - < 20% missing: KNN or group mean
   - Detection limit issues: Half minimum
4. Apply imputation
```

#### Step 3: Normalization

**Metabolomics-specific normalization:**

1. **Probabilistic Quotient Normalization (PQN)**
   - Recommended for NMR data
   - Corrects for dilution effects
   - Robust to biological variation

2. **Total Ion Current (TIC) normalization**
   - For LC-MS data
   - Corrects for instrument sensitivity
   - Applied before log transformation

3. **Internal Standard (IS) normalization**
   - Most accurate if IS available
   - Divide each metabolite by its corresponding IS
   - Corrects for extraction efficiency and instrument drift

4. **Sample weight/volume normalization**
   - Normalize to sample amount
   - Essential for absolute quantification

**Workflow:**
```
1. Select normalization method
2. For IS normalization:
   - Specify IS columns
   - Map metabolites to appropriate IS
3. Preview before/after distributions
4. Apply normalization
5. Log2 transform (recommended after normalization)
```

---

## Metabolite Identification

### Level of Identification

**MSI (Metabolomics Standards Initiative) levels:**
- **Level 1**: Confirmed with authentic standards
- **Level 2**: Putative annotation (spectral library match)
- **Level 3**: Putative characterization (class)
- **Level 4**: Unknown compounds

**How to include in analysis:**
```
1. Add "Confidence_Level" column to data
2. Filter by confidence level if needed
3. Report identification level in results
```

### Compound Name Standardization

**Recommended practices:**
- Use PubChem or HMDB names
- Be consistent (Glucose vs. D-Glucose)
- Include systematic names for lipids
- Document synonym mapping

---

## Statistical Workflows

### Workflow 1: Untargeted Metabolomics Discovery

**Objective**: Identify metabolites differing between groups

**Steps:**

1. **Load and preprocess data**
   ```
   - Load LC-MS/GC-MS data
   - Apply QC filtering (RSD < 30% in QC)
   - Impute missing values (half minimum)
   - Normalize (TIC or PQN)
   - Log2 transform
   ```

2. **Exploratory analysis**
   ```
   - Generate PCA plot
   - Check for batch effects
   - Identify outliers
   - Assess group separation
   ```

3. **Statistical testing**
   ```
   - One-way ANOVA (for 3+ groups) or t-test (2 groups)
   - Apply FDR correction (Benjamini-Hochberg)
   - Set threshold: p < 0.05, FDR < 0.05
   ```

4. **Effect size filtering**
   ```
   - Calculate fold changes
   - Filter by FC > 1.5 or 2.0
   - Calculate Cohen's d
   - Keep only large effects (d > 0.8)
   ```

5. **Visualization**
   ```
   - Volcano plot (p-value vs. fold change)
   - Heatmap of significant metabolites
   - Box plots for top metabolites
   - PCA for structure review
   - Machine Learning tab classification for supervised modeling
   ```

6. **Export results**
   ```
   - Significant metabolites table
   - Include m/z, RT, fold change, p-value, FDR
   - Ready for pathway analysis or identification
   ```

### Workflow 2: Targeted Metabolomics Quantification

**Objective**: Quantify specific metabolites across conditions

**Steps:**

1. **Data preparation**
   ```
   - Load targeted data with compound names
   - Include calibration curve data if available
   - Apply IS normalization
   ```

2. **Quality checks**
   ```
   - Check %CV for each metabolite
   - Verify linearity (if standards included)
   - Assess accuracy and precision
   ```

3. **Statistical comparison**
   ```
   - Choose appropriate test:
     * Two groups: t-test or Mann-Whitney
     * Multiple groups: ANOVA or Kruskal-Wallis
     * Time series: Repeated measures ANOVA
   - Apply multiple testing correction
   ```

4. **Post-hoc analysis**
   ```
   - Pairwise comparisons with Tukey HSD
   - Calculate effect sizes
   - Determine clinically relevant changes
   ```

5. **Visualization**
   ```
   - Box plots with significance brackets
   - Bar charts with error bars (mean ± SEM)
   - Line plots for time course
   ```

6. **Report generation**
   ```
   - Concentration tables (with units)
   - Statistical summary
   - Quality metrics
   ```

### Workflow 3: Time-Series Metabolomics

**Objective**: Analyze metabolite changes over time

**Steps:**

1. **Data structure**
   ```
   Format: Sample_ID | Subject_ID | Time | Group | Metabolites...
   Example:
   S001 | P001 | 0h  | Control | ...
   S002 | P001 | 6h  | Control | ...
   S003 | P001 | 24h | Control | ...
   ```

2. **Repeated measures analysis**
   ```
   - Two-way ANOVA with repeated measures
   - Factors: Group × Time
   - Subject as random effect
   ```

3. **Trend analysis**
   ```
   - Linear, quadratic, or cubic trends
   - Identify early, middle, or late responders
   - Calculate area under curve (AUC)
   ```

4. **Visualization**
   ```
   - Line plots: Time vs. metabolite levels
   - Separate lines per group
   - Include confidence bands (95% CI)
   ```

---

## Visualization for Metabolomics

### PCA for Metabolomics QC

**Purpose**: Check data quality and batch effects

**Best practices:**
```
1. Plot all samples including QC
2. Color by:
   - First: Sample type (QC, samples)
   - Then: Batch
   - Finally: Group

3. Expectations:
   - QC samples cluster tightly
   - No batch separation
   - Some group separation
```

**Interpretation:**
- QC spread indicates instrument variability
- Batch trends suggest batch effects
- Outliers may be failed samples

### Pathway-Level Visualization

**Approach**: Aggregate metabolites by pathway

1. **Map metabolites to pathways**
   ```
   - Use KEGG, Reactome, or custom pathway database
   - Group metabolites by biological pathway
   ```

2. **Pathway enrichment**
   ```
   - Calculate pathway scores
   - Test for enriched pathways
   - Visualize pathway activity
   ```

3. **Heatmap by pathway**
   ```
   - Group metabolites by pathway
   - Annotate pathways on heatmap
   - Show pathway boundaries
   ```

### Metabolite Correlation Networks

**Purpose**: Identify co-regulated metabolites

**Steps:**
```
1. Calculate pairwise correlations
   - Spearman correlation (robust)
   - Set threshold (|r| > 0.6)

2. Create network
   - Nodes = metabolites
   - Edges = correlations

3. Cluster analysis
   - Identify modules (communities)
   - Annotate with pathway information

4. Visualization
   - Color by metabolite class or pathway
   - Size by significance
   - Layout for clarity
```

---

## Lipidomics Analysis

### Lipid Nomenclature

**Understanding lipid names:**
```
Example: PC(16:0/18:1)
- PC: Phosphatidylcholine (lipid class)
- 16:0: Fatty acid 1 (16 carbons, 0 double bonds)
- 18:1: Fatty acid 2 (18 carbons, 1 double bond)
```

**Lipid classes:**
- **Glycerophospholipids**: PC, PE, PS, PI, PG, PA
- **Sphingolipids**: Cer, SM, HexCer, Hex2Cer
- **Glycerolipids**: TAG, DAG, MAG
- **Sterols**: Cholesterol esters, free cholesterol
- **Fatty acids**: Saturated, monounsaturated, polyunsaturated

### Lipidomics-Specific Preprocessing

**Normalization methods:**

1. **Total lipid amount**
   ```
   - Sum all lipid species
   - Normalize each to % of total
   - Useful for relative quantification
   ```

2. **Lipid class normalization**
   ```
   - Normalize within each lipid class
   - Compare PC species to each other
   - Preserves class-specific biology
   ```

3. **IS per lipid class**
   ```
   - Use class-specific internal standards
   - PC species normalized to PC-IS
   - Most accurate for lipidomics
   ```

### Lipid Class Analysis

**Steps:**

1. **Annotate lipid classes**
   ```
   - Automatically extract from lipid names
   - Group by class (PC, PE, TAG, etc.)
   ```

2. **Class-level statistics**
   ```
   - Sum lipids within each class
   - Compare class totals across groups
   - Test for class-level changes
   ```

3. **Fatty acid composition**
   ```
   - Extract chain length and unsaturation
   - Calculate saturation index
   - Analyze fatty acid distributions
   ```

4. **Visualization**
   ```
   - Stacked bar charts (class composition)
   - Saturation/unsaturation plots
   - Chain length distributions
   ```

### Lipid-Specific Plots

**Double bond distribution:**
```
1. Extract number of double bonds from each lipid
2. Create histogram by group
3. Compare saturation profiles
```

**Chain length analysis:**
```
1. Extract total carbon number
2. Plot distribution
3. Identify elongation/shortening patterns
```

**Class composition:**
```
1. Calculate % of each lipid class
2. Stacked bar chart by sample or group
3. Test for class ratio changes
```

### Example Lipidomics Workflow

**Objective**: Compare lipid profiles between disease and control

**Steps:**

1. **Load lipidomics data**
   ```
   - Data with lipid names (e.g., PC(16:0/18:1))
   - Group labels
   - QC samples included
   ```

2. **Quality control**
   ```
   - Check QC clustering
   - Filter by detection frequency
   - Remove features with high %CV in QC
   ```

3. **Normalization**
   ```
   - Apply class-specific IS normalization
   - Or use total lipid normalization
   - Log2 transform
   ```

4. **Lipid class annotation**
   ```
   - Automatic extraction from names
   - Verify annotation accuracy
   - Group by class
   ```

5. **Class-level analysis**
   ```
   - Sum lipids within classes
   - Compare class totals (t-test or ANOVA)
   - Generate class composition plots
   ```

6. **Species-level analysis**
   ```
   - Individual lipid species testing
   - FDR correction
   - Fold change filtering (FC > 1.5)
   ```

7. **Fatty acid analysis**
   ```
   - Calculate saturation indices
   - Compare chain length distributions
   - Test for compositional changes
   ```

8. **Visualization**
   ```
   - PCA for overall profiles
   - Heatmap of significant lipids
   - Volcano plot
   - Class composition stacked bars
   - Saturation/chain length plots
   ```

9. **Biological interpretation**
   ```
   - Identify enriched lipid classes
   - Relate to biological processes:
     * Membrane fluidity (saturation)
     * Energy storage (TAG levels)
     * Signaling (PA, DAG, ceramides)
   ```

---

## Pathway Integration

### Exporting for Pathway Analysis

**Prepare data for external pathway tools:**

1. **MetaboAnalyst format**
   ```
   Export:
   - Two columns: Compound name, fold change
   - Or: Compound name, t-statistic
   - Save as CSV
   ```

2. **KEGG pathway mapping**
   ```
   Export:
   - Column 1: KEGG IDs (C00001, etc.)
   - Column 2: p-values or fold changes
   - Use for KEGG Mapper
   ```

3. **HMDB integration**
   ```
   Export:
   - HMDB IDs
   - Statistical results
   - Upload to HMDB pathway tools
   ```

### Internal Pathway Scoring

**Approach**: Aggregate statistical results by pathway

1. **Map metabolites to pathways**
   ```
   - Load pathway database
   - Match compound names or IDs
   ```

2. **Calculate pathway scores**
   ```
   Methods:
   - Mean fold change of pathway members
   - GSEA (Gene Set Enrichment Analysis) for metabolites
   - Over-representation analysis (ORA)
   ```

3. **Pathway visualization**
   ```
   - Heatmap of pathway scores
   - Bar chart of top pathways
   - Network of interconnected pathways
   ```

---

## Common Metabolomics Scenarios

### Scenario 1: Biomarker Discovery

**Goal**: Find metabolites distinguishing disease from control

**Analysis plan:**
```
1. Univariate tests (ANOVA or t-test)
   - FDR < 0.05
   - Fold change > 2

2. Multivariate models (Random Forest and other ML tab models)
   - Feature selection by VIP or importance
   - Cross-validation

3. Validation
   - ROC curves
   - AUC > 0.8
   - External validation cohort

4. Panel construction
   - Combine top metabolites
   - Logistic regression
   - Optimize sensitivity/specificity
```

### Scenario 2: Drug Treatment Response

**Goal**: Identify metabolic changes due to drug

**Analysis plan:**
```
1. Paired comparison (pre vs. post treatment)
   - Paired t-test
   - Repeated measures if multiple timepoints

2. Responders vs. non-responders
   - Two-way ANOVA (Response × Time)
   - Interaction effects

3. Dose-response analysis
   - Correlation with drug concentration
   - Linear trend testing

4. Mechanism investigation
   - Pathway enrichment
   - Network analysis
   - Integration with drug targets
```

### Scenario 3: Multi-Omics Integration

**Goal**: Integrate metabolomics with transcriptomics/proteomics

**Data requirements:**
```
- Matched samples (same subjects for each omics)
- Consistent sample IDs across datasets
- Comparable statistical results
```

**Integration approaches:**
```
1. Correlation analysis
   - Metabolite-gene correlations
   - Metabolite-protein correlations

2. Pathway integration
   - Map both to same pathways
   - Concordance analysis
   - Joint pathway enrichment

3. Network integration
   - Multi-layer networks
   - Identify hub metabolites/genes
   - Causal network inference

4. Visualization
   - Combined heatmaps
   - Sankey diagrams (data flow)
   - Integrated pathway maps
```

---

## Best Practices

### Metabolomics-Specific Considerations

✅ **DO:**
- Report instrument type (LC-MS, GC-MS, NMR)
- Include QC samples (pooled samples)
- Randomize sample order
- Use appropriate internal standards
- Document sample preparation protocol
- Check for ion suppression
- Validate putative identifications
- Report metabolite identification confidence levels

❌ **DON'T:**
- Analyze without QC assessment
- Ignore batch effects
- Pool results from different platforms
- Forget to normalize
- Compare absolute intensities across batches
- Over-interpret untargeted data without validation
- Ignore isomers and isobaric compounds

### Quality Control

**Essential QC checks:**
```
1. Pooled QC samples:
   - RSD < 30% for reliable metabolites
   - Tight clustering in PCA
   - Consistent signal across run

2. Blanks:
   - Check for contamination
   - Remove background signals

3. Standard curves:
   - Verify linearity (R² > 0.99)
   - Check calibration range

4. Retention time stability:
   - Drift < 0.1 min for LC
   - Consistent across batch
```

### Sample Size Considerations

**Recommended sample sizes:**
```
Metabolomics discovery:
- Pilot: n = 10-20 per group
- Full study: n = 30-50 per group
- Validation: n = 50-100 per group

Targeted metabolomics:
- Validation: n = 50+ per group
- Clinical: n = 100+ per group

Considerations:
- Higher for small effect sizes
- Account for sample dropout (QC failure)
- Plan for independent validation cohort
```

### Data Archiving

**What to save:**
```
- Raw data files (.mzML, .raw, .d)
- Processed data (peak tables)
- Sample metadata
- Batch information
- QC summary statistics
- Processing parameters
- Statistical results
- Identification information (spectral matches)
```

**Metadata to record:**
```
- Sample collection date/time
- Storage conditions
- Extraction protocol
- Instrument parameters
- Column type (LC-MS)
- Acquisition date
- Batch assignment
- QC pass/fail
```

---

## Troubleshooting

### Common Issues

**Problem**: QC samples not clustering
```
Solution:
- Check for sample mix-up
- Verify pooling protocol
- Check instrument stability
- Consider removing problematic samples
```

**Problem**: Batch effects visible in PCA
```
Solution:
- Apply batch correction (Combat, QC-RLSC)
- Include batch as covariate in analysis
- Use batch-balanced design if possible
```

**Problem**: Too many missing values
```
Solution:
- Lower detection frequency threshold
- Check if biological (absent in some groups)
- Use appropriate imputation
- Remove unreliable features
```

**Problem**: Identification ambiguity
```
Solution:
- Use multiple databases
- Check retention time against standards
- Consider MS/MS fragmentation
- Report as putative if uncertain
- Use confidence levels (MSI standards)
```

---

## Resources

### Databases
- **HMDB**: Human Metabolome Database
- **KEGG**: Kyoto Encyclopedia of Genes and Genomes
- **MetaboAnalyst**: Online metabolomics analysis
- **LipidMaps**: Lipid structure database
- **PubChem**: Chemical compound database

### Tools
- **MS-DIAL**: Untargeted metabolomics processing
- **XCMS**: LC-MS data processing
- **MZmine**: Mass spectrometry data processing
- **Compound Discoverer**: Thermo Fisher workflow
- **Progenesis QI**: Waters workflow

### Standards
- **MSI**: Metabolomics Standards Initiative
- **MIAME**: Minimum Information About a Metabolomics Experiment
- **mzML**: Standard data format

---

**Next Steps:**
- [Custom Omics Guide](CUSTOM_OMICS_GUIDE.md) for non-metabolomics data
- [User Guide](USER_GUIDE.md) for general features
- [Troubleshooting](TROUBLESHOOTING.md) for common issues
