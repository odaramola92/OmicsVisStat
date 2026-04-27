# Machine Learning Guide

## Overview

The Machine Learning tab supports classification-focused workflows with configurable preprocessing, validation, filtering, and export.

Core UI structure:
- Steps 1-2 panel: data import/verification and groups
- Step 3 panel: ML configuration
- Step 4 panel: run actions and results log

---

## Workflow

### Step 1: Import and Verify Data

Actions:
- Import Excel input data
- Set working folder
- Verify columns

Expected outcome:
- Feature and sample columns are recognized
- Data status indicates readiness

### Step 2: Configure Groups

Actions:
- Configure group IDs and labels
- Add/remove groups
- Assign sample columns to groups
- Optional pattern-based auto-assignment

Expected outcome:
- All required sample columns are assigned to intended groups

### Step 3: Configure ML Settings

Analysis mode
- Classification workflow (with optional pairwise add-on in multiclass contexts)

Model selection
- Select one or multiple models for evaluation

Core parameters
- Test size
- CV folds
- Scaling method
- Class weight

Additional controls
- Linear-model regularization options
- Robustness testing (repeated runs, seed, feature stability)
- Optional auto-generated figure settings

Advanced validation and tuning options
- Hyperparameter tuning
- Tuning strategy (grid/random)
- Repeated stratified CV
- Nested CV
- SVM calibration
- Permutation runs
- Imputation method
- Feature selection options

Feature filters (combinable)
- Replicate filtering
- Endogenous-only filter
- HMDB-present filter
- P-value filter

Pairwise p-value mapping
- Verify pairwise p-value columns when p-value filter is enabled for pairwise ML workflows

### Step 4: Run and Review

Actions:
- Run analysis
- Test models
- Clear log

Outputs typically include:
- Performance metrics
- Model comparison summaries
- Feature ranking information
- Optional generated figures
- Exported result files in working/output folders

---

## Recommended Run Order

1. Verify columns first
2. Configure groups completely
3. Start with simple settings (single model, moderate CV)
4. Enable advanced tuning/validation after baseline run
5. Export and archive results with settings

---

## Practical Defaults

Suggested starting configuration:
- Test size around 0.2 to 0.3
- CV folds around 5
- Standard scaling unless outlier-heavy data suggests robust scaling
- Fixed seed for reproducibility

Then iterate by:
- Trying additional models
- Enabling repeated CV
- Enabling tuning/nested CV for final model selection

---

## Common Issues

Problem: Run button stays disabled
- Usually means required setup is incomplete (data verify/groups)

Problem: Unexpectedly weak performance
- Re-check group assignments
- Re-check data filtering and class balance
- Compare multiple models before concluding

Problem: P-value filter not behaving as expected
- Verify pairwise p-value column mapping
- Ensure required p-value columns were correctly assigned during verification

---

## Related Docs

- USER_GUIDE.md
- ML_QUICK_START.md
- TROUBLESHOOTING.md
- ML_IMPLEMENTATION_SUMMARY.md

---

Updated: April 27, 2026
