# Machine Learning Quick Start

This is a practical first run for the current Machine Learning tab workflow.

## Goal

Run a baseline classification analysis from import to results review.

---

## 5-Minute Walkthrough

### 1. Open Machine Learning Tab

- Launch the application.
- Open the Machine Learning tab.

### 2. Import and Verify Data

- Click Import Excel.
- Select your input file.
- Click Verify Columns and confirm sample/feature mappings.

Checkpoint:
- Data status shows loaded and verified.

### 3. Configure Groups

- Click Config Groups.
- Assign each sample column to a group label.
- Use pattern-based auto-assignment if needed.
- Save and close the group dialog.

Checkpoint:
- Group assignments are complete.

### 4. Set Baseline ML Configuration

Suggested baseline:
- Model: Random Forest
- Test size: 0.2 to 0.3
- CV folds: 5
- Scaling: standard
- Keep advanced tuning disabled for first pass

Optional:
- Set a working folder for outputs
- Keep default filters unless your dataset requires stricter filtering

### 5. Run Analysis

- Click Run Analysis.
- Wait for completion.
- Review the results log and performance summary.

Checkpoint:
- Run finishes without validation errors.

---

## Interpreting First Results

Look for:
- Reasonable test performance relative to training performance
- Stable CV behavior (not extremely volatile)
- No obvious setup warnings in the log

If training is much higher than test:
- Reduce model complexity
- Tighten feature filtering
- Increase sample size if possible

---

## Common First-Run Issues

Issue: Run stays disabled
- Usually data verification or group configuration is incomplete.

Issue: Performance is weak
- Re-check group assignment quality.
- Try multiple models.
- Review class balance.

Issue: P-value filter mismatch in pairwise mode
- Verify pairwise p-value column mapping before rerun.

---

## Suggested Next Steps

1. Enable repeated CV for stability checks.
2. Compare multiple models.
3. Enable tuning/nested CV for final model selection.
4. Keep exported outputs with run settings for reproducibility.

---

Related:
- MACHINE_LEARNING_GUIDE.md
- USER_GUIDE.md
- TROUBLESHOOTING.md

Updated: April 27, 2026
