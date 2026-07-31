"""
Machine Learning Models for Metabolomics Data Analysis

This module provides machine learning capabilities for metabolomics data including:
- Classification models (Random Forest, SVM, Gradient Boosting, etc.)
- Dimensionality reduction (PCA, LDA)
- Feature importance and selection
- Model validation and evaluation

Author: Metabolomics Statistics Tool
Date: 2025
"""

import pandas as pd
import numpy as np
import logging
import os
import time
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Machine learning imports
from sklearn.model_selection import (
    train_test_split,
    cross_val_score,
    cross_validate,
    StratifiedKFold,
    RepeatedStratifiedKFold,
    GridSearchCV,
    RandomizedSearchCV,
    permutation_test_score,
    cross_val_predict,
)
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler, RobustScaler, LabelEncoder, label_binarize
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.base import clone
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif, SelectFromModel, RFE
from sklearn.metrics import (
    classification_report, confusion_matrix, 
    roc_auc_score, roc_curve, accuracy_score,
    precision_score, recall_score, f1_score,
    balanced_accuracy_score,
)

from main_script.metabolite_statistics_analysis import apply_imputation

logger = logging.getLogger(__name__)


class MetabolomicsMLAnalysis:
    """
    Machine learning analysis for metabolomics data.
    
    This class handles the complete ML pipeline including:
    - Data preparation and preprocessing
    - Model training and validation
    - Feature importance extraction
    - Results reporting
    """
    
    # Available classification models
    CLASSIFICATION_MODELS = {
        'Random Forest': RandomForestClassifier,
        'Gradient Boosting': GradientBoostingClassifier,
        'SVM (RBF)': SVC,
        'Logistic Regression': LogisticRegression,
        'Linear Discriminant Analysis': LinearDiscriminantAnalysis
    }
    
    # Default hyperparameters for each model
    DEFAULT_PARAMS = {
        'Random Forest': {
            'n_estimators': 100,
            'max_depth': None,
            'min_samples_split': 2,
            'random_state': 42,
            'n_jobs': -1
        },
        'Gradient Boosting': {
            'n_estimators': 100,
            'learning_rate': 0.1,
            'max_depth': 3,
            'random_state': 42
        },
        'SVM (RBF)': {
            'kernel': 'rbf',
            'C': 1.0,
            'gamma': 'scale',
            'probability': True,
            'random_state': 42
        },
        'Logistic Regression': {
            'max_iter': 1000,
            'random_state': 42
        },
        'Linear Discriminant Analysis': {
            'solver': 'svd'
        }
    }
    
    def __init__(self, data_df: pd.DataFrame, group_assignments: Dict[str, str], 
                 feature_columns: List[str], metadata_columns: Optional[List[str]] = None,
                 feature_id_col: Optional[str] = None,
                 covariate_data: Optional[pd.DataFrame] = None,
                 covariate_cols: Optional[List[str]] = None):
        """
        Initialize ML analysis.
        
        Args:
            data_df: DataFrame with metabolite data (features as rows, samples as columns)
            group_assignments: Dict mapping sample names to group labels
            feature_columns: List of column names that contain sample data
            metadata_columns: Optional list of metadata column names to exclude
            feature_id_col: Optional column name containing feature IDs/names (from verify columns)
        """
        self.data_df = data_df
        self.group_assignments = group_assignments
        self.feature_columns = feature_columns
        self.metadata_columns = metadata_columns or []
        self.covariate_data = covariate_data.copy() if isinstance(covariate_data, pd.DataFrame) else None
        self.covariate_cols = list(covariate_cols or [])
        
        # Extract feature names (metabolite names) for feature importance display
        # Priority: 1) User-specified feature_id_col, 2) 'Name', 3) 'LipidIon', 4) generic Feature_N
        if feature_id_col and feature_id_col in data_df.columns:
            self.feature_names = data_df[feature_id_col].astype(str).tolist()
        elif 'Name' in data_df.columns:
            self.feature_names = data_df['Name'].tolist()
        elif 'LipidIon' in data_df.columns:
            self.feature_names = data_df['LipidIon'].tolist()
        else:
            self.feature_names = [f"Feature_{i}" for i in range(len(data_df))]
        
        self.scaler = None
        self.model = None
        self.label_encoder = LabelEncoder()
        self.results = {}
        
        logger.info(f"Initialized ML analysis with {len(feature_columns)} samples")

    def _build_covariate_design_matrix(
        self,
        sample_names: List[str],
        fit_state: Optional[Dict[str, Any]] = None
    ) -> Tuple[Optional[pd.DataFrame], Optional[Dict[str, Any]]]:
        """Build a numeric covariate matrix aligned to the requested samples."""
        if self.covariate_data is None or not self.covariate_cols:
            return None, fit_state

        cov_df = self.covariate_data.reindex(sample_names)
        available_cols = [col for col in self.covariate_cols if col in cov_df.columns]
        if not available_cols:
            return None, fit_state

        missing_token = '__missing__'

        if fit_state is None:
            numeric_cols: List[str] = []
            categorical_cols: List[str] = []
            numeric_fill_values: Dict[str, float] = {}
            prepared = pd.DataFrame(index=cov_df.index)

            for col in available_cols:
                series = cov_df[col]
                numeric_series = pd.to_numeric(series, errors='coerce')
                numeric_ratio = float(numeric_series.notna().mean()) if len(numeric_series) > 0 else 0.0

                if numeric_ratio > 0.9:
                    numeric_cols.append(col)
                    fill_value = float(numeric_series.mean()) if numeric_series.notna().any() else 0.0
                    if not np.isfinite(fill_value):
                        fill_value = 0.0
                    numeric_fill_values[col] = fill_value
                    prepared[col] = numeric_series.fillna(fill_value).astype(float)
                else:
                    categorical_cols.append(col)
                    prepared[col] = series.astype(str).replace({'nan': missing_token, 'None': missing_token, '<NA>': missing_token}).fillna(missing_token)

            encoded = pd.get_dummies(
                prepared[numeric_cols + categorical_cols],
                columns=categorical_cols,
                drop_first=True,
                dtype=float,
            )

            fit_state = {
                'numeric_cols': numeric_cols,
                'categorical_cols': categorical_cols,
                'numeric_fill_values': numeric_fill_values,
                'columns': encoded.columns.tolist(),
            }
            return encoded, fit_state

        numeric_cols = list(fit_state.get('numeric_cols', []))
        categorical_cols = list(fit_state.get('categorical_cols', []))
        numeric_fill_values = dict(fit_state.get('numeric_fill_values', {}))
        ordered_cols = list(fit_state.get('columns', []))

        prepared = pd.DataFrame(index=cov_df.index)
        for col in numeric_cols:
            if col in cov_df.columns:
                series = pd.to_numeric(cov_df[col], errors='coerce')
            else:
                series = pd.Series(np.nan, index=cov_df.index)
            prepared[col] = series.fillna(float(numeric_fill_values.get(col, 0.0))).astype(float)

        for col in categorical_cols:
            if col in cov_df.columns:
                series = cov_df[col]
            else:
                series = pd.Series(missing_token, index=cov_df.index)
            prepared[col] = series.astype(str).replace({'nan': missing_token, 'None': missing_token, '<NA>': missing_token}).fillna(missing_token)

        encoded = pd.get_dummies(
            prepared[numeric_cols + categorical_cols],
            columns=categorical_cols,
            drop_first=True,
            dtype=float,
        )
        if ordered_cols:
            encoded = encoded.reindex(columns=ordered_cols, fill_value=0.0)
        return encoded, fit_state

    @staticmethod
    def _residualize_feature_matrix(X_matrix: np.ndarray, design_matrix: Optional[pd.DataFrame]) -> np.ndarray:
        """Remove covariate effects from each feature using least-squares residualization."""
        if design_matrix is None or design_matrix.empty:
            return X_matrix

        design_values = np.asarray(design_matrix, dtype=float)
        design_with_intercept = np.column_stack([
            np.ones((design_values.shape[0], 1), dtype=float),
            design_values,
        ])

        if design_with_intercept.shape[0] <= design_with_intercept.shape[1]:
            return X_matrix

        try:
            coefficients, _, _, _ = np.linalg.lstsq(design_with_intercept, X_matrix, rcond=None)
            return X_matrix - design_with_intercept @ coefficients
        except Exception:
            return X_matrix
    
    def prepare_data(self, test_size: float = 0.3, scaling_method: str = 'standard',
                    random_state: int = 42, imputation_method: str = 'mean',
                    return_indices: bool = False) -> Tuple[np.ndarray, np.ndarray,
                                                      np.ndarray, np.ndarray]:
        """
        Prepare data for ML: transpose, split, and scale.
        
        Args:
            test_size: Proportion of data to use for testing (0.1 to 0.5)
            scaling_method: 'standard', 'robust', or 'none'
            random_state: Random seed for reproducibility
            
        Returns:
            X_train, X_test, y_train, y_test
        """
        logger.info("Preparing data for ML analysis...")
        
        # Get feature data and transpose (samples as rows, features as columns)
        X = self.data_df[self.feature_columns].T.values
        
        # Get group labels for each sample
        y_labels = [self.group_assignments.get(col, 'Unknown') for col in self.feature_columns]
        
        # Encode labels to numeric values
        y = self.label_encoder.fit_transform(y_labels)
        
        # Check for missing values and handle them when requested.
        if np.isnan(X).any() and str(imputation_method).lower() == 'mean':
            logger.warning("Missing values detected. Filling with column means.")
            col_means = np.nanmean(X, axis=0)
            inds = np.where(np.isnan(X))
            X[inds] = np.take(col_means, inds[1])
        
        # Train-test split (or use all data if test_size is 0 or None)
        sample_indices = np.arange(len(self.feature_columns))
        if test_size is None or test_size == 0.0 or test_size == 0:
            # No split - use all data for training (for PCA, LDA, etc.)
            train_indices = sample_indices
            test_indices = np.array([], dtype=int)
            X_train = X[train_indices]
            X_test = np.empty((0, X.shape[1]))
            y_train = y[train_indices]
            y_test = np.array([])
            logger.info(f"Using all {len(X_train)} samples (no train-test split)")
        else:
            # Train-test split with stratification
            try:
                train_indices, test_indices = train_test_split(
                    sample_indices,
                    test_size=test_size,
                    stratify=y,
                    random_state=random_state,
                )
                X_train = X[train_indices]
                X_test = X[test_indices]
                y_train = y[train_indices]
                y_test = y[test_indices]
                logger.info(f"Split data: {len(X_train)} train, {len(X_test)} test samples")
            except ValueError as e:
                # Fail fast instead of silently using a potentially biased random split.
                err_msg = (
                    "Stratified split failed because one or more classes are too small "
                    "for the requested test size or CV structure. Reduce the number "
                    "of classes, increase sample count, or adjust test/CV settings. "
                    f"Original error: {e}"
                )
                logger.error(err_msg)
                raise ValueError(err_msg) from e
        
        # Apply scaling
        if scaling_method == 'standard':
            self.scaler = StandardScaler()
            logger.info("Applying standard scaling")
        elif scaling_method == 'robust':
            self.scaler = RobustScaler()
            logger.info("Applying robust scaling")
        else:
            self.scaler = None
            logger.info("No scaling applied")
        
        if self.scaler:
            X_train = self.scaler.fit_transform(X_train)
            if X_test.size > 0:
                X_test = self.scaler.transform(X_test)

        if return_indices:
            return X_train, X_test, y_train, y_test, train_indices, test_indices
        
        return X_train, X_test, y_train, y_test

    def _build_cv_splitter(self, cv_folds: int, random_state: int,
                           use_repeated_cv: bool = False,
                           cv_repeats: int = 3):
        """Create CV splitter with optional repeated stratified folds."""
        if use_repeated_cv:
            return RepeatedStratifiedKFold(
                n_splits=cv_folds,
                n_repeats=max(int(cv_repeats or 1), 1),
                random_state=random_state,
            )
        return StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)

    def _build_classification_scoring(self, n_classes: int) -> Dict[str, str]:
        """Build robust multi-metric scoring dictionary for cross-validation."""
        scoring = {
            'accuracy': 'accuracy',
            'balanced_accuracy': 'balanced_accuracy',
            'f1_macro': 'f1_macro',
            'f1_weighted': 'f1_weighted',
            'precision_macro': 'precision_macro',
            'recall_macro': 'recall_macro',
        }
        if n_classes == 2:
            scoring['roc_auc'] = 'roc_auc'
        elif n_classes > 2:
            scoring['roc_auc_ovr'] = 'roc_auc_ovr_weighted'
        return scoring

    def _build_estimator(self, model_name: str, hyperparameters: Optional[Dict],
                         regularization_type: str, C: float, max_iter: int,
                         random_state: int, class_weight: Optional[str],
                         calibration_method: Optional[str] = None):
        """Build estimator instance with user/model defaults."""
        if model_name not in self.CLASSIFICATION_MODELS:
            raise ValueError(f"Unknown model_name: {model_name}")

        model_class = self.CLASSIFICATION_MODELS[model_name]

        if hyperparameters:
            params = hyperparameters.copy()
        else:
            params = self.DEFAULT_PARAMS.get(model_name, {}).copy()

            if model_name in ['Logistic Regression', 'SVM (RBF)']:
                params['C'] = C
                params['max_iter'] = max_iter
                if model_name == 'Logistic Regression':
                    params['penalty'] = regularization_type
                    if regularization_type == 'elasticnet':
                        params['solver'] = 'saga'
                        params['l1_ratio'] = 0.5
                    elif regularization_type == 'l1':
                        params['solver'] = 'saga'
                    else:
                        params['solver'] = 'lbfgs'

        if 'random_state' in params:
            params['random_state'] = random_state

        if class_weight == 'balanced' and model_name in ['Random Forest', 'SVM (RBF)', 'Logistic Regression']:
            params['class_weight'] = 'balanced'

        # Calibrated SVM option for more stable probabilities in small datasets.
        if model_name == 'SVM (RBF)' and calibration_method in {'sigmoid', 'isotonic'}:
            params['probability'] = False
            base_model = model_class(**params)
            return CalibratedClassifierCV(base_model, method=calibration_method, cv=3), params

        if model_name == 'SVM (RBF)':
            params['probability'] = True

        return model_class(**params), params

    def _build_feature_selector(self, feature_selection_method: str = 'none',
                                univariate_k: int = 50,
                                lasso_C: float = 0.1,
                                rfe_n_features: int = 50,
                                random_state: int = 42):
        """Create optional feature-selection module for training-time selection."""
        method = str(feature_selection_method or 'none').strip().lower()
        if method in {'none', ''}:
            return None
        if method == 'univariate':
            return SelectKBest(score_func=f_classif, k=max(1, int(univariate_k)))
        if method == 'lasso':
            lasso_est = LogisticRegression(
                penalty='l1',
                C=float(lasso_C),
                solver='saga',
                max_iter=5000,
                random_state=random_state,
            )
            return SelectFromModel(lasso_est)
        if method == 'rf_rfe':
            rfe_est = RandomForestClassifier(n_estimators=200, random_state=random_state, n_jobs=-1)
            return RFE(estimator=rfe_est, n_features_to_select=max(1, int(rfe_n_features)), step=0.1)
        return None

    def _extract_selected_feature_names(self, pipeline: Pipeline, feature_names: List[str]) -> List[str]:
        """Get selected feature names from a fitted pipeline selector, if present."""
        try:
            if pipeline is None or not hasattr(pipeline, 'named_steps'):
                return list(feature_names)
            selector = pipeline.named_steps.get('feature_selector')
            if selector is None or selector == 'passthrough':
                return list(feature_names)
            if hasattr(selector, 'get_support'):
                mask = selector.get_support()
                if mask is not None and len(mask) == len(feature_names):
                    return [name for name, keep in zip(feature_names, mask) if keep]
        except Exception:
            pass
        return list(feature_names)

    def _mean_std_ci95(self, values: List[float]) -> Dict[str, Optional[float]]:
        """Compute mean/std and approximate 95% CI (normal approximation)."""
        if not values:
            return {'mean': None, 'std': None, 'ci95_low': None, 'ci95_high': None}
        arr = np.asarray(values, dtype=float)
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
        se = (std / np.sqrt(arr.size)) if arr.size > 0 else 0.0
        margin = 1.96 * se
        return {
            'mean': mean,
            'std': std,
            'ci95_low': float(mean - margin),
            'ci95_high': float(mean + margin),
        }

    def _build_pipeline(self, estimator, scaling_method: str,
                        model_name: Optional[str] = None,
                        auto_skip_scaling_for_trees: bool = False,
                        feature_selector=None) -> Pipeline:
        """Build leakage-safe pipeline (scaler -> model)."""
        scaler = 'passthrough'
        if auto_skip_scaling_for_trees and model_name in {'Random Forest', 'Gradient Boosting'}:
            scaler = 'passthrough'
        elif scaling_method == 'standard':
            scaler = StandardScaler()
        elif scaling_method == 'robust':
            scaler = RobustScaler()

        selector_step = feature_selector if feature_selector is not None else 'passthrough'

        return Pipeline([
            ('scaler', scaler),
            ('feature_selector', selector_step),
            ('model', estimator),
        ])

    def _build_search_space(self, model_name: str,
                            regularization_type: str = 'l2') -> Dict[str, List[Any]]:
        """Search space for model tuning (applied on training data only)."""
        if model_name == 'Random Forest':
            return {
                'model__n_estimators': [100, 300, 500],
                'model__max_depth': [None, 5, 10, 20],
                'model__min_samples_split': [2, 5, 10],
                'model__max_features': ['sqrt', 'log2', None],
            }

        if model_name == 'SVM (RBF)':
            return {
                'model__C': [0.1, 1.0, 10.0, 100.0],
                'model__gamma': ['scale', 'auto', 0.01, 0.1, 1.0],
            }

        if model_name == 'Logistic Regression':
            grid = {'model__C': [0.01, 0.1, 1.0, 10.0]}
            if regularization_type in {'l1', 'l2', 'elasticnet'}:
                grid['model__penalty'] = [regularization_type]
            return grid

        if model_name == 'Gradient Boosting':
            return {
                'model__n_estimators': [100, 200, 400],
                'model__learning_rate': [0.01, 0.05, 0.1, 0.2],
                'model__max_depth': [2, 3, 5],
            }

        return {}
    
    def run_classification(self, model_name: str, hyperparameters: Optional[Dict] = None,
                          cv_folds: int = 5, test_size: float = 0.3,
                          scaling_method: str = 'standard',
                          regularization_type: str = 'l2', C: float = 1.0,
                          max_iter: int = 1000, random_state: int = 42,
                          class_weight: Optional[str] = None,
                          repeated_runs: int = 1,
                          tune_hyperparameters: bool = True,
                          tuning_strategy: str = 'grid',
                          tuning_iter: int = 20,
                          use_repeated_cv: bool = False,
                          cv_repeats: int = 3,
                          nested_cv: bool = False,
                          calibration_method: Optional[str] = None,
                          permutation_test_runs: int = 0,
                          imputation_method: str = 'half_min',
                          imputation_knn_neighbors: int = 5,
                          auto_skip_scaling_for_trees: bool = False,
                          feature_selection_method: str = 'none',
                          variance_percentile: float = 10.0,
                          univariate_k: int = 50,
                          lasso_C: float = 0.1,
                          rfe_n_features: int = 50,
                          stability_tracking: bool = False,
                          stability_threshold: float = 70.0,
                          stability_top_n: int = 20,
                          return_cv_predictions: bool = False) -> Dict[str, Any]:
        """
        Run classification model with cross-validation.
        
        Args:
            model_name: Name of the model from CLASSIFICATION_MODELS
            hyperparameters: Optional custom hyperparameters
            cv_folds: Number of cross-validation folds
            test_size: Proportion of data for testing
            scaling_method: Feature scaling method
            regularization_type: 'l1', 'l2', or 'elasticnet' (for linear models)
            C: Regularization strength (lower = stronger)
            max_iter: Maximum iterations for convergence
            random_state: Random seed for reproducibility
            class_weight: Class weighting strategy (None or 'balanced')
            
        Returns:
            Dictionary containing all results and metrics
        """
        logger.info(f"Running {model_name} classification...")

        repeated_runs = max(int(repeated_runs or 1), 1)

        # Full dataset (used for optional nested CV diagnostics).
        X_full = self.data_df[self.feature_columns].T.values
        y_labels_full = [self.group_assignments.get(col, 'Unknown') for col in self.feature_columns]
        y_full = self.label_encoder.fit_transform(y_labels_full)
        if np.isnan(X_full).any():
            col_means = np.nanmean(X_full, axis=0)
            inds = np.where(np.isnan(X_full))
            X_full[inds] = np.take(col_means, inds[1])

        if len(np.unique(y_full)) < 2:
            raise ValueError("At least two classes are required for classification.")

        if test_size is not None and float(test_size) > 0:
            class_counts_full = np.bincount(y_full)
            min_full = int(class_counts_full.min()) if len(class_counts_full) else 0
            min_required = max(2, int(np.ceil(1.0 / float(test_size))), int(np.ceil(1.0 / (1.0 - float(test_size)))))
            if min_full < min_required:
                raise ValueError(
                    "Insufficient per-class samples for requested train/test split. "
                    f"Minimum class size={min_full}, required at least {min_required} for test_size={test_size}."
                )

        scoring = self._build_classification_scoring(len(np.unique(y_full)))
        primary_scoring = 'roc_auc' if 'roc_auc' in scoring else 'balanced_accuracy'

        run_results: List[Dict[str, Any]] = []

        for run_idx in range(repeated_runs):
            run_seed = random_state + run_idx
            X_train, X_test, y_train, y_test, train_idx, test_idx = self.prepare_data(
                test_size=test_size,
                scaling_method='none',  # Pipeline handles scaling.
                random_state=run_seed,
                imputation_method='none',
                return_indices=True,
            )

            # Optional imputation (reuse existing Statistics implementation).
            imputation_info_train = {'applied': False, 'method': 'none'}
            imputation_info_test = {'applied': False, 'method': 'none'}
            if str(imputation_method).strip().lower() not in {'none', ''}:
                train_cols = [f"train_{i}" for i in range(X_train.shape[0])]
                test_cols = [f"test_{i}" for i in range(X_test.shape[0])]
                train_df = pd.DataFrame(X_train.T, columns=train_cols)
                test_df = pd.DataFrame(X_test.T, columns=test_cols)

                train_df_imp, imputation_info_train = apply_imputation(
                    train_df,
                    sample_cols=train_cols,
                    method=imputation_method,
                    knn_neighbors=imputation_knn_neighbors,
                )
                X_train = train_df_imp[train_cols].T.values

                if len(X_test) > 0:
                    test_df_imp, imputation_info_test = apply_imputation(
                        test_df,
                        sample_cols=test_cols,
                        method=imputation_method,
                        knn_neighbors=imputation_knn_neighbors,
                    )
                    X_test = test_df_imp[test_cols].T.values

            # Safety fallback for unresolved missing values after optional imputation.
            if np.isnan(X_train).any():
                col_means = np.nanmean(X_train, axis=0)
                inds = np.where(np.isnan(X_train))
                X_train[inds] = np.take(col_means, inds[1])
            if len(X_test) > 0 and np.isnan(X_test).any():
                col_means = np.nanmean(X_train, axis=0)
                inds = np.where(np.isnan(X_test))
                X_test[inds] = np.take(col_means, inds[1])

            covariate_adjustment_applied = False
            if self.covariate_data is not None and self.covariate_cols:
                sample_names = list(self.feature_columns)
                train_samples = [sample_names[i] for i in train_idx]
                test_samples = [sample_names[i] for i in test_idx] if len(test_idx) > 0 else []
                train_design, fit_state = self._build_covariate_design_matrix(train_samples)
                test_design = None
                if test_samples:
                    test_design, _ = self._build_covariate_design_matrix(test_samples, fit_state=fit_state)

                if train_design is not None and not train_design.empty:
                    X_train = self._residualize_feature_matrix(X_train, train_design)
                    if len(X_test) > 0 and test_design is not None and not test_design.empty:
                        X_test = self._residualize_feature_matrix(X_test, test_design)
                    covariate_adjustment_applied = True

            # Explicit check: ensure each class appears in test set (when using holdout).
            if len(X_test) > 0:
                train_classes = set(np.unique(y_train).tolist())
                test_classes = set(np.unique(y_test).tolist())
                if test_classes != train_classes:
                    raise ValueError(
                        "Hold-out split does not contain all classes. "
                        "Try larger test_size, fewer classes, or more samples."
                    )

            n_features_before_selection = X_train.shape[1]
            feature_names_current = [str(f) for f in self.feature_names]

            # Optional variance filter (training-only threshold).
            variance_info = None
            variance_mask = None
            if str(feature_selection_method).strip().lower() == 'variance':
                row_var = np.var(X_train, axis=0)
                var_pct = float(np.clip(variance_percentile, 0.0, 100.0))
                threshold = float(np.percentile(row_var, var_pct)) if len(row_var) else 0.0
                variance_mask = row_var > threshold
                if not np.any(variance_mask) and len(row_var) > 0:
                    variance_mask = row_var == np.max(row_var)

                X_train = X_train[:, variance_mask]
                if len(X_test) > 0:
                    X_test = X_test[:, variance_mask]
                feature_names_current = [name for name, keep in zip(feature_names_current, variance_mask) if keep]
                variance_info = {
                    'variance_percentile': var_pct,
                    'variance_threshold': threshold,
                    'removed_low_variance_features': int((~variance_mask).sum()) if variance_mask is not None else 0,
                }

            class_counts = np.bincount(y_train)
            min_class_count = int(class_counts.min()) if len(class_counts) > 0 else 0
            if min_class_count < cv_folds:
                raise ValueError(
                    "Stratified split failed because one or more classes are too small "
                    "for the requested test size or CV structure. Reduce the number "
                    "of classes, increase sample count, or adjust test/CV settings. "
                    f"Smallest class has {min_class_count} sample(s), but cv_folds={cv_folds}."
                )

            estimator, used_params = self._build_estimator(
                model_name=model_name,
                hyperparameters=hyperparameters,
                regularization_type=regularization_type,
                C=C,
                max_iter=max_iter,
                random_state=run_seed,
                class_weight=class_weight,
                calibration_method=calibration_method,
            )
            selector = self._build_feature_selector(
                feature_selection_method=feature_selection_method,
                univariate_k=min(max(1, int(univariate_k)), X_train.shape[1]),
                lasso_C=float(lasso_C),
                rfe_n_features=min(max(1, int(rfe_n_features)), X_train.shape[1]),
                random_state=run_seed,
            )
            pipeline = self._build_pipeline(
                estimator,
                scaling_method,
                model_name=model_name,
                auto_skip_scaling_for_trees=auto_skip_scaling_for_trees,
                feature_selector=selector,
            )

            cv_splitter = self._build_cv_splitter(
                cv_folds=cv_folds,
                random_state=run_seed,
                use_repeated_cv=use_repeated_cv,
                cv_repeats=cv_repeats,
            )

            best_params = None
            search_summary = None
            fitted_pipeline = pipeline

            if tune_hyperparameters and not hyperparameters:
                search_space = self._build_search_space(model_name, regularization_type)
                if search_space:
                    if tuning_strategy.lower() == 'random':
                        search = RandomizedSearchCV(
                            estimator=pipeline,
                            param_distributions=search_space,
                            n_iter=max(int(tuning_iter or 10), 1),
                            scoring=primary_scoring,
                            cv=cv_splitter,
                            n_jobs=-1,
                            random_state=run_seed,
                            refit=True,
                            error_score='raise',
                        )
                    else:
                        search = GridSearchCV(
                            estimator=pipeline,
                            param_grid=search_space,
                            scoring=primary_scoring,
                            cv=cv_splitter,
                            n_jobs=-1,
                            refit=True,
                            error_score='raise',
                        )

                    search.fit(X_train, y_train)
                    fitted_pipeline = search.best_estimator_
                    best_params = search.best_params_
                    search_summary = {
                        'best_score': float(search.best_score_),
                        'best_params': best_params,
                        'strategy': tuning_strategy.lower(),
                    }
                else:
                    fitted_pipeline.fit(X_train, y_train)
            else:
                fitted_pipeline.fit(X_train, y_train)

            # Cross-validation metrics on training split only.
            cv_results = cross_validate(
                estimator=clone(fitted_pipeline),
                X=X_train,
                y=y_train,
                cv=cv_splitter,
                scoring=scoring,
                n_jobs=-1,
                error_score='raise',
            )

            cv_metrics_mean = {}
            cv_metrics_std = {}
            for metric_key in scoring.keys():
                vals = np.asarray(cv_results.get(f'test_{metric_key}', []), dtype=float)
                if vals.size > 0:
                    cv_metrics_mean[metric_key] = float(np.mean(vals))
                    cv_metrics_std[metric_key] = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
                else:
                    cv_metrics_mean[metric_key] = None
                    cv_metrics_std[metric_key] = None

            cv_scores = np.asarray(cv_results.get('test_accuracy', []), dtype=float)

            # Final fit on training split.
            if not hasattr(fitted_pipeline, 'predict'):
                fitted_pipeline.fit(X_train, y_train)

            self.model = fitted_pipeline
            selected_feature_names = self._extract_selected_feature_names(fitted_pipeline, feature_names_current)
            n_features_after_selection = len(selected_feature_names)

            y_train_pred = fitted_pipeline.predict(X_train)
            train_acc = float(accuracy_score(y_train, y_train_pred))

            class_labels = self.label_encoder.classes_
            auc_score = None
            test_acc = None
            test_bal_acc = None
            test_f1_macro = None
            test_f1_weighted = None
            test_precision_macro = None
            test_recall_macro = None
            conf_matrix = None

            if len(X_test) > 0:
                y_test_pred = fitted_pipeline.predict(X_test)
                y_test_pred_proba = fitted_pipeline.predict_proba(X_test) if hasattr(fitted_pipeline, 'predict_proba') else None

                test_acc = float(accuracy_score(y_test, y_test_pred))
                test_bal_acc = float(balanced_accuracy_score(y_test, y_test_pred))
                test_f1_macro = float(f1_score(y_test, y_test_pred, average='macro', zero_division=0))
                test_f1_weighted = float(f1_score(y_test, y_test_pred, average='weighted', zero_division=0))
                test_precision_macro = float(precision_score(y_test, y_test_pred, average='macro', zero_division=0))
                test_recall_macro = float(recall_score(y_test, y_test_pred, average='macro', zero_division=0))
                conf_matrix = confusion_matrix(y_test, y_test_pred)

                if y_test_pred_proba is not None:
                    try:
                        if len(class_labels) == 2:
                            auc_score = float(roc_auc_score(y_test, y_test_pred_proba[:, 1]))
                        else:
                            auc_score = float(roc_auc_score(y_test, y_test_pred_proba, multi_class='ovr', average='macro'))
                    except Exception as e:
                        logger.warning(f"Could not compute single-model AUC: {e}")

                class_report = classification_report(
                    y_test,
                    y_test_pred,
                    target_names=class_labels,
                    output_dict=True,
                    zero_division=0,
                )
                class_report_str = classification_report(
                    y_test,
                    y_test_pred,
                    target_names=class_labels,
                    zero_division=0,
                )
            else:
                y_test_pred = np.array([])
                y_test_pred_proba = None
                class_report = {}
                class_report_str = "CV-only mode: No held-out test set. Use CV metrics as performance estimate."

            cv_predictions = None
            cv_pred_proba = None
            if return_cv_predictions:
                try:
                    cv_predictions = cross_val_predict(
                        estimator=clone(fitted_pipeline),
                        X=X_train,
                        y=y_train,
                        cv=cv_splitter,
                        method='predict',
                        n_jobs=-1,
                    )
                    if hasattr(fitted_pipeline, 'predict_proba'):
                        cv_pred_proba = cross_val_predict(
                            estimator=clone(fitted_pipeline),
                            X=X_train,
                            y=y_train,
                            cv=cv_splitter,
                            method='predict_proba',
                            n_jobs=-1,
                        )
                except Exception as e:
                    logger.warning(f"Could not compute CV predictions: {e}")

            overfitting_gap = None
            overfitting_flag = False
            if test_acc is not None:
                overfitting_gap = float(train_acc - test_acc)
                overfitting_flag = overfitting_gap > 0.10

            permutation_summary = None
            if int(permutation_test_runs or 0) > 0:
                try:
                    score, perm_scores, pvalue = permutation_test_score(
                        estimator=clone(fitted_pipeline),
                        X=X_train,
                        y=y_train,
                        scoring=primary_scoring,
                        cv=cv_splitter,
                        n_permutations=int(permutation_test_runs),
                        n_jobs=-1,
                        random_state=run_seed,
                    )
                    permutation_summary = {
                        'metric': primary_scoring,
                        'score': float(score),
                        'pvalue': float(pvalue),
                        'mean_permuted_score': float(np.mean(perm_scores)) if len(perm_scores) else None,
                        'n_permutations': int(permutation_test_runs),
                    }
                except Exception as e:
                    permutation_summary = {'error': str(e), 'n_permutations': int(permutation_test_runs)}

            nested_cv_summary = None
            if nested_cv and tune_hyperparameters and not hyperparameters:
                try:
                    X_full_nested = X_full
                    nested_covariate_adjustment_applied = False
                    if self.covariate_data is not None and self.covariate_cols:
                        full_samples = list(self.feature_columns)
                        full_design, full_fit_state = self._build_covariate_design_matrix(full_samples)
                        if full_design is not None and not full_design.empty:
                            X_full_nested = self._residualize_feature_matrix(X_full, full_design)
                            nested_covariate_adjustment_applied = True

                    search_space_nested = self._build_search_space(model_name, regularization_type)
                    if search_space_nested:
                        base_estimator_nested, _ = self._build_estimator(
                            model_name=model_name,
                            hyperparameters=None,
                            regularization_type=regularization_type,
                            C=C,
                            max_iter=max_iter,
                            random_state=run_seed,
                            class_weight=class_weight,
                            calibration_method=calibration_method,
                        )
                        nested_selector = self._build_feature_selector(
                            feature_selection_method=feature_selection_method,
                            univariate_k=min(max(1, int(univariate_k)), X_full.shape[1]),
                            lasso_C=float(lasso_C),
                            rfe_n_features=min(max(1, int(rfe_n_features)), X_full.shape[1]),
                            random_state=run_seed,
                        )
                        pipeline_nested = self._build_pipeline(
                            base_estimator_nested,
                            scaling_method,
                            model_name=model_name,
                            auto_skip_scaling_for_trees=auto_skip_scaling_for_trees,
                            feature_selector=nested_selector,
                        )
                        inner_cv = self._build_cv_splitter(cv_folds, run_seed + 101, use_repeated_cv, cv_repeats)
                        outer_cv = self._build_cv_splitter(cv_folds, run_seed + 202, use_repeated_cv, cv_repeats)

                        if tuning_strategy.lower() == 'random':
                            nested_search = RandomizedSearchCV(
                                estimator=pipeline_nested,
                                param_distributions=search_space_nested,
                                n_iter=max(int(tuning_iter or 10), 1),
                                scoring=primary_scoring,
                                cv=inner_cv,
                                n_jobs=-1,
                                random_state=run_seed,
                                refit=True,
                                error_score='raise',
                            )
                        else:
                            nested_search = GridSearchCV(
                                estimator=pipeline_nested,
                                param_grid=search_space_nested,
                                scoring=primary_scoring,
                                cv=inner_cv,
                                n_jobs=-1,
                                refit=True,
                                error_score='raise',
                            )

                        nested_scores = cross_validate(
                            estimator=nested_search,
                            X=X_full_nested,
                            y=y_full,
                            cv=outer_cv,
                            scoring=scoring,
                            n_jobs=-1,
                            error_score='raise',
                        )

                        nested_cv_summary = {'metrics': {}}
                        for metric_key in scoring.keys():
                            vals = np.asarray(nested_scores.get(f'test_{metric_key}', []), dtype=float)
                            if vals.size > 0:
                                nested_cv_summary['metrics'][metric_key] = {
                                    'mean': float(np.mean(vals)),
                                    'std': float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                                }
                        nested_cv_summary['covariate_adjustment_applied'] = nested_covariate_adjustment_applied
                except Exception as e:
                    nested_cv_summary = {'error': str(e)}

            feature_importance = self._get_feature_importance()

            run_results.append({
                'run_index': run_idx,
                'random_state': run_seed,
                'model_name': model_name,
                'hyperparameters': used_params,
                'best_params': best_params,
                'search_summary': search_summary,
                'cv_scores': cv_scores,
                'cv_metrics': {
                    'mean': cv_metrics_mean,
                    'std': cv_metrics_std,
                },
                'train_accuracy': train_acc,
                'test_accuracy': test_acc,
                'test_balanced_accuracy': test_bal_acc,
                'test_f1_macro': test_f1_macro,
                'test_f1_weighted': test_f1_weighted,
                'test_precision_macro': test_precision_macro,
                'test_recall_macro': test_recall_macro,
                'auc': auc_score,
                'confusion_matrix': conf_matrix,
                'classification_report': class_report,
                'classification_report_str': class_report_str,
                'feature_importances': feature_importance,
                'selected_features': selected_feature_names,
                'n_features_before_selection': n_features_before_selection,
                'n_features_after_selection': n_features_after_selection,
                'feature_selection_method': feature_selection_method,
                'variance_filter': variance_info,
                'class_labels': class_labels,
                'X_train': X_train,
                'X_test': X_test,
                'y_train': y_train,
                'y_test': y_test,
                'y_train_pred': y_train_pred,
                'y_test_pred': y_test_pred,
                'y_test_pred_proba': y_test_pred_proba,
                'n_features': X_train.shape[1],
                'n_train_samples': len(X_train),
                'n_test_samples': len(X_test),
                'cv_only_mode': test_acc is None,
                'class_weight': class_weight,
                'nested_cv': nested_cv_summary,
                'permutation_test': permutation_summary,
                'calibration_method': calibration_method,
                'cv_predictions': cv_predictions,
                'cv_pred_proba': cv_pred_proba,
                'overfitting_gap': overfitting_gap,
                'overfitting_flag': overfitting_flag,
                'preprocessing': {
                    'scaling_method': scaling_method,
                    'auto_skip_scaling_for_trees': bool(auto_skip_scaling_for_trees),
                    'imputation_method': imputation_method,
                    'imputation_knn_neighbors': int(imputation_knn_neighbors),
                    'imputation_train': imputation_info_train,
                    'imputation_test': imputation_info_test,
                    'split_seed': run_seed,
                    'train_class_counts': {str(k): int(v) for k, v in zip(*np.unique(y_train, return_counts=True))},
                    'test_class_counts': {str(k): int(v) for k, v in zip(*np.unique(y_test, return_counts=True))} if len(y_test) > 0 else {},
                    'n_features_before_selection': int(n_features_before_selection),
                    'n_features_after_selection': int(n_features_after_selection),
                },
            })

        if not run_results:
            raise ValueError("No successful classification run was produced.")

        # Aggregate repeated runs.
        cv_acc_means = [rr['cv_metrics']['mean'].get('accuracy') for rr in run_results if rr['cv_metrics']['mean'].get('accuracy') is not None]
        cv_bal_means = [rr['cv_metrics']['mean'].get('balanced_accuracy') for rr in run_results if rr['cv_metrics']['mean'].get('balanced_accuracy') is not None]
        cv_auc_key = 'roc_auc' if any(rr['cv_metrics']['mean'].get('roc_auc') is not None for rr in run_results) else 'roc_auc_ovr'
        cv_auc_means = [rr['cv_metrics']['mean'].get(cv_auc_key) for rr in run_results if rr['cv_metrics']['mean'].get(cv_auc_key) is not None]

        test_accs = [rr['test_accuracy'] for rr in run_results if rr['test_accuracy'] is not None]
        test_aucs = [rr['auc'] for rr in run_results if rr['auc'] is not None]
        test_bal_accs = [rr['test_balanced_accuracy'] for rr in run_results if rr['test_balanced_accuracy'] is not None]

        def _mean_std(values: List[float]) -> Tuple[Optional[float], Optional[float]]:
            if not values:
                return None, None
            arr = np.asarray(values, dtype=float)
            return float(np.mean(arr)), float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0

        cv_acc_mean, cv_acc_std = _mean_std(cv_acc_means)
        cv_bal_mean, cv_bal_std = _mean_std(cv_bal_means)
        cv_auc_mean, cv_auc_std = _mean_std(cv_auc_means)
        test_acc_mean, test_acc_std = _mean_std(test_accs)
        test_auc_mean, test_auc_std = _mean_std(test_aucs)
        test_bal_mean, test_bal_std = _mean_std(test_bal_accs)

        # Confidence intervals for repeated metrics.
        test_accuracy_ci95 = self._mean_std_ci95(test_accs)
        test_auc_ci95 = self._mean_std_ci95(test_aucs)
        test_balanced_accuracy_ci95 = self._mean_std_ci95(test_bal_accs)

        # Aggregate feature stability across repeated runs.
        stable_selected_features = []
        stable_important_features = []
        selected_frequency = {}
        important_frequency = {}
        if stability_tracking and run_results:
            n_runs = len(run_results)
            selected_counter: Dict[str, int] = {}
            important_counter: Dict[str, int] = {}
            top_n = max(1, int(stability_top_n))

            for rr in run_results:
                for feat in rr.get('selected_features', []) or []:
                    feat_s = str(feat)
                    selected_counter[feat_s] = selected_counter.get(feat_s, 0) + 1

                fi = rr.get('feature_importances') or {}
                top_features = fi.get('top_features', [])[:top_n]
                for feat, _ in top_features:
                    feat_s = str(feat)
                    important_counter[feat_s] = important_counter.get(feat_s, 0) + 1

            threshold_pct = float(stability_threshold)
            for feat, count in selected_counter.items():
                freq = 100.0 * count / n_runs
                selected_frequency[feat] = {'count': int(count), 'frequency_pct': float(freq)}
                if freq >= threshold_pct:
                    stable_selected_features.append((feat, float(freq)))

            for feat, count in important_counter.items():
                freq = 100.0 * count / n_runs
                important_frequency[feat] = {'count': int(count), 'frequency_pct': float(freq)}
                if freq >= threshold_pct:
                    stable_important_features.append((feat, float(freq)))

            stable_selected_features.sort(key=lambda x: x[1], reverse=True)
            stable_important_features.sort(key=lambda x: x[1], reverse=True)

        if test_accs:
            best_idx = int(np.argmax([rr['test_accuracy'] if rr['test_accuracy'] is not None else -np.inf for rr in run_results]))
        else:
            best_idx = int(np.argmax([rr['cv_metrics']['mean'].get('accuracy', -np.inf) for rr in run_results]))

        best_run = run_results[best_idx]
        self.model = self.model if self.model is not None else None

        self.results = best_run.copy()
        self.results.update({
            'repeated_runs': repeated_runs,
            'run_results': run_results,
            'cv_scores': best_run.get('cv_scores', np.array([])),
            'cv_mean_accuracy': cv_acc_mean,
            'cv_std_accuracy': cv_acc_std,
            'cv_mean_balanced_accuracy': cv_bal_mean,
            'cv_std_balanced_accuracy': cv_bal_std,
            'cv_mean_auc': cv_auc_mean,
            'cv_std_auc': cv_auc_std,
            'test_accuracy': test_acc_mean,
            'test_accuracy_std': test_acc_std,
            'test_balanced_accuracy': test_bal_mean,
            'test_balanced_accuracy_std': test_bal_std,
            'auc': test_auc_mean,
            'auc_std': test_auc_std,
            'test_accuracy_ci95': test_accuracy_ci95,
            'test_auc_ci95': test_auc_ci95,
            'test_balanced_accuracy_ci95': test_balanced_accuracy_ci95,
            'best_run_index': best_idx,
            'cv_metric_name_for_auc': cv_auc_key,
            'cv_folds': int(cv_folds),
            'test_size': float(test_size) if test_size is not None else None,
            'scaling_method': scaling_method,
            'split_seed': int(random_state),
            'hyperparameter_tuning': {
                'enabled': bool(tune_hyperparameters and not hyperparameters),
                'strategy': tuning_strategy.lower(),
                'iterations': int(tuning_iter),
            },
            'nested_cv_enabled': bool(nested_cv),
            'use_repeated_cv': bool(use_repeated_cv),
            'cv_repeats': int(cv_repeats),
            'calibration_method': calibration_method,
            'imputation_method': imputation_method,
            'imputation_knn_neighbors': int(imputation_knn_neighbors),
            'auto_skip_scaling_for_trees': bool(auto_skip_scaling_for_trees),
            'feature_selection_method': feature_selection_method,
            'variance_percentile': float(variance_percentile),
            'univariate_k': int(univariate_k),
            'lasso_C': float(lasso_C),
            'rfe_n_features': int(rfe_n_features),
            'stability_tracking': bool(stability_tracking),
            'stability_threshold': float(stability_threshold),
            'selected_feature_frequency': selected_frequency,
            'important_feature_frequency': important_frequency,
            'stable_selected_features': stable_selected_features,
            'stable_important_features': stable_important_features,
            'selected_features': best_run.get('selected_features', []),
            'n_features_before_selection': best_run.get('n_features_before_selection', best_run.get('n_features')),
            'n_features_after_selection': best_run.get('n_features_after_selection', best_run.get('n_features')),
            'overfitting_warning': any(bool(rr.get('overfitting_flag')) for rr in run_results),
            'overfitting_gap_mean': float(np.mean([rr['overfitting_gap'] for rr in run_results if rr.get('overfitting_gap') is not None])) if any(rr.get('overfitting_gap') is not None for rr in run_results) else None,
            'covariate_adjustment_applied': covariate_adjustment_applied,
            'covariate_columns': list(self.covariate_cols),
        })

        logger.info(
            "Classification complete. "
            f"Test accuracy(mean±std): {test_acc_mean} ± {test_acc_std}; "
            f"AUC(mean±std): {test_auc_mean} ± {test_auc_std}"
        )
        return self.results

    def run_multi_model_comparison(self, model_names: List[str], 
                                   test_size: float = 0.3, 
                                   scaling_method: str = 'standard',
                                   regularization_type: str = 'l2',
                                   C: float = 1.0,
                                   max_iter: int = 1000,
                                   random_state: int = 42,
                                   class_weight: Optional[str] = None,
                                   repeated_runs: int = 1) -> Dict[str, Any]:
        """
        Run and compare multiple classification models on the same data split.
        
        Args:
            model_names: List of model names to compare
            test_size: Proportion of data for testing (shared split)
            scaling_method: Feature scaling method
            regularization_type: 'l1', 'l2', or 'elasticnet' (for linear models)
            C: Regularization strength (lower = stronger)
            max_iter: Maximum iterations for convergence
            random_state: Random seed for reproducibility
            class_weight: Class weighting strategy (None or 'balanced')
            repeated_runs: Number of repeated comparisons to run (aggregated by mean/std)
            
        Returns:
            Dictionary containing results for all models and comparison metrics
        """
        logger.info(f"Running multi-model comparison: {', '.join(model_names)}")
        repeated_runs = max(int(repeated_runs or 1), 1)
        total_start = time.perf_counter()
        
        if test_size == 0.0 or test_size is None:
            raise ValueError("Multi-model comparison requires a held-out test set (test_size > 0)")
        
        class_labels = None
        run_results: List[Dict[str, Any]] = []

        def _build_model_params(model_name: str) -> Optional[Dict[str, Any]]:
            model_class = self.CLASSIFICATION_MODELS.get(model_name)
            if not model_class:
                return None

            params = self.DEFAULT_PARAMS.get(model_name, {}).copy()

            if model_name in ['Logistic Regression', 'SVM (RBF)']:
                params['C'] = C
                params['max_iter'] = max_iter
                if model_name == 'Logistic Regression':
                    params['penalty'] = regularization_type
                    if regularization_type == 'elasticnet':
                        params['solver'] = 'saga'
                        params['l1_ratio'] = 0.5
                    elif regularization_type == 'l1':
                        params['solver'] = 'saga'
                    else:
                        params['solver'] = 'lbfgs'

            if class_weight == 'balanced' and model_name in ['Random Forest', 'SVM (RBF)', 'Logistic Regression']:
                params['class_weight'] = 'balanced'

            return {'model_class': model_class, 'params': params}

        def _compute_metrics(y_true, y_pred, y_proba, class_labels):
            conf_matrix = confusion_matrix(y_true, y_pred)
            auc_score = None
            sensitivity_dict = {}
            specificity_dict = {}

            if len(class_labels) == 2:
                if y_proba is not None:
                    try:
                        auc_score = roc_auc_score(y_true, y_proba[:, 1])
                    except Exception as e:
                        logger.warning(f"Could not compute AUC: {e}")

                tn, fp, fn, tp = conf_matrix.ravel()
                sensitivity_dict = {str(class_labels[1]): float(tp / (tp + fn) if (tp + fn) > 0 else 0.0)}
                specificity_dict = {str(class_labels[1]): float(tn / (tn + fp) if (tn + fp) > 0 else 0.0)}
            else:
                if y_proba is not None:
                    try:
                        auc_score = roc_auc_score(y_true, y_proba, multi_class='ovr', average='macro')
                    except Exception as e:
                        logger.warning(f"Could not compute multiclass AUC: {e}")

                sensitivities = []
                specificities = []
                for i in range(len(class_labels)):
                    tp = conf_matrix[i, i]
                    fn = np.sum(conf_matrix[i, :]) - tp
                    fp = np.sum(conf_matrix[:, i]) - tp
                    tn = np.sum(conf_matrix) - tp - fn - fp
                    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
                    sensitivities.append(sens)
                    specificities.append(spec)
                    sensitivity_dict[str(class_labels[i])] = float(sens)
                    specificity_dict[str(class_labels[i])] = float(spec)

                if sensitivities:
                    sensitivity_dict['macro_avg'] = float(np.mean(sensitivities))
                if specificities:
                    specificity_dict['macro_avg'] = float(np.mean(specificities))

            return conf_matrix, auc_score, sensitivity_dict, specificity_dict

        for run_idx in range(repeated_runs):
            run_seed = random_state + run_idx
            logger.info(f"Multi-model run {run_idx + 1}/{repeated_runs} with seed={run_seed}...")
            run_start = time.perf_counter()

            X_train, X_test, y_train, y_test, train_idx, test_idx = self.prepare_data(
                test_size,
                scaling_method,
                run_seed,
                return_indices=True,
            )
            if class_labels is None:
                class_labels = self.label_encoder.classes_

            if len(X_test) == 0:
                raise ValueError("No test samples available. Cannot perform multi-model comparison.")

            if self.covariate_data is not None and self.covariate_cols:
                sample_names = list(self.feature_columns)
                train_samples = [sample_names[i] for i in train_idx]
                test_samples = [sample_names[i] for i in test_idx] if len(test_idx) > 0 else []
                train_design, fit_state = self._build_covariate_design_matrix(train_samples)
                test_design = None
                if test_samples:
                    test_design, _ = self._build_covariate_design_matrix(test_samples, fit_state=fit_state)

                if train_design is not None and not train_design.empty:
                    X_train = self._residualize_feature_matrix(X_train, train_design)
                    if len(X_test) > 0 and test_design is not None and not test_design.empty:
                        X_test = self._residualize_feature_matrix(X_test, test_design)

            run_model_results = {}
            for model_name in model_names:
                try:
                    logger.info(f"Training {model_name} (run {run_idx + 1})...")
                    model_start = time.perf_counter()

                    model_info = _build_model_params(model_name)
                    if not model_info:
                        logger.warning(f"Model {model_name} not found. Skipping.")
                        continue

                    model_class = model_info['model_class']
                    params = model_info['params']
                    if 'random_state' in params:
                        params['random_state'] = run_seed

                    model = model_class(**params)
                    model.fit(X_train, y_train)

                    feature_importance = self._extract_feature_importance_from_model(model, X_test, y_test)

                    y_test_pred = model.predict(X_test)
                    y_test_pred_proba = model.predict_proba(X_test) if hasattr(model, 'predict_proba') else None

                    test_acc = accuracy_score(y_test, y_test_pred)
                    train_acc = accuracy_score(y_train, model.predict(X_train))
                    conf_matrix, auc_score, sensitivity_dict, specificity_dict = _compute_metrics(
                        y_test, y_test_pred, y_test_pred_proba, class_labels
                    )

                    class_report = classification_report(
                        y_test, y_test_pred, target_names=class_labels, output_dict=True
                    )

                    run_model_results[model_name] = {
                        'model': model,
                        'test_accuracy': float(test_acc),
                        'train_accuracy': float(train_acc),
                        'auc': None if auc_score is None else float(auc_score),
                        'sensitivity': sensitivity_dict,
                        'specificity': specificity_dict,
                        'confusion_matrix': conf_matrix,
                        'classification_report': class_report,
                        'y_pred': y_test_pred,
                        'y_pred_proba': y_test_pred_proba,
                        'feature_importances': feature_importance,
                        'test_size': test_size,
                        'scaling_method': scaling_method,
                        'random_state': run_seed,
                    }
                    model_elapsed = time.perf_counter() - model_start
                    logger.info(f"{model_name} complete in {model_elapsed:.1f}s. Accuracy: {test_acc:.4f}, AUC: {auc_score}")

                except Exception as e:
                    logger.error(f"Error training {model_name}: {e}", exc_info=True)
                    run_model_results[model_name] = {'error': str(e)}

            run_results.append({
                'run_index': run_idx,
                'random_state': run_seed,
                'X_train': X_train,
                'X_test': X_test,
                'y_train': y_train,
                'y_test': y_test,
                'model_results': run_model_results,
                'n_train_samples': len(X_train),
                'n_test_samples': len(X_test),
            })
            logger.info(f"Run {run_idx + 1}/{repeated_runs} completed in {time.perf_counter() - run_start:.1f}s")

        # Aggregate by model across runs
        aggregate_model_results: Dict[str, Any] = {}
        comparison_data: List[Dict[str, Any]] = []
        trained_models: List[str] = []

        for model_name in model_names:
            run_entries = [
                rr['model_results'].get(model_name)
                for rr in run_results
                if rr['model_results'].get(model_name) and 'error' not in rr['model_results'].get(model_name, {})
            ]
            if not run_entries:
                aggregate_model_results[model_name] = {'error': 'All runs failed'}
                continue

            trained_models.append(model_name)
            accs = np.array([r['test_accuracy'] for r in run_entries], dtype=float)
            tr_accs = np.array([r['train_accuracy'] for r in run_entries], dtype=float)
            aucs = np.array([r['auc'] for r in run_entries if r.get('auc') is not None], dtype=float)
            importance_entries = [r.get('feature_importances') for r in run_entries if r.get('feature_importances') and r.get('feature_importances', {}).get('importances') is not None]

            class_keys = sorted({k for r in run_entries for k in r.get('sensitivity', {}).keys()})
            sensitivity_means = {}
            specificity_means = {}
            for key in class_keys:
                sens_vals = [float(r['sensitivity'][key]) for r in run_entries if key in r.get('sensitivity', {})]
                spec_vals = [float(r['specificity'][key]) for r in run_entries if key in r.get('specificity', {})]
                if sens_vals:
                    sensitivity_means[key] = float(np.mean(sens_vals))
                if spec_vals:
                    specificity_means[key] = float(np.mean(spec_vals))

            aggregated_feature_importance = None
            if importance_entries:
                importances_stack = np.vstack([np.asarray(entry['importances'], dtype=float) for entry in importance_entries])
                mean_importances = np.mean(importances_stack, axis=0)
                std_importances = np.std(importances_stack, axis=0, ddof=1) if importances_stack.shape[0] > 1 else np.zeros_like(mean_importances)
                feature_names = importance_entries[0].get('feature_names', self.feature_names)
                sorted_indices = np.argsort(mean_importances)[::-1]
                top_n = min(20, len(mean_importances))
                top_indices = sorted_indices[:top_n]
                aggregated_feature_importance = {
                    'importances': mean_importances,
                    'importances_std': std_importances,
                    'sorted_indices': sorted_indices,
                    'top_features': [(feature_names[i], float(mean_importances[i])) for i in top_indices],
                    'top_features_std': [(feature_names[i], float(std_importances[i])) for i in top_indices],
                    'feature_names': feature_names,
                    'note': f'Aggregated across {len(importance_entries)} run(s)'
                }

            aggregate_model_results[model_name] = {
                'model_name': model_name,
                'runs': run_entries,
                'n_runs': len(run_entries),
                'test_accuracy_mean': float(np.mean(accs)),
                'test_accuracy_std': float(np.std(accs, ddof=1)) if len(accs) > 1 else 0.0,
                'train_accuracy_mean': float(np.mean(tr_accs)),
                'train_accuracy_std': float(np.std(tr_accs, ddof=1)) if len(tr_accs) > 1 else 0.0,
                'auc_mean': float(np.mean(aucs)) if len(aucs) else None,
                'auc_std': float(np.std(aucs, ddof=1)) if len(aucs) > 1 else 0.0 if len(aucs) == 1 else None,
                'sensitivity': sensitivity_means,
                'specificity': specificity_means,
                'feature_importances': aggregated_feature_importance,
            }

            comparison_data.append({
                'Model': model_name,
                'Accuracy': f"{aggregate_model_results[model_name]['test_accuracy_mean']:.4f} ± {aggregate_model_results[model_name]['test_accuracy_std']:.4f}",
                'AUC': (
                    f"{aggregate_model_results[model_name]['auc_mean']:.4f} ± {aggregate_model_results[model_name]['auc_std']:.4f}"
                    if aggregate_model_results[model_name]['auc_mean'] is not None else "N/A"
                ),
                'Sensitivity': sensitivity_means,
                'Specificity': specificity_means,
            })

        first_valid = next((rr for rr in run_results if rr['model_results']), None)
        if first_valid is None:
            raise ValueError("No valid model results were produced.")
        if class_labels is None:
            class_labels = self.label_encoder.classes_ if hasattr(self.label_encoder, 'classes_') else []

        results = {
            'comparison_type': 'multi_model',
            'repeated_runs': repeated_runs,
            'random_state': random_state,
            'models_trained': trained_models,
            'model_results': aggregate_model_results,
            'comparison_table': comparison_data,
            'run_results': run_results,
            'X_train': first_valid['X_train'],
            'X_test': first_valid['X_test'],
            'y_train': first_valid['y_train'],
            'y_test': first_valid['y_test'],
            'class_labels': class_labels,
            'n_features': first_valid['X_train'].shape[1],
            'n_train_samples': first_valid['n_train_samples'],
            'n_test_samples': first_valid['n_test_samples'],
            'test_size': test_size,
            'scaling_method': scaling_method,
            'class_weight': class_weight,
        }

        logger.info(
            f"Multi-model comparison complete in {time.perf_counter() - total_start:.1f}s. "
            f"Aggregated {len(trained_models)} models across {repeated_runs} run(s)."
        )
        return results
    
    def run_pca(self, n_components: Optional[int] = None, 
                scaling_method: str = 'standard') -> Dict[str, Any]:
        """
        Run PCA for dimensionality reduction and visualization.
        
        Args:
            n_components: Number of components (None = all components)
            scaling_method: Feature scaling method
            
        Returns:
            Dictionary containing PCA results
        """
        logger.info("Running PCA analysis...")
        
        # Prepare data (no test split for PCA)
        X_train, _, y_train, _ = self.prepare_data(test_size=0.0, scaling_method=scaling_method)

        if X_train.shape[0] < 3:
            raise ValueError("PCA requires at least 3 samples.")
        if X_train.shape[1] < 2:
            raise ValueError("PCA requires at least 2 features.")
        
        # Determine number of components
        if n_components is None:
            n_components = min(X_train.shape[0], X_train.shape[1])
        
        logger.info(f"Computing {n_components} principal components...")
        
        # Run PCA
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X_train)
        
        # Get explained variance
        explained_var = pca.explained_variance_ratio_
        cumulative_var = np.cumsum(explained_var)
        
        # Get loadings (contribution of each feature to PCs)
        loadings = pca.components_
        
        # Create results dictionary
        results = {
            'pca_model': pca,
            'transformed_data': X_pca,
            'explained_variance': explained_var,
            'cumulative_variance': cumulative_var,
            'loadings': loadings,
            'groups': y_train,
            'group_labels': self.label_encoder.inverse_transform(y_train),
            'n_components': n_components,
            'feature_names': self.data_df.index.tolist() if hasattr(self.data_df.index, 'tolist') else list(range(len(loadings[0])))
        }
        
        logger.info(f"PCA complete. PC1+PC2 explain {cumulative_var[1]:.2%} of variance")
        return results
    
    def run_lda(self, scaling_method: str = 'standard') -> Dict[str, Any]:
        """
        Run Linear Discriminant Analysis for dimensionality reduction.
        
        Args:
            scaling_method: Feature scaling method
            
        Returns:
            Dictionary containing LDA results
        """
        logger.info("Running LDA analysis...")
        
        # Prepare data
        X_train, _, y_train, _ = self.prepare_data(test_size=0.0, scaling_method=scaling_method)

        if X_train.shape[0] < 3:
            raise ValueError("LDA requires at least 3 samples.")
        
        # Determine number of components (n_classes - 1)
        n_classes = len(np.unique(y_train))
        if n_classes < 2:
            raise ValueError("LDA requires at least 2 classes.")
        n_components = min(n_classes - 1, X_train.shape[1])
        
        logger.info(f"Computing {n_components} linear discriminants...")
        
        # Run LDA
        lda = LinearDiscriminantAnalysis(n_components=n_components)
        X_lda = lda.fit_transform(X_train, y_train)
        
        # Get explained variance
        explained_var = lda.explained_variance_ratio_
        
        results = {
            'lda_model': lda,
            'transformed_data': X_lda,
            'explained_variance': explained_var,
            'groups': y_train,
            'group_labels': self.label_encoder.inverse_transform(y_train),
            'n_components': n_components
        }
        
        logger.info(f"LDA complete. {n_components} discriminants computed")
        return results

    def _model_short_name(self, model_name: str) -> str:
        """Create compact folder prefix for model names."""
        mapping = {
            'Random Forest': 'RF',
            'Logistic Regression': 'LR',
            'Gradient Boosting': 'GB',
            'SVM (RBF)': 'SVM',
            'Linear Discriminant Analysis': 'LDA',
        }
        if model_name in mapping:
            return mapping[model_name]
        return ''.join(ch for ch in model_name if ch.isalnum())[:8] or 'MODEL'

    def generate_publication_figures(
        self,
        results: Dict[str, Any],
        output_root: str,
        top_n_values: Tuple[int, int] = (10, 15),
        figure_settings: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """Generate ROC, performance bar, and top-metabolite figures for single/multi model runs."""
        os.makedirs(output_root, exist_ok=True)

        def _normalize_top_n_values(values: Any) -> Tuple[int, int]:
            if not isinstance(values, (list, tuple)) or len(values) != 2:
                return (10, 15)
            try:
                cleaned = sorted({int(values[0]), int(values[1])})
            except (TypeError, ValueError):
                return (10, 15)
            if len(cleaned) != 2 or cleaned[0] <= 0:
                return (10, 15)
            return cleaned[0], cleaned[1]

        top_n_values = _normalize_top_n_values(top_n_values)

        defaults = {
            'roc_width': 8.2,
            'roc_height': 6.8,
            'roc_title_fs': 16.0,
            'roc_label_fs': 14.0,
            'roc_tick_fs': 12.0,
            'roc_legend_fs': 12.0,
            'roc_line_w': 3.0,
            'roc_axis_w': 1.6,
            'comparison_height': 6.0,
            'comparison_base_width': 3.0,
            'comparison_width_per_bar': 1.65,
            'comparison_bar_w': 0.62,
            'comparison_gap_w': 0.34,
            'comparison_error_w': 2.4,
            'comparison_title_fs': 16.0,
            'comparison_label_fs': 14.0,
            'comparison_tick_fs': 12.0,
            'comparison_value_fs': 12.0,
            'comparison_axis_w': 1.6,
            'hbar_width': 10.5,
            'hbar_base_height': 2.4,
            'hbar_height_per_feature': 0.52,
            'hbar_bar_h': 0.55,
            'hbar_gap_h': 0.28,
            'hbar_error_w': 2.2,
            'hbar_title_fs': 16.0,
            'hbar_label_fs': 14.0,
            'hbar_tick_fs': 12.0,
            'hbar_axis_w': 1.6,
        }

        fs = defaults.copy()
        if isinstance(figure_settings, dict):
            for k, v in figure_settings.items():
                if k in fs:
                    try:
                        val = float(v)
                        if val > 0:
                            fs[k] = val
                    except (TypeError, ValueError):
                        pass

        title_pad = 14

        def _style_axis(ax, tick_size: float, axis_line: float):
            ax.tick_params(axis='both', labelsize=tick_size, width=axis_line, length=6)
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontweight('bold')
            for spine in ax.spines.values():
                spine.set_linewidth(axis_line)

        def _safe_name(text: str) -> str:
            cleaned = ''.join(ch if ch.isalnum() else '_' for ch in str(text))
            while '__' in cleaned:
                cleaned = cleaned.replace('__', '_')
            return cleaned.strip('_') or 'class'

        class_labels = list(results.get('class_labels', []))
        if not class_labels and hasattr(self.label_encoder, 'classes_'):
            class_labels = list(self.label_encoder.classes_)

        # Normalize results into {model_name: [run_entries...]}
        model_runs: Dict[str, List[Dict[str, Any]]] = {}
        model_stats: Dict[str, Dict[str, Any]] = {}

        if results.get('comparison_type') == 'multi_model':
            run_results = results.get('run_results', [])
            for model_name in results.get('models_trained', []):
                entries = []
                for run in run_results:
                    mr = run.get('model_results', {}).get(model_name)
                    if mr and 'error' not in mr:
                        entries.append({
                            'y_test': run.get('y_test'),
                            'y_pred_proba': mr.get('y_pred_proba'),
                            'feature_importances': mr.get('feature_importances')
                        })
                if entries:
                    model_runs[model_name] = entries

                agg = results.get('model_results', {}).get(model_name, {})
                model_stats[model_name] = {
                    'auc_mean': agg.get('auc_mean'),
                    'auc_std': agg.get('auc_std', 0.0),
                    'acc_mean': agg.get('test_accuracy_mean'),
                    'acc_std': agg.get('test_accuracy_std', 0.0),
                    'feature_importances': agg.get('feature_importances')
                }
        else:
            model_name = results.get('model_name', 'Model')
            single_runs = []
            rr = results.get('run_results') or []
            if rr:
                for run in rr:
                    single_runs.append({
                        'y_test': run.get('y_test'),
                        'y_pred_proba': run.get('y_test_pred_proba'),
                        'feature_importances': run.get('feature_importances'),
                    })
            else:
                single_runs.append({
                    'y_test': results.get('y_test'),
                    'y_pred_proba': results.get('y_test_pred_proba'),
                    'feature_importances': results.get('feature_importances')
                })
            model_runs[model_name] = single_runs
            model_stats[model_name] = {
                'auc_mean': results.get('auc'),
                'auc_std': results.get('auc_std', 0.0),
                'acc_mean': results.get('test_accuracy'),
                'acc_std': results.get('test_accuracy_std', 0.0),
                'feature_importances': results.get('feature_importances')
            }

        # Create per-model folders and figures
        model_dirs: Dict[str, str] = {}
        binary_model_roc_data: Dict[str, Dict[str, Any]] = {}
        combined_roc_path: Optional[str] = None
        for model_name in model_runs.keys():
            short = self._model_short_name(model_name)
            model_dir = os.path.join(output_root, f"{short}_figures")
            os.makedirs(model_dir, exist_ok=True)
            model_dirs[model_name] = model_dir

            # 1) ROC per model
            runs = model_runs.get(model_name, [])
            fig, ax = plt.subplots(figsize=(fs['roc_width'], fs['roc_height']))
            mean_fpr = np.linspace(0, 1, 200)
            plotted_any = False

            if len(class_labels) == 2:
                neg_label = str(class_labels[0])
                pos_label = str(class_labels[1])
                tprs = []
                aucs = []

                for run_entry in runs:
                    y_test = run_entry.get('y_test')
                    y_proba = run_entry.get('y_pred_proba')
                    if y_test is None or y_proba is None:
                        continue

                    y_test = np.asarray(y_test)
                    y_proba = np.asarray(y_proba)
                    if y_proba.ndim != 2 or y_proba.shape[1] < 2:
                        continue

                    y_bin = (y_test == 1).astype(int)
                    if len(np.unique(y_bin)) < 2:
                        continue

                    fpr, tpr, _ = roc_curve(y_bin, y_proba[:, 1])
                    interp_tpr = np.interp(mean_fpr, fpr, tpr)
                    interp_tpr[0] = 0.0
                    tprs.append(interp_tpr)
                    try:
                        aucs.append(roc_auc_score(y_bin, y_proba[:, 1]))
                    except Exception:
                        pass

                if tprs:
                    plotted_any = True
                    mean_tpr = np.mean(tprs, axis=0)
                    mean_tpr[-1] = 1.0
                    std_tpr = np.std(tprs, axis=0) if len(tprs) > 1 else np.zeros_like(mean_tpr)
                    mean_auc = float(np.mean(aucs)) if aucs else float('nan')
                    std_auc = float(np.std(aucs, ddof=1)) if len(aucs) > 1 else 0.0
                    label = f"{pos_label} vs {neg_label} (AUC={mean_auc:.2f} $\\pm$ {std_auc:.2f})"
                    ax.plot(mean_fpr, mean_tpr, lw=fs['roc_line_w'], label=label)
                    ax.fill_between(mean_fpr, np.maximum(mean_tpr - std_tpr, 0), np.minimum(mean_tpr + std_tpr, 1), alpha=0.15)
                    binary_model_roc_data[model_name] = {
                        'mean_fpr': mean_fpr.copy(),
                        'mean_tpr': mean_tpr.copy(),
                        'std_tpr': std_tpr.copy(),
                        'mean_auc': mean_auc,
                        'std_auc': std_auc,
                        'class_label': f"{pos_label} vs {neg_label}",
                    }
            else:
                for class_idx, class_name in enumerate(class_labels):
                    tprs = []
                    aucs = []
                    for run_entry in runs:
                        y_test = run_entry.get('y_test')
                        y_proba = run_entry.get('y_pred_proba')
                        if y_test is None or y_proba is None:
                            continue
                        y_test = np.asarray(y_test)
                        y_proba = np.asarray(y_proba)
                        if y_proba.ndim != 2 or class_idx >= y_proba.shape[1]:
                            continue

                        y_bin = (y_test == class_idx).astype(int)
                        if len(np.unique(y_bin)) < 2:
                            continue

                        fpr, tpr, _ = roc_curve(y_bin, y_proba[:, class_idx])
                        interp_tpr = np.interp(mean_fpr, fpr, tpr)
                        interp_tpr[0] = 0.0
                        tprs.append(interp_tpr)
                        try:
                            aucs.append(roc_auc_score(y_bin, y_proba[:, class_idx]))
                        except Exception:
                            pass

                    if tprs:
                        plotted_any = True
                        mean_tpr = np.mean(tprs, axis=0)
                        mean_tpr[-1] = 1.0
                        std_tpr = np.std(tprs, axis=0) if len(tprs) > 1 else np.zeros_like(mean_tpr)
                        mean_auc = float(np.mean(aucs)) if aucs else float('nan')
                        std_auc = float(np.std(aucs, ddof=1)) if len(aucs) > 1 else 0.0
                        label = f"{class_name} vs rest (AUC={mean_auc:.2f} $\\pm$ {std_auc:.2f})"
                        ax.plot(mean_fpr, mean_tpr, lw=fs['roc_line_w'], label=label)
                        ax.fill_between(mean_fpr, np.maximum(mean_tpr - std_tpr, 0), np.minimum(mean_tpr + std_tpr, 1), alpha=0.15)

                        # Individual class ROC figure (multiclass only)
                        fig_cls, ax_cls = plt.subplots(figsize=(fs['roc_width'], fs['roc_height']))
                        ax_cls.plot(mean_fpr, mean_tpr, lw=fs['roc_line_w'], label=label)
                        ax_cls.fill_between(mean_fpr, np.maximum(mean_tpr - std_tpr, 0), np.minimum(mean_tpr + std_tpr, 1), alpha=0.15)
                        ax_cls.plot([0, 1], [0, 1], linestyle='--', color='gray', lw=max(1.4, fs['roc_axis_w']))
                        ax_cls.set_title(f"ROC (One-vs-Rest) - {model_name} - {class_name}", fontsize=fs['roc_title_fs'], fontweight='bold', pad=title_pad)
                        ax_cls.set_xlabel('False Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
                        ax_cls.set_ylabel('True Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
                        ax_cls.set_xlim(0.0, 1.0)
                        ax_cls.set_ylim(0.0, 1.0)
                        ax_cls.margins(x=0.0, y=0.0)
                        ax_cls.grid(alpha=0.25)
                        _style_axis(ax_cls, fs['roc_tick_fs'], fs['roc_axis_w'])
                        leg_cls = ax_cls.legend(loc='lower right', fontsize=fs['roc_legend_fs'], framealpha=0.95)
                        for txt in leg_cls.get_texts():
                            txt.set_fontweight('bold')
                        fig_cls.tight_layout()
                        class_file = f"roc_ovr_{_safe_name(class_name)}.png"
                        fig_cls.savefig(os.path.join(model_dir, class_file), dpi=300)
                        plt.close(fig_cls)

            ax.plot([0, 1], [0, 1], linestyle='--', color='gray', lw=max(1.4, fs['roc_axis_w']))
            if len(class_labels) == 2:
                ax.set_title(f"ROC - {model_name}", fontsize=fs['roc_title_fs'], fontweight='bold', pad=title_pad)
            else:
                ax.set_title(f"ROC (One-vs-Rest) - {model_name}", fontsize=fs['roc_title_fs'], fontweight='bold', pad=title_pad)
            ax.set_xlabel('False Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
            ax.set_ylabel('True Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            ax.margins(x=0.0, y=0.0)
            ax.grid(alpha=0.25)
            _style_axis(ax, fs['roc_tick_fs'], fs['roc_axis_w'])
            if plotted_any:
                leg = ax.legend(loc='lower right', fontsize=fs['roc_legend_fs'], framealpha=0.95)
                for txt in leg.get_texts():
                    txt.set_fontweight('bold')
            else:
                ax.text(0.5, 0.5, 'ROC unavailable (no probability outputs)', ha='center', va='center', fontsize=fs['roc_label_fs'], fontweight='bold')
            fig.tight_layout()
            fig.savefig(os.path.join(model_dir, 'roc_ovr_per_class.png'), dpi=300)
            plt.close(fig)

            # 3) Horizontal bars for configured top metabolites
            feat_imp = model_stats.get(model_name, {}).get('feature_importances')
            if feat_imp and feat_imp.get('top_features'):
                std_map = {feat: std for feat, std in feat_imp.get('top_features_std', [])}
                for top_n in top_n_values:
                    rows = feat_imp['top_features'][:top_n]
                    if not rows:
                        continue
                    features = [str(r[0]) for r in rows][::-1]
                    vals = [float(r[1]) for r in rows][::-1]
                    errs = [float(std_map.get(str(r[0]), 0.0)) for r in rows][::-1]

                    # Keep bar height and inter-bar spacing fixed; scale canvas to fit count.
                    bar_h = fs['hbar_bar_h']
                    gap_h = fs['hbar_gap_h']
                    unit_h = bar_h + gap_h
                    y_pos = np.arange(len(features), dtype=float) * unit_h
                    fig_h = max(6.0, fs['hbar_base_height'] + len(features) * fs['hbar_height_per_feature'])

                    fig2, ax2 = plt.subplots(figsize=(fs['hbar_width'], fig_h))
                    ax2.barh(
                        y_pos,
                        vals,
                        height=bar_h,
                        xerr=errs,
                        color='#4C78A8',
                        alpha=0.9,
                        ecolor='black',
                        capsize=5,
                        error_kw={'elinewidth': fs['hbar_error_w'], 'capthick': max(1.6, fs['hbar_error_w'] - 0.2)}
                    )
                    ax2.set_yticks(y_pos)
                    ax2.set_yticklabels(features, fontsize=fs['hbar_tick_fs'], fontweight='bold')
                    ax2.set_xlabel('Importance score', fontsize=fs['hbar_label_fs'], fontweight='bold')
                    ax2.set_title(f"Top {top_n} Metabolites - {model_name}", fontsize=fs['hbar_title_fs'], fontweight='bold', pad=title_pad)
                    ax2.set_ylim(-unit_h * 0.55, y_pos[-1] + unit_h * 0.75)
                    ax2.grid(axis='x', alpha=0.25)
                    _style_axis(ax2, fs['hbar_tick_fs'], fs['hbar_axis_w'])
                    fig2.tight_layout()
                    fig2.savefig(os.path.join(model_dir, f'top_{top_n}_metabolites.png'), dpi=300)
                    plt.close(fig2)

        # Additional combined ROC plot for binary/pairwise runs across models.
        if len(class_labels) == 2 and len(binary_model_roc_data) >= 2:
            fig_all, ax_all = plt.subplots(figsize=(fs['roc_width'], fs['roc_height']))
            color_cycle = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#72B7B2', '#B279A2']
            plotted_count = 0

            for idx, model_name in enumerate(model_runs.keys()):
                roc_info = binary_model_roc_data.get(model_name)
                if not roc_info:
                    continue

                c = color_cycle[idx % len(color_cycle)]
                mean_fpr = roc_info['mean_fpr']
                mean_tpr = roc_info['mean_tpr']
                std_tpr = roc_info['std_tpr']
                mean_auc = float(roc_info.get('mean_auc', float('nan')))
                std_auc = float(roc_info.get('std_auc', 0.0))
                label = f"{model_name} (AUC={mean_auc:.2f} $\\pm$ {std_auc:.2f})"

                ax_all.plot(mean_fpr, mean_tpr, lw=fs['roc_line_w'], color=c, label=label)
                ax_all.fill_between(
                    mean_fpr,
                    np.maximum(mean_tpr - std_tpr, 0),
                    np.minimum(mean_tpr + std_tpr, 1),
                    color=c,
                    alpha=0.10,
                )
                plotted_count += 1

            ax_all.plot([0, 1], [0, 1], linestyle='--', color='gray', lw=max(1.4, fs['roc_axis_w']))
            comparison_title = 'ROC Comparison - All Models'
            if len(class_labels) == 2:
                comparison_title = f"ROC Comparison - {class_labels[0]} vs {class_labels[1]}"
            ax_all.set_title(comparison_title, fontsize=fs['roc_title_fs'], fontweight='bold', pad=title_pad)
            ax_all.set_xlabel('False Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
            ax_all.set_ylabel('True Positive Rate', fontsize=fs['roc_label_fs'], fontweight='bold')
            ax_all.set_xlim(0.0, 1.0)
            ax_all.set_ylim(0.0, 1.0)
            ax_all.margins(x=0.0, y=0.0)
            ax_all.grid(alpha=0.25)
            _style_axis(ax_all, fs['roc_tick_fs'], fs['roc_axis_w'])

            if plotted_count > 0:
                leg_all = ax_all.legend(loc='lower right', fontsize=fs['roc_legend_fs'], framealpha=0.95)
                for txt in leg_all.get_texts():
                    txt.set_fontweight('bold')
            else:
                ax_all.text(0.5, 0.5, 'ROC unavailable (no probability outputs)', ha='center', va='center', fontsize=fs['roc_label_fs'], fontweight='bold')

            fig_all.tight_layout()
            roc_comp_name = 'model_comparison_roc.png'
            combined_roc_path = os.path.join(output_root, roc_comp_name)
            fig_all.savefig(combined_roc_path, dpi=300)
            for model_name in model_dirs:
                fig_all.savefig(os.path.join(model_dirs[model_name], roc_comp_name), dpi=300)
            plt.close(fig_all)

        # 2) Model comparison bar plots (AUC and Accuracy)
        model_names = list(model_stats.keys())

        def _plot_metric_comparison(metric_key: str, std_key: str, metric_name: str) -> None:
            means = [model_stats[m].get(metric_key) for m in model_names]
            stds = [model_stats[m].get(std_key, 0.0) for m in model_names]

            if all(v is None for v in means):
                return

            n_models = len(model_names)
            bar_w = fs['comparison_bar_w']
            gap_w = fs['comparison_gap_w']
            x = np.arange(n_models, dtype=float) * (bar_w + gap_w)
            fig_w = max(5.6, min(12.5, fs['comparison_base_width'] + n_models * fs['comparison_width_per_bar']))
            fig3, ax3 = plt.subplots(figsize=(fig_w, fs['comparison_height']))
            colors = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#B279A2']

            safe_means = [float(v) if v is not None else 0.0 for v in means]
            safe_stds = []
            for s in stds:
                try:
                    safe_stds.append(float(s) if s is not None else 0.0)
                except Exception:
                    safe_stds.append(0.0)

            bars = ax3.bar(
                x,
                safe_means,
                width=bar_w,
                yerr=safe_stds,
                capsize=6,
                error_kw={'elinewidth': fs['comparison_error_w'], 'capthick': max(1.8, fs['comparison_error_w'] - 0.2)},
                color=[colors[i % len(colors)] for i in range(len(model_names))],
                edgecolor='black',
                linewidth=1.3,
                alpha=0.9
            )

            ax3.set_xticks(x)
            label_rotation = 0 if n_models <= 3 else 12
            ax3.set_xticklabels(model_names, rotation=label_rotation, fontsize=fs['comparison_tick_fs'], fontweight='bold')
            ax3.set_ylabel(metric_name, fontsize=fs['comparison_label_fs'], fontweight='bold')
            ax3.set_title(f'Model Comparison by {metric_name} (mean $\\pm$ std)', fontsize=fs['comparison_title_fs'], fontweight='bold', pad=title_pad)
            if n_models > 0:
                ax3.set_xlim(x[0] - (bar_w * 0.7), x[-1] + (bar_w * 0.7))
            ax3.grid(axis='y', alpha=0.25)
            _style_axis(ax3, fs['comparison_tick_fs'], fs['comparison_axis_w'])

            # Put value labels above error bars and expand ylim to avoid overlap/clipping.
            top_candidates = []
            for i, b in enumerate(bars):
                m = safe_means[i]
                err = abs(safe_stds[i])
                top = m + err
                top_candidates.append(top)
                label_y = top + 0.02
                if means[i] is not None:
                    ax3.text(
                        b.get_x() + b.get_width() / 2,
                        label_y,
                        f"{m:.3f}",
                        ha='center',
                        va='bottom',
                        fontsize=fs['comparison_value_fs'],
                        fontweight='bold'
                    )

            max_top = max(top_candidates) if top_candidates else 1.0
            y_upper = max(1.0, max_top + 0.09)
            ax3.set_ylim(0.0, y_upper)

            fig3.tight_layout()
            comparison_name = f"model_comparison_{metric_name.lower()}.png"
            fig3.savefig(os.path.join(output_root, comparison_name), dpi=300)
            for model_name in model_dirs:
                fig3.savefig(os.path.join(model_dirs[model_name], comparison_name), dpi=300)
            plt.close(fig3)

        if model_names:
            _plot_metric_comparison('auc_mean', 'auc_std', 'AUC')
            _plot_metric_comparison('acc_mean', 'acc_std', 'Accuracy')

        # 4) Multi-model comparison Venn for top metabolites + Excel detail sheets
        if len(model_names) in (2, 3):
            venn_excel_path = os.path.join(output_root, 'comparison_venn_metabolites.xlsx')
            venn_sheet_written = False
            with pd.ExcelWriter(venn_excel_path, engine='openpyxl') as writer:
                for top_n in top_n_values:
                    per_model_top = {}
                    per_model_rank = {}
                    for model_name in model_names:
                        fi = model_stats.get(model_name, {}).get('feature_importances') or {}
                        top_rows = fi.get('top_features', [])[:top_n]
                        names = [str(r[0]) for r in top_rows]
                        per_model_top[model_name] = set(names)
                        per_model_rank[model_name] = {
                            str(r[0]): (idx + 1, float(r[1]))
                            for idx, r in enumerate(top_rows)
                        }

                    if not all(per_model_top.values()):
                        continue

                    # Build detailed membership table
                    all_metabs = sorted(set().union(*per_model_top.values()))
                    rows = []
                    common_all = set.intersection(*per_model_top.values())
                    for metab in all_metabs:
                        present_models = [m for m in model_names if metab in per_model_top[m]]
                        count_present = len(present_models)
                        if count_present == len(model_names):
                            category = 'Common to all models'
                        elif count_present == 1:
                            category = f"Unique to {present_models[0]}"
                        else:
                            category = 'Shared (partial overlap)'

                        row = {
                            'Metabolite': metab,
                            'Models_Count': count_present,
                            'Models_List': '; '.join(present_models),
                            'Category': category,
                        }
                        for model_name in model_names:
                            row[f'In_{model_name}'] = 1 if metab in per_model_top[model_name] else 0
                            rank_info = per_model_rank[model_name].get(metab)
                            if rank_info is None:
                                row[f'Rank_{model_name}'] = ''
                                row[f'Importance_{model_name}'] = ''
                            else:
                                row[f'Rank_{model_name}'] = rank_info[0]
                                row[f'Importance_{model_name}'] = rank_info[1]
                        rows.append(row)

                    details_df = pd.DataFrame(rows).sort_values(['Models_Count', 'Metabolite'], ascending=[False, True])
                    details_df.to_excel(writer, sheet_name=f'Top{top_n}_Details', index=False)
                    venn_sheet_written = True

                    # Draw comparison venn
                    fig_v, ax_v = plt.subplots(figsize=(7.0, 6.2))
                    ax_v.set_aspect('equal')
                    ax_v.axis('off')

                    if len(model_names) == 2:
                        m1, m2 = model_names
                        s1, s2 = per_model_top[m1], per_model_top[m2]
                        only1 = len(s1 - s2)
                        only2 = len(s2 - s1)
                        both = len(s1 & s2)

                        c1 = plt.Circle((-0.65, 0.0), 1.0, color='#4C78A8', alpha=0.35, ec='black', lw=2)
                        c2 = plt.Circle((0.65, 0.0), 1.0, color='#F58518', alpha=0.35, ec='black', lw=2)
                        ax_v.add_patch(c1)
                        ax_v.add_patch(c2)

                        ax_v.text(-1.05, 0.0, str(only1), ha='center', va='center', fontsize=16, fontweight='bold')
                        ax_v.text(1.05, 0.0, str(only2), ha='center', va='center', fontsize=16, fontweight='bold')
                        ax_v.text(0.0, 0.0, str(both), ha='center', va='center', fontsize=16, fontweight='bold')
                        ax_v.text(-0.95, 1.05, m1, ha='center', va='center', fontsize=11, fontweight='bold')
                        ax_v.text(0.95, 1.05, m2, ha='center', va='center', fontsize=11, fontweight='bold')
                        ax_v.set_xlim(-2.1, 2.1)
                        ax_v.set_ylim(-1.4, 1.6)
                    else:
                        m1, m2, m3 = model_names
                        s1, s2, s3 = per_model_top[m1], per_model_top[m2], per_model_top[m3]
                        c1 = plt.Circle((-0.78, 0.35), 1.0, color='#4C78A8', alpha=0.32, ec='black', lw=2)
                        c2 = plt.Circle((0.78, 0.35), 1.0, color='#F58518', alpha=0.32, ec='black', lw=2)
                        c3 = plt.Circle((0.0, -0.62), 1.0, color='#54A24B', alpha=0.32, ec='black', lw=2)
                        ax_v.add_patch(c1)
                        ax_v.add_patch(c2)
                        ax_v.add_patch(c3)

                        n100 = len(s1 - s2 - s3)
                        n010 = len(s2 - s1 - s3)
                        n001 = len(s3 - s1 - s2)
                        n110 = len((s1 & s2) - s3)
                        n101 = len((s1 & s3) - s2)
                        n011 = len((s2 & s3) - s1)
                        n111 = len(s1 & s2 & s3)

                        ax_v.text(-1.28, 0.55, str(n100), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(1.28, 0.55, str(n010), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(0.0, -1.25, str(n001), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(0.0, 0.72, str(n110), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(-0.55, -0.28, str(n101), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(0.55, -0.28, str(n011), ha='center', va='center', fontsize=14, fontweight='bold')
                        ax_v.text(0.0, 0.02, str(n111), ha='center', va='center', fontsize=15, fontweight='bold')

                        ax_v.text(-1.45, 1.55, m1, ha='center', va='center', fontsize=11, fontweight='bold')
                        ax_v.text(1.45, 1.55, m2, ha='center', va='center', fontsize=11, fontweight='bold')
                        ax_v.text(0.0, -1.95, m3, ha='center', va='center', fontsize=11, fontweight='bold')
                        ax_v.set_xlim(-2.1, 2.1)
                        ax_v.set_ylim(-1.9, 1.6)

                    ax_v.set_title(f'Comparison Venn (Top {top_n} Metabolites)', fontsize=14, fontweight='bold', pad=title_pad)
                    fig_v.tight_layout()
                    venn_name = f'comparison_venn_top{top_n}.png'
                    fig_v.savefig(os.path.join(output_root, venn_name), dpi=300)
                    for model_name in model_dirs:
                        fig_v.savefig(os.path.join(model_dirs[model_name], venn_name), dpi=300)
                    plt.close(fig_v)

            if not venn_sheet_written and os.path.exists(venn_excel_path):
                try:
                    os.remove(venn_excel_path)
                except Exception:
                    pass

        return {
            'output_root': output_root,
            'model_dirs': model_dirs,
            'models': list(model_runs.keys()),
            'combined_roc_path': combined_roc_path,
        }
    
    def _get_feature_importance(self) -> Optional[Dict[str, Any]]:
        """
        Extract feature importance if available from the model.
        
        Returns:
            Dictionary with importance scores and top features, or None
        """
        return self._extract_feature_importance_from_model(self.model)

    def _extract_feature_importance_from_model(self, model, X_test=None, y_test=None) -> Optional[Dict[str, Any]]:
        """Extract feature importance from any fitted model instance.
        
        Args:
            model: Fitted model instance
            X_test: Optional test features for permutation importance
            y_test: Optional test labels for permutation importance
        """
        if model is None:
            return None

        # Unwrap pipeline/calibrated wrappers to underlying estimator.
        if hasattr(model, 'named_steps') and 'model' in model.named_steps:
            model = model.named_steps['model']
        if isinstance(model, CalibratedClassifierCV):
            if hasattr(model, 'estimator'):
                model = model.estimator
            elif hasattr(model, 'base_estimator'):
                model = model.base_estimator

        feature_names = self.feature_names if hasattr(self, 'feature_names') else None

        # Check if model has feature_importances_ attribute
        if hasattr(model, 'feature_importances_'):
            importances = np.asarray(model.feature_importances_, dtype=float)
            if feature_names is None:
                feature_names = [f"Feature_{i}" for i in range(len(importances))]

            indices = np.argsort(importances)[::-1]
            top_n = min(20, len(importances))
            top_indices = indices[:top_n]
            top_features = [(feature_names[i], float(importances[i])) for i in top_indices]

            return {
                'importances': importances,
                'sorted_indices': indices,
                'top_features': top_features,
                'feature_names': feature_names,
                'note': 'Tree-based feature importance'
            }

        # For models with coef_ (like Logistic Regression, LDA)
        if hasattr(model, 'coef_'):
            coefs = np.abs(np.asarray(model.coef_, dtype=float)).mean(axis=0)
            if feature_names is None:
                feature_names = [f"Feature_{i}" for i in range(len(coefs))]

            indices = np.argsort(coefs)[::-1]
            top_n = min(20, len(coefs))
            top_indices = indices[:top_n]
            top_features = [(feature_names[i], float(coefs[i])) for i in top_indices]

            return {
                'importances': coefs,
                'sorted_indices': indices,
                'top_features': top_features,
                'feature_names': feature_names,
                'note': 'Coefficients magnitude used as importance'
            }

        # For models without native feature importance (SVM, etc.), use permutation importance
        if X_test is not None and y_test is not None:
            try:
                perm_importance = permutation_importance(
                    model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=-1
                )
                importances = perm_importance.importances_mean
                if feature_names is None:
                    feature_names = [f"Feature_{i}" for i in range(len(importances))]

                indices = np.argsort(importances)[::-1]
                top_n = min(20, len(importances))
                top_indices = indices[:top_n]
                top_features = [(feature_names[i], float(importances[i])) for i in top_indices]

                return {
                    'importances': importances,
                    'sorted_indices': indices,
                    'top_features': top_features,
                    'feature_names': feature_names,
                    'note': 'Permutation importance (model-agnostic)'
                }
            except Exception as e:
                logger.warning(f"Could not compute permutation importance: {e}")

        return None
    
    def export_results_to_excel(self, output_path: str, results_type: str = 'classification'):
        """
        Export analysis results to Excel file.
        
        Args:
            output_path: Path to save Excel file
            results_type: Type of results ('classification', 'pca', 'lda')
        """
        logger.info(f"Exporting {results_type} results to {output_path}...")
        
        def _fmt(value: Any) -> str:
            try:
                if value is None:
                    return "N/A"
                return f"{float(value):.4f}"
            except Exception:
                return str(value)

        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            wrote_any_sheet = False

            if results_type == 'classification' and self.results:
                # Multi-model comparison export
                if self.results.get('comparison_type') == 'multi_model' or (
                    'comparison_table' in self.results and 'models_trained' in self.results
                ):
                    summary_data = {
                        'Metric': [
                            'Analysis Type', 'Models Trained', 'Training Samples', 'Test Samples',
                            'Number of Features', 'Class Weight', 'Test Size', 'Scaling Method'
                        ],
                        'Value': [
                            'Multi-Model Comparison',
                            ", ".join(self.results.get('models_trained', [])),
                            self.results.get('n_train_samples', 'N/A'),
                            self.results.get('n_test_samples', 'N/A'),
                            self.results.get('n_features', 'N/A'),
                            self.results.get('class_weight') or 'none',
                            self.results.get('test_size', 'N/A'),
                            self.results.get('scaling_method', 'N/A')
                        ]
                    }
                    pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
                    wrote_any_sheet = True

                    comparison_df = pd.DataFrame(self.results.get('comparison_table', []))
                    if not comparison_df.empty:
                        comparison_df.to_excel(writer, sheet_name='Model_Comparison', index=False)

                    # Per-model detailed metrics
                    detail_rows = []
                    for model_name in self.results.get('models_trained', []):
                        model_result = self.results.get('model_results', {}).get(model_name, {})
                        if 'error' in model_result:
                            detail_rows.append({
                                'Model': model_name,
                                'Status': 'Error',
                                'Error': model_result.get('error', 'Unknown error')
                            })
                            continue

                        detail_rows.append({
                            'Model': model_name,
                            'Status': 'OK',
                            'Test Accuracy': _fmt(model_result.get('test_accuracy')),
                            'Train Accuracy': _fmt(model_result.get('train_accuracy')),
                            'AUC': _fmt(model_result.get('auc')),
                            'Sensitivity': str(model_result.get('sensitivity', {})),
                            'Specificity': str(model_result.get('specificity', {})),
                        })

                    if detail_rows:
                        pd.DataFrame(detail_rows).to_excel(writer, sheet_name='Model_Details', index=False)

                    # Aggregated feature importance / top metabolites
                    feature_rows = []
                    for model_name in self.results.get('models_trained', []):
                        model_result = self.results.get('model_results', {}).get(model_name, {})
                        feat_imp = model_result.get('feature_importances')
                        if not feat_imp or not feat_imp.get('top_features'):
                            continue

                        std_map = {feat: std for feat, std in feat_imp.get('top_features_std', [])}
                        for rank, (feature, importance) in enumerate(feat_imp['top_features'], 1):
                            feature_rows.append({
                                'Model': model_name,
                                'Rank': rank,
                                'Feature': feature,
                                'Mean Importance': float(importance),
                                'Std Importance': float(std_map.get(feature, 0.0)),
                            })

                    if feature_rows:
                        pd.DataFrame(feature_rows).to_excel(writer, sheet_name='Feature_Importance', index=False)

                # Single-model classification export
                else:
                    summary_data = {
                        'Metric': [
                            'Model', 'Repeated Runs',
                            'CV Mean Accuracy', 'CV Std Accuracy',
                            'CV Mean Balanced Accuracy', 'CV Std Balanced Accuracy',
                            'CV Mean AUC', 'CV Std AUC',
                            'Training Accuracy', 'Test Accuracy', 'Test Accuracy Std',
                            'Test Balanced Accuracy', 'Test Balanced Accuracy Std',
                            'Test Macro F1', 'Test Weighted F1',
                            'Test Macro Precision', 'Test Macro Recall',
                            'Test AUC', 'Test AUC Std',
                            'Number of Features', 'Training Samples', 'Test Samples'
                        ],
                        'Value': [
                            self.results.get('model_name', 'N/A'),
                            self.results.get('repeated_runs', 1),
                            _fmt(self.results.get('cv_mean_accuracy')),
                            _fmt(self.results.get('cv_std_accuracy')),
                            _fmt(self.results.get('cv_mean_balanced_accuracy')),
                            _fmt(self.results.get('cv_std_balanced_accuracy')),
                            _fmt(self.results.get('cv_mean_auc')),
                            _fmt(self.results.get('cv_std_auc')),
                            _fmt(self.results.get('train_accuracy')),
                            (_fmt(self.results.get('test_accuracy')) if self.results.get('test_accuracy') is not None else "N/A (CV-only)"),
                            _fmt(self.results.get('test_accuracy_std')),
                            _fmt(self.results.get('test_balanced_accuracy')),
                            _fmt(self.results.get('test_balanced_accuracy_std')),
                            _fmt(self.results.get('test_f1_macro')),
                            _fmt(self.results.get('test_f1_weighted')),
                            _fmt(self.results.get('test_precision_macro')),
                            _fmt(self.results.get('test_recall_macro')),
                            _fmt(self.results.get('auc')),
                            _fmt(self.results.get('auc_std')),
                            self.results.get('n_features', 'N/A'),
                            self.results.get('n_train_samples', 'N/A'),
                            self.results.get('n_test_samples', 'N/A')
                        ]
                    }
                    pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
                    wrote_any_sheet = True

                    # Confusion matrix
                    if self.results.get('confusion_matrix') is not None and self.results.get('class_labels') is not None:
                        conf_matrix_df = pd.DataFrame(
                            self.results['confusion_matrix'],
                            index=self.results['class_labels'],
                            columns=self.results['class_labels']
                        )
                        conf_matrix_df.to_excel(writer, sheet_name='Confusion_Matrix')

                    # Classification report
                    class_report = self.results.get('classification_report')
                    if class_report and isinstance(class_report, dict) and class_report:
                        report_df = pd.DataFrame(class_report).T
                        report_df.to_excel(writer, sheet_name='Classification_Report')

                    # Feature importance
                    feat_imp = self.results.get('feature_importances')
                    if feat_imp and isinstance(feat_imp, dict) and feat_imp.get('top_features'):
                        importance_data = feat_imp['top_features']
                        importance_df = pd.DataFrame(importance_data, columns=['Feature', 'Importance'])
                        importance_df.to_excel(writer, sheet_name='Feature_Importance', index=False)

            # Safety net: never write an empty workbook
            if not wrote_any_sheet:
                pd.DataFrame([
                    {
                        'Status': 'No exportable results',
                        'Message': 'No valid analysis results were found for export.'
                    }
                ]).to_excel(writer, sheet_name='Summary', index=False)
        
        logger.info(f"Results exported successfully to {output_path}")


def format_classification_summary(results: Dict[str, Any]) -> str:
    """
    Format classification results as readable text summary.
    
    Args:
        results: Results dictionary from run_classification
        
    Returns:
        Formatted string summary
    """
    def _fmt_float(value: Any) -> str:
        if value is None:
            return "N/A"
        try:
            if isinstance(value, float) and np.isnan(value):
                return "N/A"
        except Exception:
            pass
        try:
            return f"{float(value):.4f}"
        except Exception:
            return str(value)

    test_acc = results.get('test_accuracy')
    test_acc_text = _fmt_float(test_acc)
    if test_acc is None and results.get('cv_only_mode'):
        test_acc_text = "N/A (CV-only)"

    repeated_runs = int(results.get('repeated_runs', 1) or 1)
    test_acc_std = results.get('test_accuracy_std')
    auc_std = results.get('auc_std')

    if test_acc is not None and test_acc_std is not None:
        test_acc_text = f"{_fmt_float(test_acc)} ± {_fmt_float(test_acc_std)}"

    auc_text = _fmt_float(results.get('auc'))
    if results.get('auc') is not None and auc_std is not None:
        auc_text = f"{_fmt_float(results.get('auc'))} ± {_fmt_float(auc_std)}"

    # If these are multi-model comparison results, delegate to the dedicated formatter
    if 'model_name' not in results and 'comparison_table' in results and 'models_trained' in results:
        return format_multi_model_comparison(results)

    summary = f"""
╔══════════════════════════════════════════════════════════════╗
║         MACHINE LEARNING CLASSIFICATION RESULTS              ║
╚══════════════════════════════════════════════════════════════╝

Model: {results.get('model_name', 'N/A')}
Class Weight: {results.get('class_weight') or 'none'}
Repeated Runs: {repeated_runs}

📊 Performance Metrics:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    CV Accuracy (Mean ± Std):         {_fmt_float(results.get('cv_mean_accuracy'))} ± {_fmt_float(results.get('cv_std_accuracy'))}
    CV Balanced Accuracy (Mean ± Std): {_fmt_float(results.get('cv_mean_balanced_accuracy'))} ± {_fmt_float(results.get('cv_std_balanced_accuracy'))}
    CV AUC (Mean ± Std):              {_fmt_float(results.get('cv_mean_auc'))} ± {_fmt_float(results.get('cv_std_auc'))}
    Training Accuracy:                {_fmt_float(results.get('train_accuracy'))}
    Test Accuracy:                    {test_acc_text}
    Test Balanced Accuracy:           {_fmt_float(results.get('test_balanced_accuracy'))}
    Test Macro F1:                    {_fmt_float(results.get('test_f1_macro'))}
    Test Weighted F1:                 {_fmt_float(results.get('test_f1_weighted'))}
    Test Macro Precision:             {_fmt_float(results.get('test_precision_macro'))}
    Test Macro Recall:                {_fmt_float(results.get('test_recall_macro'))}
    Test AUC:                         {auc_text}

📈 Dataset Information:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Number of Features:              {results['n_features']}
  Training Samples:                {results['n_train_samples']}
  Test Samples:                    {results['n_test_samples']}
  Classes:                         {', '.join(results['class_labels'])}

🎯 Confusion Matrix:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    
    # Add confusion matrix
    conf_matrix = results.get('confusion_matrix')
    class_labels_raw = results.get('class_labels')
    # Handle numpy arrays and None safely
    if class_labels_raw is None:
        class_labels = []
    elif hasattr(class_labels_raw, 'tolist'):
        class_labels = class_labels_raw.tolist()
    else:
        class_labels = list(class_labels_raw) if class_labels_raw is not None else []

    if conf_matrix is None or len(class_labels) == 0:
        summary += "  (Not available in CV-only mode)\n"
    else:
        # Header
        summary += "         Predicted →\n"
        summary += "Actual ↓ " + "  ".join([f"{label:>8}" for label in class_labels]) + "\n"
        
        # Rows
        for i, label in enumerate(class_labels):
            row = f"{label:>8} " + "  ".join([f"{conf_matrix[i][j]:>8}" for j in range(len(class_labels))])
            summary += row + "\n"
    
    summary += "\n📋 Detailed Classification Report:\n"
    summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    summary += results.get('classification_report_str') or "(Not available)"

    if results.get('overfitting_warning'):
        summary += "\n\n⚠️ Overfitting Warning:\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        summary += f"  Mean train-test gap: {_fmt_float(results.get('overfitting_gap_mean'))}\n"
        summary += "  Consider stronger regularization, more samples, or stricter feature selection.\n"

    # Optional method diagnostics.
    hp_tuning = results.get('hyperparameter_tuning', {})
    if hp_tuning.get('enabled'):
        summary += "\n\n⚙️ Hyperparameter Tuning:\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        summary += f"  Strategy: {hp_tuning.get('strategy', 'grid')}\n"
        summary += f"  Iterations: {hp_tuning.get('iterations', 'N/A')}\n"
        best_params = results.get('best_params')
        if best_params:
            summary += f"  Best Params (selected run): {best_params}\n"

    if results.get('nested_cv_enabled'):
        summary += "\n🔁 Nested CV: Enabled\n"

    perm = results.get('permutation_test')
    if isinstance(perm, dict):
        summary += "\n🧪 Permutation Test:\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        if perm.get('error'):
            summary += f"  Error: {perm.get('error')}\n"
        else:
            summary += f"  Metric: {perm.get('metric', 'N/A')}\n"
            summary += f"  Score: {_fmt_float(perm.get('score'))}\n"
            summary += f"  Mean Permuted Score: {_fmt_float(perm.get('mean_permuted_score'))}\n"
            summary += f"  p-value: {_fmt_float(perm.get('pvalue'))}\n"
    
    # Add feature importance if available
    feature_importances = results.get('feature_importances')
    if feature_importances is not None and isinstance(feature_importances, dict) and feature_importances.get('top_features'):
        summary += "\n\n⭐ Top Important Features (up to 20):\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        top_features = feature_importances['top_features'][:20]
        for rank, (feature, importance) in enumerate(top_features, 1):
            feature_str = str(feature)  # Convert to string to handle int indices
            summary += f"  {rank:2d}. {feature_str:50s} {importance:.6f}\n"

    stable_selected = results.get('stable_selected_features') or []
    if stable_selected:
        summary += "\n\n🧬 Stable Selected Features Across Runs:\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        for rank, (feature, freq) in enumerate(stable_selected[:20], 1):
            summary += f"  {rank:2d}. {str(feature):50s} {freq:6.1f}%\n"

    stable_important = results.get('stable_important_features') or []
    if stable_important:
        summary += "\n\n⭐ Stable Important Features Across Runs:\n"
        summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        for rank, (feature, freq) in enumerate(stable_important[:20], 1):
            summary += f"  {rank:2d}. {str(feature):50s} {freq:6.1f}%\n"
    
    return summary


def format_pca_summary(results: Dict[str, Any]) -> str:
    """
    Format PCA results as readable text summary.
    
    Args:
        results: Results dictionary from run_pca
        
    Returns:
        Formatted string summary
    """
    n_components = min(10, results['n_components'])
    
    summary = f"""
╔══════════════════════════════════════════════════════════════╗
║      PRINCIPAL COMPONENT ANALYSIS (PCA) RESULTS              ║
╚══════════════════════════════════════════════════════════════╝

📊 Variance Explained:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    
    for i in range(n_components):
        variance = results['explained_variance'][i]
        cumulative = results['cumulative_variance'][i]
        bar_length = int(variance * 50)
        bar = "█" * bar_length + "░" * (50 - bar_length)
        summary += f"  PC{i+1:2d}:  {bar}  {variance:6.2%} (Cumulative: {cumulative:6.2%})\n"
    
    summary += f"\n✅ First 2 PCs explain {results['cumulative_variance'][1]:.2%} of total variance\n"
    summary += f"✅ First 3 PCs explain {results['cumulative_variance'][2]:.2%} of total variance\n" if results['n_components'] >= 3 else ""
    
    return summary


def format_multi_model_comparison(results: Dict[str, Any]) -> str:
    """
    Format multi-model comparison results as readable text summary.
    
    Args:
        results: Results dictionary from run_multi_model_comparison
        
    Returns:
        Formatted string summary
    """
    def _format_metric_dict(metric_dict: Dict[str, Any]) -> str:
        """Format metric dictionary into readable class=value pairs."""
        if not metric_dict:
            return "N/A"

        parts = []
        for k, v in metric_dict.items():
            key = str(k)
            try:
                value = float(v)
                parts.append(f"{key}={value:.4f}")
            except Exception:
                parts.append(f"{key}={v}")

        return ", ".join(parts)

    repeated_runs = int(results.get('repeated_runs', 1) or 1)
    random_state = results.get('random_state', 42)

    summary = f"""
╔══════════════════════════════════════════════════════════════╗
║      MULTI-MODEL COMPARISON RESULTS                          ║
║      (Fair Evaluation - Same Data Split)                     ║
╚══════════════════════════════════════════════════════════════╝

📊 Experiment Setup:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    Repeated Runs:         {repeated_runs}
    Base Seed:            {random_state}
  Training Samples:       {results['n_train_samples']}
  Test Samples:          {results['n_test_samples']}
  Features:             {results['n_features']}
  Classes:              {', '.join(results['class_labels'])}
  Test Size:            {results['test_size']:.1%}
  Scaling Method:       {results['scaling_method']}
  Class Weight:         {'balanced' if results.get('class_weight') else 'none'}

🏆 Model Comparison:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    
    # Add comparison table header
    summary += f"{'Model':<30} {'Accuracy (mean±std)':<22} {'AUC (mean±std)':<22}\n"
    summary += "─" * 54 + "\n"
    
    # Add rows from comparison_table
    for row in results['comparison_table']:
        model_name = row['Model'][:28]
        accuracy = row['Accuracy']
        auc = row['AUC']
        summary += f"{model_name:<30} {accuracy:<22} {auc:<22}\n"
    
    summary += "\n"
    
    # Add detailed metrics per model
    summary += "📋 Detailed Metrics (Test Set):\n"
    summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
    
    for model_name in results['models_trained']:
        model_result = results['model_results'].get(model_name, {})
        
        if 'error' in model_result:
            summary += f"❌ {model_name}: {model_result['error']}\n\n"
            continue
        
        summary += f"🎯 {model_name}\n"
        summary += f"   ├─ Accuracy:    {model_result.get('test_accuracy_mean', model_result.get('test_accuracy', 0)):.4f}"
        if model_result.get('test_accuracy_std') is not None:
            summary += f" ± {model_result.get('test_accuracy_std', 0):.4f}\n"
        else:
            summary += "\n"
        
        if model_result.get('auc_mean') is not None:
            auc_std = model_result.get('auc_std')
            summary += f"   ├─ AUC:         {model_result['auc_mean']:.4f}"
            if auc_std is not None:
                summary += f" ± {auc_std:.4f}\n"
            else:
                summary += "\n"
        
        # Sensitivity and Specificity
        sens = model_result.get('sensitivity', {})
        spec = model_result.get('specificity', {})
        
        if sens:
            summary += f"   ├─ Sensitivity: {_format_metric_dict(sens)}\n"
        if spec:
            summary += f"   └─ Specificity: {_format_metric_dict(spec)}\n"
        
        summary += "\n"

    # Add feature ranking summary across models
    summary += "⭐ Top Ranked Metabolites (Mean Importance Across Runs):\n"
    summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"

    any_features = False
    for model_name in results.get('models_trained', []):
        model_result = results.get('model_results', {}).get(model_name, {})
        feat_imp = model_result.get('feature_importances')
        if not feat_imp or not feat_imp.get('top_features'):
            continue

        any_features = True
        summary += f"🎯 {model_name}\n"
        summary += f"   Rank  Feature                                           Mean     Std\n"
        summary += f"   ─────────────────────────────────────────────────────────────────────\n"
        std_map = {feat: std for feat, std in feat_imp.get('top_features_std', [])}
        for rank, (feature, importance) in enumerate(feat_imp['top_features'][:10], 1):
            feature_str = str(feature)[:45]
            summary += f"   {rank:>2}.  {feature_str:<45} {float(importance):>7.4f}  {float(std_map.get(feature, 0.0)):>7.4f}\n"
        summary += "\n"

    if not any_features:
        summary += "  (No feature importance available for the selected models)\n\n"
    
    summary += "💡 Interpretation Guide:\n"
    summary += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    summary += "  • Accuracy: Overall correctness (may be misleading with imbalanced classes)\n"
    summary += "  • AUC: Area Under ROC Curve (0.5 = random, 1.0 = perfect)\n"
    summary += "  • Sensitivity: True Positive Rate (catches positive cases)\n"
    summary += "  • Specificity: True Negative Rate (avoids false alarms)\n"
    
    return summary
