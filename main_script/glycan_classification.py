"""
Glycan Classification Utility

Classifies glycans based on their structural features encoded in the feature ID.
Feature IDs are typically 4-5 digit codes (e.g., 4501 or 45010) representing:
  Position 1: HexNAc
  Position 2: Hex
  Position 3: Fucose (Fuc)
  Position 4: NeuAc (sialic acid N-acetyl)
  Position 5: NeuGc (sialic acid N-glycolyl, if present)

Classifications:
  1. Sialylated (NeuAc) - contains NeuAc (4th digit > 0) and no NeuGc
  2. Sialylated (NeuGc) - contains NeuGc (5th digit > 0) - EXCLUSIVE, no other classifications
  3. Fucosylated - contains Fuc (3rd digit > 0) and is NOT sialylated
  4. Sialofucosylated - contains both NeuAc and Fuc, but not NeuGc
  5. High mannose - starts with 25, 26, 27, 28, or 29 (can combine with other modifications)
  6. Other - doesn't fit any of the above
"""

import logging
from typing import List, Optional, Tuple
import pandas as pd
import re

logger = logging.getLogger(__name__)


def parse_glycan_id(feature_id: str) -> Optional[Tuple[int, int, int, int, int]]:
    """
    Parse a glycan feature ID into individual position values.
    
    Handles formats like:
    - "4501" → (4, 5, 0, 1, 0)
    - "4-5-0-1" → (4, 5, 0, 1, 0)
    - "45010" → (4, 5, 0, 1, 0)
    - "4-5-0-1-0" → (4, 5, 0, 1, 0)
    - "2(10)00" → (2, 10, 0, 0, 0)
    
    Returns:
        Tuple of (HexNAc, Hex, Fuc, NeuAc, NeuGc) or None if parsing fails
    """
    if not isinstance(feature_id, str):
        return None
    
    feature_id = feature_id.strip()
    
    try:
        values = []
        
        # Handle parenthetical notation like 2(10)00 → [2, 10, 0, 0]
        if '(' in feature_id and ')' in feature_id:
            # Replace parentheses with dashes and split
            parsed = re.sub(r'[()]', '-', feature_id)
            parts = parsed.split('-')
            values = [int(p) for p in parts if p]  # Remove empty strings and convert to int
        elif '-' in feature_id:
            # Dash-separated format: "4-5-0-1" or "4-5-0-1-0"
            parts = feature_id.split('-')
            values = [int(p) for p in parts if p]
        else:
            # Contiguous digits without separators
            # Assume: first 4 digits are individual positions, but need to handle multi-digit hex value
            # For now, try as single digits first
            clean_id = ''.join(c for c in feature_id if c.isdigit())
            
            if len(clean_id) == 4:
                # Standard 4-digit format: HHFN where each letter is a digit
                values = [int(d) for d in clean_id]
            elif len(clean_id) == 5:
                # 5-digit format: HHFNX
                values = [int(d) for d in clean_id]
            else:
                logger.warning(f"Unexpected glycan ID length: {feature_id} → {clean_id}")
                return None
        
        # Normalize to 5 positions: (HexNAc, Hex, Fuc, NeuAc, NeuGc)
        if len(values) == 3:
            # Add NeuAc and NeuGc (both 0)
            values.extend([0, 0])
        elif len(values) == 4:
            # Add missing 5th position (NeuGc)
            values.append(0)
        elif len(values) == 5:
            pass
        else:
            logger.warning(f"Unexpected glycan ID format: {feature_id} → {values}")
            return None
        
        return tuple(values[:5])
    except (ValueError, IndexError) as e:
        logger.warning(f"Failed to parse glycan ID '{feature_id}': {e}")
        return None


def classify_glycan(feature_id: str) -> Optional[str]:
    """
    Classify a single glycan based on its feature ID.
    
    Returns:
        Classification string or None if parsing fails
    """
    parsed = parse_glycan_id(feature_id)
    if parsed is None:
        return None
    
    hex_nac, hex_val, fuc, neu_ac, neu_gc = parsed
    
    # Get first two digits as string for high mannose check
    first_digit = hex_nac
    second_digit = hex_val
    first_two = int(f"{first_digit}{second_digit}")
    
    is_high_mannose = first_two in [25, 26, 27, 28, 29]
    has_neu_ac = neu_ac > 0
    has_neu_gc = neu_gc > 0
    has_fuc = fuc > 0
    
    # Classification logic with priority:
    # 1. NeuGc is EXCLUSIVE - if present, only return Sialylated (NeuGc)
    if has_neu_gc:
        return "Sialylated (NeuGc)"
    
    # 2. High mannose can have secondary modifications
    if is_high_mannose:
        if has_neu_ac and has_fuc:
            return "High mannose + Sialofucosylated"
        elif has_neu_ac:
            return "High mannose + Sialylated (NeuAc)"
        elif has_fuc:
            return "High mannose + Fucosylated"
        else:
            return "High mannose"
    
    # 3. Sialylated (NeuAc) without Fuc
    if has_neu_ac and not has_fuc:
        return "Sialylated (NeuAc)"
    
    # 4. Sialofucosylated (both NeuAc and Fuc)
    if has_neu_ac and has_fuc:
        return "Sialofucosylated"
    
    # 5. Fucosylated (only Fuc, no sialylation)
    if has_fuc:
        return "Fucosylated"
    
    # 6. Other
    return "Other"


def classify_glycans_batch(feature_ids: List[str]) -> List[Optional[str]]:
    """
    Classify multiple glycans at once.
    
    Args:
        feature_ids: List of glycan feature ID strings
        
    Returns:
        List of classifications (in same order as input)
    """
    return [classify_glycan(fid) for fid in feature_ids]


def process_glycan_dataframe(
    df: pd.DataFrame,
    feature_id_col: str,
    sample_cols: List[str],
    drop_na_classes: bool = True
) -> pd.DataFrame:
    """
    Process a DataFrame containing glycan data and group by classification.
    
    Steps:
    1. Add classification column
    2. Group by classification
    3. Sum all sample columns by classification group
    4. Reset index
    
    Args:
        df: Input DataFrame with feature IDs and sample columns
        feature_id_col: Name of column containing glycan feature IDs
        sample_cols: List of column names to sum (sample abundance columns)
        drop_na_classes: If True, drop rows where classification is None
        
    Returns:
        DataFrame with classifications and summed abundances
        
    Raises:
        ValueError: If required columns are missing or sample_cols is empty
    """
    if feature_id_col not in df.columns:
        raise ValueError(f"Feature ID column '{feature_id_col}' not found in DataFrame")
    
    missing_cols = [col for col in sample_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Sample columns not found: {missing_cols}")
    
    if not sample_cols:
        raise ValueError("No sample columns provided")
    
    # Create a working copy
    work_df = df.copy()
    
    # Add classification column
    work_df['Glycan_Class'] = work_df[feature_id_col].apply(classify_glycan)
    
    # Drop rows where classification is None if requested
    if drop_na_classes:
        work_df = work_df[work_df['Glycan_Class'].notna()].copy()
    
    # Group by classification and sum sample columns
    grouped = work_df.groupby('Glycan_Class', as_index=False)[sample_cols].sum()
    
    return grouped


def get_classification_summary(df: pd.DataFrame, feature_id_col: str) -> pd.DataFrame:
    """
    Get a summary of glycan classifications in a DataFrame.
    
    Returns:
        DataFrame with columns: Classification, Count, Percentage
    """
    classifications = df[feature_id_col].apply(classify_glycan)
    summary = classifications.value_counts().reset_index()
    summary.columns = ['Classification', 'Count']
    summary['Percentage'] = (summary['Count'] / summary['Count'].sum() * 100).round(2)
    
    return summary


# Test cases for validation
if __name__ == "__main__":
    print("Glycan Classification Test Cases\n" + "="*50)
    
    test_cases = [
        ("4501", "Sialylated (NeuAc)"),
        ("4511", "Sialofucosylated"),
        ("42100", "Fucosylated"),
        ("45001", "Sialylated (NeuGc)"),
        ("25000", "High mannose"),
        ("27001", "Sialylated (NeuGc)"),  # NeuGc exclusive
        ("28100", "High mannose + Fucosylated"),
        ("2300", "Other"),
        ("43000", "Other"),
        ("44000", "Other"),
        ("56111", "Sialylated (NeuGc)"),
        ("29000", "High mannose"),
        ("2(10)00", "Other"),  # 2-10-0-0
    ]
    
    print("Running classifications:")
    all_passed = True
    for feature_id, expected in test_cases:
        result = classify_glycan(feature_id)
        status = "✓" if result == expected else "✗"
        if result != expected:
            all_passed = False
        print(f"  {status} {feature_id:10} → {result:30} (expected: {expected})")
    
    print("\n" + "="*50)
    print(f"Result: {'All tests passed!' if all_passed else 'Some tests failed!'}")
