#!/usr/bin/env python3
"""
Collate the 'references' section from all summary.json files in the out directory
into a pandas DataFrame and save to file.
"""

import json
import os
import re
from pathlib import Path
import pandas as pd

EXPERIMENT_MAP = {
    'M0': 'TriComm KAN no drug',
    'M1': 'TriComm KAN 16x MIC', 
    'M2': 'TriComm KAN+PVI no drug',
    'M3': 'TriComm KAN+PVI 16x MIC',
    'M4': 'TriComm PVI no drug',
    'M5': 'TriComm PVI 16x MIC',
    'M6': 'TriComm PVB no drug',
    'M7': 'TriComm PVB 16x MIC',
}


def find_summary_files(out_dir: str) -> list[str]:
    """Find all summary.json files in the out directory, excluding old_run."""
    summary_files = []
    for root, dirs, files in os.walk(out_dir):
        # Skip old_run directory
        if "old_run" in root:
            continue
        if "summary.json" in files:
            summary_files.append(os.path.join(root, "summary.json"))
    return sorted(summary_files)


def extract_references_data(file_path: str) -> list[dict]:
    """Extract references data from a summary.json file."""
    with open(file_path, "r") as f:
        data = json.load(f)
    
    references = data.get("references", {})
    reference_dict = references.get("reference", {})
    
    # Extract path info for context
    # Path structure: out/{replicate}/{sample}/{data_type}/summary.json
    # e.g., out/pOXA48_rep1/M7_LM_clone1/output/summary.json
    # e.g., out/ancestrals/LM1/data/summary.json
    path_parts = Path(file_path).parts
    
    # Find the index of 'out' in the path to properly extract replicate and sample
    try:
        out_idx = path_parts.index("out")
        replicate = path_parts[out_idx + 1] if len(path_parts) > out_idx + 1 else ""
        sample_name = path_parts[out_idx + 2] if len(path_parts) > out_idx + 2 else ""
        data_type = path_parts[out_idx + 3] if len(path_parts) > out_idx + 3 else ""
    except ValueError:
        replicate = path_parts[-4] if len(path_parts) >= 4 else ""
        sample_name = path_parts[-3] if len(path_parts) >= 3 else ""
        data_type = path_parts[-2] if len(path_parts) >= 2 else ""
    
    results = []
    for ref_name, ref_data in reference_dict.items():
        row = {
            "file_path": file_path,
            "replicate": replicate,
            "sample_name": sample_name,
            "data_type": data_type,
            "reference_name": ref_name,
        }
        # Add all reference data fields
        row.update(ref_data)
        results.append(row)
    
    return results


def extract_experiment_name(sample_name: str) -> str | None:
    """Extract experiment name from sample name using vial number (M0, M1, etc.).
    
    Examples:
        M7_LM_clone1 -> TriComm PVB 16x MIC
        M5_198h -> TriComm PVI 16x MIC
        LM1 -> None (ancestral, not a pOXA48 sample)
    """
    # Match M0, M1, M2, etc. at the start of the sample name
    match = re.match(r'^(M\d+)_', sample_name)
    if match:
        vial_number = match.group(1)
        return EXPERIMENT_MAP.get(vial_number)
    return None


def main():
    # Define paths
    script_dir = Path(__file__).parent
    out_dir = script_dir / "out"
    output_file = script_dir / "references_collated.parquet"
    output_csv = script_dir / "references_collated.csv"
    
    print(f"Searching for summary.json files in: {out_dir}")
    
    # Find all summary.json files
    summary_files = find_summary_files(str(out_dir))
    print(f"Found {len(summary_files)} summary.json files")
    
    # Extract data from all files
    all_data = []
    for file_path in summary_files:
        try:
            data = extract_references_data(file_path)
            all_data.extend(data)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    print(f"Extracted {len(all_data)} reference entries")
    
    # Create DataFrame
    df = pd.DataFrame(all_data)
    
    # Add experiment name column for pOXA48 samples
    df['experiment_name'] = df['sample_name'].apply(extract_experiment_name)
    
    # Display summary
    print(f"\nDataFrame shape: {df.shape}")
    print(f"\nColumns: {list(df.columns)}")
    print(f"\nUnique references: {df['reference_name'].nunique()}")
    print(f"\nReference names: {sorted(df['reference_name'].unique())}")
    print(f"\nUnique sample names: {df['sample_name'].nunique()}")
    print(f"\nSample names: {sorted(df['sample_name'].unique())}")
    print(f"\nUnique experiment names: {df['experiment_name'].nunique()}")
    print(f"\nExperiment names: {sorted(df['experiment_name'].dropna().unique())}")
    
    # Also save to CSV for easy viewing
    df.to_csv(output_csv, index=False)
    print(f"Saved to: {output_csv}")
    
    # Show sample data
    print("\n--- Sample entries (first 5 rows) ---")
    print(df.head().to_string())
    
    return df


if __name__ == "__main__":
    main()
