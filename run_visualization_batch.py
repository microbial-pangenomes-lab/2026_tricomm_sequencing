#!/usr/bin/env python3
"""
Batch script to run breseq visualization on multiple samples.

This script finds all annotated.gd files in the out/ directory and creates
visualizations for each sample.
"""

import os
import subprocess
import sys

def find_gd_files(base_dir='out'):
    """
    Find all annotated.gd files in the output directory structure.
    """
    gd_files = []
    
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'annotated.gd':
                full_path = os.path.join(root, file)
                gd_files.append(full_path)
    
    return gd_files

def create_output_prefix(gd_file_path):
    """
    Create a nice output prefix from the GD file path.
    """
    # Remove the 'out/' prefix and 'data/annotated.gd' suffix
    relative_path = gd_file_path.replace('out/', '')
    sample_path = relative_path.replace('/data/annotated.gd', '')
    
    # Replace slashes with underscores for filename safety
    output_prefix = sample_path.replace('/', '_')
    
    return output_prefix

def main():
    """
    Main function to run batch visualization.
    """
    print("Finding breseq GD files...")
    gd_files = find_gd_files()

    if not gd_files:
        print("No annotated.gd files found in the out/ directory.")
        return

    print(f"Found {len(gd_files)} samples to process.")

    # Place all output files in a dedicated plots directory, creating it if
    # it does not already exist. The visualization subprocess runs with this
    # directory as its working directory so the output_prefix (which is also
    # used to build plot titles) stays clean.
    plots_dir = os.path.abspath('plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Absolute paths so they resolve regardless of the subprocess cwd.
    visualize_script = os.path.abspath('visualize_breseq.py')

    for i, gd_file in enumerate(gd_files, 1):
        print(f"\nProcessing sample {i}/{len(gd_files)}: {gd_file}")

        # Create output prefix
        output_prefix = create_output_prefix(gd_file)

        # Run the visualization script
        try:
            cmd = [sys.executable, visualize_script, os.path.abspath(gd_file), output_prefix]
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=plots_dir)

            if result.returncode == 0:
                print(f"  ✓ Success: plots/{output_prefix}.png and plots/{output_prefix}.svg created")
            else:
                print(f"  ✗ Error processing {gd_file}")
                print(f"  Error: {result.stderr}")

        except Exception as e:
            print(f"  ✗ Exception processing {gd_file}: {str(e)}")

    # Create a single legend file for the whole batch
    print("\nCreating legend file...")
    try:
        cmd = [sys.executable, visualize_script, '--create-legend']
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=plots_dir)
        if result.returncode == 0:
            print("  ✓ Legend file created: plots/mutation_legend.svg")
        else:
            print(f"  ✗ Error creating legend: {result.stderr}")
    except Exception as e:
        print(f"  ✗ Exception creating legend: {str(e)}")

    print("\nBatch processing complete!")

if __name__ == "__main__":
    main()
