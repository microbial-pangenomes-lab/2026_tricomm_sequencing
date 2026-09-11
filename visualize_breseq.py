#!/usr/bin/env python3
"""
Visualize breseq GD output files.

This script parses a breseq GD file and creates a linear genome visualization
with mutations marked as pairs of dots. Different mutation types are colored
differently based on their likely biological impact.
"""

import sys
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Define user-friendly replicon names and their actual lengths
REPLICON_NAMES = {
    'NZ_CP009273': {'name': 'Chromosome', 'length': 4631469},  # E. coli BW25113 chromosome
    'pS2313M_msfGFP': {'name': 'msfGFP plasmid', 'length': 3861},    # From GenBank file
    'pS2313RS_mScarle': {'name': 'mScarlet plasmid', 'length': 3846}, # From GenBank file
    'pS2313RS_mScarlet': {'name': 'mScarlet plasmid', 'length': 3846}, # Handle both spellings
    'NZ_MT441554': {'name': 'pOXA48 plasmid', 'length': 65499},     # From GenBank file
    'MG557852': {'name': 'PN23 plasmid', 'length': 33858}            # From GenBank file
}

# Define regions to highlight on the chromosome (from regions.py)
HIGHLIGHT_REGIONS = {
    'lysA': (2_970_996, 2_972_258),
    'pheA': (2_731_104, 2_732_264),
    'metA': (4_204_208, 4_205_137)
}

def parse_gd_file(gd_file_path):
    """
    Parse a breseq GD file and extract mutation information.

    Args:
        gd_file_path (str): Path to the GD file

    Returns:
        dict: Dictionary containing mutation data organized by replicon
    """
    mutations_by_replicon = defaultdict(list)

    # Define mutation type categories and their colors based on snippet.py
    # This maps mutation categories to the new color scheme
    mutation_type_to_index = {
        'snp_nonsynonymous': 0,
        'snp_nonsense': 1,
        'small_indel': 2,
        'snp_noncoding': 3,
        'snp_intergenic': 4,
        'snp_pseudogene': 4,
        'snp_synonymous': 5,
        'large_insertion': 6,
        'large_deletion': 6,
        'large_substitution': 6,
        'mobile_element_insertion': 7,
    }

    # Corresponding colors from snippet.py
    colors = [
        'xkcd:pale red',      # 0: snp_nonsynonymous
        'xkcd:bright red',    # 1: snp_nonsense
        'xkcd:rust',          # 2: small_indel
        'xkcd:aqua green',    # 3: snp_noncoding
        'xkcd:grey',          # 4: snp_intergenic/snp_pseudogene
        'xkcd:steel',         # 5: snp_synonymous
        'xkcd:neon purple',   # 6: large_deletion/large_insertion/large_substitution
        'xkcd:violet',        # 7: mobile_element_insertion
        'xkcd:violet red',    # 8: mobile_element_insertion_JC
    ]

    # Create colormap
    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    # Default color for unknown categories
    default_color = 'xkcd:light grey'

    with open(gd_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            parts = line.split('\t')
            if len(parts) < 6:
                continue

            mutation_type = parts[0]

            # Check if this is a special 2-letter code (MC, JC) that should be included in highlighted regions
            is_special_mutation = mutation_type in ['MC', 'JC']

            # Only consider mutations with 3-letter codes that represent confirmed mutations
            # Exclude RA (read anomalies) and UN (uncovered regions) as they are not confirmed variants
            # But allow MC and JC mutations in highlighted regions
            if len(mutation_type) != 3 and not is_special_mutation:
                continue

            # Parse the line to extract key information
            mutation_data = {}
            mutation_data['type'] = mutation_type
            mutation_data['id'] = parts[1]
            mutation_data['evidence'] = parts[2]
            mutation_data['replicon'] = parts[3]

            # Extract position information
            position_info = {}
            for field in parts[4:]:
                if '=' in field:
                    key, value = field.split('=', 1)
                    position_info[key] = value

            # Get position start and end
            position_start = None
            position_end = None

            if 'position_start' in position_info:
                position_start = int(position_info['position_start'])
            elif len(parts) >= 5:
                # For some mutation types, position is directly in column 5
                try:
                    position_start = int(parts[4])
                    position_end = position_start
                except ValueError:
                    # Handle ranges like "279,294"
                    if ',' in parts[4]:
                        pos_parts = parts[4].split(',')
                        if len(pos_parts) >= 2:
                            position_start = int(pos_parts[0])
                            position_end = int(pos_parts[1])
                    else:
                        position_start = int(parts[4])
                        position_end = position_start

            # For MC/JC mutations, position_end is in column 6
            if is_special_mutation and len(parts) >= 6:
                try:
                    position_end = int(parts[5])
                except ValueError:
                    pass  # Keep the position_end from above

            if 'position_end' in position_info:
                position_end = int(position_info['position_end'])
            elif position_end is None:
                position_end = position_start

            mutation_data['position_start'] = position_start
            mutation_data['position_end'] = position_end

            # Get mutation category
            mutation_category = 'unknown'
            if 'mutation_category' in position_info:
                mutation_category = position_info['mutation_category']
            elif mutation_type in ['RA', 'UN', 'MOB', 'INS', 'AMP', 'CON']:
                # Skip these types as they are not confirmed variants
                continue
            elif is_special_mutation:
                # MC and JC mutations get special handling - treat like large deletions
                mutation_category = 'large_deletion'
            else:
                # For mutation types without explicit category, use the type itself
                mutation_category = mutation_type.lower()

            mutation_data['category'] = mutation_category
            
            # Store gene_name if available
            if 'gene_name' in position_info:
                mutation_data['gene_name'] = position_info['gene_name']

            # For MC/JC mutations, check if they are in highlighted regions OR on pOXA48 plasmid
            if is_special_mutation:
                mutation_in_region = False

                # Include MC mutations on pOXA48 plasmid (whole length), but NOT JC
                if mutation_data['replicon'] == 'NZ_MT441554' and mutation_type == 'MC':  # pOXA48 plasmid, MC only
                    mutation_in_region = True
                # For chromosome, check if mutation overlaps any highlighted region
                elif mutation_data['replicon'] == 'NZ_CP009273':  # Chromosome
                    for region_start, region_end in HIGHLIGHT_REGIONS.values():
                        # Check if mutation overlaps with this region
                        if not (position_end < region_start or position_start > region_end):
                            mutation_in_region = True
                            break

                # Only include MC/JC mutations if they meet the criteria
                if not mutation_in_region:
                    continue

            # Assign color based on category using the new scheme
            if mutation_category in mutation_type_to_index:
                color_index = mutation_type_to_index[mutation_category]
                mutation_data['color'] = colors[color_index]
                mutation_data['category_index'] = color_index
            else:
                # For categories not in the main mapping, use judgment
                if mutation_category.startswith('snp_'):
                    # Unknown SNP type - use intergenic color
                    mutation_data['color'] = colors[4]  # xkcd:grey
                    mutation_data['category_index'] = 4
                elif mutation_category in ['ra', 'un']:
                    # Read anomalies and uncovered regions - use a neutral color
                    mutation_data['color'] = 'xkcd:light blue'
                    mutation_data['category_index'] = 99  # Custom index
                elif mutation_category in ['mob', 'ins', 'amp', 'con']:
                    # Mobile elements and other insertions - use mobile element color
                    mutation_data['color'] = colors[7]  # xkcd:violet
                    mutation_data['category_index'] = 7
                else:
                    # Unknown category - use default grey
                    mutation_data['color'] = default_color
                    mutation_data['category_index'] = 999

            mutations_by_replicon[mutation_data['replicon']].append(mutation_data)

    return mutations_by_replicon

def get_replicon_info(mutations_by_replicon):
    """
    Extract replicon information and sort them in the specified order.
    """
    # Define the desired order of replicons
    replicon_order = [
        'NZ_CP009273',  # Chromosome (BW25113)
        'NZ_MT441554',  # pOXA48 plasmid
        'MG557852',     # PN23 plasmid
        'pS2313M_msfGFP',  # Plasmid
        'pS2313RS_mScarlet',  # Plasmid
    ]

    # Get all replicons and their lengths
    replicons = []
    for replicon_id in mutations_by_replicon.keys():
        if replicon_id not in replicons:
            replicons.append(replicon_id)

    # Sort replicons according to the desired order, then alphabetically
    sorted_replicons = []
    for desired_replicon in replicon_order:
        if desired_replicon in replicons:
            sorted_replicons.append(desired_replicon)
            replicons.remove(desired_replicon)

    # Add any remaining replicons in alphabetical order
    sorted_replicons.extend(sorted(replicons))

    # Calculate replicon lengths
    replicon_lengths = {}
    for replicon_id in sorted_replicons:
        max_position = 0
        for mutation in mutations_by_replicon[replicon_id]:
            if mutation['position_end'] > max_position:
                max_position = mutation['position_end']
        replicon_lengths[replicon_id] = max_position

    return sorted_replicons, replicon_lengths

def create_visualization(mutations_by_replicon, output_prefix):
    """
    Create a visualization of mutations on the genome with separate panels for each replicon.
    """
    # Redefine colors for use in this function (since parse_gd_file scope is different)
    colors = [
        'xkcd:pale red',      # 0: snp_nonsynonymous
        'xkcd:bright red',    # 1: snp_nonsense
        'xkcd:rust',          # 2: small_indel
        'xkcd:aqua green',    # 3: snp_noncoding
        'xkcd:grey',          # 4: snp_intergenic/snp_pseudogene
        'xkcd:steel',         # 5: snp_synonymous
        'xkcd:neon purple',   # 6: large_deletion/large_insertion/large_substitution
        'xkcd:violet',        # 7: mobile_element_insertion
    ]

    # Get replicon information
    sorted_replicons, replicon_lengths = get_replicon_info(mutations_by_replicon)

    # Parse vial information from output_prefix for title
    vial_title = "Unknown"
    if "ancestrals" in output_prefix:
        # Format: "Ancestral clone, LM, replicate 1" or "Ancestral clone, LM, PVB donor"
        sample = output_prefix.split('_', 1)[1] if '_' in output_prefix else ''
        strain = sample[:2].upper()
        suffix = sample[2:]
        if strain in ('LM', 'PL', 'PM'):
            if suffix.isdigit():
                vial_title = f"Ancestral clone, {strain}, replicate {suffix}"
            elif suffix.startswith('_'):
                vial_title = f"Ancestral clone, {strain}, {suffix[1:]} donor"
            else:
                vial_title = f"Ancestral clone, {strain}"
    elif not output_prefix.startswith("old_run") and not output_prefix.startswith("test"):
        # Format: "pOXA48, rep1, LM clone1, TriComm KAN no drug" or "PN23, rep1, 309h, ..."
        parts = output_prefix.split('_')

        # Determine strain name from the first part of the prefix
        strain_name = parts[0]  # e.g. 'pOXA48' or 'PN23'

        rep_info = ""
        strain_info = ""
        time_info = ""
        vial_code = ""

        for i, part in enumerate(parts):
            if part.startswith('M') and part[1:].isdigit():
                vial_code = part
            elif part == 'rep1' or part == 'rep2':
                rep_info = part.replace('rep', 'rep')
            elif part == 'clone1' or part == 'clone2':
                time_info = "clone" + part[-1]
            elif part in ['LM', 'PL', 'PM']:
                strain_info = part
            elif re.fullmatch(r'\d+h', part):
                time_info = part

        # Get vial description - use different maps for pOXA48 and PN23
        poxa48_vial_descriptions = {
            'M0': 'TriComm KAN no drug',
            'M1': 'TriComm KAN 16x MIC',
            'M2': 'TriComm KAN+PVI no drug',
            'M3': 'TriComm KAN+PVI 16x MIC',
            'M4': 'TriComm PVI no drug',
            'M5': 'TriComm PVI 16x MIC',
            'M6': 'TriComm PVB no drug',
            'M7': 'TriComm PVB 16x MIC',
        }
        # The two PN23 replicates were run with different vial layouts, so the
        # map has to be chosen per replicate (see collate_references.py and
        # data/PN23/parsing.ipynb in the 2026_tricomm repository).
        pn23_vial_descriptions = {
            'rep1': {
                'M1': 'TriComm KAN+PN23 2xMIC',
                'M2': 'TriComm KAN no drug',
                'M3': 'TriComm KAN 1xMIC',
                'M4': 'TriComm KAN 2xMIC',
                'M5': 'TriComm PN23 no drug',
                'M6': 'TriComm PN23 1xMIC',
                'M7': 'TriComm PN23 2xMIC',
            },
            'rep2': {
                'M0': 'TriComm KAN no drug',
                'M1': 'TriComm PN23 no drug',
                'M2': 'TriComm KAN 1xMIC',
                'M3': 'TriComm PN23 1xMIC',
                'M4': 'TriComm KAN 2xMIC',
                'M5': 'TriComm PN23 2xMIC',
            },
        }

        if strain_name == 'PN23':
            if rep_info not in pn23_vial_descriptions:
                raise ValueError(
                    f"cannot label {output_prefix}: the PN23 vial layout differs "
                    f"between replicates and no replicate was found in the name")
            vial_desc = pn23_vial_descriptions[rep_info].get(vial_code, "Unknown condition")
        else:
            vial_desc = poxa48_vial_descriptions.get(vial_code, "Unknown condition")

        # Build title components
        title_components = [strain_name]
        if rep_info:
            title_components.append(rep_info)
        if strain_info:
            title_components.append(strain_info)
        if time_info:
            title_components.append(time_info)
        title_components.append(vial_desc)

        vial_title = ", ".join(title_components)

    # Determine the number of panels needed
    num_replicons = len(sorted_replicons)

    # Create a figure with subplots - one panel per replicon
    # Use smaller height per panel (2 inches instead of 4)
    fig, axes = plt.subplots(num_replicons, 1, figsize=(15, 2 * num_replicons))

    # If there's only one replicon, axes will be a single Axes object, not an array
    if num_replicons == 1:
        axes = [axes]

    # Find the maximum genome length for consistent x-axis scaling
    max_genome_length = max(replicon_lengths.values()) if replicon_lengths else 1000000

    # Plot each replicon in its own panel
    for i, (replicon_id, ax) in enumerate(zip(sorted_replicons, axes)):
        replicon_length = replicon_lengths[replicon_id]

        # Use user-friendly replicon name and actual length
        replicon_info = REPLICON_NAMES.get(replicon_id, {'name': replicon_id, 'length': replicon_length})
        friendly_name = replicon_info['name']
        actual_length = replicon_info.get('length', replicon_length)

        # Set panel title (without "Replicon:" prefix)
        ax.set_title(friendly_name, fontsize=12, fontweight='bold')

        # Add highlighted regions for chromosome (before drawing mutations)
        if replicon_id == 'NZ_CP009273':  # Chromosome
            for gene_name, (region_start, region_end) in HIGHLIGHT_REGIONS.items():
                # Highlight region with very light grey
                ax.fill_between([region_start, region_end],
                               [-0.5, -0.5], [0.5, 0.5],
                               color='xkcd:light grey', alpha=0.3)

                # Add gene name label in the middle of the region
                region_mid = (region_start + region_end) / 2
                ax.text(region_mid, 0.3, gene_name,
                       fontsize=8, ha='center', va='center',
                       bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

        # Plot mutations
        for mutation in mutations_by_replicon[replicon_id]:
            start_pos = mutation['position_start']
            end_pos = mutation['position_end']
            color = mutation['color']

            # Draw dots at start and end positions
            ax.plot([start_pos, end_pos], [0, 0],
                   'o-', markersize=6, linewidth=1,
                   color=color, alpha=0.8)

            # Add gene name annotation for 3-letter mutations if gene_name is available
            if 'gene_name' in mutation and mutation['gene_name'] and mutation['gene_name'] != '[]':
                gene_name = mutation['gene_name']
                # Clean up gene name (remove brackets if present)
                if gene_name.startswith('[') and gene_name.endswith(']'):
                    gene_name = gene_name[1:-1]
                # Add gene name below the mutation (larger font)
                ax.text((start_pos + end_pos) / 2, -0.3, gene_name, 
                       fontsize=8, ha='center', va='top', alpha=0.7)

        # Set x-axis limits using the actual replicon length with padding
        # Use smaller padding for plasmids, larger padding for chromosome
        if replicon_id == 'NZ_CP009273':  # Chromosome
            x_padding = max(5000, actual_length * 0.01)  # At least 5kb or 1%
        else:  # Plasmids - use much smaller padding
            x_padding = min(100, actual_length * 0.01)  # At most 100bp or 1%
        ax.set_xlim(-x_padding, actual_length + x_padding)
        ax.set_ylim(-0.5, 0.5)

        # Customize the panel
        ax.set_xlabel('Position', fontsize=10)
        ax.grid(True, alpha=0.3)

        # Remove y-ticks since we're only showing mutations
        ax.set_yticks([])
        ax.set_yticklabels([])

        # Remove y-axis label to save space
        ax.set_ylabel('')

        # Add some padding
        ax.margins(y=0.2)
        plt.tight_layout()

    # Add overall title with vial information (without prefix)
    fig.suptitle(vial_title, fontsize=16, y=1.02)

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    png_path = f'{output_prefix}.png'
    svg_path = f'{output_prefix}.svg'

    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, format='svg', bbox_inches='tight')

    print(f'Visualization saved to: {png_path}')
    print(f'Visualization saved to: {svg_path}')

    plt.close(fig)

def create_legend_file(output_path='mutation_legend.svg'):
    """
    Create a separate legend file for mutation types.
    """
    # Create legend elements for the color scheme
    legend_elements = []
    mutation_type_names = [
        'snp_nonsynonymous', 'snp_nonsense', 'small_indel', 'snp_noncoding',
        'snp_intergenic', 'snp_synonymous', 'large_deletion',
        'mobile_element_insertion'
    ]

    colors = [
        'xkcd:pale red',      # 0: snp_nonsynonymous
        'xkcd:bright red',    # 1: snp_nonsense
        'xkcd:rust',          # 2: small_indel
        'xkcd:aqua green',    # 3: snp_noncoding
        'xkcd:grey',          # 4: snp_intergenic/snp_pseudogene
        'xkcd:steel',         # 5: snp_synonymous
        'xkcd:neon purple',   # 6: large_deletion/large_insertion/large_substitution
        'xkcd:violet',        # 7: mobile_element_insertion
    ]

    for i, mut_type in enumerate(mutation_type_names):
        if i < len(colors):
            legend_elements.append(patches.Patch(facecolor=colors[i], edgecolor='black',
                                                label=mut_type.replace('_', ' ').title()))

    # Create a figure just for the legend
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.axis('off')  # Hide axes

    # Add legend
    fig.legend(handles=legend_elements, loc='center',
               title='Mutation Type Legend', fontsize=12, ncol=4)

    # Save the legend
    fig.savefig(output_path, format='svg', bbox_inches='tight', dpi=300)
    plt.close(fig)

    print(f"Legend saved to: {output_path}")

def main():
    """
    Main function to run the visualization.
    """
    if len(sys.argv) != 3:
        print("Usage: python visualize_breseq.py <input_gd_file> <output_prefix>")
        print("Example: python visualize_breseq.py out/sample/data/annotated.gd sample_mutations")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    print(f"Parsing GD file: {input_file}")
    mutations_by_replicon = parse_gd_file(input_file)

    if not mutations_by_replicon:
        print("No mutations found in the GD file.")
        sys.exit(0)

    print(f"Found mutations in {len(mutations_by_replicon)} replicons:")
    for replicon_id, mutations in mutations_by_replicon.items():
        print(f"  {replicon_id}: {len(mutations)} mutations")

    print(f"Creating visualization...")
    create_visualization(mutations_by_replicon, output_prefix)

    # Create annotated version with gene names
    print(f"Creating annotated visualization...")
    create_visualization(mutations_by_replicon, output_prefix + "_annotated")

    print("Done!")

if __name__ == "__main__":
    # Check if we should create a legend file
    if len(sys.argv) == 2 and sys.argv[1] == '--create-legend':
        create_legend_file()
        sys.exit(0)
    else:
        main()
