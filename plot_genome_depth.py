#!/usr/bin/env python
# mamba activate python3 # pandas numpy matplotlib
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
import numpy as np
import sys

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a chromosome-style depth plot from mosdepth region data.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input regions file from mosdepth (gzipped or plain text)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output image file name. Format is determined by extension (e.g., plot.svg, plot.png)."
    )
    parser.add_argument(
        "--breaks",
        nargs=2,
        type=float,
        default=[75.0, 100.0],
        metavar=('LOW', 'HIGH'),
        help="Two depth breakpoints for coloring. Default: 75 100."
    )
    parser.add_argument(
        "--highlight",
        type=str,
        help="Optional region to highlight. Format: 'ChromosomeName:start-end'\nExample: 'hap1Chr25:0-1000000'"
    )
    parser.add_argument(
        "--highlight_legend",
        type=str,
        default="Highlight",
        help="Legend label for the highlighted region. Default: 'Highlight'."
    )
    return parser.parse_args()

def parse_chrom_name(name):
    """
    Parses a chromosome name to find its grouping factor ('ChrXX')
    and a numerical key for sorting within that group.
    """
    # Find the main chromosome group (e.g., 'Chr25'), case-insensitive
    chrom_group_match = re.search(r'(Chr\d+)', name, re.IGNORECASE)
    if not chrom_group_match:
        return None, None, name # Cannot determine group, return name as label

    chrom_group = chrom_group_match.group(1)

    # Find a sorting key by looking for numbers outside the main group part
    remaining_part = name.replace(chrom_group, '')
    sort_key_match = re.search(r'(\d+)', remaining_part)
    
    # If a number is found, use it for sorting. Otherwise, default to 0.
    sort_key = int(sort_key_match.group(1)) if sort_key_match else 0
    
    return chrom_group, sort_key, name

def main():
    """Main function to generate the plot."""
    args = parse_arguments()
    
    # --- Data Loading and Preparation ---
    try:
        # Automatically handle gzipped or plain text files
        compression = 'gzip' if args.input.endswith('.gz') else 'infer'
        df = pd.read_csv(args.input, sep='\t', header=None, compression=compression)
        df.columns = ['chromosome', 'start', 'end', 'depth']
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading or parsing input file: {e}", file=sys.stderr)
        sys.exit(1)


    # Apply chromosome name parsing
    parsed_names = df['chromosome'].apply(parse_chrom_name)
    df[['chrom_group', 'sort_key', 'label']] = pd.DataFrame(parsed_names.tolist(), index=df.index)
    df = df.dropna(subset=['chrom_group', 'sort_key'])

    # Define color mapping based on command-line breaks
    break1, break2 = args.breaks
    def get_color(depth):
        if depth < break1:
            return 'green'
        elif break1 <= depth <= break2:
            return 'yellow'
        else:
            return 'red'
    df['color'] = df['depth'].apply(get_color)

    # Get sorted unique chromosome groups
    chrom_groups = sorted(df['chrom_group'].unique(), key=lambda x: int(re.search(r'\d+', x).group()))

    # --- Plotting ---
    n_groups = len(chrom_groups)
    if n_groups == 0:
        print("Error: No valid chromosome groups found. Check chromosome name format.", file=sys.stderr)
        sys.exit(1)
        
    n_cols = 3
    n_rows = int(np.ceil(n_groups / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 2.5), squeeze=False)

    # Plotting with top-to-bottom, then left-to-right order
    for i, chrom_group_name in enumerate(chrom_groups):
        col = i // n_rows
        row = i % n_rows
        ax = axes[row, col]

        group_df = df[df['chrom_group'] == chrom_group_name].copy()
        
        # Sort contigs within the group by the parsed sort_key
        contigs = sorted(group_df['label'].unique(), key=lambda c: group_df[group_df['label'] == c]['sort_key'].iloc[0])
        contig_y_pos = {contig: j for j, contig in enumerate(contigs)}
        
        for _, r in group_df.iterrows():
            y_pos = contig_y_pos[r['label']]
            ax.add_patch(patches.Rectangle((r['start'], y_pos - 0.4), r['end'] - r['start'], 0.8, facecolor=r['color']))

        # Highlight region if specified
        if args.highlight:
            try:
                hl_chrom, hl_range = args.highlight.split(':')
                hl_start, hl_end = map(int, hl_range.split('-'))
                if hl_chrom in contig_y_pos:
                    y_highlight = contig_y_pos[hl_chrom]
                    ax.add_patch(patches.Rectangle((hl_start, y_highlight - 0.4), hl_end - hl_start, 0.8, 
                                                   facecolor='none', edgecolor='blue', linewidth=2))
            except (ValueError, IndexError):
                print(f"Warning: Could not parse highlight string '{args.highlight}'. Skipping.", file=sys.stderr)

        # --- Axis Formatting ---
        ax.set_yticks(list(contig_y_pos.values()))
        ax.set_yticklabels(contigs)
        ax.set_title(chrom_group_name)
        ax.tick_params(axis='x', rotation=45)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f'{x/1e6:.0f}Mb'))
        ax.grid(True, axis='x', linestyle='--', alpha=0.6)
        ax.set_ylim(-0.5, len(contigs) - 0.5)
        ax.set_xlim(0, group_df['end'].max())
        ax.invert_yaxis()

    # --- Legend and Final Touches ---
    for i in range(len(axes.flatten())):
        col = i // n_rows
        row = i % n_rows
        if i >= n_groups:
            axes[row, col].set_visible(False)

    legend_patches = [
        patches.Patch(color='green', label=f'< {break1}x'),
        patches.Patch(color='yellow', label=f'{break1}x-{break2}x'),
        patches.Patch(color='red', label=f'> {break2}x'),
    ]
    if args.highlight:
        legend_patches.append(patches.Patch(facecolor='none', edgecolor='blue', linewidth=2, label=args.highlight_legend))

    # Place legend in the first available empty subplot
    if n_groups < n_rows * n_cols:
        legend_ax = axes.flatten()[n_groups]
        legend_ax.set_visible(True)
        legend_ax.axis('off')
        legend_ax.legend(handles=legend_patches, loc='center', ncol=1)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved as {args.output}")

if __name__ == '__main__':
    main()