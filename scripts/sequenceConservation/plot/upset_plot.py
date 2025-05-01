#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import plot as upsetplot
import numpy as np
import time
import argparse
from collections import Counter

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate UpSet plot for Origin groups')
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('--output', default='origin_upset_plot.png', help='Output file path')
    parser.add_argument('--min-subset-size', type=int, default=0,
                        help='Minimum size for a subset to be included in the plot')
    parser.add_argument('--sample', type=int, default=None,
                        help='Number of rows to sample (for testing on smaller dataset)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output figure')
    parser.add_argument('--figsize', type=str, default='12,10', help='Figure size (width,height)')
    parser.add_argument('--filter-lncrna', action='store_true', help='Filter for lncRNA transcripts only')
    return parser.parse_args()

def process_origin_data(file_path, sample_size=None, filter_lncrna=False):
    """
    Process the Origin column into a format suitable for UpSet plot

    Args:
        file_path: Path to the TSV file
        sample_size: Number of rows to sample (optional)
        filter_lncrna: Whether to filter for lncRNA transcripts only

    Returns:
        pandas.DataFrame: Data formatted for UpSet plot
    """
    print(f"Reading data from {file_path}...")
    start_time = time.time()

    # Define columns to read
    usecols = ['Origin']
    if filter_lncrna:
        usecols.append('Transcript Function')
        print("Filtering for lncRNA transcripts only")

    # Use chunksize for efficient processing of large files
    chunks = pd.read_csv(file_path, sep='\t', usecols=usecols, chunksize=100000)

    # Count occurrences of each Origin value
    origin_counter = Counter()

    for i, chunk in enumerate(chunks):
        if i % 10 == 0:
            print(f"Processing chunk {i}...")

        # Filter for lncRNA if requested
        if filter_lncrna:
            chunk = chunk[chunk['Transcript Function'] == 'lncRNA']
            if chunk.empty:
                continue

        # Update counter with origins in this chunk
        origin_counter.update(chunk['Origin'].value_counts().to_dict())

        # If sampling, break after reaching sample size
        if sample_size and sum(origin_counter.values()) >= sample_size:
            # Adjust the counter to match the sample size
            excess = sum(origin_counter.values()) - sample_size
            if excess > 0:
                # Remove excess counts proportionally
                for key in list(origin_counter.keys()):
                    reduction = int(excess * origin_counter[key] / sum(origin_counter.values()))
                    origin_counter[key] = max(0, origin_counter[key] - reduction)
                    excess -= reduction
                    if excess <= 0:
                        break
            break

    elapsed = time.time() - start_time
    print(f"Data processed in {elapsed:.2f} seconds")
    print(f"Found {len(origin_counter)} unique Origin values with {sum(origin_counter.values())} total entries")

    # Remove the header row if present
    if 'Origin' in origin_counter:
        print("Removing header row...")
        del origin_counter['Origin']

    # Create proper format for upsetplot
    # For upsetplot, we need a pandas Series with MultiIndex
    data_for_upset = {}
    
    # Map origin values to set memberships
    for origin, count in origin_counter.items():
        if origin == 'temp':
            # Skip UNK as it doesn't belong to any set (mudei pra temp só pra pular essa parte)
            continue
        elif origin == 'SBAR':
            data_for_upset[(True, False, False, False)] = count  # S_barberi only
        elif origin == 'SOFF':
            data_for_upset[(False, True, False, False)] = count  # S_officinarum only
        elif origin == 'SSPO':
            data_for_upset[(False, False, True, False)] = count  # S_spontaneum only
        elif origin == 'CommonSOFF_SBAR':
            data_for_upset[(True, True, False, False)] = count  # S_barberi and S_officinarum
        elif origin == 'CommonSSPO_SBAR':
            data_for_upset[(True, False, True, False)] = count  # S_barberi and S_spontaneum
        elif origin == 'CommonSSPO_SOFF':
            data_for_upset[(False, True, True, False)] = count  # S_officinarum and S_spontaneum
        elif origin == 'Common':
            data_for_upset[(True, True, True, False)] = count  # All three
        elif origin == 'UNK':
            data_for_upset[(False, False, False, True)] = count # Unknown only

    # Convert to Series with proper index
    plot_data = pd.Series(data_for_upset)
    # Convert the tuple keys to MultiIndex
    plot_data.index = pd.MultiIndex.from_tuples(
        plot_data.index, names=['S. barberi', 'S. officinarum', 'S. spontaneum', 'Unknown']
    )
    
    return plot_data

def create_upset_plot(plot_data, output_file='origin_upset_plot.png', min_subset_size=0,
                     figsize=(12, 10), dpi=300, filter_lncrna=False):
    """
    Create UpSet plot visualization

    Args:
        plot_data: Series with MultiIndex for upsetplot
        output_file: Path to save the output figure
        min_subset_size: Minimum size for a subset to be included in the plot
        figsize: Figure size as (width, height) tuple
        dpi: DPI for the output figure
        filter_lncrna: Whether the data is filtered for lncRNA (for title)
    """
    print("Creating UpSet plot...")

    # Filter by minimum subset size if specified
    if min_subset_size > 0:
        plot_data = plot_data[plot_data >= min_subset_size]
        print(f"Filtered to {len(plot_data)} sets with min_subset_size={min_subset_size}")

    # Create figure
    plt.figure(figsize=figsize)

    # Generate the upset plot
    upsetplot(plot_data, sort_by='cardinality',  show_counts=True, show_percentages=True, element_size=None)
    #upset = upsetplot(plot_data, sort_by='cardinality', show_counts=False, show_percentages=False)
    
    # Remove grid from all subplots
    for ax in plt.gcf().get_axes():
        ax.grid(False)

    # Calculate the total number of unique transcripts
    total_transcripts = sum(plot_data)
 
    # Add title
    title = 'UpSet Plot of Origin Groups'
    if filter_lncrna:
        title += f'({total_transcripts} lncRNAs)'
    plt.suptitle(title, fontsize=18)
    plt.tight_layout()

    # Save figure
    print(f"Saving plot to {output_file}...")
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Plot saved successfully!")

def main():
    args = parse_arguments()

    # Parse figsize from string
    width, height = map(float, args.figsize.split(','))
    figsize = (width, height)

    # Process data
    plot_data = process_origin_data(
        args.input_file, 
        sample_size=args.sample,
        filter_lncrna=args.filter_lncrna
    )

    # Create and save plot
    create_upset_plot(
        plot_data,
        output_file=args.output,
        min_subset_size=args.min_subset_size,
        figsize=figsize,
        dpi=args.dpi,
        filter_lncrna=args.filter_lncrna
    )

    print("Done!")

if __name__ == "__main__":
    main()
