import pandas as pd
from plotnine import *
import gzip
import numpy as np
import json

def parse_compute_matrix(matrix_file):
    """
    Parses the output of deeptools computeMatrix.

    Args:
        matrix_file (str): Path to the computeMatrix output file (gzipped).

    Returns:
        pandas.DataFrame: A tidy DataFrame with columns:
                            'group', 'feature', 'bin', 'value'.
    """
    with gzip.open(matrix_file, 'rt') as f:
        header = f.readline().strip()
        assert header.startswith('@'), "Invalid computeMatrix file format"
        header = json.loads(header[1:])

        # Extract information from the header
        sample_labels = header['sample_labels']
        group_labels = header['group_labels']
        
        # Read the rest of the file into a DataFrame
        df = pd.read_csv(f, sep='\t', header=None)
        
    # Drop the metadata columns to get the matrix
    matrix = df.drop(columns=[0, 1, 2, 3, 4, 5]).to_numpy()

    # Create a tidy DataFrame
    n_samples = len(sample_labels)
    n_bins_per_sample = matrix.shape[1] // n_samples
    
    tidy_data = []
    for feature_idx, row in enumerate(matrix):
        for sample_idx, sample in enumerate(sample_labels):
            for bin_idx in range(n_bins_per_sample):
                value = row[sample_idx * n_bins_per_sample + bin_idx]
                group_idx = 0
                for i in range(len(header['group_boundaries']) - 1):
                    if header['group_boundaries'][i] <= feature_idx < header['group_boundaries'][i+1]:
                        group_idx = i
                        break

                tidy_data.append({
                    'group': group_labels[group_idx],
                    'feature': df.iloc[feature_idx, 3],
                    'sample': sample,
                    'bin': bin_idx,
                    'value': value
                })

    return pd.DataFrame(tidy_data), header

def plot_profile(df, header, color='sample'):
    """
    Generates a profile plot from a tidy DataFrame.

    Args:
        df (pandas.DataFrame): Tidy DataFrame from parse_compute_matrix.
        color (str): Column to group by for color.

    Returns:
        plotnine.ggplot: A ggplot object.
    """
    # Convert bins to base pairs
    df['bp'] = df['bin'] * header['bin size'][0]

    summary_df = df.groupby(['bp', color])['value'].mean().reset_index()
    
    body_start = header['upstream'][0]
    body_end = body_start + header['body'][0]

    p = (
        ggplot(summary_df, aes(x='bp', y='value', color=color)) +
        geom_line() +
        geom_vline(xintercept=[body_start, body_end], linetype='dashed', color='black', alpha=0.4) +
        labs(x="Genomic Region (bp)", y="Mean Signal") +
        theme_minimal()
    )
    return p

def plot_heatmap(df, header, facet_by='sample', region_labels=('TSS', 'TES')):
    """
    Generates a heatmap from a tidy DataFrame, ordering features by mean signal.

    Args:
        df (pandas.DataFrame): Tidy DataFrame from parse_compute_matrix.
        header (dict): Header dictionary from parse_compute_matrix.
        facet_by (str): Column to facet by.

    Returns:
        plotnine.ggplot: A ggplot object.
    """
    # Convert bins to base pairs
    df['bp'] = df['bin'] * header['bin size'][0]

    body_start = int(header['upstream'][0])
    body_end = int(body_start + header['body'][0])
    
    # Calculate the mean value for each feature within the body
    feature_mean = (
        df[(df['bp'] >= body_start) & (df['bp'] <= body_end)]
        .groupby('feature')['value']
        .mean()
        .sort_values(ascending=False)
        .index
    )

    # Order the features in the DataFrame
    df['feature'] = pd.Categorical(df['feature'], categories=reversed(feature_mean), ordered=True)
    df['bp'] = pd.Categorical(df['bp'])
    start = list(df['bp'].cat.categories).index(body_start) + 1
    end = list(df['bp'].cat.categories).index(body_end) + 1
    
    p = (
        ggplot(df, aes(x='bp', y='feature', fill='value')) +
        geom_tile() +
        facet_wrap(f'~{facet_by}') +
        geom_vline(xintercept=[start, end], linetype='dashed', color='white') +
        labs(x="Genomic Region (bp)", y="Feature") +
        scale_x_discrete(
            breaks=(body_start, body_end),
            labels=region_labels
        ) +
        theme_minimal() +
        theme(
            axis_text_y=element_blank(),
            #axis_text_x=element_blank(),
            axis_ticks_major_y=element_blank()
        )
    )
    return p

if __name__ == '__main__':
    # This is an example of how to use the functions.
    # You would need a computeMatrix output file to run this.
    #
    # Example usage:
    #
    # df, header = parse_compute_matrix('path/to/your/matrix.gz')
    #
    # # Generate profile plot
    # profile_plot = plot_profile(df, header)
    # profile_plot.save("profile_plot.png")
    #
    # # Generate heatmap
    # heatmap_plot = plot_heatmap(df, header, facet_by='sample')
    # heatmap_plot.save("heatmap_plot.png")
    pass