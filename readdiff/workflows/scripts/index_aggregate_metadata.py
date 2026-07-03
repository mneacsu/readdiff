#!/usr/bin/env python3

import pandas as pd

def aggregate_csv_files(input_files, output_file):

    all_dfs = []
    for i, csv_file in enumerate(input_files):
        df = pd.read_csv(csv_file)
        df['index_id'] = i  # Ajouter l'ID de l'index
        all_dfs.append(df)

    # Concaténer tous les DataFrames
    combined_df = pd.concat(all_dfs, ignore_index=True).rename(columns = {"sample": "sample_id"})

    # Trier par sample et index_id
    combined_df = combined_df.sort_values(['sample_id', 'index_id'])

    # Sauvegarder le résultat
    combined_df.to_csv(output_file, index=False)

def main():

    window_counts = list(snakemake.input.window_counts)
    read_counts = list(snakemake.input.read_counts)
    metadata_files = list(snakemake.input.metadata)
    global_window_count = snakemake.output.global_window_count
    global_read_count = snakemake.output.global_read_count
    global_metadata = snakemake.output.global_metadata
    index_metadata = snakemake.output.index_metadata
    n_indices = snakemake.params.n_indices
    paired_end = snakemake.params.paired_end
    compress = snakemake.params.compress

    # Aggregate window counts
    aggregate_csv_files(window_counts, global_window_count)

    # Aggregate read counts
    aggregate_csv_files(read_counts, global_read_count)

    # Aggregate metadata
    aggregate_csv_files(metadata_files, global_metadata)

    # Create index metadata file
    with open(index_metadata, 'w') as f:
        f.write("n_indices,n_samples,paired_end,compressed\n")
        
        # Count samples from global_window_count
        with open(global_window_count, 'r') as gcf:
            n_samples = len(gcf.readlines()) - 1  # Subtract header
        
        f.write(f"{n_indices},{n_samples},{paired_end},{compress}\n")

if __name__ == "__main__":
    main()