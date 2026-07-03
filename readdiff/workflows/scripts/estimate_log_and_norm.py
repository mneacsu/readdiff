#!/usr/bin/env python3

import anndata
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import shutil

def main():

    raw = snakemake.input.raw
    counts = snakemake.params.counts
    normalize = snakemake.params.normalize
    log10 = snakemake.params.log10
    transformed = snakemake.output.transformed


    if normalize and (counts is None):
        raise ValueError("Window counts are required for normalization. Please provide them.")

    # Read the .h5ad file
    print(f"Reading input AnnData from {raw}")
    adata = anndata.read_h5ad(raw)

    # Ensure the matrix is in CSR format for efficient row-based slicing
    # (Only necessary if we want fast row-based operations on a sparse matrix.)
    if not sp.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()

    if normalize:
        # Read the CSV file containing counts
        print(f"Reading window counts from {counts}")
        df_window_counts = pd.read_csv(counts, index_col=0)

        if "window_count" in df_window_counts.columns:
            scaling_factor = 1e11
            df_window_counts = df_window_counts.rename(columns={"window_count": "count"})
        elif "read_count" in df_window_counts.columns:
            scaling_factor = 1e9
            df_window_counts = df_window_counts.rename(columns={"read_count": "count"})
        else:
            raise ValueError(f"File {counts} does not contain window_count / read_count column.")

        # Check if all samples in the .h5ad file are present in the counts file
        missing = set(adata.obs_names) - set(df_window_counts.index)
        if missing:
            raise ValueError(
                f"Sample(s) {missing} were not found in {counts}."
            )
        
        print("Performing window-count normalization...")
        counts = df_window_counts.loc[adata.obs_names, "count"]
        if any(counts == 0):
            raise ValueError(f"Window counts for sample(s) '{[adata.obs_names[i] for i, v in enumerate(counts) if v == 0]}' are 0.")
        factors = sp.diags(list(scaling_factor/counts))
        adata.X = factors @ adata.X
        print("Normalization complete.")

    # Apply logtransformation if specified
    if log10:
        adata.X.data = np.log10(adata.X.data + 1)
        print("Logtransformed the values using log10(1+x)")

    # Round estimates to 2 decimal places
    adata.X.data = adata.X.data.round(2)

    # Write the resulting AnnData to the output file
    adata.write_h5ad(transformed)

    # Copy to estimates/estimates_transformed.h5ad
    estimates_dir = os.path.dirname(transformed)
    shutil.copy(transformed, os.path.join(estimates_dir, f"estimates_transformed.h5ad"))

    
if __name__ == "__main__":
    main()
