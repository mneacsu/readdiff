#!/usr/bin/env python3
import os
import time
import anndata
import math

def main():

    sparse = snakemake.input.sparse
    subsets_dir = snakemake.output.subsets_dir
    n_subsets = snakemake.params.n_subsets

    # Create output directory
    os.makedirs(subsets_dir, exist_ok=True)

    # Remove existing files in the directory
    for f in os.listdir(subsets_dir):
        file_path = os.path.join(subsets_dir, f)
        if os.path.isfile(file_path):
            os.remove(file_path)

    start_time = time.time()

    print(f"Reading input file: {sparse}")
    adata = anndata.read_h5ad(sparse)
    print(f"Input shape: {adata.shape}")

    # Get feature names (in shuffled order)
    all_features = adata.var_names.tolist()
    total_features = len(all_features)

    nb_features = math.ceil(total_features/n_subsets)

    print(f"Total features: {total_features}")
    print(f"Features per subset: {nb_features}")

    for subset in range(n_subsets):

        # Calculate batch boundaries
        start_idx = subset * nb_features
        end_idx = start_idx + nb_features


        # Select features for this batch (may get fewer if at end)
        selected_features = all_features[start_idx:min(end_idx, total_features)]

        print(f"Selecting features [{start_idx}:{min(end_idx, total_features)}] -> {len(selected_features)} features")

        # Subset the AnnData
        adata_subset = adata[:, selected_features].copy()

        # Write output
        adata_subset.write_h5ad(f"{subsets_dir}/subset_{subset}.h5ad")

    elapsed_time = time.time() - start_time
    print(f"Execution time: {elapsed_time:.2f} seconds")



if __name__ == "__main__":
    main()