#!/usr/bin/env python3
import os
import shutil
import gc
import time
import anndata
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, issparse
import glob

def _write_sparse_incremental(output_path, h5_group_path, input_paths, n_obs,
                              total_vars, total_nnz, layer_key=None):
    """Write a sparse matrix incrementally to h5ad file using h5py.

    Builds a CSC matrix on disk by reading one input file at a time,
    converting its matrix to CSC, and appending data/indices/indptr.
    """
    # Detect dtypes from first input
    first = anndata.read_h5ad(input_paths[0])
    X = first.layers[layer_key] if layer_key else first.X
    X_csc = X.tocsc() if issparse(X) else csc_matrix(X)
    data_dtype = X_csc.data.dtype
    indices_dtype = X_csc.indices.dtype
    del first, X, X_csc
    gc.collect()

    with h5py.File(output_path, "r+") as f:
        if h5_group_path in f:
            del f[h5_group_path]

        grp = f.create_group(h5_group_path)
        grp.attrs["encoding-type"] = "csc_matrix"
        grp.attrs["encoding-version"] = "0.1.0"
        grp.attrs["shape"] = np.array([n_obs, total_vars], dtype=np.int64)

        data_ds = grp.create_dataset("data", shape=(total_nnz,), dtype=data_dtype)
        indices_ds = grp.create_dataset("indices", shape=(total_nnz,), dtype=indices_dtype)
        indptr_ds = grp.create_dataset("indptr", shape=(total_vars + 1,), dtype=np.int64)

        data_offset = 0
        var_offset = 0

        for i, path in enumerate(input_paths):
            print(f"    [{i+1}/{len(input_paths)}] {path}")
            adata = anndata.read_h5ad(path)

            X = adata.layers[layer_key] if layer_key else adata.X
            X_csc = X.tocsc() if issparse(X) else csc_matrix(X)
            nnz = X_csc.nnz
            n_vars_i = X_csc.shape[1]

            if nnz > 0:
                data_ds[data_offset:data_offset + nnz] = X_csc.data
                indices_ds[data_offset:data_offset + nnz] = X_csc.indices

            # Write indptr (overlapping boundary values are identical)
            indptr_ds[var_offset:var_offset + n_vars_i + 1] = X_csc.indptr.astype(np.int64) + data_offset

            data_offset += nnz
            var_offset += n_vars_i

            del adata, X, X_csc
            gc.collect()


def _write_dense_incremental(output_path, h5_group_path, input_paths, n_obs,
                             total_vars, layer_key=None):
    """Write a dense matrix incrementally to h5ad file using h5py."""
    # Detect dtype from first input
    first = anndata.read_h5ad(input_paths[0])
    X = first.layers[layer_key] if layer_key else first.X
    data_dtype = X.dtype if not issparse(X) else X.toarray().dtype
    del first, X
    gc.collect()

    with h5py.File(output_path, "r+") as f:
        if h5_group_path in f:
            del f[h5_group_path]

        ds = f.create_dataset(
            h5_group_path,
            shape=(n_obs, total_vars),
            dtype=data_dtype,
            chunks=(min(n_obs, 1000), min(total_vars, 10000)),
        )
        ds.attrs["encoding-type"] = "array"
        ds.attrs["encoding-version"] = "0.2.0"

        var_offset = 0
        for i, path in enumerate(input_paths):
            print(f"    [{i+1}/{len(input_paths)}] {path}")
            adata = anndata.read_h5ad(path)

            X = adata.layers[layer_key] if layer_key else adata.X
            data = X.toarray() if issparse(X) else np.asarray(X)
            n_vars_i = data.shape[1]

            ds[:, var_offset:var_offset + n_vars_i] = data
            var_offset += n_vars_i

            del adata, X, data
            gc.collect()

def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment, NaN-preserving."""
    p = np.asarray(pvalues, dtype=float)
    mask = ~np.isnan(p)
    n = int(mask.sum())
    if n == 0:
        return p.copy()
    pv = p[mask]
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(1, n + 1)
    adj = pv * n / ranks
    sorted_adj = adj[order]
    sorted_adj = np.minimum.accumulate(sorted_adj[::-1])[::-1]
    sorted_adj = np.minimum(sorted_adj, 1.0)
    out_valid = np.empty(n)
    out_valid[order] = sorted_adj
    out = p.copy()
    out[mask] = out_valid
    return out

def main():

    subset_counts = list(snakemake.input.subset_counts)
    query_dir = snakemake.params.query_dir
    n_subsets = snakemake.params.n_subsets
    counts = snakemake.output.counts
    results = snakemake.output.results

    start_time = time.time()

    print("Merging DESeq2 normalized counts...")

    # === First pass: collect metadata and compute sizes ===
    var_dfs = []
    obs = None
    n_obs = None
    total_vars = 0
    total_nnz_x = 0
    x_is_sparse = True

    for i, path in enumerate(subset_counts):
        print(f"  [{i+1}/{n_subsets}] Scanning: {path}")
        adata = anndata.read_h5ad(path)

        if i == 0:
            obs = adata.obs.copy()
            uns = {k: v for k, v in adata.uns.items()}
            n_obs = adata.n_obs
            x_is_sparse = issparse(adata.X)
        else:
            if set(adata.obs_names.tolist()) != set(obs.index.tolist()):
                print(f"  Warning: Subset {i+1} has different observations")
            if len(adata.obs_names) != n_obs:
                print(f"  Warning: Subset {i+1} has {len(adata.obs_names)} obs, expected {n_obs}")

        var_dfs.append(adata.var.copy())
        total_vars += adata.n_vars
        print(f"    Shape: {adata.shape}")

        if issparse(adata.X):
            total_nnz_x += adata.X.nnz
        else:
            total_nnz_x += np.count_nonzero(adata.X)

        del adata
        gc.collect()

    # === Merge var DataFrames ===
    print("\nConcatenating var metadata...")
    merged_var = pd.concat(var_dfs, axis=0)
    del var_dfs
    gc.collect()

    # === Build minimal AnnData and write ===
    print(f"Building merged AnnData ({n_obs} x {total_vars})...")
    empty_x = csc_matrix((n_obs, total_vars))
    merged = anndata.AnnData(X=empty_x, obs=obs, var=merged_var)

    print(f"Writing base h5ad to: {counts}")
    merged.write_h5ad(counts)
    del merged, empty_x, obs, merged_var
    gc.collect()

    # === Second pass: fill X incrementally via h5py ===
    if x_is_sparse:
        print(f"\nWriting X incrementally (sparse, nnz={total_nnz_x:,})...")
        _write_sparse_incremental(
           counts, "X", subset_counts, n_obs, total_vars, total_nnz_x,
        )
    else:
        print("\nWriting X incrementally (dense)...")
        _write_dense_incremental(
            counts, "X", subset_counts, n_obs, total_vars,
        )

    print("Merging DESeq2 results...")

    factors = [d for d in os.listdir(os.path.join(query_dir, "subset_0")) if os.path.isdir(os.path.join(query_dir, "subset_0", d))]
    
    dfs = []
    for factor in factors:

        level_pairs = [os.path.basename(f).removesuffix('.tsv') for f in glob.glob(f"{query_dir}/subset_0/{factor}/*_vs_*.tsv")]

        for pair in level_pairs:

            files = [f"{query_dir}/subset_{subset}/{factor}/{pair}.tsv" for subset in range(int(n_subsets))]
            df = pd.concat(
                    [pd.read_csv(f, sep="\t", index_col=0) for f in files],
                )
            
            if "pvalue" not in df.columns:
                raise ValueError(f"No pvalue column found")
        
            pvals = df["pvalue"].to_numpy()
            padj = bh_fdr(pvals)

            df["padj"] = padj

            n_valid = int((~np.isnan(pvals)).sum())
            n_sig = int(np.nansum(padj < 0.05))

            print(f"{factor} {pair.replace('_vs_', ' vs ')}: {n_valid} non-NaN p-values, {n_sig} significant at FDR 0.05")

            df["factor"] = factor
            df["test"] = pair.replace('_vs_', ' vs ')
            df.index.name = "feature"

            dfs.append(df.reset_index())
    
    pd.concat(dfs, ignore_index=True).to_csv(os.path.join(query_dir, "deseq2_merged_results.tsv"), index=False, sep='\t')

    # === Summary ===
    print(f"\n{'='*60}")
    print("MERGE SUMMARY")
    print(f"{'='*60}")
    print(f"Input files: {n_subsets}")
    print(f"Total features: {total_vars}")
    print(f"Total samples: {n_obs}")

    print(f"Execution time: {time.time() - start_time:.2f} seconds")
    print(f"{'='*60}")

    # Copy results to parent directory
    results_parent = os.path.join(os.path.dirname(query_dir), "deseq2_merged_results.tsv")
    shutil.copy(results, results_parent)

if __name__ == "__main__":
    main()