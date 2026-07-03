#!/usr/bin/env python3
import os
import pandas as pd
import anndata
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from itertools import combinations
import shutil

def prepare_counts_dataframe(adata):
    """
    Prepare counts DataFrame for PyDESeq2.
    PyDESeq2 expects samples as columns and features as rows.
    """
    if hasattr(adata.X, "toarray"):
        counts = adata.X.toarray()
    else:
        counts = adata.X

    counts_df = pd.DataFrame(counts, columns=adata.var_names, index=adata.obs_names)
    counts_df = counts_df.round().astype(int)
    return counts_df

def build_design(factor, batches):
    terms = []
    if batches:
        terms.extend(batches)
    terms.append(factor)
    return " + ".join(terms)

def main():
    sparse = snakemake.input.sparse
    test_factors = snakemake.params.test_factors
    batch_factors = snakemake.params.batch_factors
    threads = snakemake.threads
    norm_counts_file = snakemake.output.norm_counts

    out_dir = os.path.dirname(norm_counts_file) 
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    print(f"Reading input file: {sparse}")
    adata = anndata.read_h5ad(sparse)
    print(f"Input shape: {adata.shape} (obs x var)")

    # Verify batch columns exist in the metadata
    if batch_factors:
        if set(batch_factors) - set(adata.obs.columns):
            raise ValueError(
                f"Batch column(s) '{list(set(batch_factors) - set(adata.obs.columns))}' not found in obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        print(f"Using batch correction with column(s): {batch_factors}")

    # Prepare data
    counts_df = prepare_counts_dataframe(adata)

    metadata = adata.obs.astype(str).copy()
    metadata.index = adata.obs_names           

    # Run DESeq2 for each condition
    for factor in test_factors:

        column = factor.get("column")
        test_mode = factor.get("test_mode")
        control = factor.get("control")

        # Verify test column exists in the metadata
        if column not in metadata.columns:
            raise ValueError(
                f"Condition column(s) '{column}' not found in obs. "
                f"Available columns: {list(metadata.columns)}"
            )
        
        # Prepare output path
        os.makedirs(os.path.join(out_dir, column), exist_ok=True)
        
        conditions = metadata[column].unique().tolist()
        
        if len(conditions) < 2:
            raise ValueError(f"Column {column} does not have at least 2 levels. Differential expression testing is not possible.")
        
        design = build_design(column, batch_factors)

        print(f"\nRunning DESeq2 for factor: {column}")
        print("Design:", design)

        if test_mode == "one-vs-control":

            if control not in conditions:
                raise ValueError(f"Control condition '{control}' not found in column '{column}'.")

            dds = DeseqDataSet(
                counts=counts_df, 
                metadata=metadata,
                design=design,
                refit_cooks=True,
                n_cpus = threads
            )
            dds.deseq2()
        
            for cond in conditions:
                if cond != control:
                    stat_res = DeseqStats(dds, contrast=(column, cond, control))
                    stat_res.summary()
                    stat_res.results_df.to_csv(os.path.join(out_dir, column, f"{cond}_vs_{control}.tsv"), sep="\t")

        elif (test_mode == "pairwise") or (len(conditions)==2):

            print(conditions)

            if conditions == ["False", "True"]:
                conditions = ["True", "False"]

            dds = DeseqDataSet(
                counts=counts_df, 
                metadata=metadata,
                design=design,
                refit_cooks=True,
                n_cpus = threads
            )
            dds.deseq2()
        
            for cond1, cond2 in list(combinations(conditions, 2)):
                stat_res = DeseqStats(dds, contrast=(column, cond1, cond2))
                stat_res.summary()
                stat_res.results_df.to_csv(os.path.join(out_dir, column, f"{cond1}_vs_{cond2}.tsv"), sep="\t")
        
        elif test_mode == "one-vs-rest":

            for cond in conditions:
                
                meta_tmp = metadata.copy()
                meta_tmp[column] =  meta_tmp[column].apply(
                            lambda x: cond if x == cond else "rest"
                        )

                dds = DeseqDataSet(
                    counts=counts_df, 
                    metadata=meta_tmp,
                    design=design,
                    refit_cooks=True,
                    n_cpus = threads
                )
                dds.deseq2()

                stat_res = DeseqStats(dds, contrast=(column, cond, "rest"))
                stat_res.summary()
                stat_res.results_df.to_csv(os.path.join(out_dir, column, f"{cond}_vs_rest.tsv"), sep="\t")
        else:
            raise ValueError("Unknown test mode")
        
    # Save normalized counts
    norm_counts = anndata.AnnData(
        X=dds.layers["normed_counts"], obs=adata.obs.copy(), var=adata.var.copy()
    )
    norm_counts.write_h5ad(norm_counts_file)


if __name__ == "__main__":
    main()