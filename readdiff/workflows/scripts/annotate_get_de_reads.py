#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import numpy as np

def main():

    deseq2 = snakemake.input.deseq2
    features = snakemake.input.features
    padj_threshold = snakemake.params.padj_threshold if snakemake.params.padj_threshold else 0.05
    log2fc_threshold = snakemake.params.log2fc_threshold if snakemake.params.log2fc_threshold else 0
    de_features = snakemake.output.de_features

    print(f"Loading DESeq2 results: {deseq2}")
    results = pd.read_csv(deseq2, sep='\t')
        
    significant = results[(results['padj'] < padj_threshold) & (np.abs(results['log2FoldChange']) > log2fc_threshold)]["feature"].unique()
    
    if not len(significant):
        print(f"No feature found with padj < {padj_threshold} and |log2FoldChange| > {log2fc_threshold}")
        return
    
    print(f"Significant features detected: {len(significant)}")

    found=[]
    with open(de_features, "w") as out:
        for record in SeqIO.parse(features, "fasta"):
            if record.id in significant:
                SeqIO.write(record, out, "fasta")
                found.append(record.id)
    
    if set(significant) - set(found):
        raise ValueError(f"Significant features not found in input fasta: {set(significant) - set(found)}")

if __name__ == "__main__":
    main()