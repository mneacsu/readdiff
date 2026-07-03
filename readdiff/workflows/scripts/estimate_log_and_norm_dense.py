#!/usr/bin/env python3
import pandas as pd
import numpy as np

def main():
    estimate = snakemake.input.estimate
    log_norm_estimate = snakemake.output.log_norm_estimate
    counts = snakemake.params.counts
    normalize = snakemake.params.norm_estimates
    log_transform = snakemake.params.log_estimates

    if normalize and (counts is None):
        raise ValueError("Counts are required for normalization. Please provide them with --counts")

    estimates = pd.read_csv(estimate, sep="\t", index_col=0)

    if normalize:
        counts = pd.read_csv(counts, index_col=0)

        if "window_count" in counts.columns:
            scaling_factor = 1e11
            counts = counts.rename(columns={"window_count": "count"})
        elif "read_count" in counts.columns:
            scaling_factor = 1e9
            counts = counts.rename(columns={"read_count": "count"})
        else:
            raise ValueError("Counts file does not contain window_count / read_count column.")

        missing_obs = set(estimates.columns) - set(counts.index)
        if missing_obs:
            raise ValueError(f"Missing counts for: {missing_obs}")
        
        counts = counts.loc[estimates.columns, "count"]
        if any(counts == 0):
            raise ValueError(f"Counts for sample(s) '{[estimates.columns[i] for i, v in enumerate(counts) if v == 0]}' are 0.")

        estimates = estimates.mul(scaling_factor/counts)

    if log_transform:
        estimates = np.log10(estimates+1)

    estimates.round(2).to_csv(log_norm_estimate, sep="\t")


if __name__ == "__main__":
    main()