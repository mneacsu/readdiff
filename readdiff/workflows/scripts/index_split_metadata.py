#!/usr/bin/env python3
import os
import pandas as pd

def main():

    split_metadata = snakemake.output.split_metadata
    metadata = snakemake.params.metadata
    samples = list(snakemake.params.samples
)
    if metadata is None:
        metadata_subindex=pd.DataFrame({"sample_id": samples})
    else:
        # Read metadata
        metadata = pd.read_csv(metadata, index_col=0).rename_axis('sample_id')
        missing_obs = set(samples) - set(metadata.index)
        if missing_obs:
            raise ValueError(
                f"Samples missing from metadata.csv: {missing_obs}"
            )
        metadata_subindex = metadata[metadata.index.isin(samples)]

    os.makedirs(os.path.dirname(split_metadata), exist_ok=True)
    metadata_subindex.to_csv(split_metadata)

if __name__ == "__main__":
    main()  

