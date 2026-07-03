# Documentation readdiff 0.1.0

## Commands

### `run`

```
readdiff run index   [WORKDIR] [-v 0|1|2] [--dry-run] [--unlock]
readdiff run features    [WORKDIR] [-v 0|1|2] [--dry-run] [--unlock]
readdiff run estimate   [WORKDIR] [-v 0|1|2] [--dry-run] [--unlock]
readdiff run diff   [WORKDIR] [-v 0|1|2] [--dry-run] [--unlock]
readdiff run annotate   [WORKDIR] [-v 0|1|2] [--dry-run] [--unlock]
readdiff run explorer   [WORKDIR]
```

| Pipeline | Description                                          | Outputs                                                 |
|----------|------------------------------------------------|---------------------------------------------------------|
| index    | Builds Needle index from samples | `index/` |
| features | Extracts features from reads| `features/features.fasta`             |
| estimate | Estimates the expression levels of extracted features on the Needle index | `estimates/estimates_raw.h5ad` |
| diff| Performs differential expression testing with DESeq2| `differential_testing/deseq2_results.tsv`, `differential_testing/deseq2_normalized_counts.h5ad` |
| annotate | Annotates significant features by aligning to a reference genome | `annotations/annotation_table.tsv`|
| explorer | Runs explorer app ||

Snakemake resolves dependencies natively, so `readdiff run diff` will trigger
`index`, `features`, and `estimate` if their outputs are missing. However, to run the explorer app, the differential testing results and, optionally, the annotations must already be present.

```
      index  ─┐
              |→ estimate → diff → (annotate) → explorer
   features  ─┘
```


Verbosity: `0` silent, `1` progress bar, `2` progress bar + scrolling logs.

### Advanced features
**Subindices**

By default, readdiff places all samples in a single index. If you want to split samples across multiple indices, you can specify this in `config.yaml`:
```yaml
n_indices: 4 
```
Increasing the number of indices will reduce the size of each index, thereby decreasing the RAM usage during estimation.

**Custom query files**

You can provide your own query file to be used instead of randomly sampled features.
```yaml
query_file: "path/to/your/query.fasta" 
```

**Normalization and log transformation**

If you would like normalized (by window or read counts) and/or log transformed Needle estimates:
```yaml
 normalize: "windows" # or "reads"
 log10: True
```
In this case, the `estimate` module will produce an additional `estimates/estimates_transformed.h5ad` file. However, keep in mind that the `diff` module will still use the raw estimates as input, since DESeq2 does its own normalization.


### `status`

Renders a 6-column table — pipeline name, state, rule rerun
counts, output counts, dependency list, and age of the most recent output.

```
readdiff status [WORKDIR] 
```

Output:

```
  workdir: /…/_data/test_ws

  pipeline    state          rules        outputs    deps           last run
  alpha       ● done         3/3 ok       2/2        —              5s ago
  beta        ● stale        2/2 rerun    2/2        —              7s ago
  gamma       ● stale        3/3 rerun    3/3        alpha, beta    6s ago
  delta       ○ pending      2/2 queued   0/2        gamma          —

```

The glyph carries two signals:

- **shape** — `●` filled means the output exists, `○` empty means it doesn't.
- **color** — green = OK, orange = needs redo, default = your turn, dim = waiting.

| State     | Glyph       | Meaning                                                                         |
|-----------|-------------|---------------------------------------------------------------------------------|
| `done`    | `●` green   | All outputs present, dependencies done, Snakemake reports nothing to redo.      |
| `stale`   | `●` orange  | Outputs present but Snakemake would rerun something here (input newer, etc.).   |
| `missing` | `○` default | At least one output absent; dependencies are done so this pipeline can run now. |
| `pending` | `○` dim     | At least one dependency not done; this pipeline can't run yet.                  |

Stale detection runs `snakemake -n delta` once and parses which rules would
fire — it matches exactly what an actual run would do.

Topology is communicated via the `deps` column; we deliberately don't draw an
ASCII DAG (tried, the table won on every axis: scannability, density, no
dependence on terminal width).

Exit code is `0` iff every pipeline is `done` (useful in CI).

### `clean`

Dry-run by default. Pass `--yes` to actually delete.

```bash
readdiff clean [WORKDIR] --pipeline [PIPELINE]          # preview
readdiff clean [WORKDIR] --pipeline [PIPELINE] --yes    # delete outputs of specific pipeline
readdiff clean [WORKDIR] --all --yes                    # delete every pipeline + .readdiff/
```

### `config`

The CLI loads parameters from three sources, with the following priority (highest to lowest) in case of conflicts:

1. `WORKDIR/config.yaml` — *project* scope
2. `~/.config/readdiff/config.yaml` — *global* scope
3. embedded `readdiff/default_config.yaml`


```bash
readdiff config show                              # active scope, path, candidates, values
readdiff config get [PARAM]
readdiff config set [PARAM] [VALUE]               # writes to WORKDIR/config.yaml
readdiff config set [PARAM] [VALUE] --global      # writes to ~/.config/readdiff/config.yaml
readdiff config unset [PARAM]
EDITOR=nvim readdiff config edit                  # auto-creates if absent
readdiff config validate                          # checks unknown keys / type mismatches
```

There is no `init` step: `set` and `edit` create the target file on first write.

Write commands target the project scope by default and the global scope with
`--global`.

## Parameters

| Key               | Description | Default            |
|-------------------|--------------------|-----------|
| sample_dir | Directory containing sample files | samples |
| paired | Paired-end library layout | True |
| metadata | Sample metadata file (tabular) | None |
| kmer_size | Needle k-mer size | 18 |
| window_size | Needle window size | 24|
| n_bins | Needle no. of expression thresholds | 20 |
| cutoff | Needle cutoff | 2 |
| n_indices | No. of subindices to build (samples will be split between subindices) | 1 |
| random_seed | Random seed used for feature extraction | 42 |
| n_features | No. of features to extract | 100000 |
| feature_len | Feature length (in bp) | 100 |
| balance | Metadata columns to use for balanced sampling upon feature extraction | None | 
| query_file | Custom query file to use instead of randomly sampled features | None |
| normalize | Normalize estimates by window counts | False |
| log10 |Apply log10 transformation to estimates | False |
| test_factors | Metadata columns to test for differences | None |
| batch_factors | Metadata columns to use for batch correction | None |
| padj_threshold | Significance (p-value) threshold for subsetting features for downstream analysis| 0.05 |
| log2fc_threshold | Log2 fold change threshold for subsetting features for downstream analysis| 0 |
| star_index | STAR index to use for annotation | None |
| gtf_annotation | GTF annotation file (VCF) | None |
| homer_genome | Homer reference genome | hg38 |
| volcano_config | Config file for volcano explorer | volcano_explorer/volcano_config.yaml |
| max_ram | RAM available (MB) | 64000 |
| snakemake_cores | Number of cores to use | 4 |
| snakemake_profile | Snakemake profile | None |

