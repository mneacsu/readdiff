"""Pipeline registry — single source of truth for the CLI.

Adding a pipeline requires:
  1. An entry here (deps, outputs, rules).
  2. A target rule and its sub-rules in `workflows/Snakefile`.

`outputs` are paths relative to WORKDIR. They define when a pipeline is "done":
all listed outputs must exist.

`rules` is the list of Snakemake rule names that produce this pipeline's files.
It's used by `status.py` to map "snakemake -n would rerun X" back to a
pipeline-level state.
"""

from __future__ import annotations

PIPELINES: dict[str, dict] = {
    "index": {
        "deps": [],
        "outputs": [
            "index/window_count.csv",
        ],
        "rules": [
            "index_merge_paired_end",
            "index_build_subindex",
            "index_split_metadata",
            "index_aggregate_metadata"
        ],
    },
    "features": {
        "deps": [],
        "outputs": [
            "features/features.fasta"
        ],
        "rules": [
            "features_subsample"
        ],
    },
    "estimate": {
        "deps": ["index", "features"],
        "outputs": [
            "estimates/estimates_raw.h5ad"
        ],
        "rules": [
            "estimate_estimate",
            "estimate_merge",
            "estimate_log_and_norm",
            "estimate_to_sparse",
            "estimate_transformed_to_sparse"
        ],
    },
    "diff": {
        "deps": ["estimate"],
        "outputs": [
            "differential_testing/deseq2_merged_results.tsv"
        ],
        "rules": [
            "diff_filter",
            "diff_batch",
            "diff_deseq2",
            "diff_merge_deseq2"
        ],
    },
    "annotate": {
        "deps": ["diff"],
        "outputs": [
            "annotations/analysis_table.tsv"
        ],
        "rules": [
            "annotate_get_de_reads",
            "annotate_star_align",
            "annotate_concat_sam",
            "annotate_sam2bam",
            "annotate_bam2bed",
            "annotate_annotate",
            "annotate_homer",
            "annotate_merge_annotations",
            "annotate_table",
        ],
    },
}


def pipeline_names() -> list[str]:
    return list(PIPELINES.keys())


def transitive_deps(name: str) -> list[str]:
    """Return all upstream pipelines (deps of deps...), in topological order."""
    seen: list[str] = []
    def visit(n: str) -> None:
        for d in PIPELINES[n]["deps"]:
            if d not in seen:
                visit(d)
                seen.append(d)
    visit(name)
    return seen
