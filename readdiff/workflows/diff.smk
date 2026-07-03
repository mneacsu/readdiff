# =============================================================================
# RULE BATCH - Deterministic batch splitting
# =============================================================================


checkpoint diff_batch:
    input:
        sparse = f"{WORKDIR}/estimates/{{query}}_raw.h5ad",
    output:
        subsets_dir = directory(f"{WORKDIR}/differential_testing/subsets/{{query}}"),
    params:
        n_subsets = N_BATCHES,
    resources:
        mem_mb = lambda wildcards, attempt, input: max(2000, int(get_size_mb(input.sparse)**1.1 * 3)) * attempt**2,
        runtime = lambda wildcards, attempt, input: max(10, int(get_size_mb(input.sparse) / 1000 * 5)) * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/diff_batch/{{query}}.benchmark.txt",
    script:
        "scripts/diff_batch.py"

def get_subset(wildcards):
    ckpt = checkpoints.diff_batch.get(query=wildcards.query)
    return f"{WORKDIR}/differential_testing/subsets/{wildcards.query}/subset_{wildcards.subset}.h5ad"

# =============================================================================
# RULE DESEQ2 - Run PyDESeq2 for all conditions
# =============================================================================
rule diff_deseq2:
    input:
        sparse = get_subset,
    output:
        norm_counts = temp(f"{WORKDIR}/differential_testing/{{query}}/subset_{{subset}}/deseq2_normalized_counts.h5ad"),
    params:
        test_factors = TEST_FACTORS,
        batch_factors = BATCH_FACTORS,
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(5000, int(get_size_mb(input.sparse)**0.5 * 10)) * 2 ** attempt),
        runtime = lambda wildcards, attempt, input: max(15000, int(get_size_mb(input.sparse) * 50)) * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/diff_deseq2/{{query}}.subset_{{subset}}.benchmark.txt",
    script:
        "scripts/diff_deseq2.py"


# =============================================================================
# RULE MERGE_DESEQ2 - Merge results from all subsets
# =============================================================================

rule diff_merge_deseq2:
    input:
        subset_counts = lambda wildcards: [
            f"{WORKDIR}/differential_testing/{{query}}/subset_{subset}/deseq2_normalized_counts.h5ad" 
            for subset in range(N_BATCHES)
            ],
    output:
        counts = f"{WORKDIR}/differential_testing/{{query}}/deseq2_norm_counts.h5ad",
        results = f"{WORKDIR}/differential_testing/{{query}}/deseq2_merged_results.tsv" 
    params:
        query_dir = f"{WORKDIR}/differential_testing/{{query}}",
        n_subsets = N_BATCHES,
    resources:
        mem_mb = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(2000, max([get_size_mb(str(result)) for result in input.subset_counts]) * 5)) * attempt,
        runtime = lambda wildcards, attempt, input: max(10, int(get_size_mb(input.subset_counts, path_list=True) / 1000 * 5)) * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/diff_merge_deseq2/{{query}}.benchmark.txt",
    script:
        "scripts/diff_merge_deseq2.py"