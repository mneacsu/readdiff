import os
from typing import Literal

class Query_estimate:
    def __init__(
                self,
                index,
                query_file,
                log_estimates=False,
                normalize:Literal["windows", "reads", True, False]="windows",
            ):
        self.log_estimates = log_estimates
        self.normalize = "windows" if (normalize==True) else normalize
        self.query_file = query_file
        self.index = index
        self.n_subindices = INDEX.n_indices

rule estimate_estimate:
    input:
        IBF_Level_0 = lambda wildcards: os.path.join(WORKDIR, "index", f"index_{wildcards.idx}", "IBF_Level_0"),
        IBF_FPRs = lambda wildcards: os.path.join(WORKDIR, "index", f"index_{wildcards.idx}", "IBF_FPRs.fprs"),
        IBF_Levels = lambda wildcards: os.path.join(WORKDIR, "index", f"index_{wildcards.idx}", "IBF_Levels.levels"),
        stored_files = lambda wildcards: os.path.join(WORKDIR, "index", f"index_{wildcards.idx}", "Stored_Files.txt"),
        query = QUERY,
    output:
        estimate = temp(f"{WORKDIR}/estimates/{{query}}.index_{{idx}}.estimate.out"),
    params:
        index = lambda wildcards: os.path.join(WORKDIR, "index", f"index_{wildcards.idx}"),
        tmp_dir = f"{WORKDIR}/tmp/estimate/index_{{idx}}",
        chunk_length = 100000000,
        verbose = True,
    resources:
        mem_mb = lambda wildcards, input, attempt: min(MAX_RAM_MB, max(1000, get_size_mb(input.IBF_Level_0) * 1.5 + get_size_mb(input.query) * 10)* attempt), 
        runtime = lambda wildcards, input, attempt: max(600, 15 * attempt),
        attempt = lambda wildcards, attempt: attempt,
    threads: 2,
    benchmark: f"{WORKDIR}/benchmark/estimate_estimate/{{query}}.index_{{idx}}.benchmark.txt",
    script:
        "scripts/estimate_estimate.py"


rule estimate_merge:
    input:
        sub_estimates = lambda wildcards: [f"{WORKDIR}/estimates/{wildcards.query}.index_{idx}.estimate.out" 
                            for idx in range(QUERY_ESTIMATE.n_subindices)]
    output:
        merged_estimate = temp(f"{WORKDIR}/estimates/{{query}}.merged_estimate.out"),
    resources:
        mem_mb = lambda wildcards, attempt: min(MAX_RAM_MB, 1000 * attempt),
        runtime = lambda wildcards, attempt: 10 * attempt,
        attempt = lambda wildcards, attempt: attempt,
    params:
    benchmark: f"{WORKDIR}/benchmark/estimate_merge/{{query}}.benchmark.txt"
    script:
        "scripts/estimate_merge.py"

rule estimate_to_sparse:
    input:
        estimate = f"{WORKDIR}/estimates/{{query}}.merged_estimate.out",
        metadata = f"{WORKDIR}/index/metadata.csv"
    output:
        sparse = f"{WORKDIR}/estimates/{{query}}_raw.h5ad",
    params:
        sparse_tmp1 = lambda wildcards: f"{WORKDIR}/estimates/{wildcards.query}_raw.sparse.tmp1.h5ad",
        sparse_tmp2 = lambda wildcards: f"{WORKDIR}/estimates/{wildcards.query}_raw.sparse.tmp2.h5ad",
        verbose = True,
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(500, get_size_mb(input.estimate) * 15) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.estimate) / 500 * 5) * attempt, 
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/estimate_to_sparse/{{query}}.benchmark.txt",
    script:
        "scripts/estimate_to_sparse.py"

rule estimate_log_and_norm:
    input:
        raw = f"{WORKDIR}/estimates/{{query}}_raw.h5ad",
    output:
        transformed = f"{WORKDIR}/estimates/{{query}}_transformed.h5ad",
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(500, get_size_mb(input.raw) * 15) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.raw) / 500 * 5) * attempt, 
        attempt = lambda wildcards, attempt: attempt,
    params:
        counts = lambda wildcards: (
            f"{WORKDIR}/index/window_count.csv" if QUERY_ESTIMATE.normalize == "windows" else
            f"{WORKDIR}/index/read_count.csv" if QUERY_ESTIMATE.normalize == "reads" else None
            ), 
        normalize = lambda wildcards: QUERY_ESTIMATE.normalize,
        log10 = lambda wildcards: QUERY_ESTIMATE.log_estimates,
    benchmark: f"{WORKDIR}/benchmark/estimate_log_norm/{{query}}.benchmark.txt",
    script:
        "scripts/estimate_log_and_norm.py"
        
