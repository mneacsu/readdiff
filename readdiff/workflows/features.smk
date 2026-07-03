rule features_subsample:
    input:
        fastqs = SAMPLES,
        metadata = METADATA,
    output:
        features = f"{WORKDIR}/features/random_{{n_features}}features_seed{{random_seed}}.fasta",
    params:
        balance_cols = BALANCE,
        chunk_length = FEATURE_LEN,
        paired = PAIRED,
        shuffle = True,
        check_unique = True,
        verbose = True
    threads: lambda wildcards, input: min(len(input.fastqs), 8)
    resources:
        mem_mb = lambda wildcards, attempt: min(MAX_RAM_MB, 4000 * attempt),
        runtime = lambda wildcards, attempt, input: 120 * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/features_subsample/{{n_features}}.{{random_seed}}.benchmark.txt"
    script: "scripts/features_subsample.py"
