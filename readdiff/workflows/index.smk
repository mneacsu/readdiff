import os
import glob
import math
from pathlib import Path
import heapq

wildcard_constraints:
    index_id = "[0-9]+",
    prefix = "[A-Za-zΑ-Ωα-ω0-9-_.]+",

class SplitIndex:
    def __init__(
                self,
                samples,
                metadata=None,
                paired_end=True,
                compress=True,
                test_index=False,
                # Options de découpage
                n_indices=None,      # Nombre d'indices voulus
                samples_per_index=None,  # Nombre de samples par index
                max_size_gb=None,    # Taille max en Go par index
                # Paramètres needle
                w=24,
                k=18,
                cutoff=2,
                l=20,
            ):
        self.samples = samples
        self.samples.sort()
        self.paired_end = paired_end
        self.compress = compress
        self.test_index = test_index

        if metadata is None or os.path.exists(metadata):
            self.metadata = metadata
        else:
            raise ValueError(f"File {metadata} does not exist")

        # Paramètres needle
        self.w = w
        self.k = k
        self.cutoff = cutoff
        self.l = l
        if l < 0 or l > 255:
            raise ValueError(f"l must be between 0 and 255, got {l}")

        # Détermination du nombre d'indices
        self.n_indices = self._compute_n_indices(n_indices, samples_per_index, max_size_gb)

        # Répartition des samples dans les indices
        self.index_samples = self._distribute_samples(samples_per_index=samples_per_index)

        # Préparation des chemins pour les fichiers mergés/trimmés
        self._prepare_sample_paths()

    def _compute_n_indices(self, n_indices, samples_per_index, max_size_gb):
        """Calcule le nombre d'indices selon les paramètres fournis"""

        if n_indices is not None:
            return n_indices
        elif samples_per_index is not None:
            # Si paired_end, on compte les paires
            n_samples = len(self.samples) // 2 if self.paired_end else len(self.samples)
            return math.ceil(n_samples / samples_per_index)
        elif max_size_gb is not None:
            # Estimation basée sur la taille totale
            total_size_gb = self._estimate_total_size_gb()
            return max(1, math.ceil(total_size_gb / max_size_gb))
        else:
            # Par défaut: 1 index pour tous les samples
            return 1

    def _estimate_total_size_gb(self):
        """Estime la taille totale des fichiers fastq en Go"""
        total_size = 0
        for sample in self.samples:
            if os.path.exists(sample):
                total_size += os.path.getsize(sample)
        return total_size / (1024 ** 3)  # Conversion en Go
    
    def _get_sample_sizes_gb(self):
        # Tolerate samples not existing yet at parse-time: they may be produced
        # by an upstream rule. The real sizes are read by lambdas in the rules.
        def _size(p):
            return os.path.getsize(p) if os.path.exists(p) else 0
        return (
            [(_size(self.samples[2*i]) + _size(self.samples[2*i+1])) / (1024 ** 3)
                for i in range(int(len(self.samples)/2))] if self.paired_end else 
            [_size(self.samples[i])/ (1024 ** 3)
                for i in range(len(self.samples))]
        )
    def _distribute_samples(self, samples_per_index=None):
        """
        Distributes samples into subindices so that total file sizes are approximately equal.
        """

        indexed_sample_sizes = list(enumerate(self._get_sample_sizes_gb()))
        indexed_sample_sizes.sort(key=lambda x: x[1], reverse=True)

        heap = [(0, i) for i in range(self.n_indices)]
        heapq.heapify(heap)

        groups = [[] for _ in range(self.n_indices)]

        for idx, value in indexed_sample_sizes:
            curr_sum, group_id = heapq.heappop(heap)
            groups[group_id].append(idx)
            if (samples_per_index is None) or (len(groups[group_id]) < samples_per_index):
                heapq.heappush(heap, (curr_sum + value, group_id))

        groups = [g for g in groups if g]

        index_samples = {}
        for i in range(len(groups)):
            index_samples[i] = ([ 
                s
                for idx in groups[i]
                for s in self.samples[2 * idx : 2 * idx + 2]
            ] if self.paired_end else [ 
                self.samples[idx]
                for idx in groups[i]
            ]

            )
        self.n_indices = len(groups)

        return index_samples

    
    def _prepare_sample_paths(self):
        """Prépare les chemins pour les fichiers mergés (paired-end uniquement)"""
        self.merged_samples = {}
        self.merged_from = {}

        for idx, samples in self.index_samples.items():
            self.merged_samples[idx] = {}
            self.merged_from[idx] = {}

            if self.paired_end:
                # Pour paired-end: créer les chemins de merge
                for i in range(0, len(samples), 2):
                    if i + 1 < len(samples):
                        sample1, sample2 = samples[i], samples[i + 1]

                        # Extraire le préfixe commun
                        base1 = os.path.basename(sample1)
                        base2 = os.path.basename(sample2)
                        common_prefix = self._get_common_prefix(base1, base2)

                        if '.' in common_prefix:
                            raise ValueError(f"Common prefix of {base1} and {base2} contains '.'")

                        # Chemins mergés
                        merged_path = f"tmp/merged/index_{idx}/{common_prefix}.fastq"
                        self.merged_samples[idx][common_prefix] = merged_path

                        # Sources pour le merge (toujours les fichiers originaux)
                        self.merged_from[idx][common_prefix] = [sample1, sample2]
            else:
                for i in range(len(samples)):
                    base = os.path.basename(samples[i]).removesuffix(".gz").removesuffix(".fastq").removesuffix(".fq")
                    if '.' in base:
                        raise ValueError(f"Basename of {base} contains '.'")
                    self.merged_samples[idx][base] = samples[i]
                    self.merged_from[idx][base] = samples[i]

    def _get_common_prefix(self, s1, s2):
        """Trouve le préfixe commun de deux chaînes"""
        prefix = ""
        for c1, c2 in zip(s1, s2):
            if c1 != c2:
                break
            prefix += c1
        # Nettoyer le préfixe
        prefix = prefix.rstrip("R").rstrip("._")
        return prefix

    def get_samples_for_index(self, idx):
        """Retourne les samples à utiliser pour créer un index spécifique"""
        return list(self.merged_samples[idx].values())
    
    def get_sample_names_for_index(self, idx):
        return list(self.merged_samples[idx].keys())
    
    @property
    def path(self):
        """Chemin du splitindex (sans slash final)"""
        return f"index"

    @property
    def index_paths(self):
        """Liste des chemins des indices individuels"""
        return [f"{self.path}/index_{idx}/" for idx in range(self.n_indices)]

# MERGE PAIRED-END FILES =======================================================
rule index_merge_paired_end:
    input:
        R1 = lambda wildcards: INDEX.merged_from[int(wildcards.index_id)][wildcards.prefix][0],
        R2 = lambda wildcards: INDEX.merged_from[int(wildcards.index_id)][wildcards.prefix][1],
    params:
        tmp_file = "tmp/merged/index_{index_id}/{prefix}.tmp.fastq",
        cat_command = lambda wildcards, input: "zcat" if input.R1.endswith(".gz") else "cat"
    output:
        merged = temp("tmp/merged/index_{index_id}/{prefix}.fastq")
    resources:
        mem_mb = lambda wildcards, attempt: min(MAX_RAM_MB, 1000 * attempt),
        runtime = lambda wildcards, input, attempt: max(10, 10 * get_size_mb(input.R1)/500) * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/index_merge_paired_end/index_{{index_id}}.{{prefix}}.benchmark.txt"
    shell:
        """
        mkdir -p $(dirname {output.merged})

        # Merge paired-end files using simple concatenation
        # R1 and R2 reads are concatenated into a single file
        {params.cat_command} {input.R1} {input.R2} > {params.tmp_file}
        mv {params.tmp_file} {output.merged}
        """


# CREATE INDEX ==================================================================
rule index_build_subindex:
    input:
        samples = lambda wildcards: INDEX.get_samples_for_index(int(wildcards.index_id)),
    output:
        IBF_Data =  f"{WORKDIR}/index/index_{{index_id}}/IBF_Data",
        IBF_FPRs =  f"{WORKDIR}/index/index_{{index_id}}/IBF_FPRs.fprs",
        IBF_Level_0 =  f"{WORKDIR}/index/index_{{index_id}}/IBF_Level_0",
        IBF_Levels =  f"{WORKDIR}/index/index_{{index_id}}/IBF_Levels.levels",
        stored_files =  f"{WORKDIR}/index/index_{{index_id}}/Stored_Files.txt",
        window_count =  f"{WORKDIR}/index/index_{{index_id}}/window_count.csv",
        read_count =  f"{WORKDIR}/index/index_{{index_id}}/read_count.csv",
    params:
        index_folder =  f"{WORKDIR}/index/index_{{index_id}}/",
        old_stored_files =  f"{WORKDIR}/index/index_{{index_id}}/Stored_Files_old.txt",
        compress_option = lambda wildcards: INDEX.compress,
        param_w = lambda wildcards: INDEX.w,
        param_k = lambda wildcards: INDEX.k,
        param_cutoff = lambda wildcards: INDEX.cutoff,
        param_l = lambda wildcards: INDEX.l,
    resources:
        mem_mb =lambda wildcards, input, attempt: min(MAX_RAM_MB, max(1000, max([get_size_mb(file) for file in input.samples]) * 20 * attempt)),
        runtime =lambda wildcards, input, attempt: max(10, max([get_size_mb(file) for file in input.samples]) / 1000 * len(input.samples) * 7) * attempt,
        attempt = lambda wildcards, attempt: attempt,
    threads: 8
    benchmark: f"{WORKDIR}/benchmark/index_build_subindex/index_{{index_id}}.benchmark.txt"
    script:
        "scripts/index_build_subindex.py"

# SAMPLE METADATA SPLITTING ====================================================
rule index_split_metadata:
    output:
        split_metadata =  f"{WORKDIR}/index/index_{{idx}}/metadata.csv"
    params:
        metadata = lambda wildcards: INDEX.metadata,
        samples = lambda wildcards: INDEX.get_sample_names_for_index(int(wildcards.idx)),
    benchmark: f"{WORKDIR}/benchmark/index_split_metadata/index_{{idx}}.benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt, input: min(MAX_RAM_MB, 500 * attempt),
        runtime = lambda wildcards, attempt, input: 5 * attempt,
        attempt = lambda wildcards, attempt: attempt,
    script:
        "scripts/index_split_metadata.py"

# INDEX METADATA AGGREGATION ==========================================================
rule index_aggregate_metadata:
    input:
        metadata = lambda wildcards: [f"{WORKDIR}/index/index_{idx}/metadata.csv"
                                         for idx in range(INDEX.n_indices)],
        window_counts = lambda wildcards: [f"{WORKDIR}/index/index_{idx}/window_count.csv"
                                         for idx in range(INDEX.n_indices)],
        read_counts = lambda wildcards: [f"{WORKDIR}/index/index_{idx}/read_count.csv"
                                         for idx in range(INDEX.n_indices)],
    output:
        global_window_count =  f"{WORKDIR}/index/window_count.csv",
        global_read_count =  f"{WORKDIR}/index/read_count.csv",
        global_metadata =  f"{WORKDIR}/index/metadata.csv",
        index_metadata =  f"{WORKDIR}/index/index_metadata.csv",
    params:
        n_indices = lambda wildcards: INDEX.n_indices,
        paired_end = lambda wildcards: INDEX.paired_end,
        compress = lambda wildcards: INDEX.compress,
    priority: 15
    resources:
        mem_mb = lambda wildcards, input, attempt: min(MAX_RAM_MB, 500 * attempt),
        runtime = lambda wildcards, input, attempt: len(input.window_counts) * attempt,
    benchmark: f"{WORKDIR}/benchmark/index_aggregate_metadata/benchmark.txt"
    script:
        "scripts/index_aggregate_metadata.py"