
rule annotate_get_de_reads:
    input:
        deseq2 = f"{WORKDIR}/differential_testing/{{query}}/deseq2_merged_results.tsv",
        features = f"{WORKDIR}/features/{{query}}.fasta"
    output:
        de_features = temp(f"{WORKDIR}/annotations/de_features/de_{{query}}.fasta")
    params:
        padj_threshold = PADJ_THRESHOLD,
        log2fc_threshold = LOG2FC_THRESHOLD,
    resources:
        mem_mb = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(4000, int(get_size_mb(input.features) * 2))) * attempt,
        runtime = lambda wildcards, attempt: 10 * attempt,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_get_de_reads/{{query}}.benchmark.txt"  
    script:
        "scripts/annotate_get_de_reads.py"


rule annotate_star_align:
    input:
        fasta = f"{WORKDIR}/annotations/de_features/de_{{query}}.fasta",
        star_index = STAR_INDEX,
    output:
        sam_main = temp(f"{WORKDIR}/annotations/alignments/{{query}}/Aligned.out.sam"),
        sam_chim = temp(f"{WORKDIR}/annotations/alignments/{{query}}/Chimeric.out.sam"),
    params:
        tmp_dir = f"{WORKDIR}/annotations/alignments/{{query}}", 
    threads: 8,
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(10000, get_size_mb(input.fasta) + 30000 * 1.5) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.fasta) / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_star_align/{{query}}.benchmark.txt",
    shell:
        """
        STAR --twopassMode Basic \
            --runThreadN {threads} \
            --genomeDir {input.star_index} \
            --readFilesIn {input.fasta} \
            --outFileNamePrefix {params.tmp_dir}/ \
            --outFilterMultimapNmax 1000 \
			--winAnchorMultimapNmax 1000 \
            --outSAMtype SAM \
            --outSAMunmapped Within \
            --outSAMattributes NH HI NM MD AS nM jM jI \
            --chimOutType SeparateSAMold \
            --bamRemoveDuplicatesType UniqueIdentical \
			--outFilterMismatchNoverLmax 0.04 \
			--outMultimapperOrder Random \
			--chimSegmentMin 10 \
			--chimJunctionOverhangMin 10 \
            --limitOutSAMoneReadBytes 1000000
        """

rule annotate_concat_sam:
    input:
        sam_main = f"{WORKDIR}/annotations/alignments/{{query}}/Aligned.out.sam",
        sam_chim = f"{WORKDIR}/annotations/alignments/{{query}}/Chimeric.out.sam",
    output:
        sam_concat = temp(f"{WORKDIR}/annotations/alignments/{{query}}/merged.sam"),
    params:
        tmp_sam_concat = f"{WORKDIR}/annotations/alignments/{{query}}.tmp.sam",
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(1000, get_size_mb(input.sam_main)) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.sam_chim) / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_concat_sam/{{query}}.benchmark.txt",
    shell:
        """
        samtools merge -f {params.tmp_sam_concat} {input.sam_main} {input.sam_chim}

        mv {params.tmp_sam_concat} {output.sam_concat}
        """

rule annotate_sam2bam:
    input:
        sam = f"{WORKDIR}/annotations/alignments/{{query}}/merged.sam",
    output:
        bam = temp(f"{WORKDIR}/annotations/alignments/{{query}}.bam"),
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(1000, get_size_mb(input.sam)) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.sam) / 100 * 5 / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_sam2bam/{{query}}.benchmark.txt",
    shell:
        """
        samtools view -bS {input.sam} > {output.bam}
        """

rule annotate_bam2bed:
    input:
        bam = f"{WORKDIR}/annotations/alignments/{{query}}.bam"
    output:
        bed = temp(f"{WORKDIR}/annotations/alignments/{{query}}.bed")
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(1000, get_size_mb(input.bam)) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.bam) / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_bam2bed/{{query}}.benchmark.txt",
    shell:
        """   
        bedtools bamtobed -split -i {input.bam} > {output.bed}
        """

rule annotate_annotate:
    input:
        bed = f"{WORKDIR}/annotations/alignments/{{query}}.bed",
        gtf = GTF_ANNOTATION or [],
    output:
        genes_bed = temp(f"{WORKDIR}/annotations/annotations/{{query}}.genes.bed"),
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(70000, get_size_mb(input.bed) + 3000) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.bed) / 100 * 5 + 3000 / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_annotate/{{query}}.benchmark.txt",
    shell:
        """   
        bedtools intersect -wao -a {input.bed} -b {input.gtf} > {output.genes_bed}
        """

rule annotate_homer:
    input:
        bed = f"{WORKDIR}/annotations/alignments/{{query}}.bed",
    output:
        homer_annot = temp(f"{WORKDIR}/annotations/annotations/{{query}}.homer.txt"),
    params:
        homer_genome = HOMER_GENOME,
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(50000, get_size_mb(input.bed)) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.bed) / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_homer/{{query}}.benchmark.txt",
    shell:
        """
        annotatePeaks.pl {input.bed} {params.homer_genome} > {output.homer_annot}
        """


rule annotate_merge_annotations:
    input:
        sam = f"{WORKDIR}/annotations/alignments/{{query}}/merged.sam",
        genes_bed = f"{WORKDIR}/annotations/annotations/{{query}}.genes.bed",
        homer_annot = f"{WORKDIR}/annotations/annotations/{{query}}.homer.txt",
    output:
        annotated_sam = temp(f"{WORKDIR}/annotations/{{query}}.annotated.sam"),
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(5000, get_size_mb(input.sam) + get_size_mb(input.genes_bed) * 10) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.sam) / 100 * 5 + get_size_mb(input.genes_bed) / 100 * 5) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_merge_annotations/{{query}}.benchmark.txt",
    script:
        "scripts/annotate_merge_annotations.py"


rule annotate_table:
    input:
        sam = f"{WORKDIR}/annotations/{{query}}.annotated.sam"
        #kraken_taxonomy = lambda wildcards: get_fasta_analysis(wildcards.project).kraken_taxonomy,
    output:
        table = f"{WORKDIR}/annotations/{{query}}.analysis_table.tsv",
    params:
        verbose = True,
        sortbyname = False
    resources:
        mem_mb  = lambda wildcards, attempt, input: min(MAX_RAM_MB, max(1000, get_size_mb(input.sam) ) * attempt**2),
        runtime = lambda wildcards, attempt, input: max(10, get_size_mb(input.sam) / 100 * 5 ) * attempt**2,
        attempt = lambda wildcards, attempt: attempt,
    benchmark: f"{WORKDIR}/benchmark/annotate_table/{{query}}.benchmark.txt",
    script:
        "scripts/annotate_table.py"