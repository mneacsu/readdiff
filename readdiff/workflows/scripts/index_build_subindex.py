#!/usr/bin/env python3
"""Wrapper script for index_build_subindex rule."""
import os
import csv
import subprocess
import shutil
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed

def count_windows_in_fastq(args):
    """Compte le nombre de windows dans un fichier FASTQ.

    Le nombre de windows d'un read de longueur L avec une taille de window w
    est max(0, L - w + 1).
    """
    fastq_path, window_size = args
    total_windows = 0
    total_reads = 0
    open_func = gzip.open if fastq_path.endswith('.gz') else open

    try:
        with open_func(fastq_path, 'rt') as f:
            line_num = 0
            for line in f:
                line_num += 1
                # Les lignes de séquence sont à la position 2, 6, 10, ... (1-indexed)
                # Soit line_num % 4 == 2
                if line_num % 4 == 2:
                    read_length = len(line.strip())
                    windows_in_read = max(0, read_length - window_size + 1)
                    total_windows += windows_in_read
                    total_reads += 1
    except Exception as e:
        print(f"Error reading {fastq_path}: {e}")
        return os.path.basename(fastq_path), 0

    return os.path.basename(fastq_path), total_windows, total_reads

def process_file_paths(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)

        for line in infile:
            # Supprime les caractères de nouvelle ligne et espaces superflus
            trimmed_line = line.strip()
            # Extrait le basename du chemin du fichier
            basename = os.path.basename(trimmed_line)
            # Conserve tout ce qui précède le premier '.'
            modified_name = basename.split(".")[0]
            writer.writerow([modified_name])

def main():
    samples = list(snakemake.input.samples)
    index_folder = snakemake.params.index_folder
    old_stored_files = snakemake.params.old_stored_files
    param_w = snakemake.params.param_w
    param_k = snakemake.params.param_k
    param_cutoff = snakemake.params.param_cutoff
    param_l = snakemake.params.param_l
    threads = snakemake.threads
    stored_files = snakemake.output.stored_files
    window_counts_path = snakemake.output.window_count
    read_counts_path = snakemake.output.read_count
    compress_option = snakemake.params.compress_option

    # Create parent directory
    os.makedirs(index_folder, exist_ok=True)

    # Run needle ibf command
    needle_cmd = [
        "needle", "ibf", *samples,
        "-w", str(param_w),
        "-k", str(param_k),
        "--cutoff", str(param_cutoff),
        "-l", str(param_l),
        "-n", "5",
        "-f", "0.05",
        "-t", str(threads),
        "--experiment-names", "1",
        "-o", index_folder,
        "--version-check", "0",
        "--ram"
    ]
    if compress_option:
        needle_cmd.append("--compressed")

    subprocess.run(needle_cmd, check=True)

    # Count windows and reads

    window_counts = {}
    read_counts = {}

    # Préparer les arguments pour le traitement parallèle
    task_args = [(fastq, param_w) for fastq in samples]

    # Traitement parallèle
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(count_windows_in_fastq, arg): arg[0]
                  for arg in task_args}

        for future in as_completed(futures):
            sample_name, window_count, read_count = future.result()
            # Nettoyer le nom du sample
            sample_name = sample_name.replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '')
            window_counts[sample_name] = window_count
            read_counts[sample_name] = read_count

    # Écrire le CSV
    with open(window_counts_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample_id', 'window_count'])
        for sample, count in sorted(window_counts.items()):
            writer.writerow([sample, count])
    with open(read_counts_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample_id', 'read_count'])
        for sample, count in sorted(read_counts.items()):
            writer.writerow([sample, count])

    # Update stored files
    shutil.move(stored_files, old_stored_files)
    process_file_paths(old_stored_files, stored_files)
    os.remove(old_stored_files)

if __name__ == "__main__":
    main()