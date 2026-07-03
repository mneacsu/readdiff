#!/usr/bin/env python3

import os
import random
import string
import shutil
import subprocess
from tqdm import tqdm

# -------------------------------------------------------------------------
# HELPER: Run needle estimate
# -------------------------------------------------------------------------
def estimate(query_fasta, index_dir, output, needle_executable, threads):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    needle_cmd = needle_executable if needle_executable else "needle"

    command = (
        f"{needle_cmd} estimate {query_fasta} "
        f"-i {index_dir}/ "
        f"-o {output} "
        f"-t {threads} "
        "--version-check false"
    )

    process = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {process.returncode}:\n{process.stderr}"
        )


# -------------------------------------------------------------------------
# HELPER: Line cleaner to remove trailing tab
# -------------------------------------------------------------------------
def clean_line(line):
    """
    Strip trailing newline and trailing tab, then split on '\t'.
    Returns (list_of_columns, cleaned_line).
    """
    line = line.rstrip("\n").rstrip("\t")
    columns = line.split("\t")
    cleaned_line = "\t".join(columns)
    return columns, cleaned_line


# -------------------------------------------------------------------------
# AGGREGATION + HEADER
# -------------------------------------------------------------------------
def aggregate_estimate(estimates, aggregated_estimate_path, stored_files):
    """
    Concatenates all 'needle estimate' outputs line by line into one file.
    Removes any trailing tab from each line.

    Prepend a single header line from 'stored_files':
        0\t<sample1>\t<sample2>\t...\t<sampleN>
    """
    if not estimates:
        raise ValueError("No estimated files found for aggregation.")

    # 1) We open the final aggregated file for writing
    with open(aggregated_estimate_path, "w") as fout:

        # 2) We open the *first* estimate file to read its first line
        with open(estimates[0], "r") as ffirst:
            first_line = ffirst.readline()
            if not first_line:
                raise ValueError(f"First file {estimates[0]} is empty.")

            # Remove trailing tab, parse columns
            first_cols, cleaned_first_line = clean_line(first_line)
            n_col = len(first_cols)

            # Read the stored_files lines (sample names)
            with open(stored_files, "r") as sf:
                sample_names = [ln.strip() for ln in sf]
            if len(sample_names) != (n_col - 1):
                raise ValueError(
                    f"Mismatch between sample names in '{stored_files}' "
                    f"({len(sample_names)}) and the data columns (n_col - 1 = {n_col - 1}). "
                    f"Each needle-estimate line has {n_col} columns total (after cleaning)."
                )

            # Build the header line and write it out
            header_line = "0\t" + "\t".join(sample_names) + "\n"
            fout.write(header_line)

            # Write the first line from the first file (already cleaned)
            fout.write(cleaned_first_line + "\n")

            # Write the rest of that first file
            for rest_line in ffirst:
                _, cleaned = clean_line(rest_line)
                fout.write(cleaned + "\n")

        # 3) Append the remaining estimate files
        for est_file in estimates[1:]:
            with open(est_file, "r") as fin:
                for line in fin:
                    _, cleaned = clean_line(line)
                    fout.write(cleaned + "\n")


# -------------------------------------------------------------------------
# TEMP DIR & FASTA CHUNKING
# -------------------------------------------------------------------------
def prepare_tmp_dir(tmp_dir):
    os.makedirs(tmp_dir, exist_ok=True)
    while True:
        run_name = "run_" + "".join(
            random.choices(string.ascii_letters + string.digits, k=8)
        )
        tmp_rundir = os.path.join(tmp_dir, run_name)
        if not os.path.exists(tmp_rundir):
            os.makedirs(tmp_rundir)
            break
    return tmp_rundir


def slice_fasta(fasta_path, chunk_length, slice_dir, verbose=False):
    """
    Splits the FASTA into chunks where each chunk has a total
    length of <= chunk_length nucleotides.

    If chunk_length <= 0, the entire FASTA is placed in a single chunk file.

    We never split a single sequence across multiple chunks. If one sequence
    alone exceeds chunk_length, that sequence alone occupies its own chunk.
    """
    os.makedirs(slice_dir, exist_ok=True)

    # If chunk_length <= 0, create a single chunk with all sequences
    if chunk_length <= 0:
        single_path = os.path.join(slice_dir, "single_chunk.fasta")
        shutil.copyfile(fasta_path, single_path)
        return [single_path]

    slice_files = []
    current_file = None
    current_filename = None
    current_chunk_length = 0
    seq_length = 0
    current_seq_header = None
    current_seq_lines = []

    def start_new_chunk():
        nonlocal current_file, current_filename, current_chunk_length
        if current_file:
            current_file.close()
        new_index = len(slice_files)
        current_filename = os.path.join(slice_dir, f"query_chunk_{new_index}.fasta")
        current_file = open(current_filename, "w")
        slice_files.append(current_filename)
        current_chunk_length = 0

    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                # If we already have a sequence buffered, write it out
                if current_seq_header is not None:
                    # If no chunk started yet OR the next sequence doesn't fit
                    if current_file is None or (
                        current_chunk_length + seq_length > chunk_length
                    ):
                        start_new_chunk()
                    # Write the buffered sequence
                    current_file.write(current_seq_header)
                    for seq_line in current_seq_lines:
                        current_file.write(seq_line)
                    current_chunk_length += seq_length

                # Prepare for the new sequence
                current_seq_header = line
                current_seq_lines = []
                seq_length = 0
            else:
                seq_length += len(line.strip())  # Count nucleotides
                current_seq_lines.append(line)

        # Write the last sequence if it exists
        if current_seq_header is not None:
            if current_file is None or (
                current_chunk_length + seq_length > chunk_length
            ):
                start_new_chunk()
            current_file.write(current_seq_header)
            for seq_line in current_seq_lines:
                current_file.write(seq_line)
            current_chunk_length += seq_length

    if current_file:
        current_file.close()

    return slice_files

def main():

    query = snakemake.input.query
    index = snakemake.params.index
    stored_files = snakemake.input.stored_files
    tmp_dir = snakemake.params.tmp_dir
    chunk_length = snakemake.params.chunk_length
    threads = snakemake.threads
    output_file = snakemake.output.estimate
    verbose = snakemake.params.verbose

    rundir = prepare_tmp_dir(tmp_dir)

    # 1) Slice the FASTA
    slice_files = slice_fasta(
        fasta_path=query,
        chunk_length=chunk_length,
        slice_dir=os.path.join(rundir, "slices"),
        verbose=verbose,
    )

    # 2) Estimate each chunk (sequentially)
    estimated_files = []
    for chunk_path in tqdm(slice_files, desc="Estimating", disable=not verbose):
        out_path = os.path.join(rundir, os.path.basename(chunk_path) + ".estimate.out")
        estimate(
            query_fasta=chunk_path,
            index_dir=index,
            output=out_path,
            needle_executable="needle",
            threads=threads,
        )
        estimated_files.append(out_path)

    # 3) Aggregate results with header from --stored-files
    aggregated_path = os.path.join(rundir, "aggregated_estimate.tsv")
    aggregate_estimate(estimated_files, aggregated_path, stored_files)

    # 4) Move final aggregated file to --output
    shutil.move(aggregated_path, output_file)

    # 5) Cleanup
    shutil.rmtree(rundir)

    # Remove tmp directory if it exists
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    main()


