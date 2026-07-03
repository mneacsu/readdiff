#!/usr/bin/env python3
"""
Sélectionne aléatoirement des reads parmi plusieurs fichiers FASTQ
et produit un fichier FASTA.
Utilise des outils bas niveau (wc -l, awk, cat) pour minimiser
l'empreinte mémoire et le temps d'exécution.
"""

import gzip
import hashlib
import os
import random
import subprocess
import sys
import tempfile
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import shutil

def deterministic_hash(s: str) -> int:
    """Retourne un hash déterministe d'une chaîne (contrairement à hash() de Python)."""
    return int(hashlib.md5(s.encode()).hexdigest(), 16) % (2**32)


def is_gzipped(filename: str) -> bool:
    return filename.endswith(".gz")


def get_cat_command(filename: str) -> list[str]:
    """Retourne la commande pour lire le fichier (zcat pour .gz, cat sinon)."""
    if is_gzipped(filename):
        return ["zcat", filename]
    return ["cat", filename]


def count_lines(filename: str) -> int:
    """Compte le nombre de lignes avec wc -l (très rapide)."""
    cat_cmd = get_cat_command(filename)
    # Pipe: zcat/cat | wc -l
    cat_proc = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wc_proc = subprocess.Popen(
        ["wc", "-l"],
        stdin=cat_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    cat_proc.stdout.close()
    stdout, stderr = wc_proc.communicate()
    if wc_proc.returncode != 0:
        raise RuntimeError(f"wc -l failed for {filename}: {stderr.decode()}")
    return int(stdout.decode().strip())

def get_sample_name(filename: str) -> str:
    """Extrait le nom du sample (avant le premier . dans le basename)."""
    return os.path.basename(filename).removesuffix('.gz').removesuffix('.fastq').removesuffix('.fq')

def get_sample_name_paired(filename1: str, filename2: str) -> str:
    """Extrait le nom du sample (avant le premier . dans le basename)."""
    prefix = ""
    for c1, c2 in zip(os.path.basename(filename1), os.path.basename(filename2)):
        if c1 != c2:
            break
        prefix += c1
    # Nettoyer le préfixe
    prefix = prefix.rstrip("R").rstrip("._")
    return prefix

def extract_reads(
    filename: str,
    read_positions: list[int],
    output_file: str,
) -> int:
    """
    Extrait les reads sélectionnés d'un fichier FASTQ et les écrit en FASTA.
    Utilise awk pour une extraction efficace sans charger le fichier en mémoire.
    Retourne le nombre de reads extraits.
    """
    if not read_positions:
        # Créer un fichier vide
        Path(output_file).touch()
        return 0

    # Trier les positions pour un traitement séquentiel
    sorted_positions = sorted(read_positions)

    # Calculer les numéros de lignes header et seq pour chaque read
    # Read N (1-indexed) -> header: 4*(N-1)+1, seq: 4*(N-1)+2
    last_seq_line = 4 * sorted_positions[-1]  + 2

    # Créer un fichier temporaire avec les numéros de lignes
    # Format: H pour header, S pour sequence, suivi du numéro de ligne
    lines_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".lines")
    try:
        for pos in sorted_positions:
            header_line = 4 * pos + 1
            seq_line = 4 * pos + 2
            lines_file.write(f"H {header_line}\n")
            lines_file.write(f"S {seq_line}\n")
        lines_file.close()

        awk_script = f'''
            BEGIN {{
                last_line = {last_seq_line}
                # Lire le fichier de lignes
                while ((getline line < "{lines_file.name}") > 0) {{
                    split(line, parts, " ")
                    if (parts[1] == "H") header_set[parts[2]] = 1
                    else if (parts[1] == "S") seq_set[parts[2]] = 1
                }}
                close("{lines_file.name}")
            }}
            NR in header_set {{
                name = substr($0, 2)
                print ">" name
            }}
            NR in seq_set {{
                print
            }}
            NR >= last_line {{
                exit
            }}
            '''

        cat_cmd = get_cat_command(filename)

        # Pipeline: cat/zcat | awk > output
        with open(output_file, "w") as outf:
            cat_proc = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            awk_proc = subprocess.Popen(
                ["awk", awk_script],
                stdin=cat_proc.stdout,
                stdout=outf,
                stderr=subprocess.PIPE,
            )
            cat_proc.stdout.close()
            _, stderr = awk_proc.communicate()

            if awk_proc.returncode != 0:
                raise RuntimeError(f"awk failed for {filename}: {stderr.decode()}")

    finally:
        # Nettoyer le fichier temporaire
        if os.path.exists(lines_file.name):
            os.remove(lines_file.name)

    return len(read_positions)

def extract_features(
    filename: str,
    n_features: int,
    chunk_length: int,
    seed: int,
    sample_name: str,
    output_file: str,
) -> int:
    """
    Extrait les reads sélectionnés d'un fichier FASTQ et les écrit en FASTA.
    Utilise awk pour une extraction efficace sans charger le fichier en mémoire.
    Retourne le nombre de reads extraits.
    """
    if not n_features:
        # Créer un fichier vide
        Path(output_file).touch()
        return 0

    n_extracted = 0

    with open(output_file, "w") as out:
        for i, record in enumerate(SeqIO.parse(filename, "fasta")):
            if (n_extracted < n_features) and (len(record.seq) >= chunk_length) and ('N' not in str(record.seq)):
                local_random = random.Random(seed + deterministic_hash(record.id))
                selected_position = local_random.sample(range(len(record.seq) - chunk_length + 1), 1)[0]
                selected_window = SeqIO.SeqRecord(record.seq[selected_position : (selected_position + chunk_length)],
                                        id=f"{sample_name}_{record.id}_{selected_position+1}", description="")
                SeqIO.write(selected_window, out, "fasta")
                n_extracted += 1

    return n_extracted


def process_one_file(
    args_tuple: tuple,
) -> tuple[str, str, int]:
    """
    Traite un fichier FASTQ: compte les reads, sélectionne aléatoirement,
    extrait en FASTA.

    args_tuple = (filename, nfeatures_to_extract, seed, temp_dir)

    Retourne (filename, temp_output_path, nreads_extracted)
    """
    filename, nreads_to_extract, chunk_length, seed, temp_dir = args_tuple

    # 1. Compter les lignes
    n_lines = count_lines(filename)
    n_reads = n_lines // 4

    if n_reads < nreads_to_extract:
        raise ValueError(
            f"{filename} ne contient que {n_reads} reads, "
            f"mais on en demande {nreads_to_extract}."
        )

    # 2. Sélectionner aléatoirement les positions des reads
    # Utiliser une seed déterministe basée sur le fichier
    local_random = random.Random(seed + deterministic_hash(filename))
    selected_positions = sorted(local_random.sample(range(n_reads), math.ceil(nreads_to_extract * 1.2)))

    # 3. Extract the reads
    sample_name = get_sample_name(filename)
    temp_output = os.path.join(temp_dir, f"{sample_name}.tmp.fasta")

    _ = extract_reads(
        filename, selected_positions, temp_output
    )

    # 4. Extract the features
    temp_output2 = os.path.join(temp_dir, f"{sample_name}.tmp2.fasta")
    n_extracted = extract_features(
        temp_output, nreads_to_extract, chunk_length, seed, sample_name, temp_output2
    )
    if n_extracted < nreads_to_extract:
        raise ValueError(
            f"{nreads_to_extract} could not be extracted from {filename}. Consider reducing the feature length or filtering out reads containing Ns."
        )
    return filename, temp_output2, n_extracted


def shuffle_fasta(input_file: str, output_file: str, seed: int):
    """
    Mélange les reads d'un fichier FASTA.
    Lit les reads par paires de lignes (header + seq), mélange, réécrit.
    Pour les très gros fichiers, on utilise un fichier temporaire indexé.
    """
    # Lire tous les reads (paires header + seq)
    reads = []
    with open(input_file, "r") as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
            else:
                if header is not None:
                    reads.append((header, line))
                    header = None

    # Mélanger avec la seed
    local_random = random.Random(seed)
    local_random.shuffle(reads)

    # Écrire
    with open(output_file, "w") as f:
        for header, seq in reads:
            f.write(header + "\n")
            f.write(seq + "\n")


def ensure_unique_names(input_file: str, output_file: str, verbose: bool = False) -> int:
    """
    Vérifie que les noms de séquences dans un fichier FASTA sont uniques.
    Si des doublons sont trouvés, ajoute _unqX à la fin du nom (X = compteur).

    Retourne le nombre de doublons corrigés.
    """
    # Lire tous les reads
    reads = []
    with open(input_file, "r") as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
            else:
                if header is not None:
                    reads.append([header, line])  # Liste pour pouvoir modifier
                    header = None

    # Compter les occurrences de chaque nom
    name_counts: dict[str, int] = {}
    for read in reads:
        name = read[0]  # Header avec >
        name_counts[name] = name_counts.get(name, 0) + 1

    # Identifier les noms en doublon
    duplicates = {name for name, count in name_counts.items() if count > 1}

    if not duplicates:
        # Pas de doublons, copier le fichier tel quel
        subprocess.run(["cp", input_file, output_file], check=True)
        return 0

    # Corriger les doublons
    n_duplicates_fixed = 0
    seen_names: dict[str, int] = {}

    for read in reads:
        name = read[0]
        if name in duplicates:
            if name not in seen_names:
                # Première occurrence, on garde le nom tel quel
                seen_names[name] = 1
            else:
                # Doublon, on ajoute _unqX
                counter = seen_names[name]
                # Trouver un nom unique
                base_name = name  # >nom_original
                new_name = f"{base_name}_unq{counter}"
                while new_name in seen_names or new_name in name_counts:
                    counter += 1
                    new_name = f"{base_name}_unq{counter}"
                read[0] = new_name
                seen_names[name] = counter + 1
                seen_names[new_name] = 1  # Marquer le nouveau nom comme vu
                n_duplicates_fixed += 1
        else:
            seen_names[name] = 1

    # Écrire le fichier corrigé
    with open(output_file, "w") as f:
        for header, seq in reads:
            f.write(header + "\n")
            f.write(seq + "\n")

    if verbose and n_duplicates_fixed > 0:
        print(
            f"  {n_duplicates_fixed} duplicate sequence names corrected (added _unqX suffix)",
            file=sys.stderr,
        )

    return n_duplicates_fixed

def main():

    # Get parameters from Snakemake environment
    output_file = snakemake.output.features
    input_files = list(snakemake.input.fastqs)
    input_files.sort()
    metadata_path = snakemake.input.metadata
    paired = snakemake.params.paired
    n_features = int(snakemake.wildcards.n_features)
    seed = int(snakemake.wildcards.random_seed)
    chunk_length = snakemake.params.chunk_length
    balance_cols = snakemake.params.balance_cols
    threads = snakemake.threads
    shuffle = snakemake.params.shuffle
    check_unique = snakemake.params.check_unique
    verbose = snakemake.params.verbose

    # Create parent directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    
    # Vérifier que les fichiers existent
    for f in input_files:
        if not os.path.exists(f):
            print(f"Error: file {f} does not exist.", file=sys.stderr)
            sys.exit(1)

    if verbose:
        print(f"Seed: {seed}", file=sys.stderr)

    try:
        metadata = pd.read_csv(metadata_path, index_col=0)
    except:
        print(f"Error: File {metadata_path} does not exist or it is empty.", file=sys.stderr)
        sys.exit(1)

    if set(balance_cols) - set(metadata.columns):
        print(f"Error: Column(s) {set(balance_cols) - set(metadata.columns)} not in the metadata.", file=sys.stderr)
        sys.exit(1)

    if paired:
        samples = [get_sample_name_paired(f1, f2) for f1, f2 in zip(input_files[::2], input_files[1::2])]
    else:
        samples = [get_sample_name(f) for f in input_files]

    if set(samples) - set(metadata.index):
        print(f"Error: Sample(s) {set(samples) - set(metadata.index)} not in the metadata.", file=sys.stderr)
        sys.exit(1)
    metadata = metadata.loc[samples, :]

    # 1. Calculer combien de reads prendre dans chaque fichier

    random.seed(seed)

    if balance_cols:
        samples_per_condition = metadata.groupby(balance_cols).size().to_frame("n_samples")
        n_balance_cols = len(samples_per_condition)
        features_per_condition = pd.DataFrame(index=samples_per_condition.index)
        features_per_condition["nfeats"] = math.floor(n_features/n_balance_cols)
        remainder = n_features % n_balance_cols
        features_per_condition.loc[random.sample(list(features_per_condition.index), remainder),"nfeats"] += 1

        features_per_sample = pd.DataFrame(index=samples)
        features_per_sample["nfeats"] = pd.NA

        for conds in samples_per_condition.index:
            if len(balance_cols)==1:
                condition_samples = list(metadata.index[metadata[balance_cols[0]]==conds])
            else:
                condition_samples = list(metadata.index[metadata.apply(lambda row: tuple(row[balance_cols])==conds, axis=1)])
            features_per_sample.loc[condition_samples,"nfeats"] = math.floor(features_per_condition.at[conds, "nfeats"]/samples_per_condition.at[conds, "n_samples"])
            remainder = features_per_condition.at[conds, "nfeats"] % samples_per_condition.at[conds, "n_samples"]
            features_per_sample.loc[random.sample(condition_samples, remainder),"nfeats"] += 1
    else:
        features_per_sample = pd.DataFrame({
            "nfeats": math.floor(n_features/len(samples))
        }, index=samples)
        remainder = n_features % len(samples)
        features_per_sample.loc[random.sample(samples, remainder),"nfeats"] += 1
        
    features_per_sample = list(features_per_sample.loc[samples, "nfeats"])
    features_per_file = [x for s in features_per_sample for x in (math.ceil(s/2), math.floor(s/2))] if paired else features_per_sample

    if verbose:
        print(f"Repartition of reads by file:", file=sys.stderr)
        for f, n in zip(input_files, features_per_file):
            print(f"  {f}: {n} features", file=sys.stderr)

    # 2. Traiter les fichiers en parallèle
    temp_dir = tempfile.mkdtemp(prefix="fastqs_random_reads_")
    temp_files = []

    try:
        # Préparer les arguments pour chaque fichier
        process_args = [
            (f, n, chunk_length, seed, temp_dir)
            for f, n in zip(input_files, features_per_file)
        ]

        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Soumettre tous les jobs et garder l'ordre
            futures = [executor.submit(process_one_file, arg) for arg in process_args]

            # Attendre les résultats dans l'ordre original
            for i, future in enumerate(futures):
                filename = input_files[i]
                try:
                    fname, temp_path, n_extracted = future.result()
                    temp_files.append(temp_path)
                    if verbose:
                        print(
                            f"  {fname}: {n_extracted} features extracted -> {temp_path}",
                            file=sys.stderr,
                        )
                except Exception as e:
                    print(f"Error for {filename}: {e}", file=sys.stderr)
                    sys.exit(1)

        # 3. Merger les fichiers temporaires
        merged_file = os.path.join(temp_dir, "merged.fasta")

        # Utiliser cat pour merger (très efficace)
        with open(merged_file, "w") as outf:
            for tf in temp_files:
                subprocess.run(["cat", tf], stdout=outf, check=True)

        if verbose:
            print(f"Files merged into {merged_file}", file=sys.stderr)

        # 4. Shuffle si demandé
        if shuffle:
            if verbose:
                print("Shuffling features...", file=sys.stderr)
            shuffled_file = os.path.join(temp_dir, "shuffled.fasta")
            shuffle_fasta(merged_file, shuffled_file, seed)
            current_file = shuffled_file
        else:
            current_file = merged_file

        # 5. Vérifier et corriger les noms en doublon si demandé
        if check_unique:
            if verbose:
                print("Ensuring that the sequence names are unique...", file=sys.stderr)
            unique_file = os.path.join(temp_dir, "unique.fasta")
            ensure_unique_names(current_file, unique_file, verbose=verbose)
            final_file = unique_file
        else:
            final_file = current_file

        # 6. Écrire le fichier de sortie (compressé si .gz)
        if is_gzipped(output_file):
            with open(final_file, "rb") as inf:
                with gzip.open(output_file, "wb") as outf:
                    # Copier par blocs pour l'efficacité mémoire
                    while True:
                        chunk = inf.read(1024 * 1024)  # 1 MB
                        if not chunk:
                            break
                        outf.write(chunk)
        else:
            # Simple copie
            os.makedirs(os.path.basename(output_file), exist_ok=True)
            subprocess.run(["cp", final_file, output_file], check=True)

        if verbose:
            print(f"Output file: {output_file}", file=sys.stderr)

    
    finally:
        # Nettoyer les fichiers temporaires
        for tf in temp_files:
            if os.path.exists(tf):
                os.remove(tf)
        if os.path.exists(temp_dir):
            # Supprimer tous les fichiers restants
            for f in os.listdir(temp_dir):
                os.remove(os.path.join(temp_dir, f))
            os.rmdir(temp_dir)

    if verbose:
        print("Finished.", file=sys.stderr)


    # Copy to features.fasta
    features_dir = os.path.dirname(output_file)
    shutil.copy(output_file, os.path.join(features_dir, "features.fasta"))

if __name__ == "__main__":
    main()