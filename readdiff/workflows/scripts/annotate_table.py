#!/usr/bin/env python3
import os
import shutil
import subprocess
import pandas as pd
import re
from tqdm import tqdm
import re

def get_length(cigar, use_padded=False):
    """
    Calcule la longueur de l'alignement sur la référence à partir de la chaîne CIGAR.
    La longueur est déterminée par la somme des opérations qui consomment la référence :
      - M : match ou mismatch
      - D : deletion
      - N : région sautée (par exemple, des introns)
      - = : match exact
      - X : mismatch
    Si use_padded est True, la longueur des insertions (I) est également prise en compte.
    (Par défaut, ces opérations ne consomment pas la référence.)

    Paramètres:
      cigar (str): chaîne CIGAR (ex: "10M100N10M")
      use_padded (bool): si True, les insertions (I) sont comptabilisées. Par défaut, False.

    Retourne:
      int: la longueur de l'alignement sur la référence.
    """
    alignment_length = 0
    # Expression régulière pour extraire les groupes (nombre, opération) du CIGAR
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    for length_str, op in pattern.findall(cigar):
        num = int(length_str)
        if op in ("M", "D", "N", "=", "X"):
            alignment_length += num
        elif op == "I" and use_padded:
            alignment_length += num
        # Les opérations S, H, P (et autres) n'affectent pas la longueur de l'alignement sur la référence.
    return alignment_length

def get_readtype(gtf_content, cigar, TE=False):
    match_length = get_length(cigar)
    gtf_lines = gtf_content.split(",")

    def type_and_match(gtf_line):
        columns = [col.strip() for col in gtf_line.split("|")]
        return columns[1], int(columns[-1])

    types_and_matches = [type_and_match(line) for line in gtf_lines]
    UTRs = [m for t, m in types_and_matches if t == "UTR"]
    exons = [m for t, m in types_and_matches if t == "exon"]
    genes = [m for t, m in types_and_matches if t == "gene"]
    transcripts = [m for t, m in types_and_matches if t == "transcript"]

    # NO JUCTION ------------------
    if not "N" in cigar:
        if len(genes) + len(transcripts) == 0:
            return "intergenic" if not TE else "intergenic TE"
        if not (
            any([m >= match_length for m in genes])
            or any([m >= match_length for m in transcripts])
        ):
            return "ambiguous"
        if len(exons) == 0:
            return "intronic" if not TE else "intronic TE"
        if not any([m >= match_length for m in exons]):
            return "exon-intron"
        if len(UTRs) == 0:
            return "CDS"
        if not any([m >= match_length for m in UTRs]):
            return "UTR-CDS"
        else:
            return "UTR"

    # 1 JUNCTION ------------------
    elif cigar.count("N") == 1:
        if len(genes) + len(transcripts) == 0:
            return "intergenic junction"
        if not (
            any([m >= match_length for m in genes])
            or any([m >= match_length for m in transcripts])
        ):
            return "multigenic junction"
        if len(exons) >= 2:
            return "exon-exon junction"
        else:
            return "novel intragenic junction"

    # MORE THAN 2 JUNCTION ------------------
    else:
        return "multijunction"


def get_cigartype(cigar):
    cigar_operations = re.findall(r"[A-Z]", cigar)

    if "N" in cigar:
        return "junction"
    if "I" in cigar or "D" in cigar:
        return "indels"
    if len([M for M in cigar_operations if M == "M"]) == 1:
        return "match"
    return "other"

def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """
    Parse une chaîne CIGAR en liste de (longueur, opération).
    Ex: "3M1I2M1D2M" -> [(3,'M'), (1,'I'), (2,'M'), (1,'D'), (2,'M')]
    """
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def parse_md(md: str) -> list[str | int]:
    """
    Parse la chaîne MD en liste de tokens typés.
    Chaque token est soit :
      - un int    → nombre de matches consécutifs
      - une str   → base de référence (mismatch, 1 caractère)
      - une str   → délétion (commence par '^', ex: '^ACG')

    Ex: "2A3^TG2" -> [2, 'A', 3, '^TG', 2]
    """
    tokens: list[str | int] = []
    i = 0
    while i < len(md):
        # Nombre de matches
        m = re.match(r"\d+", md[i:])
        if m:
            tokens.append(int(m.group()))
            i += len(m.group())
            continue
        # Délétion : ^XYZ
        if md[i] == "^":
            m = re.match(r"\^[A-Za-z]+", md[i:])
            if m:
                tokens.append(m.group())          # ex: '^ACG'
                i += len(m.group())
                continue
        # Mismatch : lettre isolée
        if md[i].isalpha():
            tokens.append(md[i])
            i += 1
            continue
        raise ValueError(f"Token MD inattendu à la position {i} : '{md[i:]}'")
    return tokens


# ---------------------------------------------------------------------------
# Reconstruction
# ---------------------------------------------------------------------------

def reconstruct_reference(
    read_seq: str,
    cigar: str,
    md: str,
) -> str:
    """
    Reconstruit la séquence de référence locale alignée au read.

    Paramètres
    ----------
    read_seq : séquence du read (sans soft-clip si déjà tronquée,
               sinon la fonction gère le S dans le CIGAR)
    cigar    : chaîne CIGAR (ex: "5M2I3M1D4M")
    md       : valeur du tag MD, sans le préfixe "MD:Z:"
               (ex: "2A3^TG2")

    Retourne
    --------
    La séquence de référence locale (str), correspondant uniquement
    aux bases alignées (opérations M/=/X/D).
    """
    cigar_ops = parse_cigar(cigar)
    md_tokens = parse_md(md)

    ref_bases: list[str] = []
    read_pos = 0          # curseur dans read_seq
    md_idx = 0            # index dans md_tokens
    md_consumed = 0       # bases consommées dans le token numérique courant

    def current_md_token() -> Optional[str | int]:
        return md_tokens[md_idx] if md_idx < len(md_tokens) else None

    def advance_md():
        nonlocal md_idx, md_consumed
        md_idx += 1
        md_consumed = 0

    for length, op in cigar_ops:

        # ------------------------------------------------------------------
        # M / = / X  →  alignement (match ou mismatch)
        # ------------------------------------------------------------------
        if op in ("M", "=", "X"):
            for _ in range(length):
                tok = current_md_token()

                if tok is None:
                    raise ValueError(
                        f"MD trop court pour le CIGAR à read_pos={read_pos}"
                    )

                # Token numérique : match
                if isinstance(tok, int):
                    if md_consumed < tok:
                        ref_bases.append(read_seq[read_pos])
                        md_consumed += 1
                        if md_consumed == tok:
                            advance_md()
                    else:
                        advance_md()
                        # réessayer sur le token suivant
                        tok = current_md_token()
                        if isinstance(tok, str) and not tok.startswith("^"):
                            ref_bases.append(tok.upper())
                            advance_md()
                        else:
                            raise ValueError(
                                f"Token MD inattendu après un nombre : {tok}"
                            )

                # Token lettre : mismatch (la référence avait cette base)
                elif isinstance(tok, str) and not tok.startswith("^"):
                    ref_bases.append(tok.upper())
                    advance_md()

                else:
                    raise ValueError(
                        f"Token MD inattendu pendant M/=/X : {tok}"
                    )
                read_pos += 1

        # ------------------------------------------------------------------
        # I  →  insertion dans le read (absent de la référence)
        # ------------------------------------------------------------------
        elif op == "I":
            read_pos += length          # on avance dans le read, rien dans ref

        # ------------------------------------------------------------------
        # D  →  délétion dans le read (présent dans la référence)
        # ------------------------------------------------------------------
        elif op == "D":
            tok = current_md_token()
            if tok is None or not isinstance(tok, str) or not tok.startswith("^"):
                raise ValueError(
                    f"Attendu un token de délétion '^...' dans MD, obtenu : {tok}"
                )
            deleted_bases = tok[1:]     # retire le '^'
            if len(deleted_bases) != length:
                raise ValueError(
                    f"Longueur délétion CIGAR ({length}) ≠ MD ({len(deleted_bases)})"
                )
            ref_bases.extend(list(deleted_bases.upper()))
            advance_md()
            # pas de mouvement dans read_pos

        # ------------------------------------------------------------------
        # S  →  soft-clip : bases présentes dans le read, non alignées
        # ------------------------------------------------------------------
        elif op == "S":
            read_pos += length          # on saute ces bases du read

        # ------------------------------------------------------------------
        # H / P / N  →  ignorés ici
        # ------------------------------------------------------------------
        elif op in ("H", "P", "N"):
            pass

        else:
            raise ValueError(f"Opération CIGAR non gérée : {op}")

    return "".join(ref_bases)


def parse_cigar2(cigar):
    """
    Parse a CIGAR string and return totals for each operation.
    """
    # Match one or more digits followed by one letter (operation)
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    totals = {"N": 0, "I": 0, "D": 0, "M": 0, "H": 0, "S": 0}
    for num, op in pattern.findall(cigar):
        num = int(num)
        if op in totals:
            totals[op] += num
    return totals


def get_alignment_signature(cigar, nmismatch, md):
    """
    Classify alignment type based on CIGAR, NM, and MD.
    """
    ops = set(re.findall(r"[A-Z]", cigar))
    special = ops - {"M"}

    # Categorize special CIGAR operations
    categories = set()
    if "D" in special:
        categories.add("deletion")
    if "I" in special:
        categories.add("insertion")
    if "N" in special:
        categories.add("spliced")
    if "S" in special or "H" in special:
        categories.add("clipped")

    # Multiple event types = complex
    if len(categories) > 1:
        return "complex"

    # Single event type
    if len(categories) == 1:
        return categories.pop()

    # Only M operations: use NM to distinguish perfect/SNV/MNV
    if nmismatch is None or nmismatch == 0:
        return "perfect"
    if nmismatch == 1:
        return "SNV"
    return "MNV"


def parse_gene_annotation(gene_ann):
    """
    Extract gene_type and gene_name from a gene annotation string.
    """
    if gene_ann == "NoGeneAnnotation":
        return "NoGeneTypeAnnotation", "NoGeneNameAnnotation"
    gt = re.search(r'gene_type "([^"]+)"', gene_ann)
    gn = re.search(r'gene_name "([^"]+)"', gene_ann)
    gene_type = gt.group(1) if gt else "NoGeneTypeAnnotation"
    gene_name = gn.group(1) if gn else "NoGeneNameAnnotation"
    return gene_type, gene_name


def parse_TX(tx):
    """
    Parse taxon annotation from a TX:Z: tag.
    """
    if tx == "NoAnnotation":
        return "NoTaxonAnnotation"
    prefix = "TX:Z:"
    if tx.startswith(prefix):
        return tx[len(prefix) :]
    return "NoTaxonAnnotation"


def parse_SE(se):
    """
    Parse LCA mapping annotation from an SE:Z: tag.
    """
    if se == "NoAnnotation":
        return "NotFound"
    prefix = "SE:Z:"
    if se.startswith(prefix):
        return se[len(prefix) :]
    return "NotFound"


def parse_homer_annotation(homer_ann):
    """
    Parse HOMER annotation string (semicolon-delimited).
    Return dict of {homer_column: value}.
    """
    homer_columns = [
        "homer_Annotation",
        "homer_DetailedAnnotation",
        "homer_DistanceToTSS",
        "homer_NearestPromoterID",
        "homer_EntrezID",
        "homer_NearestUnigene",
        "homer_NearestRefseq",
        "homer_NearestEnsembl",
        "homer_GeneName",
        "homer_GeneAlias",
        "homer_GeneDescription",
        "homer_GeneType",
    ]

    if homer_ann == "NoHomerAnnotation":
        base = {col: "NoHomerAnnotation" for col in homer_columns}
        base.update(
            {
                "homer_TE_name": "NoTEAnnotation",
                "homer_TE_family": "NoTEAnnotation",
                "homer_TE_class": "NoTEAnnotation",
            }
        )
        return base

    # Split the annotation string (semicolon-delimited in previous step)
    fields = homer_ann.split(";")
    while len(fields) < len(homer_columns):
        fields.append("")

    ann_dict = {col: val.strip() for col, val in zip(homer_columns, fields)}

    # --- Post-processing ---
    # 1. homer_Annotation: keep only first token before space
    annotation_val = ann_dict["homer_Annotation"]
    if annotation_val and annotation_val != "NoHomerAnnotation":
        ann_dict["homer_Annotation"] = annotation_val.split(" ")[0]

    # 2. homer_DetailedAnnotation: split into TE_name, TE_family, TE_class if "|"
    detailed_val = ann_dict["homer_DetailedAnnotation"]
    if "|" in detailed_val:
        parts = [p.strip() for p in detailed_val.split("|")]
        while len(parts) < 3:
            parts.append("NoTEAnnotation")
        ann_dict["homer_TE_name"] = parts[0]
        ann_dict["homer_TE_class"] = parts[1]
        ann_dict["homer_TE_family"] = parts[2]
    else:
        ann_dict["homer_TE_name"] = "NoTEAnnotation"
        ann_dict["homer_TE_family"] = "NoTEAnnotation"
        ann_dict["homer_TE_class"] = "NoTEAnnotation"

    return ann_dict

def main():
    sam = snakemake.input.sam
    table = snakemake.output.table
    sortbyname = snakemake.params.sortbyname
    verbose = snakemake.params.verbose

    # If verbose, count total reads (non-header lines) using wc -l with grep -v '^@'
    total_reads = None
    nb_duplicates = 0
    if verbose:
        try:
            cmd = "grep -v '^@' {} | wc -l".format(sam)
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                universal_newlines=True,
            )
            total_reads = int(result.stdout.strip())
            print("Total reads (non-header lines):", total_reads)
        except Exception as e:
            print("Error counting reads with wc/grep:", e)

    # Initialize the BestName class to parse taxon names.
    # best_name = BestName(
    #    taxon_names,
    #    ["genbank common name", "scientific name", "authority"],
    #    custom_match={0: "unknown"},
    # )

    # Initialize a dictionary to collect the processed data.
    # We will use QNAME as the index.
    data = {}
    # Open the SAM file
    with open(sam, "r") as samfile:
        # If total_reads is known, use tqdm progress bar.
        if total_reads is not None:
            iterator = tqdm(samfile, total=total_reads, desc="Processing reads")
        else:
            iterator = samfile

        for line in iterator:
            if line.startswith("@"):
                continue
            line = line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) < 11:
                continue

            qname = fields[0]
            try:
                flag = int(fields[1])
            except ValueError:
                continue

            is_secondary = (flag & 0x100) != 0
            is_unmapped = (flag & 0x4) != 0

            if is_secondary and not is_unmapped:
                continue

            mapped = 0 if is_unmapped else 1

            rname = fields[2]
            try:
                pos = int(fields[3])
                mapq = int(fields[4])
            except ValueError:
                pos, mapq = None, None
            cigar = fields[5]
            seq = fields[9]
            seqlen = len(seq)

            strand = -1 if (flag & 0x10) != 0 else 1

            # Unmapped reads carry different SAM tags than mapped reads, so
            # parse tags by prefix rather than by index. Non-tag-formatted
            # fields are trailing annotation columns (gene, repeat, HOMER).
            tags = {}
            annotations = []
            for f in fields[11:]:
                if len(f) >= 6 and f[2] == ":" and f[4] == ":" and f[3] in "AifZHB":
                    tags[f[:2]] = f[5:]
                else:
                    annotations.append(f)

            def _opt_int(key):
                v = tags.get(key)
                if v is None:
                    return None
                try:
                    return int(v)
                except ValueError:
                    return None

            nmap = _opt_int("NH")
            nmismatch = _opt_int("NM")
            md = tags.get("MD")
            alignment_score = _opt_int("AS")
            nMismatch = _opt_int("nM")

            cigar_totals = parse_cigar2(cigar)
            N_tot = cigar_totals["N"]
            I_tot = cigar_totals["I"]
            D_tot = cigar_totals["D"]
            M_tot = cigar_totals["M"]
            H_tot = cigar_totals["H"]
            S_tot = cigar_totals["S"]

            # taxon_raw = "TX:Z:" + tags["TX"] if "TX" in tags else "NoAnnotation"
            LCA_raw = "SE:Z:" + tags["SE"] if "SE" in tags else "NoAnnotation"
            #taxon_id = parse_TX(taxon_raw)
            #taxon_name = best_name.get_best_name(taxon_id)
            LCA_mapping = parse_SE(LCA_raw)

            while len(annotations) < 2:
                annotations.append("NA")
            gene_annotation = annotations[0]
            homer_annotation = annotations[1]

            if is_unmapped:
                gene_type = gene_name = "NA"
                homer_dict = parse_homer_annotation("NoHomerAnnotation")
                for k in homer_dict:
                    homer_dict[k] = "NA"
                readtype_noTE = readtype = cigartype = "unmapped"
                alignment_signature = "unmapped"
            else:
                gene_type, gene_name = parse_gene_annotation(gene_annotation)
                homer_dict = parse_homer_annotation(homer_annotation)
                readtype_noTE = get_readtype(gene_annotation, cigar, TE=False)
                readtype = get_readtype(
                    gene_annotation,
                    cigar,
                    TE=False,
                )
                cigartype = get_cigartype(cigar)
                alignment_signature = get_alignment_signature(cigar, nmismatch, md)

            # Reconstruct reference sequence from seq, cigar, and MD tag
            if cigar != "*" and md is not None:
                try:
                    seq_reference = reconstruct_reference(seq, cigar, md)
                except Exception:
                    seq_reference = None
            else:
                seq_reference = None

            # Build the row dictionary.
            row = {
                "mapped": mapped,
                "nmap": nmap,
                "chromosome_alt": rname,
                "chromosome": rname.split("_")[0],
                "mapq": mapq,
                "strand": strand,
                "pos": pos,
                "cigar": cigar,
                "seq": seq,
                "seqlen": seqlen,
                "Nmismatch": nmismatch,
                "md": md,
                "seq_reference": seq_reference,
                "nMismatch": nMismatch,
                "aligment_score": alignment_score,
                "N_tot": N_tot,
                "I_tot": I_tot,
                "D_tot": D_tot,
                "M_tot": M_tot,
                "H_tot": H_tot,
                "S_tot": S_tot,
                #"repeat_class": repeat_class,
                #"repeat_family": repeat_family,
                "gene_type": gene_type,
                "gene_name": gene_name,
                #"taxon_id": taxon_id,
                #"taxon_name": taxon_name,
                "LCA_mapping": LCA_mapping,
                "readtype_noTE": readtype_noTE,
                "readtype": readtype,
                "cigartype": cigartype,
                "alignment_signature": alignment_signature,
                "weight": 1,
            }
            # add HOMER fields
            row.update(homer_dict)

            # Check if QNAME already exists; we will check for duplicates later.
            if qname in data:
                # We add the read again; duplicate checking is done at the end.
                # data[qname + "_DUPLICATE_{}".format(len(data))] = row
                nb_duplicates += 1
            else:
                data[qname] = row

    # Construct a pandas DataFrame. Use QNAME as the index.
    df = pd.DataFrame.from_dict(data, orient="index")

    # Check for duplicate QNAMEs in the index (ignoring the artificially created duplicate keys)
    original_qnames = [q for q in data.keys() if "_DUPLICATE_" not in q]
    if len(original_qnames) != len(set(original_qnames)):
        raise ValueError("Duplicate QNAME found in the SAM file.")

    # If --sortbyname flag is set, sort the DataFrame by index (QNAME)
    if sortbyname:
        df = df.sort_index()

    # Write the dataframe to the output TSV file.
    df.to_csv(table, sep="\t", index=True)

    if verbose:
        total_reads_read = len([1 for key in data])
        unique_reads = len(set(original_qnames))
        print("Number of reads processed (lines added to dataframe):", total_reads_read)
        print("Number of unique reads (QNAME):", unique_reads)
        print("Total number of non-header reads in file:", total_reads)
        print("Number of duplicates found:", nb_duplicates)

    # Copy to downstream directory
    table_parent = os.path.join(os.path.dirname(table), "analysis_table.tsv")
    shutil.copy(table, table_parent)

if __name__ == "__main__":
    main()
    