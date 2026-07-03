#!/usr/bin/env python3

def parse_bed_annotations(bed_file, annotation_fields_start=7):
    """
    Extracts metadata annotations from a BED file.
    Returns a dictionary mapping read names to annotation strings.
    """
    annotations = {}
    
    with open(bed_file, 'r') as bed:
        for line in bed:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < annotation_fields_start:
                continue  # Skip malformed lines
            
            read_name = fields[3]  # Read identifier
            annotation_data = fields[annotation_fields_start:]  # Extract relevant metadata
            
            annotation_string = " | ".join(annotation_data)
            
            if read_name in annotations:
                annotations[read_name].add(annotation_string)
            else:
                annotations[read_name] = {annotation_string}
    
    # Convert sets to comma-separated strings
    return {read: ", ".join(info) for read, info in annotations.items()}

def parse_homer_annotations(homer_file):
    """
    Extracts HOMER annotations from a .homer.txt file.
    Returns a dictionary mapping PeakID (read names) to annotation strings.
    """
    annotations = {}
    with open(homer_file, 'r') as hf:
        header = hf.readline().strip().split('\t')  # skip header
        for line in hf:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            peak_id = fields[0]  # PeakID = read identifier

            # You can decide how much of HOMER info to keep.
            # Here: everything except Chr/Start/End/Strand.
            annotation_data = fields[7:]  

            annotation_string = " ; ".join(annotation_data)
            
            if peak_id in annotations:
                annotations[peak_id].add(annotation_string)
            else:
                annotations[peak_id] = {annotation_string}

    return {read: ", ".join(info) for read, info in annotations.items()}


def merge_annotations(sam_file, genes_bed_file, homer_txt_file, output_sam_file):
    """
    Merges gene, repeat, and HOMER annotations into the SAM file.
    """
    gene_annotations = parse_bed_annotations(genes_bed_file)
    homer_annotations = parse_homer_annotations(homer_txt_file)

    with open(sam_file, 'r') as sf, open(output_sam_file, 'w') as of:
        for line in sf:
            if line.startswith('@'):  # Keep header lines unchanged
                of.write(line)
                continue
            
            fields = line.strip().split('\t')
            read_name = fields[0]

            # Retrieve annotations if they exist
            gene_annotation = gene_annotations.get(read_name, "NoGeneAnnotation")
            homer_annotation = homer_annotations.get(read_name, "NoHomerAnnotation")

            # Append annotations
            fields.append(gene_annotation)
            fields.append(homer_annotation)

            # Write modified line
            of.write("\t".join(fields) + "\n")

def main():
    sam = snakemake.input.sam
    genes_bed = snakemake.input.genes_bed
    homer_annot = snakemake.input.homer_annot
    annotated_sam = snakemake.output.annotated_sam

    merge_annotations(sam, genes_bed, homer_annot, annotated_sam)

if __name__ == "__main__":
    main()
