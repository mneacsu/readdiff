#!/usr/bin/env python3


def main():

    input_files = list(snakemake.input.sub_estimates)
    output_file = snakemake.output.merged_estimate

    if not input_files:
        print("No input files provided")
        with open(output_file, "w") as f:
            f.write("query\n")
        return

    # Ouvrir tous les fichiers
    input_files = [open(f, 'r') for f in input_files]

    with open(output_file, 'w') as out_f:
        # Traiter ligne par ligne
        line_num = 0
        while True:
            lines = [f.readline() for f in input_files]

            # Si le premier fichier est terminé, on arrête
            if not lines[0]:
                break

            # Prendre toute la ligne du premier fichier (sans le \n)
            merged_line = lines[0].rstrip('\n\r')

            # Ajouter les colonnes 2+ des autres fichiers
            for line in lines[1:]:
                # Enlever le \n et split sur tab
                line_clean = line.rstrip('\n\r')
                # Trouver la première tab et prendre tout après
                tab_pos = line_clean.find('\t')
                if tab_pos != -1:
                    # Ajouter tout ce qui est après le premier tab
                    merged_line += line_clean[tab_pos:]

            out_f.write(merged_line + '\n')
            line_num += 1

    # Fermer tous les fichiers
    for f in input_files:
        f.close()

    print(f"Merged {len(input_files)} files into {output_file}")
    print(f"Total lines: {line_num}")

if __name__ == "__main__":
    main()