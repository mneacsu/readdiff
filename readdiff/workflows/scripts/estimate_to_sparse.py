#!/usr/bin/env python3
import os
import shutil
import pandas as pd
import anndata
import time
import numpy as np
from tqdm import tqdm
from pandas.errors import EmptyDataError
import re
from scipy import sparse

def main():

    estimate = snakemake.input.estimate
    metadata = snakemake.input.metadata
    sparse_tmp1 = snakemake.params.sparse_tmp1
    sparse_tmp2 = snakemake.params.sparse_tmp2
    sparse_final = snakemake.output.sparse
    verbose = snakemake.params.verbose

    # Convert to .h5ad

    start_time = time.time()

    df = pd.io.parsers.read_csv(estimate, sep="\t", header=None)
    obs_names = df.iloc[0, 1:].values.astype(str)
    var_names = df.iloc[1:, 0].values.astype(str)
    numeric_data = df.iloc[1:, 1:].values.astype(np.float32)
    if verbose:
        print(f"Numeric data shape before transpose: {numeric_data.shape}")
    X_sparse = sparse.csr_matrix(numeric_data)
    X_sparse = X_sparse.transpose()
    adata = anndata.AnnData(X=X_sparse)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.write(sparse_tmp1)

    if verbose:
        print(f"Number of rows in TSV file: {df.shape[0]}")
        print(f"Number of columns in TSV file: {df.shape[1]}")
        print(f"Number of observations: {adata.shape[0]}")
        print(f"Number of variables: {adata.shape[1]}")
        print("Head of the dense matrix (first 5x5):")
        print(adata.X[:5, :5].toarray())
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.2f} seconds")

    # Add sample metadata

    # Lecture de l'objet AnnData
    if verbose:
        print(f"Lecture du fichier AnnData : {sparse_tmp1}")

    try:
        adata = anndata.read_h5ad(sparse_tmp1)
    except Exception as e:
        raise ValueError(f"Erreur lors de la lecture du fichier h5ad : {e}")

    # Lecture du CSV des métadonnées
    if verbose:
        print(f"Lecture du fichier CSV de métadonnées : {metadata}")

    try:
        df_obs = pd.read_csv(metadata, index_col=0)
    except EmptyDataError:
        # Si le fichier est vide (p. ex. 0 octet), on crée un DataFrame vide
        if verbose:
            print(
                "Le fichier CSV est vide, seule la colonne 'index' sera ajoutée à l'objet AnnData."
            )
        df_obs = pd.DataFrame()
    except Exception as e:
        raise ValueError(f"Erreur lors de la lecture du fichier CSV : {e}")

    # Ajout de la colonne 'index' et gestion des métadonnées
    if not df_obs.empty:
        # Vérification de la présence de toutes les observations de l'AnnData dans le CSV
        missing_obs = set(adata.obs_names) - set(df_obs.index)
        if missing_obs:
            raise ValueError(
                "Les observations suivantes de l'AnnData ne sont pas présentes dans "
                f"l'index du CSV :\n{missing_obs}"
            )
        
        df_obs.columns = [re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_") for c in df_obs.columns] 

        # Sauvegarde des anciens noms de colonnes pour identifier ensuite lesquelles ont été ajoutées
        old_columns = set(adata.obs.columns)

        # Jointure des métadonnées (left join pour conserver l'ordre des obs existantes)
        adata.obs = adata.obs.join(df_obs, how="left")

        # Affichage des nouveaux champs ajoutés (toujours affiché)
        new_columns = set(adata.obs.columns) - old_columns
        if new_columns:
            print("Colonnes de métadonnées ajoutées :")
            for col in sorted(new_columns):
                print(f"  - {col}")
            print(f"Total : {len(new_columns)} colonnes ajoutées")
        else:
            print("Aucune nouvelle colonne de métadonnées n'a été ajoutée.")
    else:
        # Si le CSV est vide, ajouter uniquement la colonne 'index'
        print("Le fichier CSV est vide")

    # Sauvegarde de l'objet AnnData modifié
    if verbose:
        print(f"Sauvegarde du nouvel objet AnnData dans : {sparse_tmp2}")
    try:
        adata.write_h5ad(sparse_tmp2)
    except Exception as e:
        raise ValueError(f"Erreur lors de l'écriture du fichier h5ad : {e}")

    # Remove temporary file
    if os.path.exists(sparse_tmp1):
        os.remove(sparse_tmp1)

    # Move to final output
    shutil.move(sparse_tmp2, sparse_final)

    # Copy to alternate name
    estimates_dir = os.path.dirname(sparse_final)
    extension = sparse_final.split("_")[-1]

    if extension.startswith("raw") or extension.startswith("transformed"):
        shutil.copy(sparse_final, os.path.join(estimates_dir, f"estimates_{extension}"))

if __name__ == "__main__":
    main()
