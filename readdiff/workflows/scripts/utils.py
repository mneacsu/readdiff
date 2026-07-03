import os
import re

def slugify(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")

def ensure_list(x):
    return x if isinstance(x, list) else [x,]

def get_size_mb(paths, path_list=None):
    """
    Calcule la taille totale en mégaoctets d'un fichier ou dossier,
    ou d'une liste de fichiers/dossiers.

    :param paths: str ou iterable de str
    :return: taille en mégaoctets (float)
    :raises TypeError: si paths n'est ni une str ni un iterable de str
    """
    def calculate_size(path):
        total_size = 0
        if os.path.isfile(path):
            total_size = os.path.getsize(path)
        else:
            for dirpath, dirnames, filenames in os.walk(path):
                for f in filenames:
                    fp = os.path.join(dirpath, f)
                    if not os.path.islink(fp):  # exclut les liens symboliques
                        total_size += os.path.getsize(fp)
        return total_size

    # Cas d'une seule chaîne de caractères
    if isinstance(paths, str):
        total_bytes = calculate_size(paths)

    # Cas d'un itérable (autre qu'une str)
    else:
        try:
            iterator = iter(paths)
        except TypeError:
            raise TypeError(
                "L'argument doit être une chaîne de caractères "
                "ou un itérable de chaînes de caractères."
            )
        total_bytes = 0
        for p in iterator:
            if not isinstance(p, str):
                raise TypeError("Chaque élément de l'itérable doit être une chaîne de caractères.")
            total_bytes += calculate_size(p)

    # Conversion en mégaoctets
    return total_bytes / (1024 * 1024)