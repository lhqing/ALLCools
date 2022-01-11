import re
import numpy as np
import pandas as pd
import ALLCools
from .motifs import Motif, MotifSet

PACKAGE_DIR = ALLCools.__path__[0]
DEFAULT_MOTIF_DIR = f"{PACKAGE_DIR}/motif/default_motif_set/"

N_SITES_PATTERN = re.compile(r"(?<=nsites= )\d+")
SPACE_PATTERN = re.compile(r"[ \t]")


def parse_header_lines(lines):
    _alphabet = None
    for _line in lines:
        if _line.startswith("ALPHABET"):
            _alphabet = _line.split("=")[1].strip()
    if _alphabet is None:
        raise ValueError

    _background = None
    next_line = False
    for _line in lines:
        if _line.startswith("Background letter frequencies"):
            next_line = True
            continue
        if next_line:
            values = _line.strip().split(" ")
            _background = {
                values[0]: float(values[1]),
                values[2]: float(values[3]),
                values[4]: float(values[5]),
                values[6]: float(values[7]),
            }
            break
    return _alphabet, _background


def parse_motif(lines, alphabet, background):
    motif_name = lines[0][5:].strip()

    txt = lines[1].split(":")[1].strip()
    _m = N_SITES_PATTERN.search(txt)
    n_sites = int(_m[0]) if _m is not None else np.NaN

    _data = []
    for _line in lines[2:]:
        try:
            values = SPACE_PATTERN.split(_line.strip())
            values = [v.strip() for v in values if v.strip() != ""]
            _data.append(list(map(float, values)))
        except ValueError:
            # some additional rows that are not PFM
            continue
    _data = np.array(_data)
    counts = np.round(_data * n_sites).astype(int)
    try:
        counts = {base: list(counts[:, i]) for i, base in enumerate(alphabet)}
    except IndexError:
        print(lines)
        raise

    _motif = Motif(alphabet=alphabet, counts=counts)
    _motif.background = background
    _motif.name = motif_name
    _motif.pseudocounts = 0.5
    return _motif


def parse_meme_database(meme_path, meta_path):
    total_motifs = {}
    header = True
    header_lines = []
    motif_lines = []
    alphabet, background = None, None

    with open(meme_path) as f:
        for line in f:
            if line[:5].upper() == "MOTIF":
                alphabet, background = parse_header_lines(header_lines)
                header = False
            if header:
                header_lines.append(line)
                continue

            if line.strip() == "":
                continue

            if (not header) and line[:5].upper() == "MOTIF":
                if len(motif_lines) > 0:
                    motif = parse_motif(
                        motif_lines, alphabet=alphabet, background=background
                    )
                    if motif.name in total_motifs:
                        raise ValueError
                    else:
                        total_motifs[motif.name] = motif
                    motif_lines = [line]
                else:
                    motif_lines = [line]
            else:
                motif_lines.append(line)

        # last motif
        if len(motif_lines) > 1:
            motif = parse_motif(motif_lines, alphabet=alphabet, background=background)
            if motif.name in total_motifs:
                raise ValueError
            else:
                total_motifs[motif.name] = motif

    meta_table = pd.read_csv(meta_path, index_col=0)
    motif_set = MotifSet(list(total_motifs.values()), meta_table=meta_table)
    return motif_set


def get_default_motif_set(database="three_databases"):
    if database == "three_databases":
        motif_set = parse_meme_database(
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.meme",
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.metadata.csv",
        )

        # default thresholds, fnr/fpr = 1000
        motif_set.thresholds = pd.read_csv(
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.thresholds.csv",
            header=None,
            index_col=0,
            squeeze=True,
        ).to_dict()
        return motif_set
    else:
        # TODO: allow user create motif set by providing meme and metadata (optional)
        raise NotImplementedError
