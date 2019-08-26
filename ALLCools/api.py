from ._allc_to_bigwig import allc_to_bigwig
from ._allc_to_region_count import allc_to_region_count
from ._bam_to_allc import bam_to_allc
from ._extract_allc import extract_allc
from ._merge_allc import merge_allc_files
from ._open import open_allc
from .count_matrix.mcds import generate_mcds
from .motif.cmotif import generate_cmotif_database
from .motif.motif_scan import allc_motif_scan
from .motif.ame import ame
from .utilities import profile_allc, tabix_allc, standardize_allc

_ = allc_to_bigwig.__doc__
_ = allc_to_region_count.__doc__
_ = bam_to_allc.__doc__
_ = extract_allc.__doc__
_ = merge_allc_files.__doc__
_ = open_allc.__doc__
_ = profile_allc.__doc__
_ = tabix_allc.__doc__
_ = standardize_allc.__doc__
_ = generate_mcds.__doc__
_ = generate_cmotif_database.__doc__
_ = allc_motif_scan.__doc__
_ = ame.__doc__
