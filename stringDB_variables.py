
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# Containing the variables that more than one script use
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

from config import (
    TIMESTAMP, TAXON_ORDER_7, TAXON_DICT, TAXON_DICT_GO,
)

# Backward-compatible aliases
current_time_abbrev = TIMESTAMP
taxon_list = TAXON_ORDER_7
taxon_dict = {tid: TAXON_DICT[tid] for tid in TAXON_ORDER_7}
taxon_dict_go = {tid: TAXON_DICT_GO[tid] for tid in TAXON_ORDER_7}
