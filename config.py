# Orthologs Pipeline — Central Configuration
#
# Single source of truth for taxon definitions, API config, and shared constants.
# All scripts should import from this module instead of defining their own taxon dicts.
#
# Code written by Zoltan Dul, PhD

from datetime import datetime

# Timestamp for output filenames
now = datetime.now()
TIMESTAMP = now.strftime("%Y%m%d-%H%M%S-%f")

###############################################################################
# Taxon Configuration
###############################################################################

# 10 canonical target species (TaxID → species name)
TAXON_DICT = {
    '9606':  'H. sapiens',
    '7955':  'D. rerio',
    '6239':  'C. elegans',
    '3702':  'A. thaliana',
    '7227':  'D. melanogaster',
    '4896':  'S. pombe',
    '4932':  'S. cerevisiae',
    '10090': 'M. musculus',
    '10116': 'R. norvegicus',
    '8364':  'X. tropicalis',
}

# Ordered canonical TaxID lists
TAXON_ORDER_10 = ('9606', '7955', '6239', '3702', '7227', '4896', '4932', '10090', '10116', '8364')
TAXON_ORDER_7 = ('9606', '7955', '6239', '3702', '7227', '4896', '4932')

# Alias TaxIDs → canonical TaxID (strain-level IDs mapped to species-level)
TAXON_ALIASES = {
    '284812': '4896',   # S. pombe strain 972h-
    '559292': '4932',   # S. cerevisiae S288c
}

# Extended dict including aliases (for lookup when data uses alias TaxIDs)
TAXON_DICT_WITH_ALIASES = {
    **TAXON_DICT,
    '284812': 'S. pombe',
    '559292': 'S. cerevisiae',
}

# All TaxIDs that should be recognized (canonical + aliases)
TARGET_TAXIDS = set(TAXON_DICT.keys()) | set(TAXON_ALIASES.keys())

# Taxon GO column names (for pandas DataFrames in STRING analysis)
TAXON_DICT_GO = {tid: f'{name} Hit' for tid, name in TAXON_DICT.items()}

# Extended species name dict including subspecies aliases (for QuickGO/UniProt data)
TAXON_DICT_NAMES = {
    '9606':   'H. sapiens',
    '63221':  'H. sapiens',       # H. sapiens neanderthalensis
    '741158': 'H. sapiens',       # H. sapiens ssp. Denisova
    '7955':   'D. rerio',
    '6239':   'C. elegans',
    '3702':   'A. thaliana',
    '7227':   'D. melanogaster',
    '4896':   'S. pombe',
    '284812': 'S. pombe',
    '4932':   'S. cerevisiae',
    '559292': 'S. cerevisiae',
    '10090':  'M. musculus',
    '10116':  'R. norvegicus',
    '8364':   'X. tropicalis',
}

###############################################################################
# STRING API Configuration
###############################################################################

STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "tsv-no-header"
STRING_METHOD_PPI = "ppi_enrichment"
STRING_CALLER_IDENTITY = "zdul_orthologs"

def get_string_request_url(method=STRING_METHOD_PPI):
    return "/".join([STRING_API_URL, STRING_OUTPUT_FORMAT, method])

###############################################################################
# GO SLIM Generic Terms
###############################################################################

GO_SLIM_GENERIC = {
    'GO:1901135': 'carbohydrate derivative metabolic process',
    'GO:0140053': 'mitochondrial gene expression',
    'GO:0140014': 'mitotic nuclear division',
    'GO:0140013': 'meiotic nuclear division',
    'GO:0098754': 'detoxification',
    'GO:0098542': 'defense response to other organism',
    'GO:0072659': 'protein localization to plasma membrane',
    'GO:0071941': 'nitrogen cycle metabolic process',
    'GO:0071554': 'cell wall organization or biogenesis',
    'GO:0065003': 'protein-containing complex assembly',
    'GO:0061024': 'membrane organization',
    'GO:0061007': 'hepaticobiliary system process',
    'GO:0055086': 'nucleobase-containing small molecule metabolic process',
    'GO:0055085': 'transmembrane transport',
    'GO:0055065': 'metal ion homeostasis',
    'GO:0051604': 'protein maturation',
    'GO:0050886': 'endocrine process',
    'GO:0050877': 'nervous system process',
    'GO:0048870': 'cell motility',
    'GO:0048856': 'anatomical structure development',
    'GO:0044782': 'cilium organization',
    'GO:0042254': 'ribosome biogenesis',
    'GO:0042060': 'wound healing',
    'GO:0036211': 'protein modification process',
    'GO:0034330': 'cell junction organization',
    'GO:0032200': 'telomere organization',
    'GO:0031047': 'gene silencing by RNA',
    'GO:0030198': 'extracellular matrix organization',
    'GO:0030163': 'protein catabolic process',
    'GO:0030154': 'cell differentiation',
    'GO:0023052': 'signaling',
    'GO:0022600': 'digestive system process',
    'GO:0022414': 'reproductive process',
    'GO:0016192': 'vesicle-mediated transport',
    'GO:0016073': 'snRNA metabolic process',
    'GO:0016071': 'mRNA metabolic process',
    'GO:0015979': 'photosynthesis',
    'GO:0012501': 'programmed cell death',
    'GO:0007568': 'aging',
    'GO:0007163': 'establishment or maintenance of cell polarity',
    'GO:0007155': 'cell adhesion',
    'GO:0007059': 'chromosome segregation',
    'GO:0007040': 'lysosome organization',
    'GO:0007031': 'peroxisome organization',
    'GO:0007018': 'microtubule-based movement',
    'GO:0007010': 'cytoskeleton organization',
    'GO:0007005': 'mitochondrion organization',
    'GO:0006954': 'inflammatory response',
    'GO:0006914': 'autophagy',
    'GO:0006913': 'nucleocytoplasmic transport',
    'GO:0006886': 'intracellular protein transport',
    'GO:0006790': 'sulfur compound metabolic process',
    'GO:0006766': 'vitamin metabolic process',
    'GO:0006629': 'lipid metabolic process',
    'GO:0006575': 'cellular modified amino acid metabolic process',
    'GO:0006520': 'cellular amino acid metabolic process',
    'GO:0006486': 'protein glycosylation',
    'GO:0006457': 'protein folding',
    'GO:0006399': 'tRNA metabolic process',
    'GO:0006355': 'regulation of transcription, DNA-templated',
    'GO:0006351': 'transcription, DNA-templated',
    'GO:0006325': 'chromatin organization',
    'GO:0006310': 'DNA recombination',
    'GO:0006281': 'DNA repair',
    'GO:0006260': 'DNA replication',
    'GO:0006091': 'generation of precursor metabolites and energy',
    'GO:0005975': 'carbohydrate metabolic process',
    'GO:0003016': 'respiratory system process',
    'GO:0003014': 'renal system process',
    'GO:0003013': 'circulatory system process',
    'GO:0003012': 'muscle system process',
    'GO:0002376': 'immune system process',
    'GO:0002181': 'cytoplasmic translation',
    'GO:0000910': 'cytokinesis',
    'GO:0000278': 'mitotic cell cycle',
}

GO_SELECTED_10 = {
    'GO:0007049': 'cell cycle',
    'GO:0000902': 'cell morphogenesis',
    'GO:0006259': 'DNA metabolic process',
    'GO:0008361': 'regulation of cell size',
    'GO:0051726': 'regulation of cell cycle',
    'GO:0051301': 'cell division',
    'GO:0006412': 'translation',
    'GO:0006099': 'tricarboxylic acid cycle',
    'GO:0000502': 'proteasome complex',
    'GO:0009295': 'nucleoid',
}
