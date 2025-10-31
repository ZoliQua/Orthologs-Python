# Orthologs-Python

A bioinformatics pipeline for analyzing **protein-protein orthologs** across 10 model organisms. Integrates data from three major databases — **eggNOG**, **STRING**, and **UniProt/QuickGO** — to identify ortholog groups, assess protein-protein interaction enrichment, and visualize cross-species conservation patterns.

## Target Species

| TaxID | Species | Common Name |
|-------|---------|-------------|
| 9606 | *H. sapiens* | Human |
| 10090 | *M. musculus* | Mouse |
| 10116 | *R. norvegicus* | Rat |
| 7955 | *D. rerio* | Zebrafish |
| 8364 | *X. tropicalis* | Western clawed frog |
| 7227 | *D. melanogaster* | Fruit fly |
| 6239 | *C. elegans* | Nematode |
| 3702 | *A. thaliana* | Thale cress |
| 4896 | *S. pombe* | Fission yeast |
| 4932 | *S. cerevisiae* | Baker's yeast |

## Analysis Workflow

```
1. eggNOG parsing       →  Ortholog group extraction (v5 & v6)
2. UniProt ID mapping   →  Cross-database protein identifier linking
3. QuickGO annotation   →  GO term & GO_SLIM functional annotation
4. Ortholog grouping    →  Proteins grouped by GO terms per species
5. STRING enrichment    →  PPI p-value calculation (known vs. random pools)
6. Statistical summary  →  Correlation analysis & Excel export
7. Visualization        →  Venn diagrams for species-level overlap
```

## Repository Structure

```
orthologs-python/
│
├── eggnogg5_1_parser_by_species.py         # eggNOG v5: parse 2759_members.tsv → 10-species groups
├── eggnogg6_1_filter_2759_level_and_taxons.py  # eggNOG v6 step 1: filter Eukaryota level
├── eggnogg6_2_analyzed_filtered_eggnog6.py     # eggNOG v6 step 2: per-species protein lists
│
├── stringDB_variables.py                   # Shared config: taxon lists, API URLs
├── stringDB_functions.py                   # Utility functions: UniProt conversion, API calls
├── stringDB_test-functions.py              # Function testing harness
├── stringDB_test-API-network-image.py      # Download PNG network images from STRING API
├── stringDB_networkx_create.py             # NetworkX protein interaction graph analysis
├── stringDB_p-value_basic.py               # Basic PPI enrichment p-value from STRING
├── stringDB_p-value_random_parser.py       # Random protein pool p-values
├── stringDB_p-value_go-list_bottle_parser.py           # GO-annotated pool p-values
├── stringDB_p-value_go-list_random-bottle_parser.py    # Random pool for specified GO terms
├── stringDB_p-value_go-list_half-plus-random-bottle_parser.py  # Hybrid: half orthologs + half random
├── stringDB_data_summarizer.py             # Aggregate p-value results
│
├── uniprotDB_files_parser.py               # Parse UniProt ID mapping files → filtered TSV
├── ortholog_uniprot_name_retriever.py      # Query UniProt REST API for protein metadata
├── ortholog_go_groupper.py                 # Group proteins by GO terms (basic)
├── ortholog_go_groupper_quickGO-7-species-query.py     # 7-species QuickGO query
├── ortholog_go_groupper_quickGO-full-query.py          # Full 10-species QuickGO query
├── ortholog_go_plot_data_summairizer.py    # Summarize GO data → Excel
│
├── analysis_summarizer.py                  # Summarize eggNOG/STRING results
├── analysis_correlation.py                 # P-value correlation analysis → Excel
├── qucikgo_1_export_parser_tester.py       # QuickGO export parsing test
│
├── quickgo-parser/                         # QuickGO API integration module
│   ├── quickGO_functions_container.py      #   Shared utilities: logging, CSV handling
│   ├── quickGO_get_GOslim.py              #   Retrieve latest GO_SLIM subset
│   ├── quickGO_query.py                   #   Main GO term query executor
│   ├── quickGO_query_1_children_terms.py  #   Fetch child terms of GO annotations
│   ├── quickGO_query_3_one_goterm_alltaxon.py  #   Single GO term across all taxa
│   ├── quickGO_query_4_one_goterm_stat.py      #   Statistical summary for GO term
│   └── quickGO_show_annotations_of_a_GOterm.py #   Display GO term annotations
│
├── venn_diagram/                           # Venn diagram visualization
│   ├── run_venn-first.py                  #   Initial overlap analysis
│   ├── run_venn2-hit.py                   #   Ortholog "hit" overlaps
│   ├── run_venn2-total.py                 #   Total protein overlaps
│   ├── source/                            #   SVG templates (ortholog_venn_*.svg)
│   └── output/                            #   Generated SVG diagrams
│
├── goatools_slim/                          # GO tools integration
│   └── goatools_query.py                  #   GOATools-based GO term queries
│
├── discover/                               # DISCOVER statistical package (submodule)
│
├── data/                                   # Input data
│   ├── eggnog/                            #   eggNOG v5 & v6 ortholog databases
│   ├── uniprot/                           #   UniProt ID mapping files
│   ├── go/                                #   Gene Ontology annotation TSVs
│   └── string-values/                     #   STRING p-value results
│
├── STRING_protein_info/                    # STRING protein metadata (aliases & info)
├── complexes/                              # Protein complex data (GO:0140014)
├── export/                                 # Per-species output (by TaxID)
└── output/                                 # Summary statistics & Excel files
```

## Scripts Overview

### eggNOG Database Processing

Parse eggNOG ortholog databases (v5 and v6) to extract ortholog groups at the Eukaryota (2759) taxonomic level, filtered for the 10 target species. Handles TaxID alias normalization (e.g., *S. pombe* 284812 → 4896, *S. cerevisiae* 559292 → 4932).

### STRING Protein Interaction Analysis

Query the STRING database API to calculate **PPI enrichment p-values** for ortholog groups. Three comparison strategies:
- **Known orthologs**: proteins from a specific GO term
- **Random pools**: randomly selected proteins as negative control
- **Hybrid pools**: half known orthologs + half random (sensitivity test)

Uses NetworkX for graph-based network property analysis.

### UniProt & Gene Ontology Integration

- Parse UniProt ID mapping files to link protein identifiers across STRING and eggNOG
- Query the QuickGO REST API for GO term annotations using GO_SLIM subsets
- Group orthologous proteins by shared functional annotations

### Visualization

Generate Venn diagrams (SVG) showing ortholog group overlaps across species, using BeautifulSoup for SVG template manipulation.

## Dependencies

- **pandas** — Data manipulation and Excel I/O
- **requests** — HTTP API calls (STRING, UniProt, QuickGO)
- **networkx** — Protein interaction network analysis
- **matplotlib** — Data visualization
- **beautifulsoup4** — SVG/HTML parsing
- **openpyxl** — Excel file writing
- **goatools** — Gene Ontology analysis tools

## Data Sources

- [eggNOG v5/v6](http://eggnog5.embl.de/) — Orthologous groups database
- [STRING v12](https://string-db.org/) — Protein-protein interaction network
- [UniProt](https://www.uniprot.org/) — Protein sequence and annotation database
- [QuickGO](https://www.ebi.ac.uk/QuickGO/) — Gene Ontology browser (EBI)

## Author

**Zoltán Dul**
