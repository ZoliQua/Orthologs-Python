# eggNOG v6 — STEP 2: Filtered file → 10-species groups + per-species protein lists
#
# Bemenet:  e6.og2seqs_and_species_filtered_10_species_2759_level.tsv
# Kimenetek:
#   eggnog6db_10_species_by_groups.tsv   — OG-nként soronként a 10 faj fehérjelistája
#   eggnog6db_{TaxID}_protein_list.tsv   — fajankénti egyedi fehérjelista
#
# v1: 7 species
# v2: 10 species — hozzáadva: M. musculus (10090), R. norvegicus (10116), X. tropicalis (8364)
#      + bugfix: S. pombe alias (284812→4896), S. cerevisiae alias (559292→4932)
#      + oszlopfejlécek a kimenetben
#      + per-species fájlok is mentve

import csv
import sys
import time

from config import TAXON_DICT, TAXON_ORDER_10, TAXON_ALIASES

csv.field_size_limit(sys.maxsize)

filtered_eggnog_file = "data/eggnog/eggnog6/e6.og2seqs_and_species_filtered_10_species_2759_level.tsv"
export_groups_file   = "data/eggnog/eggnog6db_10_species_by_groups.tsv"
export_protein_tmpl  = "data/eggnog/eggnog6db_{taxid}_protein_list.tsv"

taxon_dict = TAXON_DICT
canonical_taxid = TAXON_ALIASES
TAXON_ORDER = list(TAXON_ORDER_10)

# ─────────────────────────────────────────────────────────────────────────────
# Feldolgozás
# ─────────────────────────────────────────────────────────────────────────────
protein_lists = {tid: set() for tid in TAXON_ORDER}
groups_output = []
total = 0
start_time = time.time()

print(f"Bemenet: {filtered_eggnog_file}")
print("Feldolgozás...")

with open(filtered_eggnog_file, 'r') as f:
    for line in f:
        total += 1
        columns = line.strip().split('\t')

        if len(columns) < 6:
            continue

        og_name      = columns[1]
        num_species  = columns[2]  # step1 által frissítve
        members_raw  = columns[5].split(',')

        # Per-species fehérjelista összeállítása
        species_proteins = {tid: [] for tid in TAXON_ORDER}

        for member in members_raw:
            parts = member.split('.', 1)
            if len(parts) != 2:
                continue
            raw_tid, protein_id = parts[0], parts[1]
            # Alias normalizálás
            tid = canonical_taxid.get(raw_tid, raw_tid)
            if tid in species_proteins:
                species_proteins[tid].append(protein_id)
                protein_lists[tid].add(protein_id)

        # Sor összeállítása: og_name, num_species, majd 10 fehérjeoszlop
        actual_species_count = sum(1 for tid in TAXON_ORDER if species_proteins[tid])
        row = [og_name, str(actual_species_count)]
        for tid in TAXON_ORDER:
            row.append(','.join(species_proteins[tid]))

        groups_output.append(row)

elapsed = time.time() - start_time
print(f"  {total:,} OG feldolgozva ({elapsed:.1f}s)")

# ─────────────────────────────────────────────────────────────────────────────
# Kiírás: OG-csoportok fájl
# ─────────────────────────────────────────────────────────────────────────────
header = ['og_name', 'num_species'] + [taxon_dict[tid] + f' ({tid})' for tid in TAXON_ORDER]

with open(export_groups_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    for row in groups_output:
        writer.writerow(row)

print(f"\nOG-csoport fájl: {export_groups_file}  ({len(groups_output):,} OG)")

# ─────────────────────────────────────────────────────────────────────────────
# Kiírás: per-species fehérjelisták
# ─────────────────────────────────────────────────────────────────────────────
print("\nPer-species fehérjelisták:")
for tid in TAXON_ORDER:
    out_path = export_protein_tmpl.format(taxid=tid)
    proteins = sorted(protein_lists[tid])
    with open(out_path, 'w') as f:
        for p in proteins:
            f.write(p + '\n')
    print(f"  [{taxon_dict[tid]:<22}] {out_path}  ({len(proteins):,} egyedi fehérje)")

print("\nKész.")
