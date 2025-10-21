# eggNOG v5 — Parser: 2759_members.tsv → 10-species groups + per-species protein lists
#
# Bemenet:  data/eggnog/eggnog5/2759_members.tsv
# Kimenetek:
#   eggnog5db_10_species_by_groups.tsv   — OG-nként soronként a 10 faj fehérjelistája
#   eggnog5db_{TaxID}_protein_list.tsv   — fajankénti egyedi fehérjelista
#
# v1: 7 species (hs, dr, ce, at, dm, sp, sc)
# v2: 10 species — hozzáadva: M. musculus (10090), R. norvegicus (10116), X. tropicalis (8364)
#      + bugfix: S. pombe alias (284812→4896), S. cerevisiae alias (559292→4932)
#      + oszlopfejlécek a kimenetben
#
# A 2759_members.tsv oszlopai (tab-elválasztott):
#   [0] level_taxid   (2759 = Eukaryota)
#   [1] og_name
#   [2] num_species
#   [3] num_members
#   [4] members_list  (vesszővel elválasztott "TaxID.ProteinID" stringek)
#   [5] species_list  (vesszővel elválasztott TaxID-k)
#
# MEGJEGYZÉS: Ha a 3 új faj (mm, rn, xt) hiányzik a v5 adatból,
#   le kell tölteni: http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/

import csv
import sys
import time

csv.field_size_limit(sys.maxsize)

source_file          = "data/eggnog/eggnog5/2759_members.tsv"
export_groups_file   = "data/eggnog/eggnog5db_10_species_by_groups.tsv"
export_protein_tmpl  = "data/eggnog/eggnog5db_{taxid}_protein_list.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# 10 célspecies — kanonikus TaxID-k
# ─────────────────────────────────────────────────────────────────────────────
taxon_dict = {
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

# Alias → kanonikus TaxID
canonical_taxid = {
    '284812': '4896',   # S. pombe strain
    '559292': '4932',   # S. cerevisiae S288c
}

TAXON_ORDER = ['9606', '7955', '6239', '3702', '7227', '4896', '4932', '10090', '10116', '8364']
TARGET_TAXIDS = set(taxon_dict.keys()) | set(canonical_taxid.keys())

# ─────────────────────────────────────────────────────────────────────────────
# Feldolgozás
# ─────────────────────────────────────────────────────────────────────────────
protein_lists = {tid: set() for tid in TAXON_ORDER}
groups_output = []
total = 0
kept = 0
start_time = time.time()

print(f"Forrás: {source_file}")
print("Feldolgozás...")

with open(source_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    try:
        for columns in reader:
            total += 1

            if total % 500_000 == 0:
                elapsed = time.time() - start_time
                print(f"  {total:,} sor ({elapsed:.0f}s), megtartva: {kept:,}")

            if len(columns) < 5:
                continue

            og_name     = columns[1]
            members_raw = columns[4].split(',')
            specs_raw   = columns[5].split(',') if len(columns) > 5 else []

            # Gyors szűrés
            if not any(s in TARGET_TAXIDS for s in specs_raw):
                continue

            species_proteins = {tid: [] for tid in TAXON_ORDER}
            has_target = False

            for member in members_raw:
                parts = member.split('.', 1)
                if len(parts) != 2:
                    continue
                raw_tid, protein_id = parts[0], parts[1]
                tid = canonical_taxid.get(raw_tid, raw_tid)
                if tid in species_proteins:
                    species_proteins[tid].append(protein_id)
                    protein_lists[tid].add(protein_id)
                    has_target = True

            if not has_target:
                continue

            actual_species_count = sum(1 for tid in TAXON_ORDER if species_proteins[tid])
            row = [og_name, str(actual_species_count)]
            for tid in TAXON_ORDER:
                row.append(','.join(species_proteins[tid]))

            groups_output.append(row)
            kept += 1

    except csv.Error as e:
        sys.exit(f'CSV hiba, sor {reader.line_num}: {e}')

elapsed = time.time() - start_time
print(f"\n{'─'*55}")
print(f"Összes sor:      {total:>12,}")
print(f"Megtartott OG:   {kept:>12,}")
print(f"Futási idő:      {elapsed:>11.1f}s")
print(f"{'─'*55}")

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
    count_str = f"{len(proteins):,} fehérje"
    missing = " ← HIÁNYZIK a v5 adatból!" if len(proteins) == 0 else ""
    print(f"  [{taxon_dict[tid]:<22}] TaxID {tid}: {count_str}{missing}")

print("\nKész.")
