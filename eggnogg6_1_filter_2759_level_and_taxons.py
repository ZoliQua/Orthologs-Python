# eggNOG v6 — STEP 1: Filter e6.og2seqs_and_species.tsv
#
# Szűrés: csak a 2759 (Eukaryota) szintű sorok, és csak a 10 célspecies.
#
# v1: 7 species (hs, dr, ce, at, dm, sp, sc)
# v2: 10 species — hozzáadva: M. musculus (10090), R. norvegicus (10116), X. tropicalis (8364)
#
# A fájl oszlopai (tab-elválasztott):
#   [0] level_taxid   (pl. "2759" = Eukaryota)
#   [1] og_name
#   [2] num_species
#   [3] num_members
#   [4] species_list  (vesszővel elválasztott TaxID-k)
#   [5] members_list  (vesszővel elválasztott "TaxID.ProteinID" stringek)

import sys
import time

from config import TAXON_DICT_WITH_ALIASES, TAXON_ALIASES, TARGET_TAXIDS

source_eggnog_file = "data/eggnog/eggnog6/e6.og2seqs_and_species.tsv"
export_eggnog_file = "data/eggnog/eggnog6/e6.og2seqs_and_species_filtered_10_species_2759_level.tsv"

taxon_dict = TAXON_DICT_WITH_ALIASES
canonical_taxid = TAXON_ALIASES

# ─────────────────────────────────────────────────────────────────────────────
# Feldolgozás
# ─────────────────────────────────────────────────────────────────────────────
total_lines = 0
kept_lines = 0
skipped_wrong_level = 0
skipped_no_species = 0

start_time = time.time()

print(f"Forrás: {source_eggnog_file}")
print(f"Kimenet: {export_eggnog_file}")
print(f"Szűrés: level=2759, {len(set(taxon_dict.values()))} célspecies")
print("Feldolgozás...")

with open(source_eggnog_file, 'r') as original_file, open(export_eggnog_file, 'w') as filtered_file:
    for line in original_file:
        total_lines += 1

        if total_lines % 500_000 == 0:
            elapsed = time.time() - start_time
            print(f"  {total_lines:,} sor feldolgozva ({elapsed:.0f}s), megtartva: {kept_lines:,}")

        columns = line.strip().split('\t')

        # Csak a 2759 (Eukaryota) szintű sorok
        if len(columns) < 6 or columns[0].strip() != "2759":
            skipped_wrong_level += 1
            continue

        members = columns[5].split(',')
        species_in_row = columns[4].split(',')

        # Gyors ellenőrzés: van-e célspecies a sorban?
        if not any(s in TARGET_TAXIDS for s in species_in_row):
            skipped_no_species += 1
            continue

        # Szűrjük a member listát a célspeciesekre
        filtered_members = [m for m in members if m.split('.')[0] in TARGET_TAXIDS]

        if not filtered_members:
            skipped_no_species += 1
            continue

        # Kanonikus TaxID-ra normalizált fajlista (alias-mentes)
        filtered_species_canonical = set()
        for m in filtered_members:
            tid = m.split('.')[0]
            filtered_species_canonical.add(canonical_taxid.get(tid, tid))

        # Oszlopok frissítése
        columns[2] = str(len(filtered_species_canonical))
        columns[3] = str(len(filtered_members))
        columns[4] = ','.join(sorted(filtered_species_canonical))
        columns[5] = ','.join(filtered_members)

        filtered_file.write('\t'.join(columns) + '\n')
        kept_lines += 1

elapsed = time.time() - start_time
print(f"\n{'─'*55}")
print(f"Összes sor:             {total_lines:>12,}")
print(f"Más szint (nem 2759):   {skipped_wrong_level:>12,}")
print(f"Nincs célspecies:       {skipped_no_species:>12,}")
print(f"Megtartott OG-k:        {kept_lines:>12,}")
print(f"Futási idő:             {elapsed:>11.1f}s")
print(f"{'─'*55}")
print(f"Kimenet: {export_eggnog_file}")
