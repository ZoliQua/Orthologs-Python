
# EGGNOG GROUPPER v1.5
#
# What this file do?
# This file converts downloaded id_mapping files (from data/uniprot folder) to a reduced file size. Filtering out STRING and eggNOG db related lines.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import csv
import sys
import operator

from config import GO_SLIM_GENERIC, GO_SELECTED_10, TAXON_DICT_NAMES

csv.field_size_limit(sys.maxsize)
print_details = False

list_of_goSLIM_generic = GO_SLIM_GENERIC
list_of_selected_10 = GO_SELECTED_10

# Taxon specific collectors
taxon_list = ['3702', '6239', '7227', '7955', '9606', '559292', '284812']
taxon_dict = {tid: [] for tid in taxon_list}
taxon_dict_names = TAXON_DICT_NAMES

###############################
# Reading QuickGO export file #
###############################

# file for GO_SLIM generic
filename_go = "data/go/QuickGO-annotations-GOslim-generic-20220126.tsv"
# file for our selected 10
# filename_go = "data/go/QuickGO-annotations-special-10-20240227.tsv"
# list_of_goSLIM_generic = list_of_selected_10


filenames_go = [
    "data/go/QuickGO-annotations-GOslim-generic-20220126-at.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-ce.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-dm.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-dr.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-hs.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-sc.tsv",
    "data/go/QuickGO-annotations-GOslim-generic-20220126-sp.tsv"
]

filenames_go = [
    "data/go/QuickGO-annotations-GOslim-generic-20240320.tsv"
]

go_goslim = {}
for go_id in list_of_goSLIM_generic.keys():
    go_goslim[go_id] = []

go_subgo = {}
counter_for_goslim = 0

def read_go_files(filenames):
    global go_subgo
    global counter_for_goslim
    for filename_go in filenames:
        with open(filename_go, newline='') as f:
            #readerInDict = csv.DictReader(f, fieldnames=('type', 'uniprot', 'name', 'goslim', 'subgo', 'longname', 'taxid', 'gene_product_name', 'gene_product_synonyms', 'go_aspect'), delimiter='\t')
            readerInDict = csv.DictReader(f, fieldnames=('type', 'uniprot', 'name', 'qualifier', 'goslim', 'subgo', 'longname', 'eco_id', 'go_evidence', 'reference','with_from', 'taxid', 'assigned_by', 'annot_ext', 'go_aspect'), delimiter='\t')
            counter = 0

            try:
                for row in readerInDict:
                    if row['type'] != 'UniProtKB':
                        continue
                    if row['goslim'] not in list_of_goSLIM_generic.keys():
                        # print('Not in: ', row['subgo'], row['goslim'])
                        counter_for_goslim += 1
                    go_goslim.setdefault(row['goslim'], []).append(row['uniprot'])
                    go_subgo[row['uniprot']] = row['subgo']
                    counter += 1

            except csv.Error as e:
                sys.exit(f'file {filename_go}, line {readerInDict.line_num}: {e}')

        print(f"Parser have {counter} lines processed from {filename_go} file.")

    return True

read_go_files(filenames_go)

# In case our selection of GOSlim contains more origins of GOs - print out the number
if counter_for_goslim > 0:
    print(counter_for_goslim)


################################
# Reading UniProt Convert file #
################################

# Source filename for eggNOG_data (based on Uniprot database)
filename_uniprot_convert = "data/uniprot/uniprot_convert_merged.tsv"
# Declare eggNOG database output variable
eggNOG_database = {}

with open(filename_uniprot_convert, newline = '') as f:
    reader = csv.DictReader(f, fieldnames = ('uniprot', 'db', 'taxid', 'groupid'), delimiter = '\t')
    counter = 0
    try:
        for row in reader:
            counter += 1

            # Filtering eggNOG DB out
            if row['db'] == 'eggNOG':

                if row['groupid'] in eggNOG_database:
                    if row['taxid'] in eggNOG_database[row['groupid']]:
                        eggNOG_database[row['groupid']][row['taxid']].append(row['uniprot'])
                    else:
                        eggNOG_database[row['groupid']][row['taxid']] = [row['uniprot']]
                else:
                    eggNOG_database[row['groupid']] = {row['taxid']: [row['uniprot']]}

    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename_uniprot_convert, reader.line_num, e))

    print("Parser have", counter, "lines processed from UniProt's eggNOG merged file.")


################################
# Groupping #
# ################################

counter_for_slim = 0
for this_goslim in list_of_goSLIM_generic.keys():

    counter_for_slim += 1

    write_the_output = []
    write_the_output_hit = []
    write_the_output_hit_dict_4items = {}
    write_the_output_hit_dict_4order = {}

    counter = 0
    for groupid in eggNOG_database:

        this_line = ""
        this_mezok = {}
        this_hit_spec_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_hit_prot_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_total_spec_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}
        this_total_prot_count = {'3702': 0, '6239': 0, '7227': 0, '7955': 0, '9606': 0, '559292': 0, '284812': 0}

        for sor_taxid in taxon_list:
            sor_taxid2 = "H" + sor_taxid + "H"
            if sor_taxid in eggNOG_database[groupid]:
                this_mezok[sor_taxid] = []
                this_mezok[sor_taxid2] = []
                this_total_spec_count[sor_taxid] = 1
                for uniprot in eggNOG_database[groupid][sor_taxid]:

                    this_mezok[sor_taxid].append(uniprot)
                    this_total_prot_count[sor_taxid] += 1
                    if uniprot in go_goslim[this_goslim]:
                        this_hit_spec_count[sor_taxid] = 1
                        this_hit_prot_count[sor_taxid] += 1
                        this_mezok[sor_taxid2].append(uniprot)
            else:
                this_mezok[sor_taxid] = []
                this_mezok[sor_taxid2] = []

        average_hm_list = []

        hit_spec_count = 0
        hit_prot_count = 0
        total_spec_count = 0
        total_prot_count = 0
        for taxid in taxon_list:
            hit_spec_count += this_hit_spec_count[taxid]
            hit_prot_count += this_hit_prot_count[taxid]
            total_spec_count += this_total_spec_count[taxid]
            total_prot_count += this_total_prot_count[taxid]

            if this_hit_spec_count[taxid] > 0:
                average_hm_list.append(float(this_hit_prot_count[taxid] / this_total_prot_count[taxid]))

        average_hm_total = 0
        if len(average_hm_list) == 0:
            average_hm = 0
        else:
            for value in average_hm_list:
                average_hm_total += value
            average_hm = float(average_hm_total / total_spec_count)

        if hit_spec_count > 0:
            total_hm = hit_prot_count / total_prot_count
        else:
            total_hm = 0

        # Writing out the cells in the order of appearance

        this_line += groupid + "\t"
        this_line += str(float("{:.5f}".format(average_hm))) + "\t"
        this_line += str(float("{:.5f}".format(total_hm))) + "\t"
        this_line += str(hit_prot_count) + "\t"
        this_line += str(total_prot_count) + "\t"
        this_line += str(hit_spec_count) + "\t"
        this_line += str(total_spec_count) + "\t"

        if print_details:
            # Adding detailed information
            for taxid in taxon_list:
                this_line += str(this_hit_spec_count[taxid]) + "\t"
                this_line += str(this_hit_prot_count[taxid]) + "\t"
                #this_line += str(this_total_spec_count[taxid]) + "\t"
                #this_line += str(this_total_prot_count[taxid]) + "\t"

                if this_total_prot_count[taxid] > 0:
                    fraction = this_hit_prot_count[taxid] / this_total_prot_count[taxid]
                else:
                    fraction = 0

                #this_line += str(float("{:.5f}".format(fraction))) + "\t"
        else:

            for sor_taxid in taxon_list:
                this_line += ",".join(this_mezok[sor_taxid]) + "\t"
            for sor_taxid in taxon_list:
                sor_taxid2 = "H" + sor_taxid + "H"
                this_line += ",".join(this_mezok[sor_taxid2]) + "\t"

        counter += 1
        #write_the_output.append(this_line)
        if hit_spec_count > 0:
            if total_spec_count > 3:
                write_the_output_hit.append(this_line)
                write_the_output_hit_dict_4items[counter] = this_line
                write_the_output_hit_dict_4order[counter] = int(average_hm * 10000000)

    print(f"Cycle Nr {counter_for_slim} - Parser have {counter} elements processed from {this_goslim}.")

    sorted_dicdata = sorted(write_the_output_hit_dict_4order.items(), key = operator.itemgetter(1), reverse = True)

    if print_details:
        export_filename = "data/go/go-" + this_goslim.replace('GO:', '') + "-ordered-detailed-2024-02-correct.tsv"
    else:
        export_filename = "data/go/go-" + this_goslim.replace('GO:', '') + "-ordered-2024-02.tsv"

    counter = 0

    with open(export_filename, mode = 'w') as export_file:
        writer = csv.writer(export_file, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)

        this_line = "Group ID\t"
        this_line += "Average H/M\t"
        this_line += "Total H/M\t"
        this_line += "Hit Proteins\t"
        this_line += "Total Proteins\t"
        this_line += "Hit Species\t"
        this_line += "Total Species\t"

        for sor_taxid in taxon_list:
            this_line += taxon_dict_names[sor_taxid] + " Total\t"
        for sor_taxid in taxon_list:
            this_line += taxon_dict_names[sor_taxid] + " Hit\t"

        this_line = this_line.split("\t")
        writer.writerow(this_line)

        for rowid in sorted_dicdata:
            counter += 1
            this_line = write_the_output_hit_dict_4items[rowid[0]].split("\t")
            writer.writerow(this_line)

        print(f"Cycle Nr {counter_for_slim} - Parser have {counter} lines wrote in {export_filename} (merged) file.")
