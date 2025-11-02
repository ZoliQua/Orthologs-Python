
# EGGNOG GROUPPER v1.0
#
# What this file do?
# This file converts downloaded id_mapping files (from data/uniprot folder) to a reduced file size. Filtering out STRING and eggNOG db related lines.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import csv
import sys
import operator

from config import TAXON_DICT_WITH_ALIASES

csv.field_size_limit(sys.maxsize)

# list_of_go_tags = ["0006412", "0006629", "0006914", "0007049", "0007165", "0007568", "0008361", "0042254", "0051301", "0051726"]
# list_of_go_tags = ["0006259", "0006397", "0006399", "0006457"]

# list_of_go_tags = ["0008361", "0002376", "0009295", "0000902", "0006099", "0003013", "0000502", "0006399", "0000910", "0005975"]

list_of_go_tags = ["0008361", "0002376", "0009295", "0000902", "0006099", "0003013", "0000502", "0006399", "0000910", "0005975", "0006629", "0006914", "0007568", "0006259", "0051726", "0051301", "0006397", "0007049", "0006412"]

for go_tag in list_of_go_tags:
    filename_go = "data/go/GO-" + go_tag + ".tsv"
    taxon_dicter = TAXON_DICT_WITH_ALIASES

    go_listofproteins = {}
    go_subgo = {}
    with open(filename_go, newline='') as f:
        reader = csv.DictReader(f, fieldnames= ('type', 'uniprot', 'name', 'subgo', 'longname', 'taxid'), delimiter='\t')
        counter = 0
        try:
            for row in reader:
                if row['type'] != 'UniProtKB':
                    continue
                go_listofproteins[row['uniprot']] = taxon_dicter[row['taxid']]
                go_subgo[row['uniprot']] = row['subgo']
                counter += 1
                # if counter == 100:
                #     break
        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(filename_go, reader.line_num, e))

    filename = "data/uniprot_convert_merged.tsv"
    taxon_list = ['3702', '6239', '7227', '7955', '9606', '559292', '284812']
    taxon_dict = {'3702': [], '6239': [], '7227': [], '7955': [], '9606': [], '559292': [], '284812': []}
    eggnog_taxlist = ["9606", "7955", "6239", "3702", "7227", "4896", "4932"]

    eggNOG_database = {}

    with open(filename, newline='') as f:
        reader = csv.DictReader(f, fieldnames= ('uniprot', 'db', 'taxid', 'groupid'), delimiter='\t')
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

                # if counter == 10000:
                #      break

        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

        print("Parser have", counter, "lines processed from eggNOG merged file.")

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
                    if uniprot in go_listofproteins:
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
                average_hm_list.append(float(this_hit_prot_count[taxid]/this_total_prot_count[taxid]))

        average_hm_total = 0
        if len(average_hm_list) == 0:
            average_hm = 0
        else:
            for value in average_hm_list:
                average_hm_total += value
            average_hm = float( average_hm_total / total_spec_count)

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

        # adding detailed information
        for taxid in taxon_list:
            this_line += str(this_hit_spec_count[taxid]) + "\t"
            this_line += str(this_hit_prot_count[taxid]) + "\t"
            this_line += str(this_total_spec_count[taxid]) + "\t"
            this_line += str(this_total_prot_count[taxid]) + "\t"

            if this_total_prot_count[taxid] > 0:
                fraction = this_hit_prot_count[taxid] / this_total_prot_count[taxid]
            else:
                fraction = 0

            this_line += str(float("{:.5f}".format(fraction))) + "\t"

                # for sor_taxid in taxon_list:
        #     this_line += ",".join(this_mezok[sor_taxid]) + "\t"
        # for sor_taxid in taxon_list:
        #     sor_taxid2 = "H" + sor_taxid + "H"
        #     this_line += ",".join(this_mezok[sor_taxid2]) + "\t"

        counter += 1
        write_the_output.append(this_line)
        if hit_spec_count > 0:
            if total_spec_count > 3:
                write_the_output_hit.append(this_line)
                write_the_output_hit_dict_4items[counter] = this_line
                write_the_output_hit_dict_4order[counter] = int(average_hm*10000000)

    print(f"Parser have {counter} elements processed from GO-{go_tag}.")
    # print(write_the_output_hit)

    sorted_dicdata = sorted(write_the_output_hit_dict_4order.items(), key=operator.itemgetter(1), reverse=True)

    export_filename = "data/go-" + go_tag + "-ordered-detailed.tsv"
    counter = 0

    with open(export_filename, mode='w') as export_file:
        writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        this_line = "Group ID\t"
        this_line += "Average H/M\t"
        this_line += "Total H/M\t"
        this_line += "Hit Proteins\t"
        this_line += "Total Proteins\t"
        this_line += "Hit Species\t"
        this_line += "Total Species\t"

        for sor_taxid in taxon_list:
            this_line += taxon_dicter[sor_taxid] + " Total\t"
        for sor_taxid in taxon_list:
            this_line += taxon_dicter[sor_taxid] + " Hit\t"

        this_line = this_line.split("\t")
        writer.writerow(this_line)

        for rowid in sorted_dicdata:
            counter += 1
            this_line = write_the_output_hit_dict_4items[rowid[0]].split("\t")
            writer.writerow(this_line)

        # for line in write_the_output_hit:
        #
        #     counter += 1
        #     this_line = line.split("\t")
        #     writer.writerow(this_line)

        print(f"Parser have {counter} lines wrote in {export_filename} (merged) file.")
