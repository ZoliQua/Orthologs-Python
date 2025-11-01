
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# x
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com

import os
import csv
import sys
import random
import requests
import pandas as pd

from config import TAXON_ORDER_7, TAXON_DICT

csv.field_size_limit(sys.maxsize)

taxon_dict = {tid: TAXON_DICT[tid] for tid in TAXON_ORDER_7}

for taxid in taxon_dict:

    filename = "data/uniprot_convert_" + taxid + ".tsv"
    write_lines_all = {}
    convert_uniprot = {}
    name_uniprot = {}
    list_of_uniprots = []

    with open(filename, newline='') as f:
        reader = csv.DictReader(f, fieldnames= ('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
        counter = 0
        try:
            for row in reader:

                # Filtering STRING to convert Uniprot to STRING db id
                if row['db'] == 'convert':
                    if row['uniprot'] in convert_uniprot:
                        continue
                    else:
                        convert_uniprot[row['uniprot']] = row['selector']
                        list_of_uniprots.append(row['uniprot'])
                        counter += 1

                # Filtering STRING convert out
                if row['db'] == 'Gene_Name':
                    if row['uniprot'] in name_uniprot:
                        continue
                    else:
                        name_uniprot[row['uniprot']] = row['selector']

        except csv.Error as e:
            sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

        print(f"FileReader have {counter} lines load from {taxid} list/file.")

    for k in range(1, 5):
        responses_array = []
        p_values = []
        print(f"STRING p-value retriever has started cycle {k}.")

        for j in range(1,11):

            list_of_random_stringids = []
            this_line = [j, j]
            for i in range(1,6):
                this_random_selection = random.choice(list_of_uniprots)
                this_string_id = convert_uniprot[this_random_selection]

                if this_random_selection in name_uniprot:
                    this_protein_name = name_uniprot[this_random_selection]
                else:
                    this_protein_name = "n/a"

                list_of_random_stringids.append(taxid + "." + this_string_id)
                this_line.append(this_protein_name)
                this_line.append(this_string_id)


            string_api_url = "https://string-db.org/api"
            output_format = "tsv-no-header"
            method = "ppi_enrichment"

            request_url = "/".join([string_api_url, output_format, method])

            params = {
                "identifiers" : "%0d".join(list_of_random_stringids),   # list of my selected proteins in STRING correct IDs
                "species" : int(taxid),                                 # species NCBI identifier, ex. human 9606
                "caller_identity" : "zdultester"                          # my random app name
            }

            ## Calling STRING
            response = requests.post(request_url, data=params)

            ## Parse and print the respons Parse and print the responsee
            # print("The following IDs have been used: ", ','.join(map(str, list_of_random_stringids)) )

            for line in response.text.strip().split("\n"):
                pvalue = line.split("\t")[5]
                this_line[1] = pvalue
                p_values.append(float(pvalue))

                # Print P-value to the console (inactive)
                # print("P-value:", pvalue)

            responses_array.append(this_line)

        pvs = pd.Series(p_values, index=range(1,11))
        responses_array.append(["SUM", "", "Min (p-val)", pvs.min(), "Max (p-val)", pvs.max(), "Mean (p-val)", pvs.mean()])
        print(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {pvs.min()}.")

        export_filename = "data/pvalues_random_uniprotids_" + taxid + ".tsv"
        counter = 0

        with open(export_filename, mode='a') as export_file:
            writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting = csv.QUOTE_MINIMAL)

            for line in responses_array:

                counter += 1
                writer.writerow(line)

            # print(f"Parser have {counter} lines wrote in {export_filename} (merged) file.")



