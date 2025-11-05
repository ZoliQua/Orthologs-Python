
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins in a given bottle.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

# Import libraries
import requests  # python -m pip install requests
import pandas as pd
# Import local functions
from stringDB_functions import *
# Import local variables
from stringDB_variables import *

#######################
# SET FILE PARAMETERS #
#######################

isTest = True
num_cycles = 10
num_request_per_cycle = 10
num_bottles = 5
dir_export = "export/"
dir_log = "logs/"
str_goid = "go-0051726"
log_filename1 = dir_log + "pvalues_" + str_goid + "_general_" + current_time_abbrev + ".tsv"
log_filename2 = dir_log + "pvalues_" + str_goid + "_detailed_" + current_time_abbrev + ".tsv"


# START LOGGING
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename1, level=logging.DEBUG)

# Current selection for test
# taxid = taxon_list[0]

##################
# Design Bottles #
##################
bottles = {1.0: [], 0.8: [], 0.6: [], 0.4: [], 0.2: []}
if num_bottles == 5:
	bottles_list = (0.0, 0.2, 0.4, 0.6, 0.8, 1)
elif num_bottles == 10:
	bottles_list = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#################
# Iterate TaxID #
#################

for taxid in taxon_list:

	# Read file then return the count of read lines
	counter = ReadUniprotConvert(taxid)

	# Print & Log count of read lines
	print(f"FileReader have {counter} lines load from {taxid} list/file.")
	logging.info(f"FileReader have {counter} lines load from {taxid} list/file.")

	########################################################
	# Open GO file in Pandas as a DataFrame for this TaxID #
	########################################################
	selected_cols = ["Group ID", "Average H/M", "Total H/M", "Hit Proteins", "Total Proteins", "Hit Species", "Total Species"]

	# Adding this taxid col_name
	selected_cols.append(taxon_dict_go[taxid])

	go = pd.read_csv("data/" + str_goid + "-ordered.tsv", sep="\t", usecols=selected_cols)

	# Print & Log pd read
	print(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")
	logging.info(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")

	# Logging
	log_calls = []

	# Create P-value array for the whole taxid
	p_val_array = {}

	bottle_counter = 0
	for num in bottles_list:

		if bottle_counter == (len(bottles_list)-1):
			continue

		nr1 = bottles_list[bottle_counter]
		nr2 = bottles_list[bottle_counter + 1]
		go2 = go["Average H/M"].between(nr1, nr2)
		this_taxon_column_name = taxon_dict_go[taxid]

		bottle_name = f"{nr1}-{nr2}"
		p_val_array[bottle_name] = []

		export_filename = str_goid + "_pvalues_" + taxid + "_" + str(num_cycles * num_request_per_cycle) + "_bottle-" +\
				   str(nr1) + "-" + str(nr2) + current_time_abbrev + ".tsv"

		bottle_counter += 1

		# Print & Log Bottle info
		print(f"Bottle between {nr1}-{nr2} started.")
		logging.info(f"Bottle between {nr1}-{nr2} started.")

		# Retrieving a protein list from the filtered pandas dataset of this species
		list_of_bottle_proteins = []
		for prot in go[go2][this_taxon_column_name]:
			prot = str(prot)

			# When there is no hit protein for this species: SKIP
			if prot == "nan":
				continue
			# When there are more than one hit protein, randomly selection one
			if (prot.find(',') != -1):
				this_prots = prot.split(",")
				selected = random.choice(this_prots)
			else:
				selected = prot

			# Adding proteins to a list
			list_of_bottle_proteins.append(selected)

		if len(list_of_bottle_proteins) < 10:
			# Print & Log warning
			print(f"Bottle between {nr1}-{nr2} has less than 10 protiens {len(list_of_bottle_proteins)}.")
			logging.warning(f"Bottle between {nr1}-{nr2} has less than 10 protiens {len(list_of_bottle_proteins)}.")
			continue

		p_values_allcycles = []

		for k in range(1, (num_cycles+1)):
			responses_array = []
			p_values = []

			# Print & Log
			print(f"STRING p-value retriever has started cycle {k} of bottle {nr1}-{nr2}.")
			logging.info(f"STRING p-value retriever has started cycle {k} of bottle {nr1}-{nr2}.")

			for j in range(1, (num_cycles+1)):

				list_of_random_stringids = []
				this_line = [j, j]
				for i in range(1, 11):
					this_random_selection = random.choice(list_of_bottle_proteins)
					try:
						this_string_id = uniprot_2_stringid[this_random_selection]
					except KeyError:
						logging.warning(f"This {this_random_selection} is not in uniprot-id array")
						i -= 1
						continue

					if this_random_selection in uniprot_2_protname:
						this_protein_name = uniprot_2_protname[this_random_selection]
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
					"identifiers": "%0d".join(list_of_random_stringids),  # list of my selected proteins in STRING correct IDs
					"species": int(taxid),  # species NCBI identifier, ex. human 9606
					"caller_identity": "zdultester"  # my random app name
				}

				## Calling STRING

				if isTest == True:
					pvalue = "{:.4f}".format(random.uniform(0.0, 1.0))

				else:
					response = requests.post(request_url, data=params)

					for line in response.text.strip().split("\n"):
						pvalue = line.split("\t")[5]

				this_line[1] = pvalue
				p_values.append(float(pvalue))
				p_values_allcycles.append(float(pvalue))

				log_calls.append([taxid, nr1, nr2, k, j, pvalue, request_url,
								  ",".join(list_of_random_stringids)]) # , response.text.strip()

				# Print P-value to the console (inactive)
				# print("P-value:", pvalue)

				responses_array.append(this_line)

			pvs = pd.Series(p_values, index=range(1, 11))
			responses_array.append(["SUM", "", "Min (p-val)", pvs.min(), "Max (p-val)", pvs.max(), "Mean (p-val)", pvs.mean()])
			# Print & Log
			print(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
			logging.info(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

			this_export_path = dir_export + taxid + "/" + export_filename
			if isTest == False:
				counter = WriteExportFile(this_export_path, responses_array)

		pvs = pd.Series(p_values_allcycles, index=range(1, 101))
		p_val_array[bottle_name] = p_values_allcycles

		# Print & Log
		print(f"Parser summary for all cycle in bottle between {nr1}-{nr2}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
		logging.info(f"Parser summary all cycle in bottle between {nr1}-{nr2}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

	# Writing detailed log file
	log_counter = WriteExportFile(log_filename2, log_calls)
	print(f"Parser has written {log_counter}: lines in {log_filename2}.")
