
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins, randomly selecting the pool.
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

isTest = False
num_cycles = 20
num_request_per_cycle = 10
dir_export = "export/"
dir_log = "logs/"
str_goid = "go-0051726"
log_filename1 = dir_log + "pvalues_" + str_goid + "_general_" + current_time_abbrev + ".tsv"
log_filename2 = dir_log + "pvalues_" + str_goid + "_detailed_" + current_time_abbrev + ".tsv"

# START LOGGING
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename1, level=logging.DEBUG)

# Current selection for test
# taxid = taxon_list[0]

# Number of all request
num_sum_of_request = num_cycles * num_request_per_cycle

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
	col_hm_name = "Average H/M"

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

	# Create array for H/M sums
	hm_sum_array = {}
	hm_sum = 0

	go_sampled = go # go.sample(n=100)
	this_taxon_column_name = taxon_dict_go[taxid]

	bottle_name = f"bottle one"
	p_val_array[bottle_name] = []

	export_filename = str_goid + "_pvalues_" + taxid + "_" + str(num_cycles * num_request_per_cycle) + "_bottle-rnd" +\
			    "-" + current_time_abbrev + ".tsv"

	# Print & Log Bottle info
	print(f"Bottle Random started.")
	logging.info(f"Bottle Random started.")

	# Retrieving a protein list from the filtered pandas dataset of this species
	list_of_bottle_proteins = []
	protein_hm_array = {}
	for index, row_data in go_sampled.iterrows():

		prot = str(row_data[this_taxon_column_name])
		hm_val = row_data[col_hm_name]

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

		# Adding poritnes to H/M array
		protein_hm_array[selected] = hm_val

	if len(list_of_bottle_proteins) < 10:
		# Print & Log warning
		print(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		logging.warning(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		continue

	p_values_allcycles = []

	for k in range(1, (num_cycles+1)):
		responses_array = []
		p_values = []

		# Print & Log
		print(f"STRING p-value retriever has started cycle {k} of Bottle Random.")
		logging.info(f"STRING p-value retriever has started cycle {k} of Bottle Random.")

		for j in range(1, (num_cycles+1)):

			list_of_random_stringids = []
			hm_sum_of_random_stringids = 0
			this_line = [j, j, j]
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
				hm_sum_of_random_stringids += protein_hm_array[this_random_selection]
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

			hm_sum = "{:.4f}".format(hm_sum_of_random_stringids)

			this_line[1] = pvalue
			this_line[2] = hm_sum
			p_values.append(float(pvalue))
			p_values_allcycles.append(float(pvalue))

			log_calls.append([taxid, hm_sum, pvalue])

			# log_calls.append([taxid, hm_sum_of_random_stringids, k, j, pvalue, request_url,
			#				  ",".join(list_of_random_stringids)]) # , response.text.strip()

			# Print P-value to the console (inactive)
			# print("P-value:", pvalue)

			responses_array.append(this_line)

		pvs = pd.Series(p_values, index=range(1, num_cycles+1))
		responses_array.append(["SUM", "", "Min (p-val)", pvs.min(), "Max (p-val)", pvs.max(), "Mean (p-val)", pvs.mean()])
		# Print & Log
		print(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
		logging.info(f"Parser summary for cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

		this_export_path = dir_export + taxid + "/" + export_filename
		if isTest == False:
			counter = WriteExportFile(this_export_path, responses_array)

	pvs = pd.Series(p_values_allcycles, index=range(1, 401))
	p_val_array[bottle_name] = p_values_allcycles

	# Print & Log
	print(f"Parser summary for all cycle in Bottle Random: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
	logging.info(f"Parser summary all cycle in Bottle Random: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

	# Writing detailed log file
	log_counter = WriteExportFile(log_filename2, log_calls)
	print(f"Parser has written {log_counter}: lines in {log_filename2}.")
