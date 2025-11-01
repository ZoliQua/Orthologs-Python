
# EGGNOG RANDOM SELECTOR
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins, randomly selecting the pool.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

# Import libraries
import pandas as pd
# Import local functions
from stringDB_functions import *
# Import local variables
from stringDB_variables import *
from config import GO_SLIM_GENERIC, GO_SELECTED_10

#######################
# SET FILE PARAMETERS #
#######################

isTest = False
num_cycles = 40
num_request_per_cycle = 10
num_proteins = 6
dir_export = "export/"
dir_log = "logs/"
filename_ext = "-2022-01"

# list_goids = ["go-0008361", "go-0002376", "go-0009295", "go-0000902",
#               "go-0006099", "go-0003013", "go-0000502", "go-0006399",
#               "go-0000910", "go-0005975"]

# list_goids = ["go-0006629", "go-0006914", "go-0007568", "go-0006259",
#               "go-0051726", "go-0051301", "go-0006397", "go-0007049",
#               "go-0006412"]

list_of_goSLIM_generic = {
    'GO:0000910': GO_SLIM_GENERIC['GO:0000910'],
    'GO:0000278': GO_SLIM_GENERIC['GO:0000278'],
    **GO_SELECTED_10,
}

# list_of_goSLIM_generic = list_of_selected_10

list_goids = []
for go_id in list_of_goSLIM_generic.keys():
	list_goids.append(go_id.replace("GO:", "go-"))

for str_goid in list_goids:

	#str_goid = "go-0008361"
	log_filename1 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_general_" + current_time_abbrev + ".tsv"
	log_filename2 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_detailed_" + current_time_abbrev + ".tsv"

	# START LOGGING
	logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename1, level=logging.DEBUG)

	# Current selection for test
	# taxid = taxon_list[0]

	# Number of all request
	num_sum_of_request = num_cycles * num_request_per_cycle

	# bottles_list = (0.00, 0.33, 0.66, 1.00)
	bottles_list = (0.00, 1.00)

	for taxid in taxon_list:

		# Read file then return the count of read lines
		counter = ReadUniprotConvert(taxid)

		# Counter of Requests
		counter_requests = 0

		# Print & Log count of read lines
		print(f"FileReader have {counter} lines load from {taxid} list/file.")
		logging.info(f"FileReader have {counter} lines load from {taxid} list/file.")

		########################################################
		# Open GO file in Pandas as a DataFrame for this TaxID #
		########################################################
		selected_cols = [	"Group ID",
							"Average H/M",
							"Total H/M",
							"Hit Proteins",
							"Total Proteins",
							"Hit Species",
							 "Total Species" ]

		column_name_hm_value = "Average H/M"

		# Adding this taxid col_name
		selected_cols.append(taxon_dict_go[taxid])

		go = pd.read_csv("data/" + str_goid + "-ordered" + filename_ext + ".tsv", sep="\t", usecols=selected_cols)

		# Print & Log pd read
		print(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")
		logging.info(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")

		# Logging
		log_calls = []

		# Create P-value array for the whole taxid
		p_val_array = {}

		bottle_counter = 0
		for num in bottles_list:

			if bottle_counter == (len(bottles_list) - 1):
				continue

			# Create array for H/M sums
			hm_sum_array = {}
			hm_sum = 0

			# go_sampled = go # go.sample(n=100)

			nr1 = bottles_list[bottle_counter]
			nr2 = bottles_list[bottle_counter + 1]
			go_between_average_hm = go["Average H/M"].between(nr1, nr2)
			bottle_counter += 1

			this_taxon_column_name = taxon_dict_go[taxid]

			bottle_name = f"{nr1}-{nr2}"
			p_val_array[bottle_name] = []

			export_filename = str_goid + "_pvalues_" + taxid + "_" + str(num_cycles * num_request_per_cycle) + "_" \
							+ bottle_name + "_" + "-" + current_time_abbrev + ".tsv"

			# Print & Log Bottle info
			print(f"Taxid {taxid} - Bottle Between {nr1}-{nr2} started.")
			logging.info(f"Taxid {taxid} - Bottle Between {nr1}-{nr2} started.")

			p_values_allcycles = []

			for k in range(1, (num_cycles+1)):
				responses_array = []
				p_values = []

				# Print & Log
				print(f"STRING DB retriever has started cycle {k} of Bottle Between {nr1}-{nr2}.")
				logging.info(f"STRING DB retriever has started cycle {k} of Bottle Between {nr1}-{nr2}.")

				return_array = ParseGODataFrame(go[go_between_average_hm], taxon_dict_go[taxid], column_name_hm_value, num_proteins)

				if not return_array:
					continue

				protein_hm_array = return_array["protein_hm_array"]
				list_of_bottle_proteins = return_array["list_of_bottle_proteins"]

				j = 1

				while j < (num_request_per_cycle + 1):

					list_of_random_stringids = []
					hm_sum_of_random_stringids = 0
					this_line = [j, j, j]
					i = 1
					count_warning = 0

					while i < (num_proteins + 1):

						this_random_selection = random.choice(list_of_bottle_proteins)
						try:
							this_string_id = uniprot_2_stringid[this_random_selection]
						except:
							count_warning += 1
							logging.warning(f"This {this_random_selection} is not in uniprot-id array")
							if count_warning <= 100:
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

						i += 1

					string_api_url = "https://string-db.org/api"
					output_format = "tsv-no-header"
					method = "ppi_enrichment"

					request_url = "/".join([string_api_url, output_format, method])

					params = {
						# list of my selected proteins in STRING correct IDs
						"identifiers": "%0d".join(list_of_random_stringids),
						# species NCBI identifier, ex. human 9606
						"species": int(taxid),
						# my random app name
						"caller_identity": "zdul_tester"
					}

					## Calling STRING
					# If test, just generate random numbers
					if isTest == True:
						pvalue = "{:.4f}".format(random.uniform(0.0, 1.0))

					else:
						response = request_in_time(request_url, params)
						if not response:
							continue
						# response = requests.post(request_url, data=params)
						counter_requests += 1
						if (counter_requests % 100) == 0:
							time.sleep(5)

						for line in response.text.strip().split("\n"):
							length = len(line.split("\t"))
							if length < 6:
								break
							number_of_nodes = line.split("\t")[0]
							number_of_edges = line.split("\t")[1]
							average_node_degree = line.split("\t")[2]
							local_clustering_coefficient = line.split("\t")[3]
							expected_number_of_edges = line.split("\t")[4]
							pvalue = line.split("\t")[5]

						if length < 6:
							j -= 1
							print(f"Taxid {taxid} - ERROR - Bottle Between {nr1}-{nr2}: Length: {length}.")
							logging.error(f"Taxid {taxid} - ERROR - Bottle Between {nr1}-{nr2}: Length: {length}, Line {line}.")
							continue


					hm_sum = "{:.4f}".format(hm_sum_of_random_stringids)

					this_line[1] = pvalue
					this_line[2] = hm_sum
					p_values.append(float(pvalue))
					p_values_allcycles.append(float(pvalue))

					log_calls.append([taxid, nr1, nr2, hm_sum, number_of_nodes, number_of_edges, average_node_degree, local_clustering_coefficient, expected_number_of_edges, pvalue])

					# log_calls.append([taxid, hm_sum_of_random_stringids, k, j, pvalue, request_url,
					#				  ",".join(list_of_random_stringids)]) # , response.text.strip()

					# Print P-value to the console (inactive)
					# print("P-value:", pvalue)

					responses_array.append(this_line)

					j += 1

				pvs = pd.Series(p_values, index=range(1, num_request_per_cycle+1))
				responses_array.append(["SUM", "", "Min (p-val)", pvs.min(), "Max (p-val)", pvs.max(), "Mean (p-val)", pvs.mean()])
				# Print & Log
				print(f"Taxid {taxid} - Parser summary for p-value in cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
				logging.info(f"Taxid {taxid} - Parser summary for p-value in cycle {k}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

				this_export_path = dir_export + taxid + "/" + export_filename
				if isTest == False:
					counter = WriteExportFile(this_export_path, responses_array)
					responses_array = []

			if len(p_values_allcycles) >= 100:
				pvs = pd.Series(p_values_allcycles, index=range(1, ((num_cycles * num_request_per_cycle) + 1) ))
				p_val_array[bottle_name] = p_values_allcycles

				# Print & Log
				print(f"Taxid {taxid} - Parser summary for p-values in all cycle in Bottle Between {nr1}-{nr2}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")
				logging.info(f"Taxid {taxid} - Parser summary for p-values in all cycle in Bottle Between {nr1}-{nr2}: Min: {pvs.min()}, Max: {pvs.max()}, Mean: {'{:.4f}'.format(pvs.mean())}.")

				# Writing detailed log file
				log_counter = WriteExportFile(log_filename2, log_calls)
				print(f"Taxid {taxid} - Parser has written {log_counter}: lines in {log_filename2}.")

			else:
				# Print & Log
				print(f"Taxid {taxid} - Parser SKIPPED cycle in Bottle Between {nr1}-{nr2}.")
				logging.info(f"Taxid {taxid} - Parser SKIPPED cycle in Bottle Between {nr1}-{nr2}.")
