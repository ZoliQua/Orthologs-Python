
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

list_of_goSLIM_generic = {**GO_SLIM_GENERIC, **GO_SELECTED_10}

list_goids = []
for go_id in list_of_goSLIM_generic.keys():
	list_goids.append(go_id.replace("GO:", "go-"))

for str_goid in list_goids:

	#str_goid = "go-0008361"
	log_filename1 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_general_" + current_time_abbrev + ".tsv"
	log_filename2 = dir_log + "data_" + str_goid + "_" + str(num_proteins) + "together_detailed_" + current_time_abbrev + ".tsv"

	# START LOGGING
	logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename1, level=logging.DEBUG)

	# Counter of Requests
	counter_requests = 0

	########################################################
	# Open GO file in Pandas as a DataFrame for this TaxID #
	########################################################
	selected_cols = [	"Group ID",
	                     "Average H/M",
	                     "Total H/M",
	                     "Hit Proteins",
	                     "Total Proteins",
	                     "Hit Species",
	                     "Total Species",
	                     "A. thaliana Total",
	                     "C. elegans Total",
	                     "D. melanogaster Total",
	                     "D. rerio Total",
	                     "H. sapiens Total",
	                     "S. cerevisiae Total",
	                     "S. pombe Total" ]

	column_name_hm_value = "Average H/M"

	go = pd.read_csv("data/" + str_goid + "-ordered" + filename_ext + ".tsv", sep="\t", usecols=selected_cols)

	# Print & Log pd read
	print(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")
	logging.info(f"Pandas DataFrame have {len(go)} lines load from {str_goid} file.")

	# Logging
	log_calls = []

	# Create P-value array for the whole taxid
	p_val_array = {}

	num_all_groups = len(go)-1


	break
