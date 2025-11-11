
# Ortholog Parser / STRING DB Reader
#
# What this file do?
# This file get p-values from STRING DB for a randomly selected pool of proteins in a given bottle.
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#
# File for functions
# Import libraries
import csv
import sys
import random
import logging
import signal
import time
import requests
from requests.adapters import HTTPAdapter, Retry

from config import STRING_API_URL, STRING_OUTPUT_FORMAT, STRING_METHOD_PPI, STRING_CALLER_IDENTITY

# maximalise file-read size
csv.field_size_limit(sys.maxsize)

uniprot_2_stringid = {}
uniprot_2_protname = {}
list_of_uniprotids = []

#############
# FUNCTIONS #
#############

#####################
# Read UniProt file #
#####################


###################################################
# Read File Function ##############################
###################################################
#  reads the source tsv file for the given taxid ##
#  takes an str = taxid as an input ###############
#  writes 3 global variables ######################
#  returns a count of read lines ##################
###################################################


def ReadUniprotConvert(taxid: str) -> int:

	global uniprot_2_protname
	global uniprot_2_stringid
	global list_of_uniprotids

	filename = "data/uniprot_convert_" + taxid + ".tsv"

	with open(filename, newline='') as f:
		reader = csv.DictReader(f, fieldnames=('uniprot', 'db', 'taxid', 'selector'), delimiter='\t')
		counter = 0
		try:
			for row in reader:

				# Filtering STRING to convert Uniprot to STRING db id
				if row['db'] == 'convert':
					if row['uniprot'] in uniprot_2_stringid:
						continue
					else:
						uniprot_2_stringid[row['uniprot']] = row['selector']
						list_of_uniprotids.append(row['uniprot'])
						counter += 1

				# Filtering STRING convert out
				if row['db'] == 'Gene_Name':
					if row['uniprot'] in uniprot_2_protname:
						continue
					else:
						uniprot_2_protname[row['uniprot']] = row['selector']

		except csv.Error as e:
			sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

		return counter


#################################
# Write Export Function #########
#################################
#  creates a tsv separated file #
#  takes an array as an input ###
#  returns the count of lines ###
#################################


def WriteExportFile(export_filename: str, write_this_array: list) -> int:

	counter = 0

	with open(export_filename, mode='a') as export_file:
		writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

		for line in write_this_array:
			counter += 1
			writer.writerow(line)

	return counter


def ParseGODataFrame(go_dataframe, column_name_taxon: str, column_name_hm_value: str, num_proteins: int) -> dict | bool:
	
	# Protein list from the filtered pandas dataset of this species
	list_of_bottle_proteins = []
	protein_hm_array = {}
	
	for index, row_data in go_dataframe.iterrows():

		protein = str(row_data[column_name_taxon])
		hm_value = row_data[column_name_hm_value]

		# When there is no hit protein for this species: SKIP
		if protein == "nan":
			continue
		# When there are more than one hit protein, randomly selection one
		if protein.find(',') != -1:
			proteins_array = protein.split(",")
			protein_selected = random.choice(proteins_array)
		else:
			protein_selected = protein

		# Adding proteins to a list
		list_of_bottle_proteins.append(protein_selected)
		# Adding poritnes to H/M array
		protein_hm_array[protein_selected] = hm_value

	if len(list_of_bottle_proteins) < num_proteins:
		# Print & Log warning
		print(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		logging.warning(f"Bottle Random has less than 10 proteins {len(list_of_bottle_proteins)}.")
		return False

	else:
		return_array = {"protein_hm_array": protein_hm_array, "list_of_bottle_proteins": list_of_bottle_proteins}
		return return_array

# Custom timeout exception
class TimeoutError(Exception): pass

#  Call this function exceeds timeout
def handler(signum: int, frame) -> None:
    raise TimeoutError()

#  Function timeout decorator
def time_out(interval: int, doc: str):
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(interval)       #  Interval seconds to send SIGALRM signals to the process
                result = func(*args, **kwargs)
                signal.alarm(0)              #  After the function is executed after the specified time is executed, close the Alarm alarm clock
                return result
            except TimeoutError as e:
                #  Capture the timeout exception, what to do
                print("The function failed to run due to timeout, func:<%s>" % doc)
                logging.warning("The function failed to run due to timeout, func:<%s>" % doc)
        return wrapper
    return decorator

# Session with retry logic for STRING API calls
_retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
_session = requests.Session()
_session.mount("https://", HTTPAdapter(max_retries=_retries))


def string_api_request(protein_ids: list[str], taxid: str, method: str = STRING_METHOD_PPI) -> requests.Response | bool:
    """Make a STRING API request with retry logic and timeout."""
    request_url = "/".join([STRING_API_URL, STRING_OUTPUT_FORMAT, method])
    params = {
        "identifiers": "%0d".join(protein_ids),
        "species": int(taxid),
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    try:
        response = _session.post(request_url, data=params, timeout=10)
        response.raise_for_status()
        return response
    except requests.RequestException as e:
        logging.warning(f"STRING API request failed: {e}")
        return False


@time_out(1, "Function call")
def request_in_time(req: str, par: dict) -> requests.Response | bool:
    try:
        response = requests.post(req, data=par)
    except requests.RequestException as e:
        logging.warning(f"Request failed: {e}")
        response = False
    return response

