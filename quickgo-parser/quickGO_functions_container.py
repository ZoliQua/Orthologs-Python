#
# QuickGO Functions / container file /
#
# What this file does?
# Containing functions and classes that are used by more than one script in the directory
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#

import os
import csv
import random
import time
import logging
from datetime import datetime
import urllib.request as rqs
import urllib.error
import json
import pandas as pd

# Creating timestamp for output filename
now = datetime.now()
current_time_abbrev = now.strftime("%Y%m%d-%H%M%S-%f")

# LOGGING SETTINGS
dir_log = "logs/"
log_filename = dir_log + "quickgo_" + current_time_abbrev + ".tsv"

# START LOGGING
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', filename=log_filename, level=logging.INFO)


def TimeNow(str_now, name_of_script=False, start_time=False):
	now = datetime.now()
	time_abb = now.strftime("%Y-%m-%d - %H:%M:%S (%f)")

	if name_of_script:
		if str_now == "start":
			print(f"Runtime of '{name_of_script}' is at {str_now} phase at {time_abb}")
			return True

	elif not start_time:
		print(f"Runtime is phase {str_now} at {time_abb}")
		return True

	if str_now == "end" or name_of_script == False:

		runtime = time.time() - start_time
		hours = runtime // 3600
		temp = runtime - 3600 * hours
		minutes = temp // 60
		seconds = temp - 60 * minutes

		if name_of_script:
			print(f"Runtime of '{name_of_script}' was", '%d hrs %d mins %d secs' % (hours, minutes, seconds),
				"(" + str(float("{:.5f}".format(runtime))) + ")")
		else:
			print(f"Script running is phase {str_now}, since start", '%d hrs %d mins %d secs elapsed' % (hours, minutes, seconds),
				"(" + str(float("{:.5f}".format(runtime))) + ")")


	return True


call_counter = 0


def Sleep(func):
	"""Sleeps a given time and returns"""
	global call_counter

	def wrapper():
		global call_counter
		call_counter += 1
		sleep_random_time = random.randrange(5, 15)

		if call_counter % 100 == 0:
			LogAndPrint(f"### Script is going to sleep {sleep_random_time} seconds. ###", False)
			time.sleep(sleep_random_time)
			return func(1)

		return func()

	return wrapper


@Sleep
def SleepWakeUp(state=0):
	"""Returns from sleep"""
	if state == 1:
		LogAndPrint("Script woke up and continue the process.", False)


def LogAndPrint(text, is_printing=True, level="info"):
	"""Takes a text, that logs into a log file and print into the console"""
	if level == "info":
		logging.info(text)
	else:
		logging.debug(text)
	if is_printing:
		print(text)
	return True


def GOSlimRequestURL(this_go_id, this_taxid):

	# Change "_" into ":" in the GO term name
	go_id_include = this_go_id.replace("_", "%3A")

	# All seven taxids to include into taxonID place
	# 4932%2C4896%2C9606%2C7227%2C7955%2C6239%2C3702%2C284812%2C559292

	requestURL = "https://" \
	            + "www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" \
				+ "selectedFields=geneProductId&" \
				+ "selectedFields=symbol&" \
				+ "selectedFields=goId&" \
				+ "selectedFields=goName&" \
                + "selectedFields=taxonId&" \
				+ "selectedFields=synonyms&" \
				+ "selectedFields=goAspect&" \
				+ "goId=" + go_id_include + "&" \
				+ "goUsage=descendants&" \
				+ "goUsageRelationships=is_a%2Cpart_of%2Coccurs_in&" \
				+ "taxonId=" + str(this_taxid) + "&" \
				+ "taxonUsage=descendants"

	return requestURL

def GOSlimRequestURL_AllTaxon(this_go_id):

	# Change "_" into ":" in the GO term name
	go_id_include = this_go_id.replace("_", "%3A")

	requestURL = "https://" \
	            + "www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" \
				+ "selectedFields=geneProductId&" \
				+ "selectedFields=symbol&" \
				+ "selectedFields=goId&" \
				+ "selectedFields=goName&" \
                + "selectedFields=taxonId&" \
				+ "selectedFields=synonyms&" \
				+ "selectedFields=goAspect&" \
				+ "goId=" + go_id_include + "&" \
				+ "goUsage=descendants&" \
				+ "goUsageRelationships=is_a%2Cpart_of%2Coccurs_in"

	requestURL = "https://" \
	             + "www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" \
	             + "goId=" + go_id_include + "&" \
	             + "goUsage=descendants&" \
	             + "goUsageRelationships=is_a%2Cpart_of%2Coccurs_in&page=3"

	return requestURL


class Children:
	"""Gets the children GO terms of a GO id"""
	def __init__(self, goid):
		"""Gets all children terms of a GeneOntology ID"""
		self.goid = goid
		self.export_filename = "export/" + goid.replace(":", "_") + "_children.tsv"
		self.exist = False
		self.no_children = False
		self.call_counter = 0
		self.list_of_goterms = []
		self.dict_of_goterms = {}
		self.export_to_tsv = []
		self.successful_run = self.GetChildren(goid, 0)
		if self.successful_run:
			self.WriteChildren()

	def SleepCall(self):
		"""Sleeps a given time and returns"""
		# Adds one to the counter
		self.call_counter += 1
		# In every 100th request, going to sleep for a random time
		if self.call_counter % 100 == 0:
			sleep_random_time = random.randrange(5, 15)
			LogAndPrint(f"### Script is going to sleep {sleep_random_time} seconds. ###", False)
			time.sleep(sleep_random_time)
			return True
		return False

	def Logger(self, text, is_printing=True, level="info"):
		"""Takes a text string, then logs into a log file and print into the console"""
		if level == "info":
			logging.info(text)
		else:
			logging.debug(text)
		if is_printing:
			print(text)
		return True

	def FileExist(self):
		"""Checks file status, if exists skip this GO id"""
		if os.path.exists(self.export_filename):
			self.exist = True
			return True
		else:
			return False

	def GetEBIRequest(self, goid):
		"""Takes a request for QuickGO, returns page content"""
		url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + goid + "/children"
		hdr = {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
		       'User-Agent': "Magick Browser"}
		page = rqs.Request(url, headers=hdr)
		try:
			content = rqs.urlopen(page).read()
			return content
		except (urllib.error.URLError, urllib.error.HTTPError) as e:
			logging.warning(f"URL request failed: {e}")
			return False

	def GetChildren(self, goid, level, parent_goid=None):
		"""Gets the children GO terms of a GO id"""
		isSleep = self.SleepCall()
		content = self.GetEBIRequest(goid)
		if not content:
			return False
		json_content = json.loads(content)

		if level == 0:
			if self.FileExist():
				return False
			self.Logger(f'Start {goid} request from EBI.')
			parent_goid = goid

		if 'children' in json_content['results'][0]:
			for children in json_content['results'][0]['children']:
				if children['hasChildren']:
					if children['id'] not in self.list_of_goterms:
						self.Logger("\t"*level
							+ "(lv-" + str(level) + ") "
							+ 'Children of ' + children['id'] + " (" + children['name'] + "):")
						self.GetChildren(children['id'], level+1, children['id'])
						self.list_of_goterms.append(children['id'])
						self.dict_of_goterms[children['id']] = children['name']
						line = [children['id'], children['name'], children['relation'], 1, parent_goid]
						self.export_to_tsv.append(line)
					else:
						self.Logger("\t"*level + "(lv-" + str(level) + ") " + children['id'] + ' (children above)')
				else:
					if children['id'] not in self.list_of_goterms:
						self.list_of_goterms.append(children['id'])
						self.dict_of_goterms[children['id']] = children['name']
						line = [children['id'], children['name'], children['relation'], 0, parent_goid]
						self.export_to_tsv.append(line)
						self.Logger("\t"*level + "(lv-" + str(level) + ") " + children['id'] + ": " + children['name'])
					else:
						self.Logger("\t"*level + "(lv-" + str(level) + ") " + children['id'] + " (repeat)")
			return True
		else:
			self.no_children = True
			return False

	def WriteChildren(self, this_mode='a'):
		"""Call self, write into a TSV file the output of GO Children"""

		counter = 0
		with open(self.export_filename, mode=this_mode) as export_file:
			writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			for line in self.export_to_tsv:
				counter += 1
				writer.writerow(line)
		log_this = "Filename: " + self.export_filename + " successful."
		self.Logger(log_this)

		return counter


class Annotations:
	"""Gets the children GO terms of a GO id"""

	def __init__(self, goid):
		"""Gets all children terms of a GeneOntology ID"""
		self.goid = goid
		self.export_filename = "export/" + goid.replace(":", "_") + "_annotations.tsv"
		self.children_filename = "export/" + goid.replace(":", "_") + "_children.tsv"
		self.annotations_filenames = ["data/QuickGO-annotations-01-20210519.tsv", "data/QuickGO-annotations-02-20210519.tsv"]
		self.annotations_taxids = [[9606, 559292, 284812], [3702, 6239, 7955, 7227]]
		self.file_exist = self.FileExistAnnotations()
		self.has_children = self.FileExistChildren()
		self.call_counter = 0
		self.list_of_goterms = []
		self.dict_of_goterms = {}
		self.export_to_tsv = []

		self.ReadQuickGOAnnotation(0)
		# self.successful_run = self.GetChildren(goid, 0)
		# if self.successful_run:
		# 	self.WriteAnnotations()

	def Logger(self, text, is_printing=True, level="info"):
		"""Takes a text string, then logs into a log file and print into the console"""
		if level == "info":
			logging.info(text)
		else:
			logging.debug(text)
		if is_printing:
			print(text)
		return True

	def FileExistChildren(self):
		"""Checks file status of Children, returns a boolean value"""
		if os.path.exists(self.children_filename):
			return True
		else:
			return False

	def FileExistAnnotations(self):
		"""Checks file status of the Export file, returns a boolean variable"""
		if os.path.exists(self.export_filename):
			return True
		else:
			return False

	def ReadQuickGOAnnotation(self, file_num):

		df = pd.read_csv(self.annotations_filenames[file_num], sep='\t', low_memory=False)

		print(df)

		return True

	def WriteAnnotations(self, this_mode='a'):
		"""Call self, write into a TSV file the output of GO Children"""
		counter = 0
		with open(self.export_filename, mode=this_mode) as export_file:
			writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			for line in self.export_to_tsv:
				counter += 1
				writer.writerow(line)
		log_this = "Filename: " + self.export_filename + " successful."
		self.Logger(log_this)

		return counter


def WriteTSVFile(go_id, tax_id, export_filename, write_this_array, split=False, this_mode='a'):
	"""Write into a TSV file the output of GO Children"""
	counter = 0

	with open(export_filename, mode=this_mode) as export_file:
		writer = csv.writer(export_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

		for line in write_this_array:
			if split:
				line = line.split(split)

			counter += 1

			try:
				writer.writerow(line)
			except (csv.Error, TypeError) as e:
				log_this = "GO id: " + go_id + "TaxID: " + str(tax_id) + "Counter: " + str(counter)
				logging.error(f"{log_this}: {e}")

	log_this = "GO id: " + go_id + "TaxID: " + str(tax_id) + "Filename: " + export_filename + "successful."
	logging.info(log_this)

	return counter
