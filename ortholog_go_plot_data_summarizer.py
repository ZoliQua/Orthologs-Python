import pandas as pd
from ortholog_uniprot_name_retriever import *

go_plt_data = pd.read_excel("data/edges_10_annotation_normal_scale.xlsx")

export_header = ["GO_ID", "ALL_KOG_GROUPS",
                 "HS", "HS_KOG_HIT_GROUPS", "HS_KOG_ALL_GROUPS", "HS_DIFF_NUM", "HS_DIFF_UNIPROTS",
                 "DM", "DM_KOG_HIT_GROUPS", "DM_KOG_ALL_GROUPS", "DM_DIFF_NUM", "DM_DIFF_UNIPROTS",
                 "CE", "CE_KOG_HIT_GROUPS", "CE_KOG_ALL_GROUPS", "CE_DIFF_NUM", "CE_DIFF_UNIPROTS",
                 "DR", "DR_KOG_HIT_GROUPS", "DR_KOG_ALL_GROUPS", "DR_DIFF_NUM", "DR_DIFF_UNIPROTS",
                 "AT", "AT_KOG_HIT_GROUPS", "AT_KOG_ALL_GROUPS", "AT_DIFF_NUM", "AT_DIFF_UNIPROTS",
                 "SC", "SC_KOG_HIT_GROUPS", "SC_KOG_ALL_GROUPS", "SC_DIFF_NUM", "SC_DIFF_UNIPROTS",
                 "SP", "SP_KOG_HIT_GROUPS", "SP_KOG_ALL_GROUPS", "SP_DIFF_NUM", "SP_DIFF_UNIPROTS"]

export_header2 = ["GO_ID", "ALL_KOG_GROUPS",
                 "HS", "HS_KOG_HIT_GROUPS", "HS_KOG_ALL_GROUPS", "HS_DIFF_NUM", "HS_DIFF_UNIPROTS", "HS_DIFF_KOGS",
                 "DM", "DM_KOG_HIT_GROUPS", "DM_KOG_ALL_GROUPS", "DM_DIFF_NUM", "DM_DIFF_UNIPROTS", "DM_DIFF_KOGS",
                 "CE", "CE_KOG_HIT_GROUPS", "CE_KOG_ALL_GROUPS", "CE_DIFF_NUM", "CE_DIFF_UNIPROTS", "CE_DIFF_KOGS",
                 "DR", "DR_KOG_HIT_GROUPS", "DR_KOG_ALL_GROUPS", "DR_DIFF_NUM", "DR_DIFF_UNIPROTS", "DR_DIFF_KOGS",
                 "AT", "AT_KOG_HIT_GROUPS", "AT_KOG_ALL_GROUPS", "AT_DIFF_NUM", "AT_DIFF_UNIPROTS", "AT_DIFF_KOGS",
                 "SC", "SC_KOG_HIT_GROUPS", "SC_KOG_ALL_GROUPS", "SC_DIFF_NUM", "SC_DIFF_UNIPROTS", "SC_DIFF_KOGS",
                 "SP", "SP_KOG_HIT_GROUPS", "SP_KOG_ALL_GROUPS", "SP_DIFF_NUM", "SP_DIFF_UNIPROTS", "SP_DIFF_KOGS",
                 "KOG_DIFF_2", "KOG_DIFF_3", "KOG_DIFF_4", "KOG_DIFF_5", "KOG_DIFF_6", "KOG_DIFF_7"]

export_file = []
export_file2 = []

diff_arr_container = []

for go_id in go_plt_data["SPECIES"]:

	file_to_read_protlist = "data/" + go_id.replace('GO_', 'go-') + "-ordered-2022-01.tsv"
	this_go2 = pd.read_csv(file_to_read_protlist, sep = '\t', skiprows = 1, header = None)
	filtered_go2 = this_go2[1].between(0.9, 1.0)

	my_range = range(7, 14)
	what_kog_arr = {}

	for n in my_range:

		k = n + 7

		# Iterate for all uniprot IDs
		arr_container = []
		for index, val in this_go2[filtered_go2].iterrows():
			if type(val[n]) != str:
				continue
			if val[n].find(",") != -1:
				arr = val[n].split(",")
				arr_container.extend(arr)
				for uniprot in arr:
					what_kog_arr[uniprot] = val[0]
			else:
				arr_container.append(val[n])
				what_kog_arr[val[k]] = val[0]

		# for uniprot in arr_container:

		# Iterate for hit uniprot IDs
		arr_hit_container = []
		for index, val in this_go2[filtered_go2].iterrows():
			if type(val[k]) != str:
				continue
			if val[k].find(",") != -1:
				arr = val[k].split(",")
				arr_hit_container.extend(arr)
				for uniprot in arr:
					what_kog_arr[uniprot] = val[0]
			else:
				arr_hit_container.append(val[k])
				what_kog_arr[val[k]] = val[0]

		diff = list(set(arr_container) - set(arr_hit_container))
		diff_arr_container.extend(diff)


# print(diff_arr_container)
#
# tt = ["P05067", "P12345", "P87176", "G2TRL1"]
# tt = diff_arr_container[0:100]
#
# job_id = submit_id_mapping(
# 	from_db="UniProtKB_AC-ID", to_db="Gene_Name", ids=tt
# 	)
#
# if check_id_mapping_results_ready(job_id):
# 	link = get_id_mapping_results_link(job_id)
# 	results = get_id_mapping_results_search(link)
#
# print(results["results"][0])

counter = 0
for go_id in go_plt_data["SPECIES"]:
	file_to_read = "data/" + go_id.replace('GO_', 'go-') + "-ordered-detailed-2022-01.tsv"

	this_go = pd.read_csv(file_to_read, sep = '\t', skiprows = 1, header = None)

	filtered_go = this_go[1].between(0.9, 1.0)

	at_hit = this_go[filtered_go][8].sum()
	at_all = this_go[filtered_go][10].sum()
	at_kog_hit = (this_go[filtered_go][8] != 0).sum()
	at_kog = (this_go[filtered_go][10] != 0).sum()
	ce_hit = this_go[filtered_go][13].sum()
	ce_all = this_go[filtered_go][15].sum()
	ce_kog_hit = (this_go[filtered_go][13] != 0).sum()
	ce_kog = (this_go[filtered_go][15] != 0).sum()
	dm_hit = this_go[filtered_go][18].sum()
	dm_all = this_go[filtered_go][20].sum()
	dm_kog_hit = (this_go[filtered_go][18] != 0).sum()
	dm_kog = (this_go[filtered_go][20] != 0).sum()
	dr_hit = this_go[filtered_go][23].sum()
	dr_all = this_go[filtered_go][25].sum()
	dr_kog_hit = (this_go[filtered_go][23] != 0).sum()
	dr_kog = (this_go[filtered_go][25] != 0).sum()
	hs_hit = this_go[filtered_go][28].sum()
	hs_all = this_go[filtered_go][30].sum()
	hs_kog_hit = (this_go[filtered_go][28] != 0).sum()
	hs_kog = (this_go[filtered_go][30] != 0).sum()
	sc_hit = this_go[filtered_go][33].sum()
	sc_all = this_go[filtered_go][35].sum()
	sc_kog_hit = (this_go[filtered_go][33] != 0).sum()
	sc_kog = (this_go[filtered_go][35] != 0).sum()
	sp_hit = this_go[filtered_go][38].sum()
	sp_all = this_go[filtered_go][40].sum()
	sp_kog_hit = (this_go[filtered_go][38] != 0).sum()
	sp_kog = (this_go[filtered_go][40] != 0).sum()

	export_row = [go_id, len(this_go.index),
	              go_plt_data["HS"][counter], hs_kog_hit, hs_kog, hs_hit, hs_all,
	              go_plt_data["DM"][counter], dm_kog_hit, dm_kog, dm_hit, dm_all,
	              go_plt_data["CE"][counter], ce_kog_hit, ce_kog, ce_hit, ce_all,
	              go_plt_data["DR"][counter], dr_kog_hit, dr_kog, dr_hit, dr_all,
	              go_plt_data["AT"][counter], at_kog_hit, at_kog, at_hit, at_all,
	              go_plt_data["SC"][counter], sc_kog_hit, sc_kog, sc_hit, sc_all,
	              go_plt_data["SP"][counter], sp_kog_hit, sp_kog, sp_hit, sp_all]



	export_file.append(export_row)

	file_to_read_protlist = "data/" + go_id.replace('GO_', 'go-') + "-ordered-2022-01.tsv"
	this_go2 = pd.read_csv(file_to_read_protlist, sep = '\t', skiprows = 1, header = None)
	filtered_go2 = this_go2[1].between(0.9, 1.0)

	diff_arr = {}
	diff_kog_arr = {}
	my_range = range(7, 14)
	what_kog_arr = {}
	what_uniprot_in_kog_arr = {}

	for n in my_range:

		k = n + 7

		# Iterate for all uniprot IDs
		arr_container = []
		for index, val in this_go2[filtered_go2].iterrows():
			if type(val[n]) != str:
				continue
			if val[n].find(",") != -1:
				arr = val[n].split(",")
				arr_container.extend(arr)
				for uniprot in arr:
					what_kog_arr[uniprot] = val[0]
			else:
				arr_container.append(val[n])
				what_kog_arr[val[k]] = val[0]

		# for uniprot in arr_container:

		# Iterate for hit uniprot IDs
		arr_hit_container = []
		for index, val in this_go2[filtered_go2].iterrows():
			if type(val[k]) != str:
				continue
			if val[k].find(",") != -1:
				arr = val[k].split(",")
				arr_hit_container.extend(arr)
				for uniprot in arr:
					what_kog_arr[uniprot] = val[0]
			else:
				arr_hit_container.append(val[k])
				what_kog_arr[val[k]] = val[0]

		diff = list(set(arr_container) - set(arr_hit_container))

		this_kog_diff = []
		for uniprot in diff:
			this_kog_diff.append(what_kog_arr[uniprot])
			if what_kog_arr[uniprot] not in what_uniprot_in_kog_arr:
				what_uniprot_in_kog_arr[what_kog_arr[uniprot]] = [uniprot]
			else:
				what_uniprot_in_kog_arr[what_kog_arr[uniprot]].append(uniprot)

		# job_id = submit_id_mapping(
		# 	from_db="UniProtKB_AC-ID", to_db="Gene_Name", ids=diff
		# 	)
		#
		# if check_id_mapping_results_ready(job_id):
		#     link = get_id_mapping_results_link(job_id)
		#     results = get_id_mapping_results_search(link)
		#
		# print(results["results"][0])
		#
		# break

		diff_arr_container.extend(diff)
		diff_arr[n] = diff
		diff_kog_arr[n] = this_kog_diff

	kog = {}
	species_arr = {7: "AT", 8: "CE", 9: "DM", 10: "DR", 11: "HS", 12: "SC", 13: "SP"}
	for nr in diff_kog_arr:
		for uniprot_element in diff_kog_arr[nr]:

			this_spec = species_arr[nr]

			if uniprot_element in kog.keys():
				if this_spec not in kog[uniprot_element]:
					kog[uniprot_element].append(this_spec)
			else:
				kog[uniprot_element] = [this_spec]

	kog_diff_2 = []
	kog_diff_3 = []
	kog_diff_4 = []
	kog_diff_5 = []
	kog_diff_6 = []
	kog_diff_7 = []

	for kog_group in kog:
		if len(kog[kog_group]) > 1:
			this_line = kog_group + ": " + ",".join(kog[kog_group]) + " - " + " ".join(what_uniprot_in_kog_arr[kog_group])
			if len(kog[kog_group]) == 2:
				kog_diff_2.append(this_line)
			elif len(kog[kog_group]) == 3:
				kog_diff_3.append(this_line)
			elif len(kog[kog_group]) == 4:
				kog_diff_4.append(this_line)
			elif len(kog[kog_group]) == 5:
				kog_diff_5.append(this_line)
			elif len(kog[kog_group]) == 6:
				kog_diff_6.append(this_line)
			elif len(kog[kog_group]) == 7:
				kog_diff_7.append(this_line)


	export_row = [go_id, len(this_go2.index),
		        go_plt_data["HS"][counter], hs_kog_hit, hs_kog, (hs_all - hs_hit), ",".join(diff_arr[11]), ",".join(diff_kog_arr[11]),
		        go_plt_data["DM"][counter], dm_kog_hit, dm_kog, (dm_all - dm_hit), ",".join(diff_arr[9]), ",".join(diff_kog_arr[9]),
		        go_plt_data["CE"][counter], ce_kog_hit, ce_kog, (ce_all - ce_hit), ",".join(diff_arr[8]), ",".join(diff_kog_arr[8]),
		        go_plt_data["DR"][counter], dr_kog_hit, dr_kog, (dr_all - dr_hit), ",".join(diff_arr[10]), ",".join(diff_kog_arr[10]),
		        go_plt_data["AT"][counter], at_kog_hit, at_kog, (at_all - at_hit), ",".join(diff_arr[7]), ",".join(diff_kog_arr[7]),
		        go_plt_data["SC"][counter], sc_kog_hit, sc_kog, (sc_all - sc_hit), ",".join(diff_arr[12]), ",".join(diff_kog_arr[12]),
		        go_plt_data["SP"][counter], sp_kog_hit, sp_kog, (sp_all - sp_hit), ",".join(diff_arr[13]), ",".join(diff_kog_arr[13]),
	            ",".join(kog_diff_2), ",".join(kog_diff_3), ",".join(kog_diff_4), ",".join(kog_diff_5), ",".join(kog_diff_6), ",".join(kog_diff_7)]

	export_file2.append(export_row)

	counter += 1

df = pd.DataFrame(export_file, columns=export_header)
df2 = df.set_index("GO_ID")

df1b = pd.DataFrame(export_file2, columns = export_header2)
df2b = df1b.set_index("GO_ID")

#df2 = df.T

#df2.to_excel("data/edges_10_annotation_normal_scale_cutoff_09_numk.xlsx")
df2b.to_excel("data/edges_10_annotation_normal_scale_cutoff_09_uniprotsk.xlsx")
