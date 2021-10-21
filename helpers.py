import pandas as pd
import math


######################################################################
## helper functions for dealing with neuron names

##	clean the post neurons list
##	in SI-3-Synapse-lists-male.csv,
##	post synaptic neurons are given as a comma-separated list
##	of neurons
##	entries need to be cleaned/removed, and sorted
##	e.g. HOA,[AVG] -> AVG,HOA
##	note if everything is removed, then returns empty string
def clean_post(post):
	return apply_fn_post(clean_neuron_name, post)


##apply function fn : str -> str to each entry of comma-separated list post
##removes entries that return empty string under fn
##sorts also
def apply_fn_post(fn, post):
	post = post.replace(" ","")
	post_list = post.split(",")
	post_list = map(fn, post_list)
	post_list = filter(lambda x : x != "" , post_list)
	post_list = sorted(post_list)
	post = ",".join(post_list)
	return post


def clean_neuron_name(neuron):
	##remove spaces - sneaky errors!
	neuron = neuron.replace(" ","")
	##remove square brackets (treat unsure as sure)
	neuron = neuron.replace("[","")
	neuron = neuron.replace("]","")
	##some names are "old", use name in SI-4 cell list
	neuron = convert_name(neuron)
	if good_cellname(neuron):
		return neuron
	else:
		return ""


##found manually
convert_dict = {
	"PVS" : "PVPR",
	"PVU" : "PVPL",
	"adp" : "mu_anal",
	"sph" : "mu_sphincter",
	"intL" : "mu_intL",
	"intR" : "mu_intR"
}
## convert name of cell that appears in SI 3 Synapse List
## to that appearing in SI 4 Celllist
def convert_name(name):
	global convert_dict
	if name in convert_dict:
		return convert_dict[name]
	else:
		return name


##get all male cells from SI-4,
##which were manually stored in file male_celllist.csv
#(but also compare https://www.wormatlas.org/celllistsulston.htm)
#seems some cells like HSNL/R, VC01(VC1) are there, but not in celllist
celllist = pd.read_csv("male_celllist.csv", header=None, comment='#')
#remove unexpected spaces...
celllist = celllist[0].apply(lambda x : x.replace(" ",""))
celllist = celllist.tolist()


## check if cellname is in cell list SI4
## expect name to have been cleaned and converted
def good_cellname(name):
	global celllist
	return name in celllist



######################################################################
##left-right pairs stuff

cellLR = pd.read_csv('male_celllist_LR.csv',index_col=False)

##e.g. 'PVPL' -> 'PVPR'
def flipLR(name):
	global cellLR
	ind = cellLR.index[cellLR.L == name].tolist()
	if len(ind) > 0:
		if len(ind) > 1:
			print('cell name appears twice in cellLR!')
		return cellLR.iloc[ind[0]].R
	ind = cellLR.index[cellLR.R == name].tolist()
	if len(ind) > 0:
		if len(ind) > 1:
			print('cell name appears twice in cellLR!')
		return cellLR.iloc[ind[0]].L
	return name

##get name without L/R
def ridLR(name):
	global cellLR
	ind = cellLR.index[cellLR.L == name].tolist() + cellLR.index[cellLR.R == name].tolist()
	if len(ind) > 0:
		if len(ind) > 1:
			print('cell name appears twice in cellLR!')
		return cellLR.iloc[ind[0]].name
	return name

##if name comes from L/R pair, get it as a list
##if not, return name as one element list
##e.g. 'PVP' -> ['PVPL','PVPR']
##e.g. 'HOA' -> ['HOA']
def uncombineLR(name):
	global cellLR
	ind = cellLR.index[cellLR.name == name].tolist()
	if len(ind) > 0:
		if len(ind) > 1:
			print('cell name appears twice in cellLR!')
		return [cellLR.iloc[ind[0]].L, cellLR.iloc[ind[0]].R]
	return [name]


######################################################################
##for computing similarity

##f(x,y) from pg 4 of Jarrell et al Supplementary Materials
##use less generic name instead of f
C_1 = 0.5
C_2 = 1
def similarity_score_single(x, y):
	return min(x,y) - C_1 * max(x,y) * (math.e ** ( - C_2 * min(x,y) ))

##takes in two pandas.Series (with same column names)
def similarity_score(x, y):
	df = pd.DataFrame()
	df = df.append(x)
	df = df.append(y)
	df.index = [0,1]
	return sum(df.apply(lambda c : similarity_score_single(c[0],c[1]),axis=0))


######################################################################
##some old stuff, used for testing


##	some cells would appear in the SI-3 but not in SI-4,
##	later found that some names were updated,
##	e.g. 'PVS' -> 'PVPR'
##	weird_names records such entries
weird_names = []
def weird_cellname(name):
	global weird_names
	global celllist
	if good_cellname_old(name) != (name in celllist):
		print('Add weird name: ' + name + '; in SI4? : ' + str(name in celllist))
		weird_names.append(name)

##detecting bad names by manually eliminating
##only for diagnosing problem
def good_cellname_old(name):
	return name == "" and name[0:3] != "unk" and name[0:3] != "obj" and name[0:6] != "contin" and name[0:4] != "frag"


#### test flipLR
##df_test = pd.DataFrame()
##df_test['flip_to_L'] = cellLR['R'].apply(flipLR)
##df_test['flip_to_R'] = cellLR['L'].apply(flipLR)
#### expect 0:
##sum(df_test.flip_to_L != cellLR.L)
##sum(df_test.flip_to_R != cellLR.R)
