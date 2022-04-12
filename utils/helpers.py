import pandas as pd
import math
import os

######################################################################
##folder where data is stored
data_source = os.path.join(os.path.dirname(__file__), '../data-source/')
data_aux = os.path.join(os.path.dirname(__file__), '../data-aux/')

######################################################################
## helper functions for cleaning neuron names

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
celllist = pd.read_csv(data_aux+'male_celllist.csv', header=None, comment='#')
#remove unexpected spaces...
celllist = celllist[0].apply(lambda x : x.replace(" ",""))
celllist = celllist.tolist()


## check if cellname is in cell list SI4
## expect name to have been cleaned and converted
def good_cellname(name):
	global celllist
	return name in celllist


######################################################################
## helper functions for post neuron names

##gets the n-th cell in post
##returns empty string if n is out of bounds
def extract_post(n, post):
    arr = post.split(',')
    if len(arr) < n + 1:
        return ''
    else:
        return arr[n]

##returns True if cell is in post (as a comma-sep string)
##when used in df.apply, it's passed as a Series?
##need to convert to string..
def in_post(cell, post):
    post = str(post)
    return cell in post.split(',')


######################################################################
##left-right pairs stuff

cellLR = pd.read_csv(data_aux+'male_celllist_LR.csv',index_col=False)

## L <-> R
## applies to one name and (comma-sep) list of names
## 'PVPL,AVAR' -> 'PVPR,AVAL'
## ['PVPL','AVAR'] -> ['PVPR','AVAL']
def flipLR(name):
	if type(name) == list:
		return map(flipLR_one,name)
	if ',' in name:
		name_list = name.split(',')
		name_list_flipped = map(flipLR_one,name_list)
		return ",".join(name_list_flipped)
	return flipLR_one(name)

##e.g. 'PVPL' -> 'PVPR'
def flipLR_one(name):
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
		return cellLR.iloc[ind[0]]['name']
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


def isLR(name):
	global cellLR
	return (name in list(cellLR.L)) or (name in list(cellLR.R))

def drop_noLR(df_edgelist):
	return df_edgelist.loc[(df_edgelist['pre'].apply(isLR)) | (df_edgelist['post'].apply(isLR))]


######################################################################
##for computing similarity

##f(x,y) from pg 4 of Jarrell et al Supplementary Materials
##use less generic name instead of f
C_1 = 0.5
C_2 = 1
def similarity_score_single(x, y):
	return min(x,y) - C_1 * max(x,y) * (math.e ** ( - C_2 * min(x,y) ))

def simscore(x,y):
	return total_similarity_score(similarity_score_single,x,y)

##L1
def L1_diff_single(a,b):
	return abs(a - b)

def L1_dist(x,y):
	return total_similarity_score(L1_diff_single,x,y)


##takes in two pandas.Series x,y (with same column names)
##and fn of two variables
##sums fn(x[i],y[i]) over columns i
def total_similarity_score(fn, x, y):
	df = pd.DataFrame()
	df = df.append(x)
	df = df.append(y)
	df.index = [0,1]
	return sum(df.apply(lambda c : fn(c[0],c[1]),axis=0))

##convert fn_single into a function on two Series x,y
##optional new function name
def total_fn(fn_single, newfnname=''):
	def totaling_fn(x, y):
		return total_similarity_score(fn_single,x,y)

	if len(newfnname) > 0:
		totaling_fn.__name__ = newfnname
	else:
		totaling_fn.__name__ = 'total_' + fn_single.__name__
	return totaling_fn


##L/R symmetry compute
##expects:
##adj : DataFrame adjacency matrix (pre against post) (may not be symmetric)
##fn : the similarity function to use to compare two rows (Series)
##polyad : True if post is given as polyad (default = True)
##combineLR : True if combine the L/R versions of post
##returns the similarity score for each neuron (LR merged)

def LR_symmetry(adj, fn, polyad=False, combineLR=False):
	#get the pre,post names
	rows = list(adj.index)
	cols = list(adj.columns)

	#apply the flipping to names
	adj_flipped = adj.copy()
	newrows = map(flipLR,rows)
	newcols = map(flipLR,cols)
	adj_flipped.index = newrows
	adj_flipped.columns = newcols

	#sometimes L/R don't both show up; add zero rows/cols
	for r in set(newrows).difference(set(rows)):
		adj.loc[r] = 0

	for c in set(newcols).difference(set(cols)):
		adj[c] = 0

	for r in set(rows).difference(set(newrows)):
		adj_flipped.loc[r] = 0

	for c in set(cols).difference(set(newcols)):
		adj_flipped[c] = 0

	#reorder the rows and cols to match
	combine_rows = list(set(rows).union(set(newrows)))
	combine_cols = list(set(cols).union(set(newcols)))
	adj = adj[combine_cols]
	adj = adj.reindex(combine_rows)
	adj_flipped = adj_flipped[combine_cols]
	adj_flipped = adj_flipped.reindex(combine_rows)

	##expect combine_rows,combine_cols to be L/R symmetric
	##only meant for testing
	#for x in combine_rows:
	#	if flipLR(x) not in combine_rows:
	#		print('combine_rows not sym: ', x)
	#for x in combine_cols:
	#	if flipLR(x) not in combine_cols:
	#		print('combine_cols not sym: ', x)

	#compute the similarity for each row
	df_out = pd.DataFrame((x, fn(adj.loc[x],adj_flipped.loc[x])) for x in combine_rows)
	df_out.columns = ['pre',fn.__name__]
	df_out['pre'] = df_out['pre'].apply(ridLR)
	df_out = df_out.groupby('pre').mean().reset_index()

	return df_out

	#for x in combine_rows:
	#	if df_out.loc[x][1] != df_out.loc[flipLR(x)][1]:
	#		print(x) ##expect empty



##LR_symmetry_pairwise
##compares the number of sections (of contact/synapse/..)
##from cell X to Y
##with X' to Y'
##X',Y' is LR mirror of X,Y
##(if X is not a LR homolog pair, then X' = X)
##(for example, HOA -> PHCL would be compared with HOA -> PHCR,
##so in a sense it computes the asymmetry of HOA)
##expect:
##	df: edgelist, with columns: pre, post, sections,
##	fn: function fn to compare two numbers
##	simscore: (optional) name of similarity score
##								(defaults to name of fn)
##output:
##
##summary of algorithm:
##start with df:
#AVAL, PVCL, 1
#AVAL, PVCR, 1
#AVAR, PVCR, 2
#(AVAR, PVCL, 0) (this row is not actually present in data because 0)
##produces df_all:
#AVAL, PVCL, 1, 2
#AVAL, PVCR, 1, 0
#AVAR, PVCR, 2, 1
#AVAR, PVCL, 0, 1
##applies the fn to compare the two numbers in each row
def LR_symmetry_pairwise(df, fn, simscore_name=''):
	df_flipped = df.copy()
	df_flipped['pre'] = df_flipped['pre'].apply(flipLR)
	df_flipped['post'] = df_flipped['post'].apply(flipLR)
	df_flipped.columns = ['pre','post','sections_flipped']
#
	df_copy = df.copy().sort_values(['pre','post'])
	df_flipped = df_flipped.sort_values(['pre','post'])
#
	#sneaky trick to get the number of sections for the flipped
	df_copy['sections_flipped'] = 0
	df_flipped['sections'] = 0
	df_all = df_copy.append(df_flipped)
	df_all = df_all.groupby(['pre','post']).sum().reset_index()
	#has columns pre, post, sections, sections_flipped
#
	#get function name, serves as new column name
	fn_name = fn.__name__
	if len(simscore_name) > 0:
		fn_name = simscore_name
#
	#apply sim score to each row
	df_all[fn_name] = df_all.apply(lambda row : fn(row.sections, row.sections_flipped), axis=1)
#
	#reorder columns
	df_all = df_all[['pre','post',fn_name,'sections','sections_flipped']]
#
	return df_all

######################################################################
##stuff?

##turn polyads into usual synapses, i.e. forget polyad data
##should assume that names have been cleaned
##input df should have columns 'pre','post','sections'
##outputs dataframe with those columns
##note that this also sums up the number of sections over
##(so there is unique row for each pre,post
##all synapses between cells
def polyad_to_monad(df):
	##max number of post in a single polyad
	max_post = max(df.apply(lambda row : len(row.post.split(',')), axis=1))
	##df to hold result
	df_separate = pd.DataFrame(columns=['pre','post','sections'])
	for n in range(max_post):
		##extract the n-th neuron in the polyad
		##some polyads don't have that many; just remove row
		df_tmp = df.copy()
		df_tmp.post = df.post.apply(lambda post : extract_post(n,post))
		df_tmp = df_tmp[df_tmp.post != '']
		df_separate  = df_separate.append(df_tmp)
	return df_separate.groupby(['pre','post']).sum().reset_index()



######################################################################

##convert edge list to adj matrix
##expect columns of edge list to be 'pre', 'post', and 'sections'
def edge_to_adj(df_edges):
	return pd.crosstab(index=df_edges['pre'],columns=df_edges['post'],values=df_edges['sections'],aggfunc='sum').fillna(0)


######################################################################
##loading data
##all from N2Y

#get synapse list
#columns: pre,post,sections
#cleans the data in the process
def get_synapse_list():
	syn = pd.read_csv(data_source+'SI-3-Synapse-lists-male.csv')
	syn = syn.query('EM_series=="N2Y" & type=="chemical"')
	syn = syn[['pre','post','sections']]

	#clean the pre neurons
	syn['pre'] = syn['pre'].apply(clean_neuron_name)
	#clean the post neurons
	syn['post'] = syn['post'].apply(clean_post)

	#after cleaning, empty string means unknown, filter them out
	syn = syn[(syn.pre != '') & (syn.post != '')]

	return syn


#get adjacency matrix, entries are amount of contact by pixel
#returns dataframe, rows and columns are indexed by cell names
def get_contact_adj():
	contact_adj = pd.read_csv(data_source+'N2Y-PAG-contact-matrix.csv',index_col=0)
	contact_adj = contact_adj.fillna(0)
	return contact_adj


#contact matrix to edge list
#note that the edgelist will have two rows for each edge,
#one for each direction
#clears the zero entries by default
##TODO bad name; should be contact_list_from_adj...
##keeping this for backwards compatibility
def contact_list_from_edge(contact_adj, clear_zeros=True):
	contact_edgelist = contact_adj.stack().reset_index()
	contact_edgelist.columns = ['pre','post','pixels']
	if clear_zeros:
		contact_edgelist = contact_edgelist[contact_edgelist['pixels'] != 0]

	return contact_edgelist

##better version of contact_list_from_edge
def edgelist_from_adj(contact_adj, row_col_val=['row','col','val'], clear_zeros=True):
	edgelist = contact_adj.stack().reset_index()
	edgelist.columns = row_col_val.copy()
	val = row_col_val[2]
	if clear_zeros:
		edgelist = edgelist[edgelist[val] != 0]
#
	return edgelist


#get contact as edge list
#returns dataframe with columns 'pre','post','pixels'
def get_contact_list():
	contact_adj = get_contact_adj()
	return contact_list_from_edg(contact_adj)


#get gene expressions of cells
#returns dataframe, row = gene, col = cells
#it seems all entries are integers
#they're presented as floats though, check entries integer with
#exp.loc[gene][cell].is_integer()
def get_gene_expression():
	exp = pd.read_csv(data_source+'Expression-matrix-Jan-2020.csv',index_col=0)
	exp = exp.fillna(0) ##empty entries interpret as 0

	##drop the last 5 columns ('Unnamed: 23*'); they're all 0
	exp = exp[exp.columns[:-5]]

	return exp



######################################################################

#apply dist to every pair of rows
#outputs an adj matrix, output_df[cell1][cell2] = fn(cell1,cell2)
def pairwise_dist(df, fn):
	rows = df.index
	output_df = pd.DataFrame(index=rows)
	for r in rows:
		print(r)
		output_df[r] = df.apply(lambda rr : fn(df.loc[r],rr), axis=1)
#
	return output_df

def pairwise_dist_col(df, fn):
	cols = df.columns
	output_df = pd.DataFrame(index=cols)
	for c in cols:
		print(c)
		output_df[c] = df.apply(lambda cc : fn(df[c],cc), axis=0)
#
	return output_df




###################################################################
##attempt to add a more convenient method to get rows of dataframe
#e.g. que(syn, {'pre' : 'AVAL'})
def que(self, search_dict):
	if len(search_dict) == 0:
		return self
#
	query_list = [f"(self['{k}'] == {repr(v)})" for (k,v) in search_dict.items()]
	query = " & ".join(query_list)
	print(query)
	return self[eval(query)]



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
