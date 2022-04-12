##3rd attempt: use deep learning
##each edge i->j is a data point
##input/features will be the gene expressions of i and j
##output is contact or synapse or ratio


import pandas as pd
import networkx as nx
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import math
import pickle #for saving models

##importing local package helpers
##need the sys stuff because it's located in a different folder
##append parent folder to sys.path so it knows to look there for utils folder
import sys
sys.path.append('../')
import utils.helpers as helpers


from sklearn.linear_model import LinearRegression
from scipy.special import expit #the logistic sigmoid function
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from math import sqrt
from sklearn.metrics import r2_score

######################################################################
# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime('%d-%m-%Y-%H%M%S')
print('now = ' + dt_string)



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

###################################################################
##load data
##get the synapses, forget polyad structure
syn = helpers.get_synapse_list()
syn = helpers.polyad_to_monad(syn)


##get gene expressions
##rows=genes, col=cells
exp = helpers.get_gene_expression()

##get contact matrix, then edgelist
#also need all the cell names in order to compare with syn
contact_adj = helpers.get_contact_adj()
contact_cells = set(contact_adj.index)
#contact_edgelist = helpers.contact_list_from_edge(contact_adj)
##put False to keep the zeroes; for predicting contact only
contact_edgelist = helpers.contact_list_from_edge(contact_adj, False)

###################################################################
##restrict to cells in common in syn and contact
syn_cells = set(syn['pre']).union(set(syn['post']))

common_cells = syn_cells.intersection(contact_cells)
common_cells = common_cells.intersection(set(exp.columns))
##for contact only
#common_cells = contact_cells.intersection(set(exp.columns))

common_cells_list = list(common_cells)

###################################################################
##restrict syn, contact, gene to common cells
syn = syn[(syn['pre'].isin(common_cells)) & (syn['post'].isin(common_cells))]
contact = contact_edgelist[(contact_edgelist['pre'].isin(common_cells)) & (contact_edgelist['post'].isin(common_cells))]
exp = exp[common_cells]

###################################################################
##combine the syn and contact
##strategy: stack the two dataframes,
##then groupby pre & post
syn['pixels'] = 0
contact['sections'] = 0
contact = contact[syn.columns] ##rearrange columns

syn_contact = syn.append(contact)
syn_contact = syn_contact.groupby(['pre','post']).sum().reset_index()

##add column of ratio = sections/pixels
##apparently some synapses form without contact?!?
syn_contact = syn_contact[syn_contact.pixels > 0]
syn_contact['ratio'] = syn_contact.apply(lambda row : row['sections'] / row['pixels'] * 10000.0, axis=1)


###################################################################
##in preparation for trying to predict synapse from exp

syn_adj = helpers.edge_to_adj(syn)
syn_filled = helpers.contact_list_from_edge(syn_adj, False)
syn_filled.columns = ['pre', 'post', 'sections']

###################################################################
##selecting subset of genes;
##not enough data points -> only about 6000
num_genes = len(exp.index)  ##set this to smaller if want to subset
genes_ind = np.random.choice(len(exp.index),num_genes)

##list of selected gene,
##also dictionary to get ind from gene
#gene_of = list(exp.index)
gene_of = [exp.index[i] for i in genes_ind]
ind_of_gene = { gene : ind for ind,gene in enumerate(gene_of) }

##likewise index the cells
cell_of = list(exp.columns)
ind_of_cell =  { cell : ind for ind,cell in enumerate(cell_of) }
num_cells = len(cell_of)

##vector of 1s, useful
ones_row = np.full((1,num_genes), 1.0)
ones_col = np.full((num_genes,1), 1.0)


##a bit more convenient to use the raw numpy array
##each row is a cell
##extract rows corresponding to selected genes
exp_subset = exp.loc[gene_of]
expT = exp_subset.values.transpose()


##the data points; populate with append_data below
##each row in syn_contact, i.e. each directed edge,
##gives one data point
x = [] ##later turn into np array
y = []

##append data point corresponding to edge i->j between cells i,j
##in our use case, E should be expT
##adds experiment y_val (the observed ratio, r['ratio'])
##adds 1d np array to x consisting of E[i] concat E[j]
def append_data(E, i, j, y_val):
	global x, y
	y.append(y_val)
	x.append(np.concatenate((E[i], E[j])))


##append data; predicting ratio
#syn_contact.apply(lambda row :
#	append_data(
#		expT,
#		ind_of_cell[row['pre']],
#		ind_of_cell[row['post']],
#		row['pixels'])
#	, axis=1)

###append data: predicting synapses
#syn_filled.apply(lambda row :
#	append_data(
#		expT,
#		ind_of_cell[row['pre']],
#		ind_of_cell[row['post']],
#		row['sections'])
#	, axis=1)

##append data: predicting contact
contact.apply(lambda row :
	append_data(
		expT,
		ind_of_cell[row['pre']],
		ind_of_cell[row['post']],
		row['pixels'])
	, axis=1)

###################################################################
##data loaded, now run!

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.30,
        random_state=40)


##this takes a long time,
##and scores are bad..
##-0.095 for the ratio, -0.32 for synapses, -1.0 for contact..

mlp_contact = MLPRegressor(
	hidden_layer_sizes=(100,50,25),
	activation='relu',
	solver='lbfgs',
	max_iter=np.inf)
mlp_contact.fit(x_train, y_train)
mlp_contact.score(x_test, y_test)



##save model
model_filename = 'model-100-50-20-contact.pkl'
with open(model_filename, 'wb') as file:
	pickle.dump(mlp_contact, file)


##load
#with open(model_filename, 'rb') as file:
#	loaded_model = pickle.load(file)


###################################################################
##save in csv, with gene names

##extract model coefficients
wstar = model.coef_[0:num_genes**2]
bstar = model.coef_[num_genes**2:]

##reshape
wstar = wstar.reshape(num_genes,num_genes)
bstar = bstar.reshape(num_genes,num_genes)

dfw = pd.DataFrame(wstar)
dfw.columns = gene_of
dfw.index = gene_of
dfw.to_csv('wstar-full-'+dt_string+'.csv',encoding='utf-8-sig')

dfb = pd.DataFrame(bstar)
dfb.columns = gene_of
dfb.index = gene_of
dfb.to_csv('bstar-full+'dt_string+'.csv',encoding='utf-8-sig')

###################################################################
##transform to log scale because numbers are large
##cut off for small values
##also negative numbers do -log(|x|)
epsilon = 1e-9
def abslog(x):
	if abs(x) < epsilon:
		return 0
	elif x > 0:
		return math.log(x)
	else:
		return -math.log(-x)


dfw_log = dfw.applymap(abslog)
dfb_log = dfb.applymap(abslog)

dfw_log.to_csv('wstar-log-full-'+dt_string+'.csv',encoding='utf-8-sig')
dfb_log.to_csv('bstar-log-full-'+dt_string+'.csv',encoding='utf-8-sig')

###################################################################
##TODO do polyads


##abandoning the linear regression model for now
##as there are too few data points to be useful

###################################################################
###################################################################
##second attempt:
##


###################################################################
##saving pairs that have no contact but have synapses
#bad_synapse = syn_contact[(syn_contact.pixels == 0)]
#bad_synapse = bad_synapse.sort_values('sections', ascending=False)
#bad_synapse.to_csv('bad_synapse.csv',encoding='utf-8-sig')



###################################################################
##stuff not used, old idea

##weight, bias matrix to be optimized
w = np.zeros((num_genes, num_genes))
w_rowsum = np.sum(w,1)
w_colsum = np.sum(w,0)
b = np.zeros((num_genes, num_genes))


##the function p_est is simple enough that we can precompute
##some weights stuff to speed up
##this should be called before p_est
def update_weight_sums():
	global w_rowsum, w_colsum
	w_rowsum = np.sum(w,1)
	w_colsum = np.sum(w,0)


##currently not using this
##
##the function, estimates amount of synapse formed i->j
##E = exp_ind = expression matrix, columns = cells
##(in actual use should be exp_ind, not the E from above)
##i,j = cells (numbers; cell i is the ith row in exp_ind)
##for example exp['PDA']['bam-2']
##(a bit weird because column index comes first)
##use rather simple function,
##which is summing the following over all pairs of genes A,B:
## w[A,B] * ( E[j][B] - E[i][A] ) + b[A,B] * E[j][B] * E[i][A]
##which simplifies to (. is dot product)
## E[j] . w_colsum - E[i] . w_rowsum + E[i]^T * b * E[j]
##
##before use, must have run update_weight_sums() if updated w
def p_est(E, i, j):
	##again warning E should be exp_ind
	return np.dot(E[j], w_colsum) - np.dot(E[i], w_rowsum) + E[i] @ b @ E[j]
	#return expit(..) #expit = 1 / 1 + e^-x



##TODO do contact, include also the 0 contact!


##genes from cengen
alt_genes = ['cam-1'
, 'bam-2'
, 'C24H10.1'
, 'C48E7.6'
, 'C54G4.4'
, 'casy-1'
, 'cdh-1'
, 'cdh-12'
, 'cdh-3'
, 'cdh-4'
, 'cdh-8'
, 'cdh-9'
, 'clc-1'
, 'clc-2'
, 'clc-3'
, 'clc-4'
, 'clr-1'
, 'crb-1'
, 'dgn-1'
, 'dig-1'
, 'dma-1'
, 'eat-20'
, 'efn-2'
, 'efn-4'
, 'egg-6'
, 'egl-15'
, 'fmi-1'
, 'fmil-1'
, 'fshr-1'
, 'glit-1'
, 'hic-1'
, 'him-4'
, 'hmr-1'
, 'hpo-30'
, 'igcm-1'
, 'igcm-2'
, 'igcm-3'
, 'igcm-4'
, 'igdb-1'
, 'igdb-2'
, 'igeg-1'
, 'igeg-2'
, 'iglr-1'
, 'iglr-2'
, 'iglr-3'
, 'ina-1'
, 'K10D6.2'
, 'lad-2'
, 'lat-1'
, 'lat-2'
, 'let-805'
, 'lin-17'
, 'lon-2'
, 'lron-10'
, 'lron-11'
, 'lron-12'
, 'lron-13'
, 'lron-14'
, 'lron-15'
, 'lron-3'
, 'lron-4'
, 'lron-5'
, 'lron-6'
, 'lron-7'
, 'lron-9'
, 'madd-4'
, 'mig-6'
, 'mua-3'
, 'mup-4'
, 'ncam-1'
, 'nlg-1'
, 'nlr-1'
, 'nrx-1'
, 'oig-1'
, 'oig-2'
, 'oig-4'
, 'oig-5'
, 'oig-8'
, 'pan-1'
, 'pat-2'
, 'pat-3'
, 'ptp-3'
, 'ptp-4'
, 'pxn-1'
, 'pxn-2'
, 'rig-1'
, 'rig-3'
, 'rig-4'
, 'rig-5'
, 'rig-6'
, 'sax-3'
, 'sax-7'
, 'sdn-1'
, 'slt-1'
, 'smp-1'
, 'smp-2'
, 'syg-1'
, 'syg-2'
, 'sym-1'
, 'T05E11.2'
, 'ten-1'
, 'tol-1'
, 'tsp-1'
, 'tsp-11'
, 'tsp-12'
, 'tsp-13'
, 'tsp-14'
, 'tsp-15'
, 'tsp-16'
, 'tsp-17'
, 'tsp-18'
, 'tsp-19'
, 'tsp-20'
, 'tsp-21'
, 'tsp-3'
, 'tsp-4'
, 'tsp-5'
, 'tsp-6'
, 'tsp-7'
, 'tsp-8'
, 'tsp-9'
, 'unc-40'
, 'unc-5'
, 'unc-52'
, 'vab-1'
, 'vab-2'
, 'vab-9'
, 'ver-3'
, 'ver-4'
, 'wrk-1'
, 'zig-1'
, 'zig-10'
, 'zig-11'
, 'zig-12'
, 'zig-2'
, 'zig-3'
, 'zig-4'
, 'zig-5'
, 'zig-7'
, 'zig-8'
, 'zig-9'
]
