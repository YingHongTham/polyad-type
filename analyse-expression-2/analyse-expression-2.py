import pandas as pd
import networkx as nx
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.special import expit #the logistic sigmoid function

##importing local package helpers
##need the sys stuff because it's located in a different folder
##append parent folder to sys.path so it knows to look there for utils folder
import sys
sys.path.append('../')
import utils.helpers as helpers

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
exp = helpers.get_gene_expression()

##get contact matrix, then edgelist
#also need all the cell names in order to compare with syn
contact_adj = helpers.get_contact_adj()
contact_cells = set(contact_adj.index)
contact_edgelist = helpers.contact_list_from_edge(contact_adj)

###################################################################
##restrict to cells in common in syn and contact
syn_cells = set(syn['pre']).union(set(syn['post']))

common_cells = syn_cells.intersection(contact_cells)
common_cells = common_cells.intersection(set(exp.columns))

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
##set up loss function
##in order to save time copying,
##mostly will pass expression values by reference via exp
##and refer to columns of exp

num_genes = len(exp.index)
#num_genes = 5 ##for testing

##we use the ordering of rows in exp to index genes
gene_of = list(exp.index)
ind_of_gene = { gene : ind for ind,gene in enumerate(gene_of) }

##likewise index the cells
cell_of = list(exp.columns)
ind_of_cell =  { cell : ind for ind,cell in enumerate(cell_of) }

##vector of 1s, useful
ones_row = np.full((1,num_genes), 1.0)
ones_col = np.full((num_genes,1), 1.0)


##a bit more convenient to use the raw numpy array
##each row is a cell
#exp = exp[:num_genes] ##for testing
exp_ind = exp.values.transpose()
##E is more convenient version,
##E[i] is essentially exp_ind[i],
##but exp_ind[i] is a 1D array,
##while E[i] is a 2D array (in row shape (1, num_genes))
E = exp_ind.reshape((len(exp_ind)), 1, num_genes)


##don't need this, doing regression on dual space instead..
##weight, bias matrix to be optimized
w = np.zeros((num_genes, num_genes))
w_rowsum = np.sum(w,1)
w_colsum = np.sum(w,0)
b = np.zeros((num_genes, num_genes))


##the ground truth of the ratio of syn to contact
r = syn_contact[['pre', 'post', 'ratio']]


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


##the data points
##each row in syn_contact, i.e. each directed edge,
##gives one data point
x = [] ##later turn into np array
y = []

##append data point corresponding to edge i->j between cells i,j
##E should use exp_ind
##adds experiment y_val (the observed ratio, r['ratio'])
##adds 1d np array to x, consisting of the dual TODO explain
def append_data(E, i, j, y_val):
	global x, y
	y.append(y_val)
#
	#coeffs of w's
	wdual_j = ones_col.dot(E[j])
	wdual_i = E[i].T.dot(ones_row)
	wdual = wdual_j - wdual_i
	bdual = E[i].T.dot(E[j])
#
	#flatten and append to x
	wdual_flat = wdual.flatten()
	bdual_flat = bdual.flatten()
	x.append(np.append(wdual_flat, bdual_flat))


#for row in range(len(r)):
r.apply(lambda row : append_data(E, ind_of_cell[row['pre']], ind_of_cell[row['post']],
	row['ratio']), axis=1)


xtest = x[:1000]
ytest = y[:1000]
model = LinearRegression().fit(xtest, ytest)
model.score(xtest,ytest)

xtest2 = x[1000:2000]
ytest2 = y[1000:2000]
model2 = LinearRegression().fit(xtest2, ytest2)
model2.score(xtest2,ytest2)

xtest3 = x[2000:3000]
ytest3 = y[2000:3000]
model3 = LinearRegression().fit(xtest3, ytest3)
model3.score(xtest3,ytest3)

###################################################################
###################################################################

