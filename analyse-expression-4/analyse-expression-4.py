##4th attempt: compare exp of cells that have similar neighbourhoods


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

######################################################################
# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime('%d-%m-%Y-%H%M%S')
print('now = ' + dt_string)


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
contact_edgelist = helpers.contact_list_from_edge(contact_adj)

###################################################################
##restrict to cells in common

syn_cells = set(syn['pre']).union(set(syn['post']))
contact_cells = set(contact_adj.index)
exp_cells = set(exp.columns)

#common_cells = exp_cells.intersection(syn_cells)
#common_cells = common_cells.intersection(contact_cells)

###################################################################
###restrict syn, contact, gene to common cells
#syn = syn[(syn['pre'].isin(common_cells)) & (syn['post'].isin(common_cells))]
#contact = contact_edgelist[(contact_edgelist['pre'].isin(common_cells)) & (contact_edgelist['post'].isin(common_cells))]
#exp = exp[common_cells]


###################################################################
##get the similarity/difference of neighborhoods of all pairs of cells
##

syn_adj = helpers.edge_to_adj(syn)

#get pairwise distance; this takes a long time, so load from saved
load_from_file = True
if (load_from_file):
	syn_pairwise_dist = pd.read_csv('syn_pairwise_dist_18-01-2022-171318.csv',index_col=0)
	exp_pairwise_dist = pd.read_csv('exp_pairwise_dist_18-01-2022-171318.csv',index_col=0)
else:
	syn_pairwise_dist = helpers.pairwise_dist(syn_adj, helpers.L1_dist)
	con_pairwise_dist = helpers.pairwise_dist(contact_adj, helpers.L1_dist)
	exp_pairwise_dist = helpers.pairwise_dist_col(exp, helpers.L1_dist)
	syn_pairwise_dist.to_csv('syn_pairwise_dist_'+dt_string+'.csv',encoding='utf-8-sig')
	exp_pairwise_dist.to_csv('exp_pairwise_dist_'+dt_string+'.csv',encoding='utf-8-sig')
	con_pairwise_dist.to_csv('con_pairwise_dist_'+dt_string+'.csv',encoding='utf-8-sig')

test_exp = exp[['LUAR', 'DVC']]
test_dist= helpers.pairwise_dist_col(test_exp, helpers.L1_dist)
testlist = ['DVC','LUAR']
tt = test_dist[testlist].loc[testlist]

##get common cells
#here only use pre cells in syn
syn_pre_cells = set(syn_adj.index)
common_cells = exp_cells.intersection(syn_pre_cells)
#common_cells = exp_cells.intersection(contact_cells)
common_list = list(common_cells)

syn_dist_common = syn_pairwise_dist[common_list].loc[common_list]
con_dist_common = con_pairwise_dist[common_list].loc[common_list]
exp_dist_common = exp_pairwise_dist[common_list].loc[common_list]

###################################################################
##turn into edgelist
##combine into one adj matrix first..

Z = pd.DataFrame(index=common_list, columns=common_list)

for x in common_list:
	for y in common_list:
		Z[x][y] = (syn_pairwise_dist[x][y], exp_pairwise_dist[x][y])

Zlist = helpers.edgelist_from_adj(Z,row_col_val=['i','j','dist'],clear_zeros=False)

Zlist['syn_dist'] = Zlist['dist'].apply(lambda r : r[0])
Zlist['exp_dist'] = Zlist['dist'].apply(lambda r : r[1])

df = Zlist[Zlist['syn_dist'] > 1000]
x = np.array(df['syn_dist'])
y = np.array(df['exp_dist'])

##this doesn't seem to work
#def zipUp(x,y):
#	return (x,y)
#
#zz = np.vectorize(np.vectorize(zipUp))
#
#Z = con_pairwise_dist.combine(exp_pairwise_dist, np.minimum)

###################################################################

#syn_dist_edgelist = helpers.edgelist_from_adj(syn_dist_common,row_col_val=['i','j','syn_dist'],clear_zeros=False)
con_dist_edgelist = helpers.edgelist_from_adj(con_dist_common,row_col_val=['i','j','con_dist'],clear_zeros=False)
exp_dist_edgelist = helpers.edgelist_from_adj(exp_dist_common,row_col_val=['i','j','exp_dist'],clear_zeros=False)

ttt = helpers.edgelist_from_adj(exp_pairwise_dist,row_col_val=['i','j','exp_dist'],clear_zeros=False)


dist_edgelist = syn_dist_edgelist.copy()
#dist_edgelist = con_dist_edgelist.copy()
dist_edgelist['exp_dist'] = exp_dist_edgelist['exp_dist']

#TODO for some reason exp_dist_edgelist has NaN values
dist_edgelist = dist_edgelist[~dist_edgelist['exp_dist'].isnull()]

###################################################################
##plot/fit line

#x = dist_edgelist['syn_dist']
#df = dist_edgelist
#df = df[(df['con_dist'] < 1e7) & (df['con_dist'] > 1e6)]
df = dist_edgelist.copy()
df = df[df['syn_dist'] == 0]
df = df[(df['syn_dist'] < 1e4) & (df['syn_dist'] > 5)]
x = df['syn_dist'].values
y = df['exp_dist'].values

#plot
plt.scatter(x,y)


##best fit line
from sklearn.linear_model import LinearRegression

#fit() takes x = array of arrays, y = array
xx = x.reshape(-1,1)

model = LinearRegression().fit(xx, y)
model.score(xx,y)


###################################################################
##selecting ray neurons
##starts with 'R'

def isRayNeuron(cell):
	return cell[0] == 'R'

def isR1A(cell):
	return cell[0:3] == 'R1A'

def isRayNeuronLabel(cell, neuronLabel):
	return cell[0:len(neuronLabel)] == neuronLabel

##use example:
#syn[syn.pre.apply(lambda cell: isRayNeuronLabel(cell, 'R4BR'))]

###################################################################
