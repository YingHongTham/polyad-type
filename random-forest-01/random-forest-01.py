import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import export_graphviz
import pickle
import pydot
# need to $ sudo apt-get install graphviz (check $ dot -V)
from sklearn.linear_model import LinearRegression

import sys
sys.path.append('../')
import utils.helpers as helpers

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

common_cells = exp_cells.intersection(syn_cells)
common_cells = common_cells.intersection(contact_cells)

syn = syn[(syn['pre'].isin(common_cells)) & (syn['post'].isin(common_cells))]
contact = contact_edgelist[(contact_edgelist['pre'].isin(common_cells)) & (contact_edgelist['post'].isin(common_cells))]
exp = exp[common_cells]

###################################################################
## as in analyse-expression-*,
## each edge in synapse graph gives one data point
## the features are given by
## -genes of pre
## -genes of post
## -pair of genes, value = 1 if both expressed, = 0 else

features = syn[['pre','post']]
labels = syn.sections

for gene in exp.index:
	features['pre-'+gene] = features.pre.apply(lambda cell : exp.loc[gene][cell])

for gene in exp.index:
	features['post-'+gene] = features.post.apply(lambda cell : exp.loc[gene][cell])

## this takes too long
#for pregene in exp.index:
#	for postgene in exp.index:
#		features['pre-post-'+pregene+'-'+postgene] = features.apply(lambda r : r['pre-'+pregene] * r['post-'+postgene], axis=1)


features = features.drop(['pre','post'], axis=1)

features_np = np.array(features)
labels_np = np.array(labels)

rf = RandomForestRegressor(n_estimators=1000, random_state=42)

rf.fit(features_np, labels_np)

##save
filename = 'model_02.sav'
pickle.dump(rf, open(filename, 'wb'))

##load
#rf = pickle.load(open(filename, 'rb'))
#feature_list = np.array(('pre-' + exp.index).append('post-' + exp.index))

###################################################################
## feature extraction by importance

importances = list(rf.feature_importances_)

feature_list = list(features.columns)

feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]

# list of pairs, (feature, importance)
feature_importances = sorted(feature_importances, key = lambda x : x[1], reverse = True)

for pair in feature_importances:
	print('Variable: {:20} Importance: {}'.format(*pair))

###################################################################
# get top 20/50 most important genes
top20 = [p[0] for p in feature_importances[:20]]
top50 = [p[0] for p in feature_importances[:50]]

# separates 'pre-gene' to ['pre',gene] (same with post)
def separate_prepost_gene(p_gene):
	ind = p_gene.find('-')
	return [p_gene[:ind], p_gene[ind+1:]]

top20sep = list(map(separate_prepost_gene, top20))

ind_of_p_gene = { p_gene : ind  for ind,p_gene in enumerate(top20) }
ind_of_p_gene_50 = { p_gene : ind  for ind,p_gene in enumerate(top50) }

###################################################################
## do linear regression based on top20 features

#prepare features

feat_lin = syn[['pre','post']]

#for p_gene in top20:
for p_gene in top50:
	p, gene = separate_prepost_gene(p_gene)
	feat_lin[p_gene] = feat_lin[p].apply(lambda cell : exp.loc[gene][cell])

feat_lin = feat_lin.drop(['pre','post'], axis=1)

#each row is one data point
feat_lin_np = np.array(feat_lin)

## from before
#labels = syn.sections
#labels_np = np.array(labels)


lin_model = LinearRegression().fit(feat_lin, labels_np)
lin_model.score(feat_lin, labels_np)

lin_model_all = LinearRegression().fit(features, labels_np)
lin_model_all.score(features, labels_np)

###################################################################
## visualizing a tree in the rf

tree = rf.estimators_[5]

export_graphviz(tree, out_file='tree.dot', feature_names=feature_list, rounded=True, precision=1)

(graph, ) = pydot.graph_from_dot_file('tree.dot')

graph.write_png('tree.png')

###################################################################
graph = pydot.Dot(graph_type='graph')
for i in range(3):
	edge = pydot.Edge('king', 'lord%d' % i)
	graph.add_edge(edge)

vassal_num = 0
for i in range(3):
	for j in range(2):
		edge = pydot.Edge('lord%d' % i, 'vassal%d' % vassal_num)
		graph.add_edge(edge)
		vassal_num += 1

graph.write_png('example_graph_1.png')
###################################################################
###################################################################
###################################################################

