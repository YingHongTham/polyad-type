import pandas as pd
import numpy as np
from helpers import *

######################################################################
##TODO refactor this, copied from get-polyad-types.py

#get synapse list
syn = pd.read_csv('SI-3-Synapse-lists-male.csv')
syn = syn.query('EM_series=="N2Y" & type=="chemical"')
syn = syn[['pre','post','sections']]

#clean the pre neurons
syn['pre'] = syn['pre'].apply(clean_neuron_name)
#clean the post neurons
syn['post'] = syn['post'].apply(clean_post)

#after cleaning, empty string means unknown, filter them out
syn_cleaned = syn[(syn.pre != '') & (syn.post != '')]

#sum total number of sections by pre,post
syn_grouped = syn_cleaned.groupby(['pre','post']).sum().reset_index()

##get synapses that are large ( >2 sections)
syn_robust = syn_cleaned[syn_cleaned.sections > 2]
syn_robust_grouped = syn_robust.groupby(['pre','post']).sum().reset_index()


######################################################################

syn_cleaned_flippedLR = syn_cleaned.copy()
syn_cleaned_flippedLR['pre'] = syn_cleaned_flippedLR['pre'].apply(flipLR)
syn_cleaned_flippedLR['post'] = syn_cleaned_flippedLR['post'].apply(lambda x : apply_fn_post(flipLR, x))

syn_grouped_flippedLR = syn_grouped.copy()
syn_grouped_flippedLR['pre'] = syn_grouped_flippedLR['pre'].apply(flipLR)
syn_grouped_flippedLR['post'] = syn_grouped_flippedLR['post'].apply(lambda x : apply_fn_post(flipLR, x))

##to adj matrix
adj_orig = pd.crosstab(index=syn_grouped['pre'],columns=syn_grouped['post'],values=syn_grouped['sections'],aggfunc='sum').fillna(0)
adj_flip = pd.crosstab(index=syn_grouped_flippedLR['pre'],columns=syn_grouped_flippedLR['post'],values=syn_grouped_flippedLR['sections'],aggfunc='sum').fillna(0)


##ensure they have the same set of columns
columns_orig = set(adj_orig.columns)
columns_flip = set(adj_flip.columns)

for col in columns_flip.difference(columns_orig):
	adj_orig[col] = 0

for col in columns_orig.difference(columns_flip):
	adj_flip[col] = 0

##ensure column names are in the same order
adj_orig = adj_orig[sorted(adj_orig.columns)]
adj_flip = adj_flip[sorted(adj_flip.columns)]

adj_orig_pre = set(adj_orig.index)
adj_flip_pre = set(adj_flip.index)
common_pre = adj_orig_pre.intersection(adj_flip_pre)


def L1_diff_single(a,b) : return abs(a - b)

L1_dist = pd.DataFrame((x, total_similarity_score(L1_diff_single, adj_orig.loc[x],adj_flip.loc[x])) for x in common_pre)
Adam_dist = pd.DataFrame((x, total_similarity_score(similarity_score_single, adj_orig.loc[x],adj_flip.loc[x])) for x in common_pre)

L1_dist.columns = ['pre','L1_dist']
Adam_dist.columns = ['pre','similarity']



##also do contact matrix

##same as contact_pixels_adj from synapse-per-contact
contact_adj = pd.read_csv('N2Y-PAG-contact-matrix.csv',index_col=0)
contact_adj = contact_adj.fillna(0)
cols = contact_adj.columns
rows = contact_adj.index

newcols = map(flipLR, list(cols))
newrows = map(flipLR, list(rows))


##rename with L/R flipped
contact_flipped = contact_adj.copy()
contact_flipped.columns = newcols
contact_flipped.index = newrows
##then reorder
contact_flipped = contact_flipped[cols]
contact_flipped = contact_flipped.reindex(rows)



#adj_orig = adj_orig.sort_index()
#adj_flip = adj_flip.sort_index()


##TODO next: compare rows of adj_orig, adj_flip
#sum(abs(adj_orig.iloc[0] - adj_flip.iloc[0]))
