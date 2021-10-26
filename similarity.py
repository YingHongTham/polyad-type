import pandas as pd
import numpy as np
from helpers import *
from datetime import datetime
import matplotlib.pyplot as plt

######################################################################
# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime('%d-%m-%Y-%H%M%S')
print('now = ' + dt_string)

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
#syn_robust = syn_cleaned[syn_cleaned.sections > 2]
#syn_robust_grouped = syn_robust.groupby(['pre','post']).sum().reset_index()

syn_monad = polyad_to_monad(syn_grouped)

##also get contacts
##same as contact_pixels_adj from synapse-per-contact
contact_adj = pd.read_csv('N2Y-PAG-contact-matrix.csv',index_col=0)
contact_adj = contact_adj.fillna(0)

contact_edgelist = contact_adj.stack().reset_index()
contact_edgelist.columns = ['pre','post','sections']
contact_edgelist = contact_edgelist[contact_edgelist.sections != 0]

######################################################################
#pairwise
syn_simscore = LR_symmetry_pairwise(syn_monad, similarity_score_single, 'simscore')
contact_simscore = LR_symmetry_pairwise(contact_edgelist, similarity_score_single, 'simscore')
contact_L1simscore = LR_symmetry_pairwise(contact_edgelist, L1_diff_single, 'L1_simscore')

#drop the pairs that have no L/R
syn_simscore_drop = drop_noLR(syn_simscore)
contact_simscore_drop = drop_noLR(contact_simscore)
contact_L1simscore_drop = drop_noLR(contact_L1simscore)

#do fraction of simscore out of max
def fraction_simscore(row):
	return row['simscore'] / max(row['sections'],row['sections_flipped'])

syn_simscore_drop['simscore_fraction'] = syn_simscore_drop.apply(fraction_simscore,axis=1)
syn_simscore_drop = syn_simscore_drop[['pre','post','simscore','simscore_fraction','sections','sections_flipped']]
contact_simscore_drop['simscore_fraction'] = contact_simscore_drop.apply(fraction_simscore,axis=1)
contact_simscore_drop = contact_simscore_drop[['pre','post','simscore','simscore_fraction','sections','sections_flipped']]

#also do fraction of L1 diff out of max for contact
def fraction_diff(row):
	return row['L1_simscore'] / max(row['sections'],row['sections_flipped'])

contact_L1simscore_drop['diff_fraction'] = contact_L1simscore_drop.apply(fraction_diff,axis=1)
contact_L1simscore_drop = contact_L1simscore_drop[['pre','post','L1_simscore','diff_fraction','sections','sections_flipped']]

#save to file
syn_simscore_drop.to_csv('synapse_simscore_pairwise_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
contact_simscore_drop.to_csv('contact_simscore_pairwise_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
contact_L1simscore_drop.to_csv('contact_L1simscore_pairwise_'+dt_string+'.csv',encoding='utf-8-sig',index=False)

plt.hist(contact_simscore_drop['simscore_fraction'])
plt.hist(contact_L1simscore_drop['diff_fraction'])

######################################################################
#"global diff"?

##to adj matrix
adj = edge_to_adj(syn_grouped)
adj_monad = edge_to_adj(syn_monad)

##functions used to compare two weights
def L1_diff_single(a,b):
	return abs(a - b)

fn_simscore = total_fn(similarity_score_single)
fn_L1 = total_fn(L1_diff_single)

##compute similarity for polyads
df_out = LR_symmetry(adj,fn_L1)
df_out.columns = ['pre','L1_dist']
df_out_Adam = LR_symmetry(adj,fn_simscore)
df_out_Adam.columns = ['pre','sim_score_Adam']

df_out.to_csv('L1_dist_01_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_Adam.to_csv('Adam_dist_01_'+dt_string+'.csv',encoding='utf-8-sig',index=False)

##compute similarity for monads
df_out_monad = LR_symmetry(adj_monad,fn_L1)
df_out_monad.columns = ['pre','L1_dist']
df_out_monad_Adam = LR_symmetry(adj_monad,fn_simscore)
df_out_monad_Adam.columns = ['pre','sim_score_Adam']

df_out_monad.to_csv('L1_dist_01_monad_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_monad_Adam.to_csv('Adam_dist_01_monad_'+dt_string+'.csv',encoding='utf-8-sig',index=False)


##compute similarity for contact
df_out_contact = LR_symmetry(contact_adj,fn_L1)
df_out_contact.columns = ['pre','L1_dist']
df_out_contact_Adam = LR_symmetry(contact_adj,fn_simscore)
df_out_contact_Adam.columns = ['pre','sim_score_dist']

df_out_contact.to_csv('L1_dist_01_contact_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_contact_Adam.to_csv('Adam_dist_01_contact_'+dt_string+'.csv',encoding='utf-8-sig',index=False)

