import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

import imp
helpers = imp.load_source('helpers','utils/helpers.py')

######################################################################
# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime('%d-%m-%Y-%H%M%S')
print('now = ' + dt_string)

######################################################################
##load data

syn_cleaned = helpers.get_synapse_list()

#sum total number of sections by pre,post
syn_polyad = syn_cleaned.groupby(['pre','post']).sum().reset_index()

##get synapses that are large ( >2 sections)
#syn_robust = syn_cleaned[syn_cleaned.sections > 2]
#syn_robust_grouped = syn_robust.groupby(['pre','post']).sum().reset_index()

syn_monad = helpers.polyad_to_monad(syn_polyad)

##also get contacts
##same as contact_pixels_adj from synapse-per-contact
##TODO fix that one
contact_adj = helpers.get_contact_adj()
contact_edgelist = helpers.get_contact_edgelist()

######################################################################
#pairwise (should really be called edge-wise)
syn_simscore = helpers.LR_symmetry_pairwise(syn_monad, similarity_score_single, 'simscore')
contact_simscore = helpers.LR_symmetry_pairwise(contact_edgelist, similarity_score_single, 'simscore')
contact_L1simscore = helpers.LR_symmetry_pairwise(contact_edgelist, L1_diff_single, 'L1_simscore')

#drop the cells that have no L/R
syn_simscore_drop = helpers.drop_noLR(syn_simscore)
contact_simscore_drop = helpers.drop_noLR(contact_simscore)
contact_L1simscore_drop = helpers.drop_noLR(contact_L1simscore)

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
adj = helpers.edge_to_adj(syn_polyad)
adj_monad = helpers.edge_to_adj(syn_monad)


fn_simscore = helpers.simscore
fn_L1 = helpers.L1_dist

##compute similarity for polyads
df_out = helpers.LR_symmetry(adj,fn_L1)
df_out.columns = ['pre','L1_dist']
df_out_Adam = helpers.LR_symmetry(adj,fn_simscore)
df_out_Adam.columns = ['pre','sim_score_Adam']

df_out.to_csv('L1_dist_01_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_Adam.to_csv('Adam_dist_01_'+dt_string+'.csv',encoding='utf-8-sig',index=False)

##compute similarity for monads
df_out_monad = helpers.LR_symmetry(adj_monad,fn_L1)
df_out_monad.columns = ['pre','L1_dist']
df_out_monad_Adam = helpers.LR_symmetry(adj_monad,fn_simscore)
df_out_monad_Adam.columns = ['pre','sim_score_Adam']

df_out_monad.to_csv('L1_dist_01_monad_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_monad_Adam.to_csv('Adam_dist_01_monad_'+dt_string+'.csv',encoding='utf-8-sig',index=False)


##compute similarity for contact
df_out_contact = helpers.LR_symmetry(contact_adj,fn_L1)
df_out_contact.columns = ['pre','L1_dist']
df_out_contact_Adam = helpers.LR_symmetry(contact_adj,fn_simscore)
df_out_contact_Adam.columns = ['pre','sim_score_dist']

df_out_contact.to_csv('L1_dist_01_contact_'+dt_string+'.csv',encoding='utf-8-sig',index=False)
df_out_contact_Adam.to_csv('Adam_dist_01_contact_'+dt_string+'.csv',encoding='utf-8-sig',index=False)


######################################################################
## getting all similarity scores of all pairs

syn_adj = helpers.edge_to_adj(syn_monad)

pre = syn_adj.index

fn_simscore = helpers.simscore
fn_L1 = helpers.L1_dist

pairwise_L1 = helpers.pairwise_dist(syn_adj,fn_L1)
pairwise_simscore = helpers.pairwise_dist(syn_adj,fn_simscore)

pairwise_L1.to_csv('pairwise_monad_L1_'+dt_string+'.csv',encoding='utf-8-sig')
pairwise_simscore.to_csv('pairwise_monad_simscore_'+dt_string+'.csv',encoding='utf-8-sig')
#read with
pairwise_L1 = pd.read_csv('similarity/pairwise_monad_L1_27-10-2021-161929.csv',index_col=0)
pairwise_simscore = pd.read_csv('similarity/pairwise_monad_simscore_27-10-2021-161929.csv',index_col=0)

##for synapses, also compare similarity of incoming edges
syn_adj_T = syn_adj.transpose()

pre = syn_adj_T.index

pairwise_T_L1 = helpers.pairwise_dist(syn_adj_T,fn_L1)
pairwise_T_simscore = helpers.pairwise_dist(syn_adj_T,fn_simscore)

pairwise_T_L1.to_csv('pairwise_T_monad_L1_'+dt_string+'.csv',encoding='utf-8-sig')
pairwise_T_simscore.to_csv('pairwise_T_monad_simscore_'+dt_string+'.csv',encoding='utf-8-sig')

