import pandas as pd
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
from helpers import *

######################################################################
# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime('%d-%m-%Y-%H%M%S')
print('now = ' + dt_string)

######################################################################
#TODO merge left and right neurons


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
## compare the amount of robust sections to all sections

syn_compare = syn_cleaned.copy()
syn_compare['robust_sections'] = syn_cleaned['sections'].apply(lambda x : x * (x > 2))
syn_compare_grouped = syn_compare.groupby(['pre','post']).sum().reset_index()
syn_compare_grouped['robust_ratio'] = syn_compare_grouped.robust_sections / syn_compare_grouped.sections

##plot; most are 0, but a sizeable 1, and then the rest are sprinkled
##in the middle
ratios = syn_compare_grouped['robust_ratio']
numbins = 100
fig, ax = plt.subplots(1,1)
ax.hist(ratios,bins=numbins)
#plt.show()


######################################################################
##get distribution of synapse sizes (#sections) by (pre) neuron
##see if there are synapse types
synapse_sections_dist = syn_cleaned.groupby('pre')['sections'].apply(list).reset_index(name='section_dist')
synapse_sections_dist['pvariance'] = synapse_sections_dist.section_dist.apply(statistics.pvariance)

##get cell with largest variation (it's PHAR)
max_var_ind = np.argmax(synapse_sections_dist['pvariance'])

##PHAR has largest variance
PHAR_dist = synapse_sections_dist.iloc[max_var_ind].section_dist
##equivalently: synapse_sections_dist[synapse_sections_dist.pre == 'PHAR'].section_dist.values[0]


dir_histograms = 'synapse-size-histograms-'+dt_string
#os.mkdir(dir_histograms)

def save_graph_by_row(row):
	mylist = row['section_dist']
	if len(mylist) < 60:
		return
	cellname = row['pre']
	numbins = max(mylist) - min(mylist) + 1
	fig, ax = plt.subplots(1,1)
	ax.hist(mylist,bins=numbins)
	ax.set_title('Distribution of synapse sizes of '+cellname+' (as pre)')
	ax.set_xlabel('Synapse size (=number of sections)')
	ax.set_ylabel('frequency')
	plt.savefig(dir_histograms+'/'+cellname+'.png')

synapse_sections_dist.apply(save_graph_by_row, axis=1)

##HOA has (by far) the biggest synapse
##>>> syn_cleaned[(syn_cleaned.pre == 'HOA') & (syn_cleaned.sections > 10)]
##      pre            post  sections
##927   HOA       PHCR,R8AR        15
##1021  HOA  CP02,PHCL,R8AL        14
##2548  HOA  CP02,PHCL,R8AL        17
##2562  HOA            PHCL        41


######################################################################
######################################################################
##save the tables
######################################################################
######################################################################

syn_grouped.to_csv('synapse-sections-'+dt_string+'.csv',encoding='utf-8-sig')
syn_robust_grouped.to_csv('synapse-sections-robust-'+dt_string+'.csv',encoding='utf-8-sig')

##save as adjacency table
#(see https://medium.com/@yangdustin5/quick-guide-to-pandas-pivot-table-crosstab-40798b33e367)
def save_df_as_adj(df,filename):
    adj_table = pd.crosstab(index=df['pre'],columns=df['post'],values=df['sections'],aggfunc='sum')
    #adj_table.to_csv('adj_table_03.csv')
    adj_table.to_csv(filename+'.csv',encoding='utf-8-sig')
    adj_table_transpose = adj_table.transpose()
    adj_table_transpose.to_csv(filename+'-transpose.csv',encoding='utf-8-sig')


save_df_as_adj(syn_grouped,'adj-table-'+dt_string)
save_df_as_adj(syn_grouped,'adj-table-robust-'+dt_string)
#also save the outdegrees (number of types)
syn_count.to_csv('outdegrees-'+dt_string+'.csv',encoding='utf-8-sig')
syn_count_robust.to_csv('outdegrees-robust-'+dt_string+'.csv',encoding='utf-8-sig')

######################################################################
##graphing time

##count number of post-synaptic types for each neuron
syn_count = syn_grouped.groupby(['pre']).size().reset_index(name='num_post_types')
##graph it
count_types = syn_count['num_post_types']
numbins = count_types.max() - count_types.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(count_types,bins=numbins)
ax.set_title('Distribution of synapse outdegrees (N2Y)')
ax.set_xlabel('outdegree (number of post-synaptic types)')
ax.set_ylabel('frequency')
#plt.show()
#save the plot
plt.savefig('histogram-outdegree-N2Y-'+dt_string+'.png')

##count number of post-synaptic types for each neuron
##but for robust
syn_count_robust = syn_robust_grouped.groupby(['pre']).size().reset_index(name='num_post_types')
##graph it
count_types_robust = syn_count_robust['num_post_types']
numbins = count_types_robust.max() - count_types_robust.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(count_types_robust,bins=numbins)
ax.set_title('Distribution of synapse outdegrees (N2Y) (robust)')
ax.set_xlabel('outdegree (number of post-synaptic types)')
ax.set_ylabel('frequency')
#plt.show()
#save the plot
plt.savefig('histogram-outdegree-N2Y-robust-'+dt_string+'.png')
################################################################
## experimental
################################################################

##get all the post synaptic stuff
##trying to sort out which values to keep/discard..
##there are the entries like [CP01],
##   which are unsure connections so we discard
## but there's also stuff like unk5590...
## even unk4733[hyp], contin5934[AVF][AVH][AVJ][PVN]...
#all_post = df['post'].tolist()
#all_post = [x.split(',') for x in all_post]
#all_post = sum(all_post,[])
#all_post = sorted(list(set(all_post)))
#all_pre = sorted(list(set(df['pre'].tolist())))


## way to get new dataframe;
## axis=1 means apply to rows
## 'expand' makes it 3 columns instead of 1
#df.apply(lambda row : [row.pre, sort_post(row.post), row.sections], axis=1, result_type='expand').head()


