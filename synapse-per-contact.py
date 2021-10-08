import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


######################################################################
##some helpers
##gets the n-th cell in post
##returns empty string if n is out of bounds
def extract_post(n, post):
    arr = post.split(',')
    if len(arr) < n + 1:
        return ''
    else:
        return arr[n]


######################################################################
##load data
contact_sections = pd.read_csv('contact-sections-07-10-2021-175538.csv', index_col=0)
synapse_sections = pd.read_csv('synapse-types-06-10-2021-173407/synapse-sections-06-10-2021-173407.csv', index_col=0)
##use non-robust for now; unclear how to get robust version of contact


######################################################################
##get polyads as multiple synapses (as in old adjacency matrix)
##unclear how to compare otherwise
#get max num post synaptic cells
max_post = max(synapse_sections.apply(lambda row : len(row.post.split(',')),axis = 1))
synapse_separate = pd.DataFrame(columns=['pre','post','sections'])
for n in range(max_post):
    df = synapse_sections.copy()
    df.post = df.post.apply(lambda post : extract_post(n,post))
    df = df[df.post != '']
    synapse_separate = synapse_separate.append(df)

synapse_regrouped = synapse_separate.groupby(["pre","post"]).sum().reset_index()


######################################################################
##contact_sections, flipping pre and post
df = contact_sections.copy()
df.columns = ['post','pre','sections']
df = df.reindex(contact_sections.columns, axis=1)
contact_sections_bi = contact_sections.append(df)


######################################################################
##combine the pre and post columns to make unique keys
##then turn them into dictionaries
##easier for comparison
synapse_dict = dict(zip(map(';'.join, zip(synapse_regrouped.pre, synapse_regrouped.post)), synapse_regrouped.sections))
contact_dict = dict(zip(map(';'.join, zip(contact_sections_bi.pre, contact_sections_bi.post)), contact_sections_bi.sections))
#that's a mouthful; step by step:
#   prepost = zip(synapse_regrouped.pre, synapse_regrouped.post)
#   prepost = map(';'.join, prepost)  ##combines (pre,post) into pre;post
#   synapse_dict = dict(zip(prepost,synapse_regrouped.sections))

contact_keys = set(contact_dict.keys())
synapse_keys = set(synapse_dict.keys())
contact_keys.intersection(synapse_keys)
##hmm it's very few, probably because synapse is from N2Y
##while contact mostly N2U??
