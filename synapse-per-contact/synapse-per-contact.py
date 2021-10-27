import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


######################################################################
##load data
contact_pixels_adj = pd.read_csv('N2Y-PAG-contact-matrix.csv',index_col=0)
#this is from get-contact-matrix-Brittin-N2U.py
#contact_sections = pd.read_csv('contact-sections-07-10-2021-175538.csv', index_col=0)
synapse_sections = pd.read_csv('synapse-types-06-10-2021-173407/synapse-sections-06-10-2021-173407.csv', index_col=0)
synapse_sections_robust = pd.read_csv('synapse-types-06-10-2021-173407/synapse-sections-robust-06-10-2021-173407.csv', index_col=0)

##use non-robust for now; unclear how to get robust version of contact

######################################################################
##convert contact matrix to edgelist
contact_pixels = contact_pixels_adj.stack().reset_index()
contact_pixels.columns = ['pre','post','pixels']
##it gives both directions,e.g.
##      pre  post        pixels
##0  AS10  AVAL  14998.156804
##102  AVAL  AS10  14998.156804

######################################################################
###USED FOR N2U
###contact_sections, flipping pre and post
#df = contact_sections.copy()
#df.columns = ['post','pre','sections']
#df = df.reindex(contact_sections.columns, axis=1)
#contact_sections_bi = contact_sections.append(df)


######################################################################
##treat polyads as multiple monads (as in old adjacency matrix)
##(unclear how to compare with contact matrix otherwise)
synapse_monad = polyad_to_monad(synapse_sections)


######################################################################
##combine the pre and post columns to make unique keys
##then turn them into dictionaries
##easier for comparison
synapse_dict = dict(zip(map(';'.join, zip(synapse_monad.pre, synapse_monad.post)), synapse_monad.sections))
contact_dict = dict(zip(map(';'.join, zip(contact_pixels.pre, contact_pixels.post)), contact_pixels.pixels))
##contact_dict = dict(zip(map(';'.join, zip(contact_sections_bi.pre, contact_sections_bi.post)), contact_sections_bi.sections))
#that's a mouthful; step by step:
#   prepost = zip(synapse_monad.pre, synapse_monad.post)
#   prepost = map(';'.join, prepost)  ##combines (pre,post) into pre;post
#   synapse_dict = dict(zip(prepost,synapse_monad.sections))

contact_keys = set(contact_dict.keys())
synapse_keys = set(synapse_dict.keys())
contact_keys.intersection(synapse_keys)
##lots of cells are not in contact;
##need to weed out muscle cells

######################################################################
##analyze self-contact / autapse
synapse_self = synapse_monad[synapse_monad.pre == synapse_monad.post]
contact_self = contact_pixels[contact_pixels.pre == contact_pixels.post]

synapse_self_cells = set(synapse_self['pre'])
contact_self_cells = set(contact_self['pre'])

##TODO very strange! some cells synapse but contact is not recorded!
## synapse_self_cells.difference(contact_self_cells)
##set(['SDQL', 'PVM', 'PDER', 'R7AL'])
contact_no_synapse_cells = contact_self_cells.difference(synapse_self_cells)
contact_no_synapse = contact_self[contact_self.pre.isin(contact_no_synapse_cells)]

common_cells = synapse_self_cells.intersection(contact_self_cells)
synapse_common = synapse_self[synapse_self.pre.isin(common_cells)]
contact_common = contact_self[contact_self.pre.isin(common_cells)]
