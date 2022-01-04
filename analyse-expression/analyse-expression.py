import pandas as pd
import networkx as nx
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
##bokeh visualizing stuff
from bokeh.io import show, save
from bokeh.models import Range1d, Circle, ColumnDataSource, MultiLine, PointDrawTool
from bokeh.plotting import figure
from bokeh.models.graphs import from_networkx

##importing local package helpers
##need the sys stuff because it's located in a different folder
##append parent folder to sys.path so it knows to look there for utils folder
import sys
sys.path.append('../')
import utils.helpers as helpers

###################################################################

###path to folder containing data
#filepath_datasource = '../data-source/'
#
##exp is dataframe, row = gene, col = cells
#exp = pd.read_csv(filepath_datasource+'Expression-matrix-Jan-2020.csv',index_col=0)
#exp = exp.fillna(0) ##empty entries interpret as 0
#
###drop the last 5 columns ('Unnamed: 23*'); they're all 0
#exp = exp[exp.columns[:-5]]


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
##noting the asymmetry in chemical synapse formation;
##some pairs, e.g. PCAR,PVX and PCAL,PVX
##have lots of contact, and many synapses in one direction,
##but few synapses in the other direction
##may be good place to look in gene expression TODO
##
##plot graph of asymmetry in number of synapse vs contact
#TODO
#syn_contact.sort_values('pixels', ascending=False)


###################################################################
##extract full subgraph of connect/contact graph
##that expresses a particular (pair of?) gene(s)


##return cells that express a particular genes
##return type: set
def get_cells_expressing_gene(gene, exp):
	return set(filter(lambda c : exp.loc[gene][c] != 0, exp.columns))


##return dataframe with rows from syn
##whose both pre,post cells express gene
def gene_subgraph(gene, syn):
	cells = get_cells_expressing_gene(gene, exp)
	return syn[(syn['pre'].isin(cells)) & (syn['post'].isin(cells))]


##testing
#gene = 'bam-2'
#cells = get_cells_expressing_gene(gene, exp)
#subgr = gene_subgraph(gene, syn)

###################################################################
##turn dataframe (list) of synapses into graph (networkx)

##networkx receives weighted edges as a tuple (v,w,{'attr' : val})
##turns rows of df into (list of) such tuples
def df_to_ebunch(syn):
	ebunch = syn.apply(lambda row : (row.pre, row.post, {'sections':row.sections}), axis=1)
	return list(ebunch)


syn_nx = nx.DiGraph()
syn_nx.add_edges_from(df_to_ebunch(syn))

#syn_subgr_nx = syn_nx.subgraph(cells)


##to draw with plt
nx.draw(syn_nx)
plt.show()

###################################################################
##visualizing with Bokeh
##stolen from https://melaniewalsh.github.io/Intro-Cultural-Analytics/06-Network-Analysis/02-Making-Network-Viz-with-Bokeh.html

#syn_nx is some nx graph, should be subgraph corresponding to gene
def bokeh_visualization(syn_nx, gene='any'):
	title = 'Chemical Synapses, cells expressing '+gene
#
	#Establish which categories will appear when hovering over each node
	#TODO do for edges too, show weights
	HOVER_TOOLTIPS = [('cell', '@index')]
#
	#Create a plot — set dimensions, toolbar, and title
	#only initializes
	plot = figure(tooltips = HOVER_TOOLTIPS,
								tools="pan,wheel_zoom,save,reset",
								active_scroll='wheel_zoom',
								x_range=Range1d(-10.1, 10.1),
								y_range=Range1d(-10.1, 10.1), title=title)
	#
	#Create a network graph object with spring layout
	# https://networkx.github.io/documentation/networkx-1.9/reference/generated/networkx.drawing.layout.spring_layout.html
	network_graph = from_networkx(syn_nx, nx.spring_layout, scale=10, center=(0, 0))
#
	#Set node size and color
	network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')
#
	#Set edge opacity and width
	network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
#
	##allow drag points.. TODO get it to work
	#plot.add_tools(PointDrawTool(renderers = [network_graph.node_renderer], empty_value = 'black'))
#
	#Add network graph to the plot
	plot.renderers.append(network_graph)
#
	#show(plot)
	save(plot, filename='subgraph_'+gene+'.html')


##do all the genes
for gene in exp.index:
	cells = get_cells_expressing_gene(gene, exp)
	syn_subgr_nx = syn_nx.subgraph(cells)
	bokeh_visualization(syn_subgr_nx, gene)



##now compare synapse to contact;
##make graph visualization with sliding bar for some parameter
##show edges with weight above parameter
##edge weight by ratio of synapse / contact

##TODO hmm perhaps easier to compare contact with gap junction
##since both are symmetric and only pairwise


