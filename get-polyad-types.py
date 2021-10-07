import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
##note to self: need to apt-get install python-tk for matplotlib

# datetime object containing current date and time
now = datetime.now()

print("now =", now)

# dd/mm/YY H:M:S
dt_string = now.strftime("%d-%m-%Y-%H%M%S")
###############################################
## helper functions for manipulating the rows

#clean the post neurons list
#splits comma-sep string into a python list
#perform clean_neuron_name on each entry
#sort, then recombine into comma-sep string
#e.g. HOA,[AVG] -> AVG,HOA
#note if everything is removed, then returns empty string
##also, if we have a dyad, and one post-syn is removed by clean,
##then we end up with a monad. potentially lose info
def clean_post(post):
    post = post.replace(" ","")
    post_list = post.split(",")
    post_list = map(clean_neuron_name, post_list)
    post_list = filter(lambda x : x != "" , post_list)
    post_list = sorted(post_list)
    post = ",".join(post_list)
    return post

def clean_neuron_name(neuron):
    ##remove spaces
    neuron = neuron.replace(" ","")
    ##remove square brackets
    neuron = neuron.replace("[","")
    neuron = neuron.replace("]","")
    neuron = convert_name(neuron)
    if good_cellname(neuron):
        return neuron
    else:
        return ""


convert_dict = {
    "PVS" : "PVPR",
    "PVU" : "PVPL",
    "adp" : "mu_anal",
    "sph" : "mu_sphincter",
    "intL" : "mu_intL",
    "intR" : "mu_intR"
}
## convert name of cell that appears in SI 3 Synapse List
## to that appearing in SI 4 Celllist
def convert_name(name):
    global convert_dict
    if name in convert_dict:
        return convert_dict[name]
    else:
        return name


## check if cellname is in cell list SI4
## expect name to have been cleaned and converted
##(used for sanity check: expect no weird_names after running
##weird_names are neurons that appear in the synapse table
##but not in the cell table SI 4)
weird_names = []
def good_cellname(name):
    global weird_names
    global celllist
    if good_cellname_old(name) != (name in celllist):
        print(name + "; in SI4? : " + str(name in celllist))
        weird_names.append(name)
    return name in celllist

##detecting bad names by manually eliminating
##only for diagnosing problem
def good_cellname_old(name):
    return name == "" and name[0:3] != "unk" and name[0:3] != "obj" and name[0:6] != "contin" and name[0:4] != "frag"


######################################################################

#get list of neurons (help with pruning bad data)
celllist = pd.read_csv("male_celllist.csv", header=None, comment='#')
#remove unexpected spaces...
celllist = celllist[0].apply(lambda x : x.replace(" ",""))
celllist = celllist.tolist()
##note to self: DataFrame object doesn't have tolist(),
##but a column object does..


#TODO get the contact data from the mysql
#TODO merge left and right neurons


#get synapse list
syn = pd.read_csv('SI-3-Synapse-lists-male.csv')
syn = syn.query('EM_series=="N2Y" & type=="chemical"')
syn = syn[["pre","post","sections"]]

#clean the pre neurons
syn["pre"] = syn["pre"].apply(clean_neuron_name)
#clean the post neurons
syn["post"] = syn["post"].apply(clean_post)

#after cleaning, empty string means unknown
syn_cleaned = syn[(syn.pre != "") & (syn.post != "")]

#sum total number of sections by pre,post
#(like mysql query)
#groupby returns DataFrame, has weird indexing (fix with reset_index)
syn_grouped = syn_cleaned.groupby(["pre","post"]).sum().reset_index()

##get synapses that are large ( >2 sections)
syn_robust = syn_cleaned[syn_cleaned.sections > 2]
syn_robust_grouped = syn_robust.groupby(["pre","post"]).sum().reset_index()
#syn_robust = syn_grouped[syn_grouped.sections > 2]


######################################################################
##graphing time

##count number of post-synaptic types for each neuron
syn_count = syn_grouped.groupby(["pre"]).size().reset_index(name="num_post_types")
##graph it
count_types = syn_count["num_post_types"]
numbins = count_types.max() - count_types.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(count_types,bins=numbins)
ax.set_title("Distribution of synapse outdegrees (N2Y)")
ax.set_xlabel("outdegree (number of post-synaptic types)")
ax.set_ylabel("frequency")
#plt.show()
#save the plot
plt.savefig("histogram-outdegree-N2Y-"+dt_string+".png")

##count number of post-synaptic types for each neuron
##but for robust
syn_count_robust = syn_robust_grouped.groupby(["pre"]).size().reset_index(name="num_post_types")
##graph it
count_types_robust = syn_count_robust["num_post_types"]
numbins = count_types_robust.max() - count_types_robust.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(count_types_robust,bins=numbins)
ax.set_title("Distribution of synapse outdegrees (N2Y) (robust)")
ax.set_xlabel("outdegree (number of post-synaptic types)")
ax.set_ylabel("frequency")
#plt.show()
#save the plot
plt.savefig("histogram-outdegree-N2Y-robust-"+dt_string+".png")

##turn "edge list" into adjacency table before saving
#(see https://medium.com/@yangdustin5/quick-guide-to-pandas-pivot-table-crosstab-40798b33e367)
def save_df(df,filename):
    adj_table = pd.crosstab(index=df["pre"],columns=df["post"],values=df["sections"],aggfunc="sum")
    #adj_table.to_csv('adj_table_03.csv')
    adj_table.to_csv(filename+'.csv',encoding='utf-8-sig')
    adj_table_transpose = adj_table.transpose()
    adj_table_transpose.to_csv(filename+'-transpose.csv',encoding='utf-8-sig')

######################################################################
##save the tables

save_df(syn_grouped,"adj-table-"+dt_string)
save_df(syn_grouped,"adj-table-robust-"+dt_string)
#also save the outdegrees (number of types)
syn_count.to_csv("outdegrees-"+dt_string+".csv",encoding='utf-8-sig')
syn_count_robust.to_csv("outdegrees-robust-"+dt_string+".csv",encoding='utf-8-sig')

################################################################
## experimental
################################################################

##get all the post synaptic stuff
##trying to sort out which values to keep/discard..
##there are the entries like [CP01],
##   which are unsure connections so we discard
## but there's also stuff like unk5590...
## even unk4733[hyp], contin5934[AVF][AVH][AVJ][PVN]...
#all_post = df["post"].tolist()
#all_post = [x.split(",") for x in all_post]
#all_post = sum(all_post,[])
#all_post = sorted(list(set(all_post)))
#all_pre = sorted(list(set(df["pre"].tolist())))


## way to get new dataframe;
## axis=1 means apply to rows
## 'expand' makes it 3 columns instead of 1
#df.apply(lambda row : [row.pre, sort_post(row.post), row.sections], axis=1, result_type='expand').head()


