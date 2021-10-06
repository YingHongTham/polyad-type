import pandas as pd
import matplotlib.pyplot as plt
##note to self: need to apt-get install python-tk for matplotlib

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
    post_list = map(clean_name, post_list)
    post_list = filter(lambda x : x != "" , post_list)
    post_list = sorted(post_list)
    post = ",".join(post_list)
    return post

#also removes the square brackets, treat them as sure
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
    if good_cellname_old(name) != (name in neurons):
        print(name, "; in SI4? : ", name in neurons)
        weird_names.append(name)
    return name in neurons

##detecting bad names by manually eliminating
def good_cellname_old(name):
    return name[0:3] != "unk" and name[0:3] != "obj" and name[0:6] != "contin" and name[0:4] != "frag"


######################################################################

#get list of neurons (help with pruning bad data)
neurons = pd.read_csv("male_neurons.csv", header=None, comment='#')
#remove unexpected spaces...
neurons = neurons[0].apply(lambda x : x.replace(" ",""))
neurons = neurons.tolist()

##note to self: DataFrame object doesn't have tolist(),
##but a column object does..
##doing DataFrame.values.tolist() gets a list of lists
##(one list for each row)

#TODO get the contact data from the mysql
#TODO merge left and right neurons

#get synapse list
df = pd.read_csv('SI-3-Synapse-lists-male.csv')
df = df.query('EM_series=="N2Y" & type=="chemical"')
df = df[["pre","post","sections"]]

#clean the pre neurons
df["pre"] = df["pre"].apply(clean_neuron_name)
#clean the post neurons
df["post"] = df["post"].apply(clean_post)

#after cleaning, empty string means unknown
df_pruned = df[(df.pre != "") & (df.post != "")]

##like mysql
##groupby returns DataFrame, but has weird indexing
df_grouped = df_pruned.groupby(["pre","post"]).sum().reset_index()

##get synapses that are large ( >2 sections)
df_robust = df_grouped[df_grouped.sections > 2]

##count number of post-synaptic types for each neuron
df_count = df_grouped.groupby(["pre"]).size().reset_index(name="count")
##graph it
counts = df_count["counts"]
numbins = counts.max() - counts.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(counts,bins=numbins)
ax.set_title("Distribution of synapse outdegrees (N2Y)")
ax.set_xlabel("outdegree")
ax.set_ylabel("frequency")
plt.show()

##count number of post-synaptic types for each neuron
##but for robust
df_count_robust = df_robust.groupby(["pre"]).size().reset_index(name="count")
##graph it
counts_robust = df_count_robust["count"]
numbins = counts_robust.max() - counts_robust.min() + 1
fig, ax = plt.subplots(1,1)
ax.hist(counts_robust,bins=numbins)
ax.set_title("Distribution of synapse outdegrees (N2Y) (robust)")
ax.set_xlabel("outdegree")
ax.set_ylabel("frequency")
plt.show()

##turn "edge list" into adjacency table
#(see https://medium.com/@yangdustin5/quick-guide-to-pandas-pivot-table-crosstab-40798b33e367)
def save_df(df,filename):
    adj_table = pd.crosstab(index=df["pre"],columns=df["post"],values=df["sections"],aggfunc="sum")
    #adj_table.to_csv('adj_table_03.csv')
    adj_table.to_csv(filename+'.csv')
    adj_table_transpose = adj_table.transpose()
    adj_table_transpose.to_csv(filename+'_transpose.csv')



################################################################
## experimental
################################################################

def ff(row):
    print("Contin number: ", row.contin_Number)
    ans = ""
    for i in range(1,row.partner_Number):
        print(row["post" + str(i)])
        ans += row["post" + str(i)]
    return ans

df.apply(ff, axis=1).head()


#get all the post synaptic stuff
#trying to sort out which values to keep/discard..
#there are the entries like [CP01],
#   which are unsure connections so we discard
# but there's also stuff like unk5590...
# even unk4733[hyp], contin5934[AVF][AVH][AVJ][PVN]...
all_post = df["post"].tolist()
all_post = [x.split(",") for x in all_post]
all_post = sum(all_post,[])
all_post = sorted(list(set(all_post)))

all_pre = sorted(list(set(df["pre"].tolist())))


## way to get new dataframe;
## axis=1 means apply to rows
## 'expand' makes it 3 columns instead of 1
df.apply(lambda row : [row.pre, sort_post(row.post), row.sections], axis=1, result_type='expand').head()


