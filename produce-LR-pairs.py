import pandas as pd
import csv


# get all N2Y cell names
df = pd.read_csv("male_celllist.csv", header=None, comment='#')
cells = list(df[0])
#with open('male_celllist.csv') as f:
#	reader = csv.reader(f)
#	data = list(reader)



## find the last non-numeric character in string
## returns position, letter
def get_last_letter(cellname):
	for i in range(len(cellname)-1,-1,-1):
		if not cellname[i].isdigit():
			return i, cellname[i]


## >>> split_cellname('dBWML1')
## ['dbWM1', 'dbWML1', 'dbWMR1']
def split_cellname(cellname):
	ind, c = get_last_letter(cellname)
	if c == 'L':
		front = cellname[0:ind]
		back = cellname[ind+1:len(cellname)]
		return [front+back, cellname, front+'R'+back]
	else:
		return ['','','']

cellLR = df.apply(lambda x : split_cellname(x[0]), axis=1, result_type='expand')
cellLR.columns = ['name','L','R']
cellLR = cellLR[cellLR.name != '']


cellLR.to_csv('male_celllist_LR.csv',index=False,encoding='utf-8-sig')
##read with:
#cellLR_read_from_file = pd.read_csv('male_celllist_LR.csv',index_col=False)


######################################################################
## some tests
## most cells with L,R pairs have L/R at the end
## but some, like dBWML1/dBWMR1 have number at the end
## also have to be careful about 'ALA' (there's not 'ARA'!)
## names containing L,R
#cell_L = filter(lambda x : 'L' in x, cells)
#cell_R = filter(lambda x : 'R' in x, cells)
#
#sum(((x[-1] != 'L') and (x[-1] != 'R')) for x in cell_L)
#len(cell_L) ##not the same!
#filter(lambda x : ((x[-1] != 'L') and (x[-1] != 'R')), cell_L)


####just checking if each L neuron has corresponding R neuron
##def split_cellname_R(cellname):
##	ind, c = get_last_letter(cellname)
##	if c == 'R':
##		front = cellname[0:ind]
##		back = cellname[ind+1:len(cellname)]
##		return [front+back, front+'L'+back, cellname]
##	else:
##		return ['','','']
##
##cellLR_test = df.apply(lambda x : split_cellname(x[0]), axis=1, result_type='expand')
##cellLR_test.columns = ['name','L','R']
##cellLR_test = cellLR_test[cellLR_test.name != '']
##
##names = set(cellLR['name'])
##names_test = set(cellLR_test['name'])
##print("Same? ", names == names_test)
