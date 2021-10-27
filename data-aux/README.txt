-contact-Brittin-N2U-*.csv are produced by script
	get-contact-matrix-Brittin-N2U.py
	it extracts data from a mysql table
	note there is some work to set up the sql server and stuff
	see README in that folder

-male_celllist.csv
	manually copied all male cell name from SI-4-Celllists.xlsx
	read from file (see also helpers.py):
	>>> celllist = pd.read_csv("male_celllist.csv", header=None, comment='#')

-male_celllist_LR.csv
	produced by script produce-LR-pairs
	table of names that have left/right versions,
	three columns, name without L/R, name with L, name with R
	read from file (see also helpers.py):
	>>> cellLR = pd.read_csv('male_celllist_LR.csv',index_col=False)

