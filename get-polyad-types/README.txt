get-polyad-types.py produces:
-2 lists:
	-synapse-sections-<date>.csv, synapse-sections-robust-<date>.csv
	-list of synapse types with number of sections
-4 tables:
	adj-table-(robust)-<date>-(transpose).csv
	-adj-table has row=pre, col=post types
	-transpose because LibreOffice is having trouble with too many columns
	-robust means we get rid of all synapses of size <= 2
-2 more tables:
	outdegrees-(robust).csv
	-counts number of types of post-synaptic types for each pre
-2 plots of the outdegrees tables
-<date> has date and current time, helpful for versioning

Overview of pipeline:

-start with two excel tables:
	-SI 3 Synapses
	-SI 4 Cell lists
-SI 4 has a list of all cells
	-manually copy the male cells (all sheets except hermaphrodite specific)
		from SI 4 into a file male_celllist.csv
	-load into list celllist
	-note: need to do some cleaning (some have extra white spaces)
-SI 3 rows correspond to synapses
	-load into Pandas DataFrame syn

-extract 3 columns "pre", "post", "sections" (=size) from syn
-"post" is a comma-separated list of post-synaptic cells
-perform cleaning on cell names (some are discarded, see below for more)
	->syn_cleaned
-the robust verions mean only consider synapses of #sections > 2
-total up size of each synapse type
	-output to adj-table csv's
-count number of postsynaptic types for each neuron
	-output to outdegrees csv's
-plot histogram

On the cleaning procedure:
-HOA,[AVG] means not sure if AVG is a post-synaptic cell
	-we will take it as a confirmed yes, so remove all square brackets
-remove whitespaces from names
-entries like cdlR/dgl5R/sph mean unsure about which cell; discard
-check if cell names are in the celllist
	-returns empty string if bad names
	-e.g. for post entry,
		"xxx,HOA,yyy" -> "HOA"
		"xxx,yyy" -> ""
-discard rows with an "" entry


NOTE: some of these procedures were later refactored into helpers.py
