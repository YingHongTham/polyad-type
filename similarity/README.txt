similarity.py does three different things,
giving three collections of tables
the first two parts just load the necessary data
similarity.py produces the following tables


##############################

detects how L/R-symmetric a pair of (pre,post) cells is
	if pre=X,post=Y, then we are comparing the similarity/difference
	between number_sections(X -> Y) and number_sections(X' -> Y')
	where X',Y' are the homologs of X,Y

synapse_simscore_pairwise
contact_simscore_pairwise
contact_L1simscore_pairwise
-columns are pre, post, score, diff_fraction, sections, sections_flipped
-pre,post = cells (only fed it monads)
	I also got rid of the cells without homolog
-sections = number of sections (synapse or contact)
-sections_flipped = sections for the mirrors(homologs),
	e.g. AVAL -> PVCL has 1 section,
		the homologs are AVAR -> PVCR, which has 2 sections,
		so sections_fliped for AVAL,PVCL is 2
-score: if pre=X,post=Y, then score=fn(sections(X -> Y), sections(X' -> Y'))
	X',Y' = homolog of X,Y
	fn = simscore or L1_simscore
	simscore is using Adam's similarity score, higher means more similar
	L1_simscore is actually L1 distance, in this case just the absolute difference
-diff_fraction: fraction of score out of max of sections and sections_flipped

##############################

detects how L/R-symmetric a single cell is, namely,
	if we flip the worm L/R, how similar/different are
	the connections of X in the original and flipped worms?
	e.g. if X -> A, A', B, C, ..
	and X' -> A', A, B', C', ..
	then we compare these two rows

L1_dist_01
Adam_dist_01
L1_dist_01_monad
Adam_dist_01_monad
L1_dist_01_contact
Adam_dist_01_contact
-two columns: pre,score (L1_dist or simscore)

##############################

simply gets the similarity/difference between every pair of cells,
no regard for L/R symmetry

output files:
pairwise(_T)_monad_L1
pairwise(_T)_monad_simscore
-adj matrix, where the weight indicates the similarity/L1_dist
	of the two cells, i.e.
	table[cell1][cell2] = fn(cell1,cell2)
-no flipping LR, just on the nose comparisons of the cell neighbours
-the _T version means we're considering the incoming synapses
