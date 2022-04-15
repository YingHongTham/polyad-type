import json

import sys
sys.path.append('../')
import utils.helpers as helpers

exp = helpers.get_gene_expression()

# list of cells in which a given gene is expressed
# e.g. cellsExpGene['zig-8'] = ['cdlL', 'cdlR']
cellsExpGene = { gene : list(filter(lambda x : exp.loc[gene][x] > 0, exp.loc[gene].keys())) for gene in exp.index }

# genes sorted by the list (order on set of lists is lexicographic)
geneSorted = sorted(cellsExpGene, key = lambda x : cellsExpGene[x])

# get unique lists
# uniqueGenes is list of indices of geneSorted
# of some gene representing the set of genes which have the same list of cells
uniqueGenes = []
for g in geneSorted:
	if len(uniqueGenes) == 0:
		uniqueGenes.append(g)
	elif cellsExpGene[g] != cellsExpGene[uniqueGenes[-1]]:
			uniqueGenes.append(g)

f = open("expression-matrix.json", "w")
f.write(json.dumps(cellsExpGene, indent=2))
f.close()

#from subprocess import Popen, PIPE
#p = Popen(['xsel','-pi'], stdin=PIPE)
#p.communicate(input='Hello, World')

