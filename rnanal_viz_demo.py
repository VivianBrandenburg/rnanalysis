#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: vivbra
examples for using rnanalysis package
"""



# read data
import os
from rnanalysis import rnaseqdata

example_data = 'example_data/DESeq_tables/'
f = os.listdir(example_data)

datanames = [x.split('.')[0] for x in f]
filenames = [example_data+x for x in f]

data = rnaseqdata(filenames, datanames,  padj=0.05, log2fold=1)



# compare 
from rnanalysis import compare
comp = compare(data.table_log2change)

# plot
from rnanalysis import compplot
compplot(comp.table)



# =============================================================================
# change included experiments
# =============================================================================

keep = datanames[:3]
comp = compare(data.table_log2change,keep=keep)
compplot(comp.table)





# =============================================================================
# custom ordering of groups
# =============================================================================
comp = compare(data.table_log2change,keep=keep)

## prepare a key function for reorderung the table 
neworder = [('wt',), ('mut1',), ('mut3',),  ('mut1', 'wt'), ('mut3', 'wt'),  ('mut3', 'mut1'), ('mut3', 'mut1', 'wt')]
keys  =lambda col: col.map(lambda key: neworder.index(key))

## reorder the data
reordered_data = comp.table.sort_values(by='group_labels', key=keys )
compplot(reordered_data)


# =============================================================================
# divide up-/down-regulated genes 
# =============================================================================

comp = compare(data.table_log2change)
table_updown = comp.make_table_updown()
compplot(table_updown)

#%%
# =============================================================================
# mark specific genes 
# =============================================================================
comp = compare(data.table_padj, keep=keep)
mark_genes = data.table_padj[data.table_padj['mut2'].notnull()].index
mark_genes = mark_genes.to_list()

table_marked_genes = comp.make_table_markedGenes(mark_genes)
compplot(table_marked_genes)




# =============================================================================
# change visualization 
# =============================================================================
compplot(table_marked_genes, 
 		 title='regulated genes in different colors', 
 		 bar_names={'marked':'some specific genes', 
 					'unmarked':'other regulated genes'}, 
 		 bar_colors={'marked':'teal', 
					 'unmarked':'maroon'}, 
 		 dot_color='dimgrey')
