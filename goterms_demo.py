#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:54:29 2024

@author: viv
"""

# =============================================================================
# prepare GO lists 
# =============================================================================
import rnanalysis.gotools as go

# Prepare the GOterms as a Dictionary with the genes as keys
terms_file = 'example_data/GOterms/GOterms_list.txt'
terms = go.get_GO_terms(terms_file)

# write terms for category *Molecular Function* to a file
go.write_gene2GO(terms,'Molecular Function', 'example_data/GOterms/Terms_MF.txt')



# =============================================================================
# prepare GOuniverse
# =============================================================================
import os 

# get a list of files to include in GOuniverse
tabledir = 'example_data/DESeq_tables/'
filelist = [os.path.join(tabledir, x) for x in os.listdir(tabledir)]

# collect the genes from that file
GOuniverse = go.genes_from_DESeqResults(filelist)

# write them to outfile 
with open('example_data/GOterms/GOuniverse.txt','w') as outf: 
    outf.write('\n'.join(GOuniverse))




# =============================================================================
# GOterms visualization
# =============================================================================


import os
import rnanalysis.gotools as go

# get a list of files in a nice sorting.
# This sorting will be used for plotting order
example_data = 'example_data/GOenrichments/'
f = os.listdir(example_data)
f = sorted(f)[:-1] # excluding wt for a nicer demo image 

filenames = [example_data + x for x in f]
datanames = [x.split('.')[0] for x in f]

# set measure that you want to use for choosing the data to plot 
significance_column='weightFisher' 
significance_cutoff = 0.01



# process the data with p-values and genecounts for each enriched term
# this will also produce a column with categories for p-values 
enrichedTerms = go.enrichedGO.from_filenames(filenames, datanames, significance_column, significance_cutoff)






####### plotting with jittering #######



import rnanalysis.gotools.visualization as viz
import seaborn as sns 
import matplotlib.pyplot as plt


### set the sizes of dots and rows as well as
### the horizontal order of dots. These will 
### be used to calculated the y_jitter values


# set hrow height in inches. 
row_height_in_inches = 0.7 # change this value if you want higher / narrower rows 
figsize = (10, len(enrichedTerms.table)*row_height_in_inches)

# setting categorical dot sizes 
size_dict = {'< 0.05':50,  '< 0.005':130, '< 0.001':300,}
size_order = list(size_dict.keys())

# setting colors 
colors = { 'mut1': '#E69F00', 'mut2': '#009E73', 'mut3': 'black'}

# setting the horizontal order of labels and dots (from top to bottom)
hue_order = ['mut1', 'mut2', 'mut3']


sorted_table = viz.y_jitter(figsize=figsize, 
                            df=enrichedTerms.table,
                            xcol='genecount', 
                            ycol='Term',
                            sizecol='p categories',
                            size_dict=size_dict, 
                            hue_order = hue_order)



# set plot style and figure size 
sns.set_style('whitegrid', {'axes.grid': False})
plt.figure(figsize=figsize)

# use seaborn for the actual plotting. 
sns.scatterplot(sorted_table, x='genecount', y='y_jitter', 
                size='p categories', sizes=size_dict, 
                size_order = size_order, 
                hue='variable', palette=colors,                
                zorder=3, linewidth=0, hue_order=hue_order)


# replace the y_jitter values with the original GOTerms on the y-axis
viz.set_y_ticks(enrichedTerms.table['Term'])

# set a nice grid on the y-axis above and below each GOTerm
viz.set_grid(enrichedTerms.table, x='genecount', y='Term')


plt.ylabel('')
plt.legend(loc='upper right')
plt.savefig('plots/GOenrichment/jittered.png', bbox_inches='tight')
plt.show()



''' ####### simple plotting #######

import matplotlib.pyplot as plt
import seaborn as sns
# setting colors
colors = {'mut1': '#E69F00', 'mut2': '#009E73', 'mut3': 'black'}
hue_order = [ 'mut1', 'mut2', 'mut3']

# setting dot sizes
size_dict = {'< 0.05':50,  '< 0.005':130, '< 0.001':300,}
size_order = list(size_dict.keys())

# setting plot style and size
row_height_in_inches = 0.7
figsize = (10, len(enrichedTerms.table)*row_height_in_inches)
sns.set_style('whitegrid')
plt.figure(figsize=figsize)



test = sns.scatterplot(enrichedTerms.table, x='genecount', y='Term', 
                size='p categories', sizes=size_dict, 
                size_order = size_order, 
                hue='variable', palette=colors,
                hue_order=hue_order)
# plt.savefig('plots/GOenrichment/simpleplot.png', bbox_inches='tight')
plt.show()
'''
