


# Installation

Gotools is part of my rnanalysis package. Download the package form ?here? and install it with pip: 

```bash
cd rnanal
pip install .
```



# Preparing GOterms

I downloaded the GO terms from [STRING](https://string-db.org/cgi/download?sessionId=baRWqSm9bX7G), the used terms are found in [GOterms_list.txt](example_data/GOterms/GOterms_list.txt). Extract a GOterms translation file from the list: 
```python
import rnanalysis.gotools as go

# Prepare the GOterms as a Dictionary with the genes as keys
terms_file = 'example_data/GOterms/GOterms_list.txt'
terms = go.get_GO_terms(terms_file)

# write terms for category *Molecular Function* to a file
go.write_gene2GO(terms,'Molecular Function', 'GOterms/Terms_Molecular_Function.txt')
```

This results in a file like [Terms_MF.txt](example_data/GOterms/Terms_MF.txt), which you can use for the enrichment analysis. To prepare a **GOuniverse** you can use this:
```python
import os 

# get a list of files to include in GOuniverse
tabledir = 'example_data/DESeq_tables/'
filelist = [os.path.join(tabledir, x) for x in os.listdir(tabledir)]

# collect the genes from that file
GOuniverse = go.genes_from_DESeqResults(filelist)

# write them to outfile 
with open('example_data/GOterms/GOuniverse.txt','w') as outf: 
    outf.write('\n'.join(GOuniverse))
```


## GOterm enrichment
I used [topGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf) for analysing gene enrichment with [this R code](Goenrichment.r).

# Visualize GO enrichment

Use the ***class enrichedGO***   to prepare the enrichment results for visualization. You can initialize an instance directly with a list of files, or create that filelist with a list of file names



First prepare a list of filenames and a list of data names. 
```python
# get a list of files in a nice sorting. This sorting will be used for plotting order
example_data = 'example_data/GOenrichments/'
f = os.listdir(example_data)
f = sorted(f)
f.insert(0, f.pop())

filenames = [example_data + x for x in f]
datanames = [x.split('.')[0] for x in f]
```

Get the enriched GOTerms: 
```python
import rnanalysis.gotools as go

# set measure that you want to use for choosing the data to plot 
measure='weightFisher' 
significance_cutoff = 0.01 

# process the data with p-values and genecounts for each enriched terms
enrichedTerms = go.enrichedGO.from_filenames(filenames, datanames, measure, significance_cutoff)

```

The attribute ***enrichedTerms.table*** is a pd.DataFrame with measure values and genecounts of all significantly enriched GOTerms. 



### Seaborn scatterplot without overlapping dots

I use matplotlib and seaborn to plot the enriched GOterms. But the plots were hard to read because some of the dots were overlapping. To solve this I wrote a workaround in which I added an extra collumn 'y_jitter' and changed the y-values so that the dots do not overlap in vertical orientation. The values of the x-axis are unchanged. 

Use gotools.enrichedGO.from_filenames to prepare a table  with p-values and genecounts for each enriched term. 
```python

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
```


Then use  gotools.visualization.y_jitter to prepare an extra column for moving the dots on the y-axis so they do not overlap. Set the sizes of dots and rows as well as the horizontal order of dots. These will be used to calculated the y_jitter values.

``` python
import rnanalysis.gotools.visualization as viz

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

# produce a new table with a column for non-overlapping dots
sorted_table = viz.y_jitter(figsize=figsize, 
                            df=enrichedTerms.table,
                            xcol='genecount', 
                            ycol='Term',
                            sizecol='p categories',
                            size_dict=size_dict, 
                            hue_order = hue_order)

```

Now you can use the new column for plotting. 
You can seperate the overlappng dots with the function  ***y_jitter()***. Specify the column names of columns for x-axis, y-axis and size selection. Also specify size_dict dictionary, the figsize and the dataframe (df) that will be used. Optionally, you can give a hueorder that controls the jittering on the y axis to stack the dots in a particular order.  

```python
import seaborn as sns 
import matplotlib.pyplot as plt

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
```
<img src="'plots/GOenrichment/jittered.png" width="600"/>
![[plots/GOenrichment/jittered.png | 600]]


For comparison, here is the image without jittering:

<img src="'plots/GOenrichment/simpleplot.png" width="600"/>
![[simpleplot.png|600]]


