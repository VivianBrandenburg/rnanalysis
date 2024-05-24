
def get_GO_terms(terms_infile):
    """
    Collect GOTerms for each gene for all GO categories from a STRING file

    Parameters
    ----------
    terms_infile : str, filename
        infile as downloaded from STRING

    Returns
    -------
    res : dict, geneIDs as keys
        values are dictioaries with GO categories as keys

    """
    res = {}    
    with open(terms_infile, 'r') as inf:
        for line in inf:
            if not line[0]=='#':
                tabs = line.split('\t')
                gene =tabs[0].split('.')[-1].strip()
                if gene not in res.keys():
                    res[gene]=[]
                function = tabs[1].split('(')[0].strip()
                term = tabs[2].strip()
                res[gene].append([function, term])
    
    # group categories for each gene
    for key,value in res.items():
        categories = list(set([x[0] for x in value]))
        categories = {cat:[] for cat in categories}
        for x in value:
            categories[x[0]].append(x[1])
        res[key]=categories
               
    return res




def write_gene2GO(terms, filterkey, outfile):
    with open(outfile, 'w') as outf:
        for key, value in terms.items():
            if filterkey in value.keys():    
                filterterms = value[filterkey]
                outf.write(key + '\t' + ','.join(filterterms) + '\n')
        


def genes_from_DESeqResults(filelist):
    genes = []
    for f in filelist: 
        with open(f,'r') as inf:
            inf=inf.readlines()
            inf=[x.split(',')[0].strip('"') for x in inf]
            genes+=inf
    genes=list(set([x for x in genes if x ]))
    return genes

