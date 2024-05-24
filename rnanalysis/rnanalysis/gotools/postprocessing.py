#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:23:13 2024

@author: vivbra
"""



import pandas as pd
import numpy as np
from rnanalysis.data_analysis import get_all_combinations



class enrichedGO():
    
    def __init__(self, data, measure, significance_cutoff, combinations=None):
        """
        Process enrichment data for plotting. 
        Alternative initiation with enrichedGO.from_filenames([list, of, files], measure, significance_cutoff)

        Parameters
        ----------
        data : dictionary 
            Dictionary of data with {filename:data}.
        measure : str
            columname of the measures you want to use for p-value categories in the plot.
        significance_cutoff : float
            

        Attributes
        -------
        table : pd.DataFrame (melted) 
            table with measure values and genecounts of all siognificantly enriched GOTerms
            entries are ordered according to the groups, eg  ( (e1) (e2) (e3) (e1,e2) (e1,e3) ...)
            e1, e2, e3 are sorted according to order of data (using __init__()) or filenames/data_names (using from_filenames())
        files : dictionary
            data that will be used with {filename:data}.  
        dataOI : table
            measures below the given significance_cutoff that will be used 
        countsOI : table
            genecounts belonging to dataOI that will be used

        """

        
        self.data = data
        
        # get measure values and counts for significant GOterms  
        self.dataOI, self.countsOI = self.__get_tableOI(data, measure, significance_cutoff)
        
        
        # sort the genecount table according to the groups that are involved 
        # groups: ( (e1) (e2) (e3) (e1,e2) (e1,e3) (e2,e3) (e1,e2,e3) )
        
        if type(combinations)==list: 
            self.all_combinations  = combinations
        else: 
            self.all_combinations  = list(get_all_combinations(data.keys()))
        if combinations == 'print':
            print(self.all_combinations)
        
        sorted_table = self.__get_combinations(self.all_combinations, list(data.keys()), self.countsOI)
        
        # melt the table for plotting 
        sorted_table = self.__unpivot( sorted_table)
        
        # add the measure values to the table
        sorted_table = sorted_table.rename({'value':'genecount'}, axis=1)
        sorted_table= pd.merge(sorted_table, self.__unpivot( self.dataOI), how='left', on=['Term', 'variable'])
        
        # get the categorical p-values (*,**,***) instead of continous values
        sorted_table['p categories']= [self.__get_p_categories( x) for x in sorted_table.value]
        
        # remove nan-entries
        self.table = sorted_table[sorted_table.value.notnull()]
        
    
    @classmethod
    def from_filenames(cls, filenames, data_names, measure, significance_cutoff, combinations=None):
        data = {d : pd.read_table(f, sep='\t') for d,f in zip(data_names,filenames)}
        return cls(data, measure, significance_cutoff, combinations)
        
    
    def __get_combinations(self, all_combinations, keys, df):
        
        sorted_table = pd.DataFrame()
        for comb in all_combinations:
            comb = list(comb)
            anticomb = [x for x in keys if x not in comb]
            select = df[df[anticomb].isnull().all(1)]
            select = select[select[comb].notnull().all(1)] 
            select = select.sort_values(by=comb[0])
            sorted_table = pd.concat([sorted_table, select], axis=0)
            sorted_table = sorted_table.reset_index(drop=True)
        return sorted_table
    
    
    
    def __get_p_categories(self, x):
        if x < 0.001:   return '< 0.001'
        elif x < 0.005: return '< 0.005'
        elif x < 0.05:  return '< 0.05'
        else:           return np.nan
        
            
    def __get_tableOI(self, filesOI, measure, threshold, names=None):
        # get a df that contains the measure for all GOs of interest for all experiments
        # and another dataframe that contains the corresponding gene counts
        if not names: names=list(filesOI.keys())
        l = filesOI[names[0]]
        valuesOI = pd.DataFrame({'GO.ID':l['GO.ID'],'Term':l['Term'], 
                                names[0]:l[measure]})
        countsOI = pd.DataFrame({'GO.ID':l['GO.ID'],'Term':l['Term'], 
                                names[0]:l['Significant']})
        
        datanames = list(self.data.keys())
        for k in datanames[1:]:
            r = filesOI[k]
            r_t = pd.DataFrame({'GO.ID':r['GO.ID'], k:r[measure]})
            valuesOI = pd.merge(valuesOI, r_t, how='outer', on='GO.ID')
            r_c = pd.DataFrame({'GO.ID':r['GO.ID'], k:r['Significant']})
            countsOI = pd.merge(countsOI, r_c, how='outer', on='GO.ID')
                
        # override insignificant values with nan 
        for k in datanames:
            valuesOI[k] = [x if x<threshold else np.nan for x in valuesOI[k]]
            countsOI[k] = [y if x<threshold else np.nan for x,y 
                           in zip(valuesOI[k], countsOI[k])]
                                                                       
        # remove entries where only nan 
        countsOI = countsOI[~countsOI[datanames].isnull().all(1)]
        valuesOI = valuesOI[~valuesOI[datanames].isnull().all(1)]
        
        return valuesOI, countsOI
        
    

    def __unpivot(self, df):
        plotdata = df[[x for x in df.columns if x !='GO.ID']]
        plotdata = pd.melt(plotdata, id_vars = 'Term', ignore_index=False)
        plotdata = plotdata.sort_index()
        return plotdata



