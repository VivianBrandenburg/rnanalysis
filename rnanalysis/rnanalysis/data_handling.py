"""
@author: vivbra
read in and process data from DESeq2 analysis
"""

import pandas as pd
import numpy as np
import os
from collections import OrderedDict



def read_file(filename):
    data = pd.read_csv(filename)
    data = data.rename({'Unnamed: 0':'ID'}, axis='columns')
    return data



def name_from_experiment(keyword):
    # function to get treatment/control genotype name without condition
    if keyword == 'name_treatment':
        return lambda x : x.split('_')[1].split('.')[0]
    if keyword == 'name_control':
        return lambda x : x.split('_')[0].split('.')[0]







        
        


class rnaseqdata():
    
    
    def __init__(self, experiment_files, experiment_names, padj=0.05, log2fold=1):
        
        self.filelist = experiment_files
        self.experiments = experiment_names
        self.table_padj_all = self.__get_data_tables('padj')
        self.table_log2change_all = self.__get_data_tables('log2FoldChange')
        self.table_padj, self.table_log2change = self.__apply_filter(padj, log2fold)
        self.table_updown = self.__make_table_updown(self.table_log2change)
        


    @classmethod
    def __from_path(cls, datapath, *args):
        """Initiate rnaseqdata object with path. 

        Parameters
        ----------
        datapath : str. Path to data
        *args : padjust and log2foldchange. optional 
            Default is padj=0.05, log2fold=1

        Returns
        -------
        rnaseqdata instance
        """
        # read filenames and select files
        filelist = [datapath + x for x in os.listdir(datapath)
                         if x[-8:] == '_all.csv' and x.find('L5C1') == -1]
        experiments = [x.replace('_all.csv', '').replace(datapath,'') for x in filelist]
        return cls(filelist, experiments, *args)




    def __get_data_tables(self, key):  
        names, files = self.experiments, self.filelist
        
        # initialize  
        first_data = read_file(files[0])
        select = first_data.loc[:,['ID',key]]
        select = select.rename(columns={key:names[0]})
        
        for name, file in zip(names[1:], files[1:]):
            data = read_file(file)
            select = pd.merge(select, data.loc[:,['ID', key]],
                                  how='outer', on='ID')
            select = select.rename(columns={key:name})
        
        select = select.set_index('ID', drop=True)
        return select


    def __filtering(self, tab_in, masked):
        tab_out = tab_in.where(masked.notna())
        tab_out = tab_out.dropna(axis=0, how='all')
        return tab_out
    
    
    def __apply_filter(self, padj, log2fold):
        masked = self.table_padj_all.mask(self.table_padj_all > padj)
        masked = masked.mask(self.table_log2change_all.abs() < log2fold)
        
        table_padj = self.__filtering(self.table_padj_all, masked)
        table_log2fold = self.__filtering(self.table_log2change_all, masked)
        
        return table_padj, table_log2fold
    
    def __make_table_updown(self, data):
        data = data.apply(np.sign).replace({-1: 'down', 1: 'up', 
                                            np.nan: np.nan})
        return data
        
        

