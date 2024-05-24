#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:24:51 2024

@author: vivbra
"""

import numpy as np
import matplotlib.pyplot as plt


def y_jitter(figsize,df, xcol, ycol, sizecol, size_dict, hue_order=None):
    """ 
    Adds a 'y_jitter' column to a pandas DataFrame o avoid overlaps of points in a seaborn.scatterplot plot with a categorical y-axis.. The 'y_jitter' column should be used as y-axis for seaborn.scatterplot, the original category names of the y-axis will be added automatically. 
    If hueorder is provided, the columns will automatically be re-ordered to stack the markers according to the hueorder. Thus, use the returned DataFrame as input for seaborn.scatterplot 
    
    
    Parameters
    ----------
    figsize : array 
        with (figure_height, figure_width) in inch
    df : pandas.DataFrame
        data to work with 
    xcol : column name (str)
        data  for x axis 
    ycol : column name (str)
        categorical (!) data of y axis 
    sizecol : column name (str)
        categorical (!) data for marker size 
    size_dict : dictionary
        sizecol -> markersize 
        
    hue_order : list of strings, optional. 
        order in which the dots are to be stacked. The default is None.

    Returns
    -------
    df : pandas.DataFrame with new column 'y_jitter', optionally reaaranged according to varaible'hue_order'

    """
    
    # write a dict in whcih categories are translated into numerical on y axis
    y_names = list(df[ycol].unique())[::-1]
    y_ids = {v:n+1 for n,v in enumerate(y_names)}
    

    # this shrinkfactor is used to adjust the estimated height/width of 1 y in inches. I thought the size of y in inches should be (inch_size_y_axis / inch_size_of_one_y). Turns out it is not and has to be corrected with this factor. I have no clue why this is needed, just i know this code needs it to work 
    shrinkfactor = 0.75
    
    # convert the sizes in lineheights relative to the overall plot size
    plot_width, plot_height = figsize
    y_in_inches = (plot_height/len(y_names))*shrinkfactor
    x_in_inches = (plot_width/(max(df[xcol])*1.05))*shrinkfactor
    
    # convert the square point size of a marker in its height in inches 
    height_in_inches = {k:np.sqrt(v)/72 for k,v in size_dict.items()}
    
    # convert the height of a marker in unit y 
    relative_dot_height = {k:v/y_in_inches for k,v in height_in_inches.items()}


    # if the markers are to stacked in a specific order, prepare this
    # hue_order is reversed to be consistent with the behavaviour of hue_order in seaborn.scatterplot
    if hue_order:
        variable_ranks = {v:n for n,v in enumerate(hue_order[::-1])}
    
    
    # find the overlapping markers with same y and a distance in x that is smaller than largest point in size_dict. 
    df['y_jitter'] = np.nan
    overlaps, overlaps_index = [], []
    for index, row in df.iterrows():
        _y = row[ycol]
        _x = row[xcol] 
        _x_overlap = max(height_in_inches.values())/x_in_inches
        
        overlap = df.copy()[(df[ycol] == _y) & (abs(df[xcol]-_x) < _x_overlap)]
        
        # if no overlapping markers, just write down the number corresponding to the category of this point
        if len(overlap) < 2:
            df.loc[index, 'y_jitter'] = y_ids[_y]
        
        # if there are overlaps, add the overlapping markers to a new dataframe and add this on the list of overlaps 
        else: 
            if index not in overlaps_index: 
                overlaps_index += list(overlap.index)
                overlaps.append(overlap)
    
    # go over all overlaps
    for overlap in overlaps:
        
        # sort the stacking markers according to the variable hue_order (optional)
        if hue_order:
            overlap['vranks'] = [variable_ranks[x] for x in overlap.variable]
            overlap = overlap.sort_values(by='vranks')    
        
        #  adjust the value of y to shift overlaping markers on the y axis
        size_sum = sum([relative_dot_height[x] for x in overlap[sizecol]]) 
        new_position = y_ids[overlap.iloc[0][ycol]] - size_sum/2
        for index, row in overlap.iterrows():
            own_size =  relative_dot_height[row[sizecol]]
            df.loc[index, 'y_jitter'] = new_position + own_size*0.5 
            new_position = new_position + own_size
    

    return df



def set_y_ticks(df_series):
    
    y_names = list(df_series.unique())[::-1]
    
    #change ytick names back from numerical to categories
    plt.yticks(range(1,len(y_names)+1))
    plt.gca().set_yticklabels(y_names)
    
    #adjust x axis for nicer looking figure
    plt.ylim(0.5, len(y_names)+0.5)
    
    
    

colors = {'WT.no_WT.H2O2.csv' : 'black',
          'LO.no_LO.H2O2.csv': '#56B4E9',
          'lsrB.no_lsrB.H2O2.csv':'#E69F00',
          'oxyR.no_oxyR.H2O2.csv':'#009E73',
          'WT.no_lsrB.no.csv':'#E69F00',
          'WT.no_LO.no.csv':'#56B4E9',
          'WT.no_oxyR.no.csv':'#009E73',
          'LO': '#56B4E9', 
          'lsrB': '#E69F00',
          'oxyR': '#009E73',
          'WT' : 'black',
          }







def set_grid(df, x, y):
    # set grid
    x = df[x].unique()
    y = df[y].unique()
    xgrid = [x-0.5 for x in range(2,len(y)+1)]
    for i in xgrid: 
        plt.axhline(i, color='lightgrey', zorder=0)
    
    max_x_margin = int(max(x) + 0.05*max(x)) 
    ygrid= range(0,max_x_margin,10)
    for i in ygrid:
        plt.axvline(i, color='lightgrey', zorder=1)
    

