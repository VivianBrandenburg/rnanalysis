import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


class compare():
    """
    takes a pandas Dataframe with results from Differential Expression Analysis and compares different experiments with each other. 

    Parameters
    ----------
    data : pandas DataFrame object 
        contains experiments in columns and geneIDs as index
    keep : str or list of str, optional
        Experiments to include in the comparison. The default is None.
        The order of this list will be used for the order of the experiments in the plot.
    remove : str or list of str, optional
        Experiments to remove from the comparison. The default is None.
    custom_comparison : nested list of strings, optional
        give a list of experiment groups to sort yout diagram 
        
        
    Methods
    ----------
    make_table_markedGenes([list, of, marked, genes]) : 
        make a table that distinguishes between marked and unmarked genes
    make_table_updown() :
        make a table that distinuguishes between up- and down-regulated genes


    Attributes
    ----------
    data : object
        similar to input data, but filtered for included experiments         

    table : table with number of regulated genes for each group. 
    	with columns:
    	included : which experiments are included
        excluded : which experiments are excluded
        data : pd.DataFrame with included genes 
        regulated_genes : IDs of the regulated genes

    """
    
    def __init__(self, data, keep=None, remove=None, marked_genes='', combinations=None):

        self.data = self.select_columns(data=data, remove=remove, keep=keep)
        if combinations:
            self.experiment_groups = combinations
        else:
            columnsOI = [x for x in list(self.data.columns)]
            self.experiment_groups = list(get_all_combinations(columnsOI))
        self.__compdict = self.__compare_all_combinations(remove, keep)
        self.table = self.make_table_perGroup(self.__compdict)
        
    
    def __compare_all_combinations(self, remove=None, keep=None):
        columnsOI = [x for x in list(self.data.columns)]
        results = []
        for columns in self.experiment_groups:
            res_dict = {}
            res_dict['included'] = included = columns
            res_dict['excluded'] = excluded = [
                x for x in columnsOI if x not in included]
            res_dict['data'] = dataOI = self.select_rows(data=self.data,
                                                         include=included, exclude=excluded)
            res_dict['regulated_genes'] = list(dataOI.index)
            results.append(res_dict)
        return results

    @staticmethod
    def __str_to_list(*str_or_list):
        import copy
        """checks if input is str or list. If str, returns list with one item. Else returns list unchanged. This allows for input flexibility in other functions. 
        """
        for i in str_or_list:
            if type(i) == str:
                i = [i]
            yield copy.deepcopy(i)

    def make_table_perGroup(self, comp):
        # add a table with the number of genes regulated in each group
        table = pd.DataFrame({'group_labels': [x['included'] for x in comp],
                                  'excluded': [x['excluded'] for x in comp],
                                  'regulated_genes': [len(x['data']) for x in comp],
                                  'gene_names': [list(x['data'].index) for x in comp]})
        return table

    @classmethod
    def select_rows(cls, data, include, exclude=None):
        """select genes that have non-nan entry in include and (optionally) have nan entry in exclude        
        Returns pandas.DataFrame of selected genes 
        """
        selection = data.copy(deep=True)
        include, exclude = cls.__str_to_list(include, exclude)
        selection = data[data[list(include)].notnull().all(1)]
        if exclude:
            selection = selection[selection[list(exclude)].isnull().all(1)]
        return selection

    @classmethod
    def select_columns(cls, data, remove=None, keep=None):
        """Remove all collumns of the data that are not needed in this comparison
        """
        selection = data.copy(deep=True)  # change variable types to list
        keep, remove = cls.__str_to_list(keep, remove)
        if keep:        # tidy up data table
            selection = selection[keep]
        if remove:
            selection = selection.drop(columns=remove)
        return selection

    def make_table_markedGenes(self, genelist):
        """
        produce a table that counts the genes that are / are not in the given genelist for each of the experiment groups

        Parameters
        ----------
        genelist : list of strings
            list of gene-IDs

        Returns
        -------
        table : pandas.DataFrame
            can be used for plotting with rnanal.compplot
            

        """
        # add a table with regulated genes in the list / not in the list
        comp = self.__compdict
        marked_included = [[g for g in x['data'].index if g in genelist] for x in comp]
        nonmarked_included = [[g for g in x['data'].index if g not in genelist] for x in comp]
        table = pd.DataFrame({'group_labels': [x['included'] for x in comp],
                                  'excluded': [x['excluded'] for x in comp],
                                  'gene_names': [list(x['data'].index) for x in comp],
                                  'marked': [len(x) for x in marked_included],
                                  'unmarked': [len(x) for x in nonmarked_included]}                               
                                 )
        return table
    
    def make_table_updown(self):
        """
        produce a table that counts the up-/downreguilated genes in each experiment group

        Returns
        -------
        table : pandas.DataFrame
            can be used for plotting with rnanal.compplot

        """
        table = pd.DataFrame(columns=['up', 'down', 'ambiguous', 'group_labels'])
        for n, i in enumerate(self.__compdict):
            dataOI = i['data']
            included = i['included']
            valuesOI = dataOI[list(included)].values
            if not 'up' in valuesOI and not 'down' in valuesOI:
                up = len(dataOI[dataOI[list(included)].gt(0).all(1)])
                down = len(dataOI[dataOI[list(included)].lt(0).all(1)])
            else:
                up = len(dataOI[dataOI[list(included)].eq('up').all(1)])
                down = len(dataOI[dataOI[list(included)].eq('down').all(1)])
            ambiguous = len(i['regulated_genes']) - up - down
            table.loc[len(table)] = {'up': up, 'down': down, 'ambiguous': ambiguous,
                                   'group_labels': included}
        return table
    



def get_all_combinations(list_of_conditions):
    """Takes a list and returns all combinations of items in that list.
    Example: [a,b] -> ((a,), (b,), (a,b))
    """
    import itertools
    num = len(list_of_conditions)
    for i in range(num):
        for j in itertools.combinations(list_of_conditions, i+1):
            yield j

# =============================================================================
# plotting stuff
# =============================================================================


class compplot():

    def __init__(self, table, title='', naming_function=None, 
                 bar_names = None, bar_colors = None, dot_color = 'black',
                 outfile='plot.png'):
        """
        Produces a stacked bar plot with a dotplot

        Parameters
        ----------
        table : pd.DataFrame.
        title : str, optional
            The deafult is None
        naming_function : function, optional
            Function to convert column names into plot labels. The default is None.
        bar_names : dictionary, optional
            To translate bar labels, eg {'marked':'A', 'unmarked':'B'}. The default is None.
        bar_colors : dictionary, optional
            To change the bar colors, similar to bar_names. The default is None.
        dot_color : color string, optional
            To change the color of the lower dot plot. The default is 'black'.
        outfile : filename, str. Also specifies the filetype. If None, no plot is saved
            The default is 'plot.png'

        Returns
        -------
        None. Prints a plot and saves it as outfile
        
        """
        
        self.table = table
        self.title = title
        self.naming_function = naming_function
        self.outfile = outfile
        self.bar_names = bar_names
        self.bar_colors = bar_colors  
        self.dot_color = dot_color 
        
        self.__set_styles()        
        self.__plot()

    def __set_styles(self):
        self.style_top = {"axes.spines.right": False, "axes.spines.left": False,
                     "axes.spines.top": False, "axes.spines.bottom": False,
                     "xtick.bottom": False, 'xtick.labelbottom': False,
                     "ytick.left": True, }
        self.style_bottom = {"axes.spines.right": False, "axes.spines.left": False,
                         "axes.spines.top": False, "axes.spines.bottom": False,
                         "xtick.bottom": False, 'xtick.labelbottom': False,
                         "ytick.left": False}       
    
    def __plot(self):
        
        # set styles top        
        sns.set_style(style="ticks", rc=self.style_top)
        
        fig = plt.figure(figsize=(10, 15))
        ax1 = plt.subplot2grid(shape=(4, 1), loc=(0, 0), rowspan=3)
        self.plot_stacked_bars(ax1)
    
        # set styles top        
        sns.set_style(style="ticks", rc=self.style_bottom)
        
        ax2 = plt.subplot2grid(shape=(4, 1), loc=(3, 0))
        self.plot_comparison_info_design(ax2)
        if self.outfile: 
            plt.savefig(self.outfile)

    def plot_stacked_bars(self, ax):
        
        self.table.plot(kind='bar', stacked=True, title=self.title, color=self.bar_colors, ax=ax)
        
        # adjust legend if needed 
        if self.bar_names:
            h, names =  ax.get_legend_handles_labels()
            names = [self.bar_names[x] for x in names]
            ax.legend(names)
        
        # remove unneeded axis labels
        plt.ylabel('number of regulated genes')
        plt.gca().set_xticklabels([])

    def plot_comparison_info_design(self, ax):
        """Generate plot showing which groups of experiment are included in each bar. 
        Grey dots mark every possible position. Blue dots mark the positions of experiments that are included.        
        """
        # prepare data
        black_dot_labels = list(self.table.group_labels)
        if self.naming_function:
            black_dot_labels = [[self.naming_function(y) for y in x]
                               for x in black_dot_labels]
        # 
        y_labels = sum([list(x) for x in black_dot_labels], [])
        y_lables = list(dict.fromkeys(y_labels))        
        black_dots, grey_dots = self.__get_dot_positions(
            y_labels, black_dot_labels)

        # plotting dots
        sns.scatterplot(data=grey_dots, x='x', y='y',
                        s=150, c='#E8E8E8', ax=ax)
        sns.scatterplot(data=black_dots, x='x', y='y', s=150, ax=ax, color=self.dot_color)

        # plotting vertical lines
        for n, i in enumerate(black_dot_labels):
            if len(i) > 1:
                ax.vlines(x=n, ymin=i[0], ymax=i[-1], color=self.dot_color)

        # remove unneeded axis labels
        plt.ylabel(''), plt.xlabel(''), plt.gca().set_xticklabels([])

    def __get_dot_positions(self, y_labels, black_dot_labels):
        """Generate mock data to plot positions of grey dots (all possible positions) and blue dots (marking experiments included in the group)
        """
        black_dots = [[n, j] for n, i in enumerate(black_dot_labels) for j in i]
        black_dots = pd.DataFrame(black_dots, columns=['x', 'y'])

        len_x = black_dots.x.max()
        grey_dots_x = [i for i in range(len_x+1) for j in y_labels]
        grey_dots_y = [j for i in range(len_x+1) for j in y_labels]
        grey_dots = pd.DataFrame({'x': grey_dots_x, 'y': grey_dots_y})
        return black_dots, grey_dots

