import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class DataWrangling:
    def __init__(self, file):
        self.file = file
        if isinstance(file, pd.DataFrame):
            self.df = file
        elif file.endswith('.csv'):
            self.df = pd.read_csv(file)
        
       
        self.columns = list(self.get_columns())
    
        self.label_font = 30
        self.legend_font =  20
        self.title_font = 35
        self.linewidth = 4

    def get_columns(self):
        return self.df.columns
    
    def get_levels(self, col):
        return self.df[col].unique()
    
    def get_data_by_level(self, col, level):
        if level not in self.get_levels(col):
            print('Error: level not in the column')
            return None
        else:
            return self.df[self.df[col] == level]
        
    def remove_columns(self, col_to_remove):
        self.df = self.df.drop(columns=col_to_remove)
        self.columns = list(self.get_columns())

        
    def get_averaged_data(self, col_to_average):
        columns = self.columns.copy()
        for col in col_to_average:
            if col not in columns:
                print('Error: column not in the dataframe')
                return None
            columns.remove(col)

        averaged_df = self.df.groupby(columns, as_index=False).agg({
            col_to_average[i]: 'mean' for i in range(len(col_to_average))
        })

        # get the standard deviation
        std_df = self.df.groupby(columns, as_index=False).agg({
            col_to_average[i]: 'std' for i in range(len(col_to_average))
        })

        for col in col_to_average:
            averaged_df[col + '_std'] = std_df[col]

        return averaged_df
    
    def switch_to_averaged_data(self, col_to_average):
        columns = self.columns.copy()
        for col in col_to_average:
            if col not in columns:
                print('Error: column not in the dataframe')
                return None
            columns.remove(col)
        
  
        averaged_df = self.df.groupby(columns, as_index=False).agg({
             col_to_average[i]: 'mean' for i in range(len(col_to_average))
        })

        std_df = self.df.groupby(columns, as_index=False).agg({
            col_to_average[i]: 'std' for i in range(len(col_to_average))
        })

        for col in col_to_average:
            averaged_df[col + '_std'] = std_df[col]

        self.df = averaged_df
        self.columns = list(self.get_columns())


    def plot_y_vs_x_by_cat(self, y, x, categories, ymin=None, ymax=None, std = True, alpha=0.2, xticks=False, title=None, save=False, save_path=None, startidx=0):        

        sns.set_theme(style="whitegrid", palette="bright", color_codes=True)
        plt.figure(figsize=(13, 10))
        plt.grid(alpha=0.5)
                
        if categories is None:
            print(self.df)
            plt.plot(self.df[x][startidx:], self.df[y][startidx:],  linewidth=self.linewidth)
            
            

            if std:
                # chek if the std column exists
                if y + '_std' not in self.columns:
                    print('Error: no standard deviation column')
                    return None
                else:
                    # fill the area between mean - std and mean + std
                    plt.fill_between(self.df[x][startidx:], self.df[y][startidx:] - self.df[y + '_std'], self.df[y][startidx:] + self.df[y + '_std'], alpha=alpha)

        elif len(categories) == 1:
            cat = categories[0]
            levels = self.get_levels(cat)

            for level in levels:
                data = self.get_data_by_level(cat, level)
                data = data[startidx:]
                print(data)
                plt.plot(data[x], data[y], label=level, linewidth=self.linewidth)

                if std:
                    # chek if the std column exists
                    if y + '_std' not in self.columns:
                        print('Error: no standard deviation column')
                        return None
                    else:
                        # fill the area between mean - std and mean + std
                        plt.fill_between(data[x], data[y] - data[y + '_std'], data[y] + data[y + '_std'], alpha=alpha)
            # put a title for the legend
            if cat == 'k':
                plt.legend(title='Size (k)', fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
            else:
                plt.legend(title=cat, fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)

        else: # 2 nested categories
            cat = categories[0]
            nest_cat = categories[1]
            levels = self.get_levels(cat)
            nest_levels = self.get_levels(nest_cat)
            for level in levels:
                for nest_level in nest_levels:
                    data = self.get_data_by_level(cat, level)
                    data = data[data[nest_cat] == nest_level]
                    data = data[startidx:]
                    plt.plot(data[x], data[y], label=level + ' - ' + nest_level, linewidth=self.linewidth)
                    if std:
                        # chek if the std column exists
                        if y + '_std' not in self.columns:
                            print('Error: no standard deviation column')
                            return None
                        else:
                            # fill the area between mean - std and mean + std
                            plt.fill_between(data[x], data[y] - data[y + '_std'], data[y] + data[y + '_std'], alpha=alpha)

            # put a title for the legend
            if cat == 'k':
                plt.legend(title='Size (k)', fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
            else:
                plt.legend(title=cat + ' - ' + nest_cat, fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)

        # ymin, ymax
        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        # xlabel, yabel
        if x == 'size':
            plt.xlabel('# of processes', fontsize= self.label_font)
        elif x == 'nthreads' or x == 'n_threads':
            plt.xlabel('# of threads', fontsize= self.label_font)
        elif x == 'k':
            plt.xlabel('Size (k)', fontsize= self.label_font)
        if y == 'time':
            plt.ylabel('Time (s)', fontsize= self.label_font)
        elif y == 'GFLOPS':
            plt.ylabel('GFLOPS', fontsize= self.label_font)
        
        
        # xticks
        plt.xticks(fontsize= self.label_font)
        plt.yticks(fontsize= self.label_font)
        if isinstance(xticks, list):
            plt.xticks(xticks)
        elif xticks:
            plt.xticks(self.df[x].unique())

        # title
        if title is not None:
            plt.title( title, fontsize=self.title_font )

        plt.tight_layout()
        if save:
            plt.savefig(save_path + '.png', dpi=300)
        plt.show()



    


############################################
# speedup
############################################


    def plot_speedup(self, y, x, categories, ymin=None, ymax=None, xticks=False, title=None, save=False, save_path=None):

        
        sns.set_theme(style="whitegrid", palette="bright", color_codes=True)
        plt.figure(figsize=(13, 10))
        plt.grid(alpha=0.5)
        
        if categories is None:
            plt.plot(self.df[x], (self.df[y][0] / self.df[y]), label='speedup', linewidth=self.linewidth)
            plt.plot(self.df[x], self.df[x], '--', label='ideal')
            plt.legend(fontsize=self.legend_font)
            plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)


        elif len(categories) == 1:
            cat = categories[0]
            levels = self.get_levels(cat)

            for level in levels:
                data = self.get_data_by_level(cat, level)
                data.reset_index(inplace=True)
                base = data[y][0]
                plt.plot(data[x], base/(data[y]), label=level, linewidth=self.linewidth)

            # plot the ideal speedup
            plt.plot(data[x], data[x], '--', label='ideal')

            # put a title for the legend
            if cat == 'k':
                plt.legend(title='Size (k)', fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
            else:
                plt.legend(title=cat, fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
        
        else:
            cat = categories[0]
            nest_cat = categories[1]
            levels = self.get_levels(cat)
            nest_levels = self.get_levels(nest_cat)
            for level in levels:
                for nest_level in nest_levels:
                    data = self.get_data_by_level(cat, level)
                    data = data[data[nest_cat] == nest_level]
                    data.reset_index(inplace=True)
                    base = data[y][0]
                    plt.plot(data[x], base/data[y], label=level + ' - ' + nest_level, linewidth=self.linewidth)

            plt.plot(data[x], data[x], '--', label='ideal')
                    
            plt.legend(title=cat + ' - ' + nest_cat, fontsize=self.legend_font)
            plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)

        
        # ymin, ymax
        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        # xlabel, yabel
        if x == 'size':
            plt.xlabel('# of processes', fontsize= self.label_font)
        elif x == 'nthreads' or x == 'n_threads':
            plt.xlabel('# of threads', fontsize= self.label_font)
        elif x == 'k':
            plt.xlabel('Size (k)', fontsize= self.label_font)
        plt.ylabel('Speedup', fontsize= self.label_font)
        
        
        # xticks
        plt.xticks(fontsize= self.label_font)
        plt.yticks(fontsize= self.label_font)
        if isinstance(xticks, list):
            plt.xticks(xticks)
        elif xticks:
            plt.xticks(self.df[x].unique())

        # title
        if title is not None:
            plt.title('Speedup - ' + title, fontsize=self.title_font )

        plt.tight_layout()
        if save:
            plt.savefig(save_path + '.png', dpi=300)
        plt.show()



    ############################################
    # efficiency
    ############################################

    def plot_efficiency(self, y, x, categories, ymin=None, ymax=None, xticks=False, title=None, save=False, save_path=None, startidxeff=1):
        
        sns.set_theme(style="whitegrid", palette="bright", color_codes=True)
        plt.figure(figsize=(13, 10))
        plt.grid(alpha=0.5)


        if categories is None:
            plt.plot(self.df[x][startidxeff:], (self.df[y][0] / self.df[y][startidxeff:])/self.df[x][startidxeff:], label='efficiency', linewidth=self.linewidth)
            plt.plot(self.df[x][startidxeff:], np.ones(len(self.df[x])-startidxeff), '--', label='ideal')
            plt.legend(fontsize=self.legend_font)
            plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)

          


        elif len(categories) == 1:
            cat = categories[0]
            levels = self.get_levels(cat)

            for level in levels:
                data = self.get_data_by_level(cat, level)
                data.reset_index(inplace=True)
                base = data[y][0]
                plt.plot(data[x][startidxeff:], (base/(data[y][startidxeff:]))/data[x][startidxeff:], label=level, linewidth=self.linewidth)


            # plot the ideal efficiency
            plt.plot(data[x][startidxeff:], np.ones(len(data[x])-startidxeff), '--', label='ideal')

            # put a title for the legend
            if cat == 'k':
                plt.legend(title='Size (k)', fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
            else:
                plt.legend(title=cat, fontsize=self.legend_font)
                plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
        
        else:
            cat = categories[0]
            nest_cat = categories[1]
            levels = self.get_levels(cat)
            nest_levels = self.get_levels(nest_cat)
            for level in levels:
                for nest_level in nest_levels:
                    data = self.get_data_by_level(cat, level)
                    data = data[data[nest_cat] == nest_level]
                    data.reset_index(inplace=True)
                    base = data[y][0]
                    plt.plot(data[x][startidxeff:], (base/(data[y][startidxeff:]))/data[x][startidxeff:], label=level+ ' - ' + nest_level, linewidth=self.linewidth)

            # plot the ideal efficiency
            plt.plot(data[x][startidxeff:], np.ones(len(data[x])-startidxeff), '--', label='ideal')
            
            # put a title for the legend
            plt.legend(title=cat + ' - ' + nest_cat, fontsize=self.legend_font)
            # font of title of the legend
            plt.setp(plt.gca().get_legend().get_title(), fontsize=self.legend_font)
        
        
        # ymin, ymax
        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        # xlabel, yabel
        if x == 'size':
            plt.xlabel('# of processes', fontsize= self.label_font)
        elif x == 'nthreads' or x == 'n_threads':
            plt.xlabel('# of threads', fontsize= self.label_font)
        elif x == 'k':
            plt.xlabel('Size (k)', fontsize= self.label_font)
        plt.ylabel('Efficiency', fontsize= self.label_font)
        
        
        # xticks
        plt.xticks(fontsize= self.label_font)
        plt.yticks(fontsize= self.label_font)
        if isinstance(xticks, list):
            plt.xticks(xticks)
        elif xticks:
            plt.xticks(self.df[x][1:].unique())

        # title
        if title is not None:
            plt.title('Efficiency - ' + title, fontsize=self.title_font )

        plt.tight_layout()
        if save:
            plt.savefig(save_path + '.png', dpi=300)
        plt.show()

        


def pipeline(file, 
             mode, 
             xaxis, 
             yaxis, 
             categories, 
             speedup=True, 
             efficiency=False, 
             ymin=None, ymax=None, 
             yminsp=None, ymaxsp=None, 
             ymineff=None, ymaxeff=None ,std = True, alpha=0.2, 
             xticks=False, 
             xtickssp=None, 
             xtickseff=None, 
             title=None, 
             save=False, save_path=None,
             startidx=0,
             startidxeff=1):
    df = DataWrangling(file)
    df = DataWrangling(df.get_data_by_level('mode', mode))
    df.switch_to_averaged_data([yaxis])
    df.plot_y_vs_x_by_cat(yaxis, xaxis, categories, ymin, ymax, std, alpha, xticks, title, save, save_path, startidx=startidx)
    if speedup:
        if save:
            save_pathsp = save_path + '_speedup'
        df.plot_speedup(yaxis, xaxis, categories, yminsp, ymaxsp, xtickssp, title, save, save_pathsp)
    if efficiency:
        if save:
            save_pathef = save_path + '_efficiency'   
        df.plot_efficiency(yaxis, xaxis, categories, ymineff, ymaxeff, xtickseff, title, save, save_pathef, startidxeff=startidxeff)

   

def pipeline_2(file, 
             xaxis, 
             yaxis, 
             categories, 
             mode = None, 
             speedup=True, 
             efficiency=False, 
             ymin=None, ymax=None, 
             yminsp=None, ymaxsp=None, 
             ymineff=None, ymaxeff=None ,std = True, alpha=0.2, 
             xticks=False, 
             xtickssp=None, 
             xtickseff=None, 
             title=None, 
             save=False, save_path=None):
    df = DataWrangling(file)
    if mode != None:
        df = DataWrangling(df.get_data_by_level('mode', mode))
    df.remove_columns(['iteration'])
    if yaxis == 'time':
        df.remove_columns(['GFLOPS'])
        df.switch_to_averaged_data([yaxis])
    elif yaxis == 'GFLOPS':
        df.remove_columns(['time'])
        df.switch_to_averaged_data([yaxis])
    print(df.columns)
    df.plot_y_vs_x_by_cat(yaxis, xaxis, categories, ymin, ymax, std, alpha, xticks, title, save, save_path)
    if speedup:
        if save:
            save_pathsp = save_path + '_speedup'
        df.plot_speedup(yaxis, xaxis, categories, yminsp, ymaxsp, xtickssp, title, save, save_pathsp)
    if efficiency:
        if save:
            save_pathef = save_path + '_efficiency'   
        df.plot_efficiency(yaxis, xaxis, categories, ymineff, ymaxeff, xtickseff, title, save, save_pathef)

   