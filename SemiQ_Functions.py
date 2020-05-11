import numpy as np
import pandas as pd

# %%
# this puts all the data into a good format
# would recommend to sort the dataframe used after this step
# to keep everything organized
# input is the file name with the dataset
# output is the desired file name
# path is the file path of the folder that contains the input
# will make the output in the same folder as the input
def revise_data(input='', output='', path=''):
    if input == '' or output == '' or path == '':
        print('\nplease check inputs')
        return 1
    

    # Loads in the dataset and fills in the missing values
    input_path = path + '/' + input
    data = pd.read_excel(input_path, sheet_name=None)

    # list of all the excel sheets
    sheets = list(data.keys())
    sheets.remove('Component')

    # the dataframe to store all the revised data
    revised_data = pd.DataFrame()

    # gets only the portion of the data that we are interested in and renames the columns to Concentration and compound name
    for sheet in sheets:
        names = data[sheet].iloc[3]
        # data[sheet].drop([0, 1, 2, 3], inplace=True)
        data[sheet].columns = names
        data[sheet].columns.name = None
        data[sheet].rename(columns={'Filename' : 'Concentration'}, inplace=True)
        data[sheet].dropna(subset=['Concentration'], inplace=True)
        data[sheet].reset_index(inplace=True)
        data[sheet].rename(columns = {'Area' : sheet}, inplace=True)

        # need to filter out the double blanks
        filter = data[sheet]['Concentration'].str.contains('\d*ng')
        
        # attaches this new subset of the data to the revised DataFrame
        if sheet == sheets[0]:
            revised_data = pd.concat([revised_data, data[sheet][filter][['Concentration', sheet]]])
            # print(revised_data)
        else:
            revised_data = pd.merge(revised_data, data[sheet][filter][['Concentration', sheet]],how='inner', on=['Concentration'])
        

    
    # outputs the data to specific output
    output_path = path + '/' + output
    writer = pd.ExcelWriter(output_path)
    revised_data.to_excel(writer,'Sheet1',index=False)
    writer.save()
    
    return 0


# %%
# log transforms the dataframe from the revise_data function
# would advise to do this to keep everything consistant.
# if you do skip this function make sure to melt the dataframe after preprocessing with 'concentration', 'compound', and 'area' as the the 3 columns
def log_transform(data=None):
    # melts the initial dataframe, so that we can log transform it
    melt_df = pd.melt(data, id_vars='concentration', var_name='compound', value_name='area')
    
    # log transform the concentration and area
    melt_df['concentration'] = np.log10(melt_df['concentration'])
    melt_df['area'] = np.log10(melt_df['area'])
    
    
    
    return melt_df








# %%
# make a function to get the similarity matrix
# also melts the dataframe
# !!! After you do the chemical_sim makes sure to check the dataframe.uniques()
# !!! to ensure they match the dataframe.unique() from the dataframe made in log_transform

def chemical_sim(input=''):
    from scipy.spatial.distance import cdist
    # reads in the chemical similarty fingerprints
    fingers = pd.read_csv(input, index_col=0)
    
    # maybe add in something to change column names if they don't match
    
    
    # calculates the jacccard or Tanimoto similarity
    sim_mat = pd.DataFrame(data=1. - cdist(fingers, fingers, metric='jaccard'), index=fingers.index, columns=fingers.index)

    # rename index to imol_compound
    sim_mat.index.name = 'imol_compound'
    sim_mat.reset_index(inplace=True)

    # find a list of the columns just removing the imol_compound
    sim_columns = list(sim_mat.columns)
    sim_columns.remove('imol_compound')
    # print(sim_columns)

    # melt the sim_mat data to show the imol_compound with jmol_compound and drop duplicates
    sim_mat_melt = pd.melt(sim_mat, id_vars='imol_compound', var_name='jmol_compound', value_name='similarity')
    sim_mat_melt = sim_mat_melt[sim_mat_melt.imol_compound != sim_mat_melt.jmol_compound].copy()

    # I have to make pair id to drop the duplicates
    sim_mat_melt['pair_id'] = (sim_mat_melt
                          .apply(lambda x: "|".join(sorted([str(x['imol_compound']),
                                                            str(x['jmol_compound'])]
                                                          )
                                                   ),
                                 axis=1))                                                        


    sim_mat_melt.sort_values(by='pair_id',inplace=True)
    sim_mat_melt.drop_duplicates(subset=('pair_id'),inplace=True)
    
    return sim_mat_melt

# %%
    



















# %%
# ranks the compound from the chemical_sim function
# sim_data is the melted dataframe from chemical_sim
# num is the maximum number of matches you want
# target is the particular compound you want to rank, if none than all compounds will be ranked
# threshold is the similarity threshold that you want to select
# this function assumes the data from the jaccard or tanimoto DataFrame has been made
# the DataFrame is also assumed to be melted
def rank_compounds(sim_data=None, num=5, target=None, threshold=0.85):
    # this is the DataFrame that is going to be returned
    # should have compound name, pair, and similarity score
    output = pd.DataFrame()
    
    # list of compounds in the data that need to looped through
    unique_compounds = list(sim_data['imol_compound'].unique())
    
    for compound in unique_compounds:
        # creates a mask for each unique compund
        # searches imol and jmol for said compound
        imol_mask = sim_data['imol_compound'] == compound
        jmol_mask = sim_data['jmol_compound'] == compound
        final_mask = imol_mask | jmol_mask
        # print(sim_mat_melt[final_mask])
        
        # creates a series for the top 5 similarity matches and stores them in knn
        if output.empty:
            output = sim_data[final_mask].nlargest(n=num, columns='similarity')
            output.insert(3, 'compound', compound, True)
        else:
            new_output = sim_data[final_mask].nlargest(n=num, columns='similarity')
            new_output.insert(3, 'compound', compound, True)
            output = pd.concat([output, new_output], ignore_index=True)
    
    # removes imol and jmol compounds with pair
    mask1 = output['compound'] != output['imol_compound']
    mask2 = output['compound'] != output['jmol_compound']
    output['pair'] = None
    output['pair'].where((mask1), output['jmol_compound'].values, inplace=True)
    output['pair'].where((mask2), output['imol_compound'].values, inplace=True)
    
    # sorts the output DataFrame and resets the index
    output.sort_values(by='compound', inplace=True)
    output = output[['compound', 'imol_compound', 'jmol_compound', 'similarity']]
    output.reset_index(inplace=True)
    output.drop('index', axis=1, inplace=True)
    
    # sorts out the imol and jmol compounds to just compound and the pair.
    mask1 = output['compound'] != output['imol_compound']
    mask2 = output['compound'] != output['jmol_compound']
    output['pair'] = None
    output['pair'].where((mask1), output['jmol_compound'].values, inplace=True)
    output['pair'].where((mask2), output['imol_compound'].values, inplace=True)
    output = output[['compound', 'pair', 'similarity']]
    
    # add something to set the threshold
    mask = output['similarity'] >= threshold
    output = output[mask]
    
    output.sort_values(by=['compound', 'similarity'], inplace=True, ascending=[True, False])
    

    
    if target != None:
        mask = output['compound'] == target
        return(output[mask].sort_values(by='similarity'))
    else:
        return(output)







# %%
# function to find the local line
# assumes that the data has already been log transformed and melted
# data is the dataframe from log_transform()
# sim_data is the melted dataframe from chemical_sim
# num is the maximum number of matches you want
# target is the particular compound you want to rank, if none than all compounds will be ranked
# threshold is the similarity threshold that you want to select
from scipy.stats import linregress

def find_line(target=None, data=None, sim_data=None, num=5, threshold=0.85, show_match=False, plot=True):
    
    # finds the x axis and a list of all the compounds
    x = data['concentration'].unique()

    # creates a dataframe with all the compounds filtered out by thresholds and total matches
    # use this dataframe to choose which compounds to include in the scatterplot and subsequent regression line
    rank_df = rank_compounds(sim_data=sim_data, target=target, num=num, threshold=threshold)
    
    compounds = list(rank_df['pair'].unique())
    # print(compounds)
    
    # a mask using a list of compounds in rank_df to filter out the dataframe in data
    mask = data['compound'].isin(compounds)
    
    # shows a numpy array of the compounds that were taken from rank_df and used as a filter in data
    if show_match == True:
        print(rank_df)
    
    # finds everything necessary for the regression line
    slope, intercept, r_value, p_value, std_err = linregress(data[mask]['concentration'], data[mask]['area'])
    
    # calculates y for all values of x using y = mx + b linear equation
    y = slope * x + intercept
    
    if plot ==True:
        import matplotlib.pyplot as plt
        import seaborn as sns

        # plot of all the data points
        sns.relplot(x='concentration', y='area', hue='compound', data=data[mask])

        # plots the regression line
        # in this case this is the line for c8
        l_line, = plt.plot(x, y, label='Local y = {}x + {}'.format(slope.round(6), intercept.round(4)))
        plt.legend(handles=[l_line])
        plt.show()
    
    
    
    return y, slope, intercept, r_value, p_value, std_err
    
    # return 0
    
    
    # want to get a list of compounds from sim_data


# %%
# finds all the regression information of each compound and outputs a DataFrame
# data is the dataframe from log_transform
from scipy.stats import linregress

# Now that we can make lines we want to get all the regression lines for each compound
def make_reg_df(data=None):
    # store all the compouds in a list
    reg_compounds = data['compound'].unique()
    # print(reg_compounds)
    
    # empty dataframe to hold all the necessary values
    # still need slope, intercept, r-value, p-value, std_err, and r^2
    reg_df = pd.DataFrame({'compound': [], 'slope': [], 'intercept': [], 
                            'r_value': [], 'p_value': [], 'std_err': []})
    
    # loops through each unique compound and finds the regression line
    # information for each compound individually
    # temp is going to be a list of dictionaries
    # the loop will go through each compound and create a temporary dictionary
    # this dict will act as a row and be appended to temp
    # the dataframe will be created using this list of dictionaries at the end
    temp = []
    for compound in reg_compounds:
        # finding the regression line information for 1 compound
        mask = data['compound'] == compound
        slope, intercept, r, p, std_err = linregress(data[mask]['concentration'],
                                                        data[mask]['area'])
                                                        
        # the temporary dictionary
        dict = {'compound': compound, 'slope': slope, 'y-intercept': intercept,
                    'r_value': r, 'p_value': p, 'std_err': std_err}
        
        # appends the dictionary to temp
        temp.append(dict)
    
    reg_df = pd.DataFrame(temp)
    reg_df['r^2'] = reg_df['r_value'] ** 2
    reg_df.set_index('compound', inplace=True)
    reg_df.index.name = None
    
    return reg_df

# %%
# function to calculate the percent difference for each compound
# the output will be a dictionary with the keys being the name of each
# compound and the values will be the dataframe
# if concat is true than the return is one huge dataframe with every compound's
# percent difference
# the reason I would want all the dataframe together is to make a histogram
# of the percent error
# data is the dataframe from log_transform, reg_df is the dataframe from make_reg_df, and concat joins all the
# dataframe in the dicitonary to create 1 dataframe if false. if true than returns dictionary with key being the compound
# and values as the DataFrame
def find_pdiff(data=None, reg_df=None, concat=False):
    # list of all the unique compounds in the data set
    compounds = list(data['compound'].unique())
    
    # dictionary to contain all the dataframes of each compound
    dict = {compound: pd.DataFrame() for compound in compounds}
    
    for compound in compounds:
        # mask to only look at the current compound
        mask = data['compound'] == compound
        
        # these are the acutal log values of the concentration and area1
        log_x = data[mask]['concentration']
        log_y = data[mask]['area']
        
        # using the reg_df we can find the predicted concentration value
        # keep in mind this is still log(x)
        # so 10^pred_x = x
        pred_x = (log_y - reg_df.loc[compound]['y-intercept']) / reg_df.loc[compound]['slope']
        pred_x = 10 ** pred_x
        
        # now that we have the predicted x value we should find the acutal x values
        # should be the same as what we made the concentration to be
        # sort of a double check?
        actual_x = 10 ** log_x
        
        # calculating the percent difference
        diff = ((pred_x - actual_x) / actual_x) * 100
        
        # making the final dataframe to append to the dictionary
        dict[compound] = pd.DataFrame({'actual_x': actual_x, 'pred_x':pred_x,
                                        'percent_difference': diff, 'compound': compound})
    
    if concat == True:
        error_df = pd.concat(dict.values(), ignore_index=True)
        return error_df
    
    else:
        return dict
        
        
        
        
        
        
        
        
        
        
        
        
        
# %%
# function to find the root mean square of each compound versus
# other compound's calibration curve
# rmse_dict will have the compound as the key and the value will 
# be a dataframe of all the rmse values for that particular compound or key
# data is dataframe from log_transform() and reg_df is the dataframe from make_re_df()
from math import sqrt
def find_rmse(data=None, reg_df=None):
    # dictionary to store all the dataframes
    rmse_dict = {}
    
    # list of all the compounds
    compounds = list(data['compound'].unique())
        
    
    # x stores the x axis 
    # this is a numpy array... need to turn this into a pandas series
    x = data['concentration'].unique()
    x = pd.Series(x)
    # print(type(x))
    # print(x)
    
    
    # loops through each target compound making the rmse calculations
    for target_compound in compounds:
        target_mask = data['compound'] == target_compound
        
        # we want to get the regression line info for everything
        # not the target compound
        reg_mask = reg_df.index != target_compound
        # print(reg_df[reg_mask])
        
        # temporary list to store the temporary dicitonaries that
        # are made from the 2nd for loop
        temp_list = []
        
        # These are rmse for each compound using the intensity data
        # of the target_compound
        for compound in reg_df[reg_mask].index:
            # getting the slope and intercept of the non-target compounds
            slope = reg_df[reg_mask]['slope'].loc[compound]
            intercept = reg_df[reg_mask]['y-intercept'].loc[compound]
            
            
            
            # now we want to calculate the intensity that was predicted
            # from the non-target compounds and use that for our rmse 
            # calculations for the targed compound
            # stored into log_y_pred
            log_y_pred = (slope * x) + intercept
            # print(slope)
            # print(log_y_pred)
            
            # completing the rmse calculations
            # 1st getting the area from the target compound - predicted area of non target compound
            rmse = data[target_mask]['area'].subtract(log_y_pred.values, level='int')
            rmse = rmse ** 2
            rmse = rmse.sum() / x.count()
            rmse = sqrt(rmse)
            # print(rmse)
            
            # temporary dicitonary to store all the rows that we are making
            # then at the end of the 2nd for loop, combine them into
            # one DataFrame and store that in rmse_dict
            # faster than creating temporary dataframes 
            temp_dict = {'target_compound': target_compound, 'pair_compound': compound,
                            'rmse': rmse}
        
            temp_list.append(temp_dict)
        
        # print(temp_list)    
        
        # df stores the DataFrame temporary, so that the columns
        # can be rearranged
        # appending the dataframe to rmse_dict
        df = pd.DataFrame(temp_list)
        df = df[['target_compound', 'pair_compound', 'rmse']]
        rmse_dict[target_compound] = pd.DataFrame(df) 
        
        # concats the entire dictionary into a DataFrame
        output = pd.concat(rmse_dict.values())
    return output














# %%
# function to create a pair_id column in the DataFrame
# quality of life function to allow working with the data easier
# this allows the target_compound and pair_compounds to remain identifiable
# rmse_df is the dataframe from the find_rmse function
def pair_id(rmse_df=None):
    output = rmse_df.copy()
    output['pair_id'] = (output
                          .apply(lambda x: "|".join(sorted([str(x['target_compound']),
                                                            str(x['pair_compound'])]
                                                          )
                                                   ),
                                 axis=1))
    
    return output













# %%
def hii():
    return(3)
# %%















