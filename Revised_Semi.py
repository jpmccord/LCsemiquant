# %%
import numpy as np
import pandas as pd
import re

# reads in the data that was made from format_data.py
path = 'rerun_Acids.xlsx'
data = pd.read_excel(path)

# adds in NaN values
data.replace('NF', np.nan, inplace=True)

# keep a old data copy just in case
old_data = data.copy()

# name all the columns into lower case and store them in a list
# fix the name of the c5 and c8 ethers
print(data.columns)
columns = list(data.columns.str.lower())
data.columns = columns



# going to fill Nan with 0 for now
# the pefluoro(methoxybutanoic)acid is questionable and could be removed???
# data.fillna(0, inplace=True)

# need to replace pfoa with 'c8' and fix column names
data.rename(columns={'pfoa' : 'c8'}, inplace=True)
data.rename(columns=lambda x: re.sub('_ݞ\d*', '...', x), inplace=True)

# make the concentration column to floats and sort the column
data.concentration.replace(regex='ng', value='', inplace=True)
data['concentration'] = data['concentration'].astype('float')
data.sort_values(by='concentration', inplace=True)
data.drop('pfmoaa', axis=1, inplace=True)

new_data = data.copy()
new_data.reset_index(inplace=True)
new_data.drop('index', axis=1, inplace=True)


new_data




# %%
# now lets load in the chemical fingerprints from ChemoTyper
from scipy.spatial.distance import cdist

# load in the fingerprint data
fingers = pd.read_csv('Toxprint_Acids_revised.csv', index_col=0)
fingers.rename(index=lambda x: re.sub('perfluorodioxadecanoic_acid', 'perfluoro-3,6,-dioxadecanoi…', x), inplace=True)
fingers.rename(index=lambda x: re.sub('perfluoro\(4-methoxybutanoic\)acid', 'perfluoro(4-methoxybutanoic…', x), inplace=True)


# calculates the jacccard or Tanimoto similarity
sim_mat = pd.DataFrame(data=1. - cdist(fingers, fingers, metric='jaccard'), index=fingers.index, columns=fingers.index)

# rename index to imol_compound
sim_mat.index.name = 'imol_compound'
sim_mat.reset_index(inplace=True)

# find a list of the columns just removing the imol_compound
sim_columns = list(sim_mat.columns)
sim_columns.remove('imol_compound')
print(sim_columns)

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
sim_mat_melt







# %%
# time to simplify the code to find the top5 similar compounds with a function
# this function assumes the data from the jaccard or tanimoto DataFrame has been made
# the DataFrame is also assumed to be melted
def rank_compounds(data=None, num=5, target=None):
    # this is the DataFrame that is going to be returned
    # should have compound name, pair, and similarity score
    output = pd.DataFrame()
    
    # list of compounds in the data that need to looped through
    unique_compounds = list(data['imol_compound'].unique())
    
    for compound in unique_compounds:
        # creates a mask for each unique compund
        # searches imol and jmol for said compound
        imol_mask = data['imol_compound'] == compound
        jmol_mask = data['jmol_compound'] == compound
        final_mask = imol_mask | jmol_mask
        # print(sim_mat_melt[final_mask])
        
        # creates a series for the top 5 similarity matches and stores them in knn
        if output.empty:
            output = data[final_mask].nlargest(n=num, columns='similarity')
            output.insert(3, 'compound', compound, True)
        else:
            new_output = sim_mat_melt[final_mask].nlargest(n=num, columns='similarity')
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
    
    if target != None:
        mask = output['compound'] == target
        return(output[mask].sort_values(by='similarity'))
    else:
        return(output)

# %%
# practicing to organize the DataFrame
# really would like to implement removing the imol and jmol compounds
# with a pair column
top5 = rank_compounds(sim_mat_melt, target=None)



# temp = rank_compounds(sim_mat_melt, target=None)
# mask1 = temp['compound'] != temp['imol_compound']
# mask2 = temp['compound'] != temp['jmol_compound']
# temp['pair'] = None
# temp['pair'].where((mask1), temp['jmol_compound'].values, inplace=True)
# temp['pair'].where((mask2), temp['imol_compound'].values, inplace=True)
# temp = temp[['compound', 'pair', 'similarity']]

top5














# %%
# we can start looking at new_data now that we can see all this
import matplotlib.pyplot as plt
import seaborn as sns


# going to make the melt new_data to make it easier to plot
new_data_melt = pd.melt(new_data, id_vars='concentration', var_name='compound', value_name='area')

# the increases across compounds are quite dramatic, so we need to do a logrithmic transformation
# however in the area there was a 0 made for some of the genx areas, so we need to change that to 1
# print(new_data_melt['concentration'])
new_data_melt['concentration'] = np.log10(new_data_melt['concentration'])
# new_data_melt.replace(0, 1, inplace=True)
new_data_melt['area'] = np.log10(new_data_melt['area'])

# plots out a individual scatter plot for each compound
# sns.relplot(x='concentration', y='area', hue='compound', data=new_data_melt,
            # legend=False, col='compound', col_wrap=4)
new_data_melt



# %%
# now that we have all the DataFrames stored which are the new_data and top5 dataframes
# we can load in a function to find the regression line slope, intercepts and y
# from Local_Regress.py import find_local
from scipy.stats import linregress
#from Local_Regress import find_local

def find_local(compound=None, x=np.nan, data=None):
    if compound == None:
        print('No Target Compund Inputed')
        return 0
    # example target
    target = compound
    mask = top5['compound'] == target
    # print(mask)

    # for loop to find all the pairings with the target compound
    pairings = []
    for n in range(5):
        if top5[mask]['imol_compound'].iloc[n] == target and top5[mask]['similarity'].iloc[n] >= 0.85:
            pairings.append(top5[mask]['jmol_compound'].iloc[n])
            # print(pairings)
        elif top5[mask]['similarity'].iloc[n] >= 0.85:
            pairings.append(top5[mask]['imol_compound'].iloc[n])
        else:
            pass
            # print(pairings)
    # print(top5[mask])
    # now pairings is a list containing the 5 mosst similar compounds
    # I want to take only columns of new_data that match the pairings
    # and put all the Area values into another DataFrame called pair_melt
    
    # keep the concentration index
    pairings.append('concentration')
    
    # this is list of all compounds we want the area for to use in the pd.melt()
    values = list(data.columns)
    values = values.remove('concentration')

    pair_melt = pd.melt(new_data[pairings], id_vars='concentration', value_vars=values,
                        value_name='area', var_name='compound')
    # print(pair_melt) 

    # this adds 
    pair_log = pair_melt.copy()
    pair_log['concentration'] = np.log10(pair_melt['concentration'])
    pair_log['area'] = np.log10(pair_melt['area'])
    # print(pair_log)
    
    # yay now we have the pair_log data of the top 5 nearest neighbors to our target compound
    # time to make a linear regression line
    slope, intercept, r_value, p_value, std_err = linregress(pair_log['concentration'], pair_log['area'])
    # print(slope, intercept, r_value, p_value, std_err)

    y = slope * x + intercept
    # line = plt.plot(x, y, label='Local y = {}x + {}'.format(slope, intercept))
    # print(y)
    
    return y, slope, intercept, r_value, p_value, std_err

# %%







# %%
# we want to find the global and linear regression lines and plot again
# finding the global linear regression line
slope, intercept, r_value, p_value, std_err = linregress(new_data_melt['concentration'], new_data_melt['area'])

# print(slope, intercept, r_value, p_value, std_err)

# rounding off the concentration values
# new_data_melt['concentration'] = new_data_melt['concentration'].round(2)


# now to make the global regression line
# x = np.arange(new_data_melt['concentration'].min(), new_data_melt['concentration'].max(), 0.25)
x = new_data_melt['concentration'].unique()
y = (slope * x) + intercept


# want to find a local regression line for c6
local_y, local_slope, local_intercept, local_p, local_r, local_std_error = find_local('c8', x, data=new_data_melt)


plt.figure(figsize=(20, 20))

sns.relplot(x='concentration', y='area', hue='compound', data=new_data_melt)
g_line, = plt.plot(x, y, label='Global y = {}x + {}'.format(slope.round(6), intercept.round(4)))
l_line, = plt.plot(x, local_y, label='Local y = {}x + {}'.format(local_slope.round(6), local_intercept.round(4)))
plt.legend(handles=[g_line, l_line])
plt.show()



# there is something wrong with the find_local function
# doesn't find new slope, y, or intercept

 



# %%
# In order to validate the model and figure out the best regression line
# we need to figure out the error values for each compound's regression line
# first things first find all the compund's regression line and put them on a table
reg_df = pd.DataFrame({'compound': [], 'slope': [], 'intercept': [], 
                        'r_value': [], 'p_value': [], 'std_err': []})

# lists all the columns in new_columns
new_columns = list(new_data.columns)
new_columns.remove('concentration')


# finds all the linear regression line variables for all the compounds individually
for compound in new_columns:
    mask = new_data_melt['compound'] == compound
    local_slope, local_intercept, local_r, local_p, local_std_error = linregress(new_data_melt[mask]['concentration'], new_data_melt[mask]['area'])
    temp = pd.DataFrame({'compound': [compound], 'slope': [local_slope], 'intercept': [local_intercept],
                    'r_value': [local_r], 'p_value': [local_p], 'std_err': [local_std_error]})
    reg_df = reg_df.append(temp, ignore_index=True)
    # print(new_data_melt[mask])

reg_df.set_index('compound', inplace=True)
reg_df.index.name = None

reg_df['r^2'] = reg_df['r_value'] ** 2

# r values look really really off
reg_df

# %%
# checking the r_values of all these lines
target = 'c8'
mask = new_data_melt['compound'] == target
c8_y = (reg_df['slope'].loc[target] * x) + reg_df['intercept'].loc[target]

plt.figure(figsize=(20,20))
sns.relplot(x='concentration', y='area', data=new_data_melt[mask], legend='brief')
test, = plt.plot(x, c8_y, label=target)
plt.legend(handles=[test])
plt.title(target)
plt.xlabel('Log(concentration)')




# root mean squared error
# incorporate retention time
# spearman rank coefficient


# %%
# manuelly check r_value
# pearson correlation coefficient
from numpy import cov
from scipy.stats import pearsonr

target = 'c8'
mask = new_data_melt['compound'] == target
covariance = cov(new_data_melt[mask]['concentration'], new_data_melt[mask]['area'])

p_r_value = pearsonr(new_data_melt[mask]['concentration'], new_data_melt[mask]['area'])
print(p_r_value)


















# %% 
# now I need to remake the linear equation for each regression line
# and plug in the respective area values for y and find the concentration

# store all dataframes in a dictionary
dict = {compound: pd.DataFrame() for compound in reg_df.index}

# the range of x values only having the concentrations of interest
new_x = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000]
new_x = np.log10(new_x)

for compound in reg_df.index:
    # mask to find targeted compound
    mask = new_data_melt['compound'] == compound
    
    # variables to hold the log data
    log_x = new_data_melt[mask]['concentration']
    log_y = new_data_melt[mask]['area']
    # print(log_y)
    
    # undo the log to see actual concentration and areas
    actual_x = 10 ** log_x
    actual_y = 10 ** log_y
    
    # writing a function to find log(x) when there is log(y)
    pred_x = 10 ** ((log_y - reg_df.loc[compound]['intercept']) / reg_df.loc[compound]['slope'])
    # print(x)
    
    # percent difference between pred_x and actual_x
    diff = (pred_x - actual_x) / actual_x * 100
    
    
    
    new_df = pd.DataFrame({'actual_x': actual_x, 'actual_y': actual_y, 
                            'pred_x': pred_x, 'percent_difference': diff, 'compound': compound})
    
    # stores this new dataframe inside the dictionary and calculates the actual x
    dict[compound] = dict[compound].append(new_df, ignore_index=True)
    
    
    # time to calculate the predicted x using log_y
    
dict['c7']    








# %%
# now that I have all the percent differences. 
# I can concat all the dataframes together into one big DataFrame to plot
error_df = pd.concat(dict.values(), ignore_index=True)

# plt.figure()
# sns.relplot(x='actual_x', y='percent_difference', hue='compound', data=error_df)
# plt.show()
# plt.close()

plt.figure()
sns.distplot(error_df['percent_difference'], bins=30, kde=False)
plt.show()
plt.close()





# %%
# now I want to find the root mean squared error for each compound versus
# other calibration curves
# rmse_dict is going to store each compound's rmse versus all the other compounds
# key = compound, value = DataFrame of rmse for other compounds
from math import sqrt
rmse_dict = {compound: pd.DataFrame() for compound in reg_df.index}

for target_compound in reg_df.index:
    # x axis
    new_x = pd.Series([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000])
    new_x = np.log10(new_x)
    # print(new_x)

    # list of all the compounds
    reg_columns = list(reg_df.index)
    length = len(reg_columns) - 1

    # mask to find the compound of interest
    mask = new_data_melt['compound'] == target_compound
    reg_mask = reg_df.index == target_compound

    # temporary dataframe to store the rmse values for all the other compounds
    temp = pd.DataFrame()

    # calculates the root mean squared value for the targeted compound
    # versus all the other compound's calibration curve
    # appends the dataframe onto temp
    for compound in reg_df[~reg_mask].index:
        # getting the correct slope and intercept
        slope = reg_df[~reg_mask]['slope'].loc[compound]
        intercept = reg_df[~reg_mask]['intercept'].loc[compound]
        
        # calculates the predicted y
        log_y_pred = (slope * new_x) + intercept
        
        # calculates the roots mean squared values
        rmse = new_data_melt[mask]['area'].subtract(log_y_pred.values, level='int')
        rmse = rmse ** 2
        rmse = rmse.sum() / new_x.count()
        rmse = sqrt(rmse)
        row = pd.DataFrame([[compound, rmse]], columns=['compound', 'rmse'])
        
        # appends to the temporary DataFrame
        temp = temp.append(pd.DataFrame([[compound, rmse]]), ignore_index=True)
        
    # renames and attaches the DataFrame into the correct dictionary key-value pair
    temp.rename(columns={0: 'compound', 1: 'rmse'}, inplace=True)
    temp['target_compound'] = target_compound
    temp = temp[['target_compound', 'compound', 'rmse']]

    rmse_dict[target_compound] = rmse_dict[target_compound].append(temp)

rmse_df = pd.concat(rmse_dict.values())
rmse_df
    


# %%
# how about doing some plotting for these rmse values
mask = rmse_df['compound'] == 'c8'
plt.figure(figsize=(20,20))
sns.catplot(x='target_compound', hue='compound', y='rmse',kind='bar', data=rmse_df[mask],height=7, aspect=3)







# %%
# rmse vs similarity
compounds = list(reg_df.index)
simil_df = rank_compounds(sim_mat_melt, num=9, target=None)

# time to make a new DataFrame combining both simil_df and rmse_df
# print(simil_df.columns)

# new column to identify the pairs in order to merge them
# I have to make pair to drop the duplicates
# this method works because even if the imol and jmol compounds swapped,
# the chemical similarity is the same
simil_df['pair'] = (simil_df
                      .apply(lambda x: "|".join(sorted([str(x['imol_compound']),
                                                        str(x['jmol_compound'])]
                                                      )
                                               ),
                             axis=1))    

simil_df.sort_values(by='pair',inplace=True)
simil_df.drop_duplicates(subset=('pair'),inplace=True)
simil_df.reset_index(drop=True, inplace=True)


# now to do the same on the rmse_df
rmse_df['pair'] = (rmse_df
                      .apply(lambda x: "|".join(sorted([str(x['target_compound']),
                                                        str(x['compound'])]
                                                      )
                                               ),
                             axis=1))
                                 
rmse_df.sort_values(by='pair',inplace=True)

# I don't want to drop duplicates because the target 
# and pair_compounds matter when determining rmse
# resets the index
rmse_df.reset_index(drop=True, inplace=True)

# merges the rmse_df and simil_df and picks out the important columns
# the completed df is saved as rvs_df
# didn't remove duplicates because of there is a change in rmse
# if the target-pair changes
rvs_df = rmse_df.merge(simil_df, left_on='pair', right_on='pair', suffixes=('_left', '_right')).copy()
rvs_df = rvs_df[['target_compound', 'compound_left', 'rmse', 'similarity', 'pair']].copy()
rvs_df.rename(columns={'compound_left': 'pair_compound'}, inplace=True)
rvs_df.sort_values(by='target_compound', inplace=True)
rvs_df.reset_index(drop=True, inplace=True)
rvs_df








# %%
# Time to create a scatter plot to visualize the correlation
sns.relplot(x='similarity', y='rmse', hue='target_compound', data=rvs_df)














# %%
# # want a function to find the percent different in a curve
# def find_error():
    








# %%
# now I want to incorporate retention time as a feature when choosing similar compounds
# loading in the retention time data
rt_df = pd.read_csv('./rt_df.csv')


# need to find delta t of each compounds retention time with each other
new_rt = pd.DataFrame()
temp = rvs_df.copy()
rt_df.loc[6, 'compound'] = temp.loc[2, 'pair_compound']
rt_df.loc[7, 'compound'] = temp.loc[1, 'pair_compound']
temp.drop('rmse', axis=1, inplace=True)

# finds the target compound's rt and the pair compounds rt
# then finds dt doing dt = target_rt - pair_rt
rt_df.rename(columns={'compound': 'target_compound'}, inplace=True)
temp = temp.merge(rt_df, on='target_compound')
temp = temp.merge(rt_df, left_on='pair_compound', right_on='target_compound', suffixes=('_pair', '_target'))
temp['dt'] = temp['rt_target'] - temp['rt_pair']
temp = temp.rename(columns={'target_compound_pair': 'target_compound'})

# moves the dt over to the rvs_df and makes a final_df as a copy
rvs_df = rvs_df.merge(temp[['target_compound', 'pair_compound', 'dt']], on=['target_compound', 'pair_compound'])
final_df = rvs_df.copy()
final_df


# %%
# Time to make a scatter plot
# creates a temp DataFrame to plot
# removes pair duplicates
# absolute value on dt to 
temp = final_df.copy()
# temp.drop_duplicates(subset='pair', inplace=True)
temp['dt'] = temp['dt'].apply(abs)


sns.relplot(x='similarity', y='dt', hue='target_compound', data=temp)






sns.relplot(x='rmse', y='dt', hue='target_compound', data=temp)
plt.show()

























