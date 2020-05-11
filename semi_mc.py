# %%
# this part just loads in the function from a separate file
# the hii() function is just to check that this worked
import SemiQ_Functions as sqf
dir(sqf)

del sqf

import SemiQ_Functions as sqf
import importlib
importlib.reload(sqf)

dir(sqf)
sqf.hii()








# %%
# my preprocessing that I did to the data.
# always do this after revise_data() to make sure the column names are correct and that the concentration column is correct

import numpy as np
import pandas as pd
import re

# reads in the data that was made from format_data.py
path = '../Semi_Quant/rerun_Acids.xlsx'
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

# resets the index
new_data = data.copy()
new_data.reset_index(inplace=True)
new_data.drop('index', axis=1, inplace=True)

# going to store all the compounds in a list
compounds = list(data.columns)
compounds.remove('concentration')
# print(compounds)

new_data




# %%
new_data_melt = sqf.log_transform(new_data)

# checking to ensure the columns from sim_mat_melt later on match up
new_data_melt.replace('perfluoro(4-methoxybutanoic…', 'perfluoro(4-methoxybutanoic)acid', inplace=True)
new_data_melt.replace('perfluoro-3,6,-dioxadecanoi…', 'perfluorodioxadecanoic_acid',inplace=True)
print(new_data_melt['compound'].unique())

new_data_melt






# %%
# I can load in the chemical fingerprints using the chemical_sim function

path = '../Semi_Quant/Toxprint_Acids_revised.csv'
sim_mat_melt = sqf.chemical_sim(path)

sim_mat_melt.reset_index(drop=True, inplace=True)
# print(sim_mat_melt['imol_compound'].unique())

sim_mat_melt












# %%
# time to rank the compounds. starting with top5
# when we make the monte carlo, we probably should loop through thresholds

rank_df = sqf.rank_compounds(sim_mat_melt,target='c8', threshold=.85)


rank_df





# %%
# finds the linear regression for c8
y, slope, intercept, r_value, p_value, std_err = sqf.find_line(target='c8', data=new_data_melt, sim_data=sim_mat_melt, show_match=True)

# import matplotlib.pyplot as plt
# import seaborn as sns
# 
# # plot of all the data points
# sns.relplot(x='concentration', y='area', hue='compound', data=new_data_melt)
# 
# # all the unique concentrations for the x-axis
# x = new_data_melt['concentration'].unique()
# 
# # plots the local regression line
# # in this case this is the line for c8
# l_line, = plt.plot(x, y, label='Local y = {}x + {}'.format(slope.round(6), intercept.round(4)))
# plt.legend(handles=[l_line])
# plt.show()











# %%
# makes the regression dataframe containing all the regression line
# information for each compound
reg_df = sqf.make_reg_df(new_data_melt)
reg_df.loc['c4']['p_value']
reg_df





# %%
p_df = sqf.find_pdiff(data=new_data_melt, reg_df=reg_df, concat=True)
print(p_df)

# example of the histogram of percent difference
# right tailing
sns.distplot(p_df['percent_difference'], bins=30, kde=False)







# %%
rmse_df = sqf.find_rmse(data=new_data_melt, reg_df=reg_df)


rmse_df

# %%
real_rmse = sqf.pair_id(rmse_df=rmse_df)

real_rmse













# %%




















