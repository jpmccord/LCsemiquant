# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:49:07 2020

@author: James
"""
import pandas as pd

#%% Function to accept filepath to targeted Orbitrap reults as input and return a formatted pandas dataframe
def orbi_reformat(file):
# Loads in the dataset as a dict
    sheet_to_df_map = pd.read_excel(file, sheet_name = None, skiprows = 4, na_values = ["NA","NF","nan","INF"])

# Concatenate sheets and drop extra trailing rows from the Orbitrap, assuming Exp Method is always filled for real data
# also removes empty columns
    df_raw = pd.concat(sheet_to_df_map).dropna(subset=["Exp Method"]).dropna(axis = 1, how = "all").reset_index()

# Get list of compounds from cells because we apparently can't trust sheet names 
    cmp_list = pd.read_excel(file, sheet_name = None, skiprows = 1, nrows =1)
    cmp_df = pd.concat(cmp_list).reset_index().dropna(axis = 1, how = "all").drop(axis = 1, labels = "level_1")

# Merge in compound list on sheetnames then drop and rename some columns
    df_rename = pd.merge(df_raw,cmp_df, how = "left", on ="level_0").set_index("Component Name").drop(axis =1, labels = "level_0").reset_index().rename(columns = {"level_1":"Run Order"})

    return df_rename

#%% Test block
"""
file = "Semi_Quant_Acids_Revised.xls"

test_def = orbi_reformat(file)
"""
#%% 

