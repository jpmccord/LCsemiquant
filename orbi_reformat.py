# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:49:07 2020

@author: James
"""

import pandas as pd

file = "NJ_targeted_test2.xlsx"

#%% Function to accept filepath as input and return a formatted pandas dataframe
def orbi_reformat(file)
# Loads in the dataset as a dict
    sheet_to_df_map = pd.read_excel(file, sheet_name = None, skiprows = 4)

# Concatenate sheets and drop extra trailing rows from the Orbitrap, assuming Exp Method is always filled for real data
# also removes empty columns
    df_raw = pd.concat(sheet_to_df_map).dropna(subset=["Exp Method"]).dropna(axis = 1, how = "all").reset_index()

    df_renamed = df_raw.rename(columns = {"level_0":"Compound", "level_1":"Run Order"})

    return df_renamed
