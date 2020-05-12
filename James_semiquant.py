# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:50:35 2020

@author: James
"""
#%%
import numpy as np
import pandas as pd
from orbi_reformat import orbi_reformat

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import cdist

#%% Import from file
file = "Semi_Quant_Acids_Revised.xls"

raw_df = orbi_reformat(file)

#%% Turn filenames into meaningful concentration values
raw_df['concentration'] = pd.to_numeric(raw_df['Filename'].replace(to_replace="ng",value =""), errors = "coerce")

raw_df.dropna(subset=['concentration'], inplace = True)

# Rename to short names
raw_df["Component Name"].replace(to_replace = ["PFOA","GenX","Perfluoro(4-methoxybutanoic)acid","Perfluoro-3,6,-dioxadecanoic_Acid"],
                                 value = ["C8","HFPO-DA","PFMOBA","PFPE-5"],
                                 inplace = True)

#%% Basic plotting to check values

raw_df['log10_conc'] = np.log10(raw_df["concentration"])
raw_df['log10_area'] = np.log10(raw_df["Area"])

sns.scatterplot(x = 'concentration', y = "Area", hue= 'Component Name', data =raw_df).legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
sns.scatterplot(x = 'log10_conc', y = "log10_area", hue= 'Component Name', data =raw_df).legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)


#%% now lets load in the chemical fingerprints from ChemoTyper

# load in the fingerprint data
fingers = pd.read_csv('../Semi_Quant/Toxprint_Acids_revised.csv', index_col=0)
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

