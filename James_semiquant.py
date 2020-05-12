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

from scipy import stats
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from sklearn.compose import TransformedTargetRegressor
from sklearn.metrics import mean_squared_error
from math import sqrt

from scipy.spatial.distance import cdist

#%% Import from file
file = "Semi_Quant_Acids_Revised.xls"

raw_df = orbi_reformat(file)

#%% Turn filenames into meaningful concentration values
raw_df['concentration'] = raw_df['Filename'].replace(to_replace = "ng", value ="", regex = True)
raw_df['concentration'] = pd.to_numeric(raw_df['concentration'], errors = "coerce")

raw_df.dropna(subset=['concentration'], inplace = True)

# Rename to short names
raw_df["Component Name"].replace(to_replace = ["PFOA","GenX","Perfluoro(4-methoxybutanoic)acid","Perfluoro-3,6,-dioxadecanoic_Acid"],
                                 value = ["C8","HFPO-DA","PFMOBA","PFPE-5"],
                                 inplace = True)

raw_df['log10_conc'] = np.log10(raw_df["concentration"])
raw_df['log10_area'] = np.log10(raw_df["Area"])

#%% Basic plotting to check values

sns.scatterplot(x = 'concentration', y = "Area", hue= 'Component Name', data =raw_df).legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
sns.scatterplot(x = 'log10_conc', y = "log10_area", hue= 'Component Name', data =raw_df).legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

#%% Generate regression lines for compound data

cmp = "C8"

#Placeholder to create reference data for regression until we have a chooser function
ref_list = raw_df["Component Name"].unique()
ref_list = ref_list[ref_list!=cmp]
ref_set = raw_df[raw_df["Component Name"].isin(ref_list)]

#subset to just the test compound
test_cmp = raw_df[raw_df["Component Name"] == cmp]
X_test = test_cmp['log10_area'].values[:,np.newaxis]
Y_test = test_cmp['log10_conc'].values

# LinearRegression expects an array of shape (n, 1) for the "Training data"
X = ref_set['log10_area'].values[:,np.newaxis]
# target data is array of shape (n,) 
Y = ref_set['log10_conc'].values

model = LinearRegression().fit(X,Y)

rmse = sqrt(mean_squared_error(Y_test, model.predict(X_test)))

plt.scatter(X_test, Y_test,color='g')
plt.plot(X_test, model.predict(X_test),color='k')
plt.show()


#%% Chemotyper Similarity
fingers = pd.read_csv('Toxprint_Acids_revised.csv', index_col=0)

# calculates the jacccard or Tanimoto similarity
sim_mat = pd.DataFrame(data=1. - cdist(fingers, fingers, metric='jaccard'), index=fingers.index, columns=fingers.index)

sim_mat = sim_mat.reset_index()
sim_log = pd.melt(sim_mat, id_vars='index', value_name='similarity')
#%% Monte Carlo Attempt
all_cmps = raw_df["Component Name"].unique().tolist()

def model_fit_goodness(cmp,similarity):
    if cmp not in sim_log["variable"].values:
        return [cmp, similarity, None]
        
    #subset to just the test compound
    test_cmp = raw_df[raw_df["Component Name"] == cmp]
    X_test = test_cmp['log10_area'].values[:,np.newaxis]
    Y_test = test_cmp['log10_conc'].values
    
    #subset to reference list for modeling
    ref_list = sim_log[sim_log["index"] == cmp]
    ref_list = ref_list[ref_list["variable"] != cmp ]
    ref_list = ref_list[ref_list['similarity'] > similarity]
    ref_set = raw_df[raw_df["Component Name"].isin(ref_list['variable'])]
    
    if len(ref_set) == 0:
        return [cmp, similarity, None]
    
    X = ref_set['log10_area'].values[:,np.newaxis]
    Y = ref_set['log10_conc'].values
    
    model = LinearRegression().fit(X,Y)
    rms = sqrt(mean_squared_error(Y_test, model.predict(X_test)))
        
    return [cmp, similarity, rms]
#%% Test all chemotyper similarity levels for decent regression
results = [model_fit_goodness(cmp, similarity) for cmp in all_cmps for similarity in np.arange(0,1,0.1)]

results_df = pd.DataFrame(results, columns = ["cmp","similarity","rmse"])
    

#%% Plot results from similarity scaling

sns.scatterplot(x = 'similarity',
                y = "rmse",
                hue= 'cmp',
                data =results_df).legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)



#%% Transformation within LinearRegression maybe
# LinearRegression expects an array of shape (n, 1) for the "Training data"
X = test_cmp['log10_area'].values[:,np.newaxis]
# target data is array of shape (n,) 
Y = test_cmp['concentration'].values
regr_trans = TransformedTargetRegressor(regressor=LinearRegression(),
                                        func=np.log1p,
                                        inverse_func=np.expm1)

model2 = regr_trans.fit(X,Y)
y_pred = np.sort(model2.predict(X))

X_pred = np.sort(X, axis = None)

plt.scatter(X, Y,color='g')
plt.plot(X_pred, y_pred,color='k')
plt.show()

