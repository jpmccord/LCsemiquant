import numpy as np
import pandas as pd

# %%
path = 'Semi_Quant_Acids_Revised.xls'

data = pd.read_excel(path, sheet_name=None)
temp_data = data.copy()

# list of all the excel sheets
sheets = list(temp_data.keys())
sheets.remove('Component')

# the dataframe to store all the revised data
revised_data = pd.DataFrame()
revised_data['compound'] = sheets
revised_data['rt'] = np.nan

# revised_data['compound'].drop(, inplace=True)
# revised_data

# %%

# now I need to loop data to access each sheet and get the retention time
# then I need to take that retention time and add it to revised_data
# have to figure out what to use retention time for because
# they seem to be all over the place with no similarities


for sheet in sheets:
    # formats the sheet to get the necessary columns
    names = temp_data[sheet].iloc[3]
    # data[sheet].drop([0, 1, 2, 3], inplace=True)
    temp_data[sheet].columns = names
    temp_data[sheet].columns.name = None
    temp_data[sheet].dropna(subset=['Filename'], inplace=True)
    temp_data[sheet].rename(columns={'Filename' : 'Concentration'}, inplace=True)
    temp_data[sheet].reset_index(inplace=True)
    temp_data[sheet].rename(columns = {'Area' : sheet}, inplace=True)
    
    # takes the rt value in temp_data and drops all the Nan values
    # also finds the average of all the rt
    # need to drop the double blanks
    pattern = '\d*ng'
    filter = temp_data[sheet]['Concentration'].str.contains(pattern)
    # print(filter)
    rt = temp_data[sheet][filter]['RT']
    rt.replace('NF', np.nan, inplace=True)
    # rt.dropna(inplace=True)
    # rt.drop(index=2,axis=0, inplace=True)
    # print(temp_data[sheet][filter])
    rt = rt.mean()
    revised_data.loc[revised_data['compound'] == sheet, 'rt'] = rt

revised_data['compound'] = revised_data['compound'].str.lower()
mask = revised_data['compound'] == 'pfmoaa'
revised_data.drop(list(revised_data.loc[mask].index), inplace=True)
revised_data.loc[revised_data['compound'] == 'pfoa', 'compound'] = 'c8'
revised_data


# %%
# outputs the data to specific output
output = '../Semi_Quant/rt_df.csv'
revised_data.to_csv(output, index=False)
# writer.save()























