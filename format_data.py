import numpy as np
import pandas as pd

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
# data[sheets[0]][filter][['Filename', 'Area']]

if __name__ == '__main__':
    revise_data('Semi_Quant_Acids_Revised.xls', 'rerun_Acids.xlsx', '../Semi_Quant')


