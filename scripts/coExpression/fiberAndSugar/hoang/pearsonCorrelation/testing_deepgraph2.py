import pandas as pd
import numpy as np
import os

X_df = pd.read_csv("/home/renato/tmp/Svi_matrix_tpm_subsample.txt", sep="\t")
Y_df = X_df.set_index("Name")
#print(Y_df.head())

#print(f'{np.corrcoef(Y_df.loc[["Sevir.3G341410.1"]].to_numpy()[0],Y_df.loc[["Sevir.3G341210.1"]].to_numpy()[0])}, expected {-0.095264}')
#print(f'{np.corrcoef(Y_df.loc[["Sevir.3G341400.1"]].to_numpy()[0],Y_df.loc[["Sevir.3G341300.1"]].to_numpy()[0])} expected: {0.004159}')

files = os.listdir('tmp/correlations/')
files.sort()
store = pd.HDFStore('e.h5', mode='w')
for f in files:
    et = pd.read_pickle('tmp/correlations/{}'.format(f))
    store.append('e', et, format='t', data_columns=True, index=False)
store.close()

# load correlation table
e = pd.read_hdf('e.h5')

e_reset = e.reset_index()
e_reset.replace({"s": X_df.to_dict()['Name']}, inplace=True)
e_reset.replace({"t": X_df.to_dict()['Name']}, inplace=True)
#print(e_reset.head())
#print(e_reset.tail())

e_reset_dict = {}

for index, row in e_reset.iterrows():
    if row['s'] in e_reset_dict.keys():
        e_reset_dict[row['s']][row['t']] = row['corr']
    else:
        e_reset_dict[row['s']] = {}
        e_reset_dict[row['s']][row['t']] = row['corr']

for gene in e_reset_dict.keys():
    print(f'{gene}:', end="")
    for gene2 in e_reset_dict[gene]:
        print(f'\t{gene2}({e_reset_dict[gene][gene2]})', end="")
    print("")
