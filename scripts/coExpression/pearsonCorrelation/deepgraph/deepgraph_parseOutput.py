#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description='Compute large correlation matrices in parallel with DeepGraph', add_help=True)
parser.add_argument('-i', '--input', dest='input', metavar='file', help='Quantification matrix', required=True)

args = parser.parse_args()
infile = args.input

# outfile = input without ".txt" + _newtwork.txt
outfile = infile[:-4] + '_newtwork.txt'

X_df = pd.read_csv(infile, sep="\t")
Y_df = X_df.set_index("Name")

# store correlation values in an hdf file
files = os.listdir('tmp/correlations/')
files.sort()
store = pd.HDFStore('e.h5', mode='w')

for f in files:
    et = pd.read_pickle('tmp/correlations/{}'.format(f))
    store.append('e', et, format='t', data_columns=True, index=False)
store.close()

# load correlation table
print('Loading correlation table ...')
e = pd.read_hdf('e.h5')

e_reset = e.reset_index()
e_reset.replace({"s": X_df.to_dict()['Name']}, inplace=True)
e_reset.replace({"t": X_df.to_dict()['Name']}, inplace=True)

e_reset_dict = {}

for index, row in e_reset.iterrows():
    if row['s'] in e_reset_dict.keys():
        e_reset_dict[row['s']][row['t']] = row['corr']
    else:
        e_reset_dict[row['s']] = {}
        e_reset_dict[row['s']][row['t']] = row['corr']

# save output
print('Saving correlation as pcc format')
with open(outfile, 'w') as out:
    for gene in e_reset_dict.keys():
        out.write(f'{gene}:')
        for gene2 in e_reset_dict[gene]:
            out.write(f'\t{gene2}({e_reset_dict[gene][gene2]}) ')
        out.write("\n")
