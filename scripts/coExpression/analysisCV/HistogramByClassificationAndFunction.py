#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="plot function and panRNAome groups of filtered expression matrix")
parser.add_argument("-cv", dest="classification_file", metavar="panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv", required=True, help="Path to the classification file")
parser.add_argument("-oc", dest="out_class", metavar="classification_CV.png", required=True, help="classification histogram")
parser.add_argument("-of", dest="out_function", metavar="function_CV.png", required=True, help="function histogram")
parser.add_argument("-b", dest="bins", type=int, default=10, help="Number of bins for the histogram")

args = parser.parse_args()

df = pd.read_csv(args.classification_file, sep='\t')

# Find min and max
cv_min = df['CV'].min()
cv_max = df['CV'].max()

# Filter and plot histogram by each classification
classifications = df['Classification'].unique()
for classification in classifications:
    plt.figure()
    df[df['Classification'] == classification]['CV'].hist(bins=args.bins)
    plt.title(f'Coefficient of variation - {classification}')
    plt.xlabel('Coefficient of variation')
    plt.ylabel('Frequency')
    plt.xlim(cv_min, cv_max)  # Definindo os limites do eixo X
    plt.savefig(f"{args.out_class}_{classification}.png")

# Filter and plot histogram by each function
functions = df['Function'].unique()
for function in functions:
    plt.figure()
    df[df['Function'] == function]['CV'].hist(bins=args.bins)
    plt.title(f'Coefficient of variation - {function}')
    plt.xlabel('Coefficient of variation')
    plt.ylabel('Frequency')
    plt.xlim(cv_min, cv_max)  # Definindo os limites do eixo X
    plt.savefig(f"{args.out_function}_{function}.png")

# Filter and plot histogram by classification
plt.figure()
for classification, group in df.groupby('Classification'):
    group['CV'].hist(alpha=0.5, label=classification, bins=args.bins)
plt.title('Coefficient of variation by classification')
plt.xlabel('Coefficient of variation')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(args.out_class)

# Filter and plot histogram by function
plt.figure()
for function, group in df.groupby('Function'):
    group['CV'].hist(alpha=0.5, label=function, bins=args.bins)
plt.title('Coefficient of variation by function')
plt.xlabel('Coefficient of variation')
plt.ylabel('Frequency')
plt.legend()
plt.savefig(args.out_function)
