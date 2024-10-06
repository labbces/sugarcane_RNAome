#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(prog="lncRNAsCVDistribution.py", description="Generate a histogram showing CV distribution of annotated lncRNAs", add_help=True)
parser.add_argument("-lncRNAs", metavar="Hoang_lncRNAs.tsv", dest="lncRNAs", help="dataset lncRNAs genes", type=str, required=True)
parser.add_argument("-CV", metavar="Hoang2017_counts_filters_VST_CV.txt", dest="CV", help="dataset CV distribution", type=str, required=True)
parser.add_argument("-out", metavar="Hoang_lncRNAs_CV_Distribution.png", dest="output", help="dataset lncRNAs CV distribution plot", type=str, required=True)
args = parser.parse_args()

lncRNAs = args.lncRNAs
CV = args.CV
output = args.output

genes = pd.read_csv(lncRNAs, header=None, names=["Gene"])
cv_file = pd.read_csv(CV, sep="\t")

merged_data = pd.merge(genes, cv_file, on="Gene")

cv_values = merged_data['CV']

plt.figure(figsize=(10, 6))
sns.histplot(cv_values, color='skyblue')
plt.title('Histograma dos valores de CV para os genes lncRNA', fontsize=16)
plt.xlabel('Coeficiente de Variação (CV)', fontsize=14)
plt.ylabel('Frequência', fontsize=14)
plt.tight_layout()

plt.savefig(output, dpi=300)

mean_cv = cv_values.mean()
median_cv = cv_values.median()
std_cv = cv_values.std()
min_cv = cv_values.min()
max_cv = cv_values.max()
quantiles_cv = cv_values.quantile([0.25, 0.5, 0.75])

print(f"Média do CV: {mean_cv:.4f}")
print(f"Mediana do CV: {median_cv:.4f}")
print(f"Desvio padrão do CV: {std_cv:.4f}")
print(f"Valor mínimo do CV: {min_cv:.4f}")
print(f"Valor máximo do CV: {max_cv:.4f}")
print(f"Quartis do CV: \n{quantiles_cv}")
