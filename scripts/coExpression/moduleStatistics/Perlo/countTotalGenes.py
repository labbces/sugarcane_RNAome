#!/usr/bin/env python

import pandas as pd

df = pd.read_csv('moduleCNCwithlncRNAs.tsv', sep='\t')
total = df[['Genes', 'Genes_coding', 'Genes_ncRNA', 'Genes_lncRNA']].sum()
total['Modules']=len(df)
result_df = pd.DataFrame(total).transpose()
result_df.to_csv('totalGenesInModulesCNCwithlncRNAs.tsv',sep='\t', index=False)
