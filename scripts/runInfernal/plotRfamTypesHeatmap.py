import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def extract_csv_gen_plot(directory_path):
    csv_files = [f for f in os.listdir(directory_path) if f.endswith('.tblout')]

    data_dict = {}

    for csv_file in csv_files:
        genotype = csv_file.split('.')[0]  # extract the first part of the filename as the genotype
        file_path = os.path.join(directory_path, csv_file)
        df = pd.read_csv(file_path, index_col=0, delimiter=',')
        data_dict[genotype] = df['Count']

    #print(data_dict)
    
    # create combined df with all genotypes and RNA types (fill NA with zeros)
    combined_df = pd.DataFrame(data_dict).fillna(0).astype(int)

    #print(combined_df)
    
    # normalize each column to scale the colors properly 
    df_norm_col = (combined_df.T-combined_df.T.mean())/combined_df.T.std()

    print(df_norm_col)

    labels = combined_df.T.values # labels for annotation
    
    df_norm_col.index.names = ['Genotypes']
    df_norm_col.columns.names = ['RNA family types']

    plt.figure(figsize=(16, 12))
    g = sns.heatmap(df_norm_col, square=False, linewidth=0, xticklabels=True, yticklabels=True,
                    annot=labels, fmt='', annot_kws={"size": 6, "weight": "light"},
                    cmap="YlGnBu", cbar_kws={'label': 'Normalized Count'})

    g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=8)
    g.set_xticklabels(g.get_xticklabels(), rotation=70, fontsize=8)
    g.set_title('Infernal Rfam family types')

    plt.tight_layout()
    plt.show()
    plt.savefig('panRNAomeRfamFamilies.png', dpi=300)

# Usage
directory_path = '/home/felipe/Documents/sugarcane_RNAome/scripts/runInfernal/individualGenotypes' 
extract_csv_gen_plot(directory_path)