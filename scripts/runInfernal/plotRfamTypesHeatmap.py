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

    #print(df_norm_col)

    # Reorder df columns -  using transposed combined_df to keep numbers in order (lncRNAs on the left)
    combined_df = combined_df.T[[' CD-box', ' leader', ' miRNA', ' riboswitch', ' tRNA', ' thermoregulator', '6S', 'Arthropod_7SK', 
                               'Bacteria_small_SRP', 'EBv-sisRNA-2', 'Fungi_SRP', 'Intron', 'Metazoa_SRP', 'OLE', 'Plant_SRP', 
                               'Protozoa_SRP', 'VA', 'enod40', 'tmRNA', ' lncRNA', ' rRNA', ' splicing', ' sRNA', ' ribozyme', 
                               'Cis-reg', ' HACA-box', ' antisense']]
    
    labels = combined_df.values # labels for annotation
    print(combined_df)

    # Reorder df columns - now I am just keeping the new column order - from normalized df (lncRNAs on the left)
    df_norm_col = df_norm_col[[' CD-box', ' leader', ' miRNA', ' riboswitch', ' tRNA', ' thermoregulator', '6S', 'Arthropod_7SK', 
                               'Bacteria_small_SRP', 'EBv-sisRNA-2', 'Fungi_SRP', 'Intron', 'Metazoa_SRP', 'OLE', 'Plant_SRP', 
                               'Protozoa_SRP', 'VA', 'enod40', 'tmRNA', ' lncRNA', ' rRNA', ' splicing', ' sRNA', ' ribozyme', 
                               'Cis-reg', ' HACA-box', ' antisense']]
    
    df_norm_col.index.names = ['Genótipos']
    df_norm_col.columns.names = ['Famílias de RNA']
    

    plt.figure(figsize=(16, 12))
    g = sns.heatmap(df_norm_col, square=False, linewidth=0, xticklabels=True, yticklabels=True,
                    annot=labels, fmt='', annot_kws={"size": 8, "weight": "light"},
                    cmap="YlGnBu", cbar_kws={'label': 'Contagens normalizadas'})

    g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=8)
    g.set_xticklabels(g.get_xticklabels(), rotation=70, fontsize=8)
    g.set_title('Famílias de RNA do Rfam')

    plt.tight_layout()
    plt.show()
    #plt.savefig('panRNAomeRfamFamilies.png', dpi=300)#, bbox_inches='tight')

# Usage
directory_path = 'C:/Users/PC/Desktop/temp/Rfam/individualGenotypes' 
#directory_path = '/home/felipe/Documents/sugarcane_RNAome/scripts/runInfernal/individualGenotypes' 

extract_csv_gen_plot(directory_path)
