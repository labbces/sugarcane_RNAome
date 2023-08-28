#!/home/felipevzps/miniconda3/envs/snakemake_conekt/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def extract_csv_gen_plot(csv_path):
    data = pd.read_csv(csv_path, index_col=0, delimiter=',')
    df = pd.DataFrame(data)

    df_norm_col=(df-df.mean())/df.std()

    labels = np.array([[959, 28, 180, 9, 16, 1, 14],
                       [885, 29, 267, 9, 10, 0, 8],
                       [797, 103, 308, 11, 17, 0, 3],
                       [635, 101, 290, 18, 29, 0, 9],
                       [585, 56, 328, 45, 31, 3, 15],
                       [2127, 138, 635, 64, 67, 6, 12],
                       [746, 40, 198, 33, 35, 0, 9],
                       [528, 88, 141, 60, 52, 3, 16],
                       [1166, 130, 448, 28, 44, 5, 10],
                       [1493, 139, 595, 28, 49, 7, 14],
                       [1077, 80, 391, 5, 18, 3, 25],
                       [81, 45, 74, 31, 15, 1, 2],
                       [49, 40, 84, 45, 0, 0, 3],
                       [1052, 67, 148, 47, 4, 1, 9],
                       [823, 62, 205, 35, 8, 4, 3],
                       [1196, 153, 303, 39, 69, 0, 7],
                       [665, 85, 96, 26, 41, 1, 3],
                       [905, 99, 159, 31, 43, 0, 5],
                       [949, 91, 220, 36, 52, 0, 3],
                       [649, 95, 78, 26, 45, 0, 2],
                       [860, 86, 164, 36, 52, 0, 3],
                       [631, 70, 74, 37, 18, 0, 3],
                       [877, 104, 171, 43, 48, 0, 6],
                       [897, 91, 159, 30, 31, 0, 3],
                       [2866, 361, 1326, 57, 96, 4, 9],
                       [746, 43, 225, 30, 26, 1, 13],
                       [2440, 371, 1214, 71, 142, 0, 4],
                       [1018, 137, 286, 44, 66, 0, 3],
                       [852, 100, 143, 29, 54, 0, 7],
                       [1084, 131, 289, 50, 179, 0, 5],
                       [1348, 162, 203, 45, 68, 0, 6],
                       [906, 117, 147, 35, 52, 4, 2],
                       [812, 123, 164, 33, 121, 0, 3],
                       [984, 100, 278, 34, 44, 1, 5],
                       [1099, 143, 290, 45, 75, 0, 6],
                       [1085, 162, 193, 47, 60, 1, 2],
                       [806, 103, 119, 30, 46, 0, 6],
                       [1147, 140, 165, 24, 59, 0, 5],
                       [877, 59, 211, 16, 31, 0, 3],
                       [506, 40, 128, 11, 12, 0, 5],
                       [2569, 277, 1156, 50, 122, 1, 3],
                       [2639, 329, 1202, 72, 102, 0, 7],
                       [1121, 18, 122, 7, 34, 1, 9],
                       [703, 143, 242, 24, 18, 2, 5],
                       [3415, 394, 1452, 92, 172, 0, 15],
                       [717, 12, 120, 20, 16, 0, 5],
                       [651, 27, 185, 12, 9, 0, 4],
                       [663, 61, 177, 26, 24, 0, 4]         
                      ])

    data.index.names = ['genotypes']
    g = sns.heatmap(df_norm_col, square=False, linewidth=0, xticklabels=True, yticklabels=True,
                    annot=labels, fmt='', annot_kws={"size": 6, "weight": "light"},
                    cmap="YlGnBu")#,fmt="d", cmap="YlGnBu", annot=True)
    
    g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=5)
    g.set_xticklabels(g.get_xticklabels(), rotation=45, fontsize=5)
    g.set_title('Infernal Rfam families')
    
    plt.tight_layout()
    plt.savefig("infernal_cpclncrnas_families_48_normalized_annotation.png", dpi=300) 
    plt.show()

extract_csv_gen_plot("heatmap_infernal_48_rfam_cpcncrnas.csv")

