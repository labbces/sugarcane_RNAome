import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

#PanGeneCount = "panRNAome/panRNAomeTrajectoryTable_0.8.tsv"                                                   # panRNAome
PanGeneCount = "panTranscriptome/panTranscriptomeTrajectoryTable_I2.8.tsv"                                     # pan-transcriptome
df = pd.read_csv(f"{PanGeneCount}", sep='\t', header=0)

# Filtrar e agrupar os dados
df_pan = df[df['Classification'] == 'pan'].groupby('NumberGenotypes').agg({'NumberGenes': 'mean', 'NumberGroups': 'mean'}).reset_index()

# Adicionar nova linha com total de genótipos, genes e grupos
#new_row = pd.DataFrame({'NumberGenotypes': [48], 'NumberGenes': [8392174], 'NumberGroups': [3407188]})        # panRNAome
new_row = pd.DataFrame({'NumberGenotypes': [50], 'NumberGenes': [5089777], 'NumberGroups': [600931]})          # pan-transcriptome

df_pan = pd.concat([df_pan, new_row], ignore_index=True)

# Ordenar os dados pelo número de genótipos (garantir a ordem correta)
df_pan = df_pan.sort_values(by='NumberGenotypes').reset_index(drop=True)

# Calcular dx e dy
dGenotypes = np.diff(df_pan['NumberGenotypes'])
print(f"dGenotypes:\n{dGenotypes}")

dGroups = np.diff(df_pan['NumberGroups'])
print(f"dGroups:\n{dGroups}")

# Calcular frações (dx/dy)
dGroups_dGenotypes = pd.DataFrame({'dGroups': dGroups / dGenotypes, 'NumberGenotypes': df_pan['NumberGenotypes'][:-1]})
print(f"dGroups_dGenotypes:\n{dGroups_dGenotypes}")

# Ajustar o modelo Locally Weighted Scatterplot Smoothing (LOWESS): https://www.statsmodels.org/dev/examples/notebooks/generated/lowess.html
# usage: lowess is for adding a smooth curve to a scatterplot
lowess_result = lowess(dGroups_dGenotypes['dGroups'], dGroups_dGenotypes['NumberGenotypes'], frac=0.5)
#print(lowess_result)

lowess_x, lowess_y = lowess_result[:, 0], lowess_result[:, 1]

# Calcular resíduos e desvio padrão
residuals = dGroups_dGenotypes['dGroups'] - np.interp(dGroups_dGenotypes['NumberGenotypes'], lowess_x, lowess_y)
std_dev = np.std(residuals)

# Plotar a primeira derivada (dx/dy) com a área do desvio padrão da LOWESS
plt.figure(figsize=(14, 7))
plt.plot(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], 'o')#, label='Data Points')
plt.plot(lowess_x, lowess_y, 'r-', label='LOWESS')

# Preencher a área do desvio padrão
plt.fill_between(lowess_x, lowess_y - std_dev, lowess_y + std_dev, color='r', alpha=0.2)#, label='LOWESS')

# OBS: A primeira derivada representa o coeficiente angular da tangente em cada um dos pontos na curva.
plt.title('Primeira derivada: Grupos do pan-transcriptoma')
plt.xlabel('Genótipos')
plt.ylabel('dx/dy')
plt.legend()
plt.grid(True)
plt.xticks(np.arange(dGroups_dGenotypes['NumberGenotypes'].min(), dGroups_dGenotypes['NumberGenotypes'].max() + 1, 1))
plt.tight_layout()

#plt.savefig("panRNAomeGroupsTrajectory_0.8_firstDerivative.png", dpi=300) 
plt.savefig("panTranscriptomeGroupsTrajectory_I2.8_firstDerivative.png", dpi=300)
plt.show()