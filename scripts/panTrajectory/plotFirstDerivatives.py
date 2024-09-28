import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm

PanGeneCount = "panTrajectory/panRNAome/panRNAomeTrajectoryTable_0.8.tsv"                                                   # panRNAome
#PanGeneCount = "panTrajectory/panTranscriptome/panTranscriptomeTrajectoryTable_I2.8.tsv"                                     # pan-transcriptome
df = pd.read_csv(f"{PanGeneCount}", sep='\t', header=0)

# Filtrar e agrupar os dados
df_pan = df[df['Classification'] == 'pan'].groupby('NumberGenotypes').agg({'NumberGenes': 'mean', 'NumberGroups': 'mean'}).reset_index()

# Adicionar nova linha com total de genótipos, genes e grupos
new_row = pd.DataFrame({'NumberGenotypes': [48], 'NumberGenes': [8392174], 'NumberGroups': [3407188]})        # panRNAome
#new_row = pd.DataFrame({'NumberGenotypes': [50], 'NumberGenes': [5089777], 'NumberGroups': [600931]})          # pan-transcriptome

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

# https://www.statsmodels.org/dev/examples/notebooks/generated/lowess.html
def lowess_with_confidence_bounds(x, y, eval_x, N=200, conf_interval=0.95, lowess_kw=None):
    """
    Perform Lowess regression and determine a confidence interval by bootstrap resampling
    """
    # Lowess smoothing
    smoothed = sm.nonparametric.lowess(exog=x, endog=y, xvals=eval_x, **lowess_kw)

    # Perform bootstrap resamplings of the data
    # and  evaluate the smoothing at a fixed set of points
    smoothed_values = np.empty((N, len(eval_x)))
    for i in range(N):
        sample = np.random.choice(len(x), len(x), replace=True)
        sampled_x = x[sample]
        sampled_y = y[sample]

        smoothed_values[i] = sm.nonparametric.lowess(
            exog=sampled_x, endog=sampled_y, xvals=eval_x, **lowess_kw
        )

    # Get the confidence interval
    sorted_values = np.sort(smoothed_values, axis=0)
    bound = int(N * (1 - conf_interval) / 2)
    bottom = sorted_values[bound - 1]
    top = sorted_values[-bound]

    return smoothed, bottom, top

eval_x = np.linspace(dGroups_dGenotypes['NumberGenotypes'].min(), dGroups_dGenotypes['NumberGenotypes'].max(), 100)

# Computar a suavização e os intervalos de confiança
smoothed, bottom, top = lowess_with_confidence_bounds(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], eval_x, lowess_kw={"frac": 0.5}, conf_interval=0.95)
smoothed90, bottom90, top90 = lowess_with_confidence_bounds(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], eval_x, lowess_kw={"frac": 0.5}, N=48, conf_interval=0.95)
smoothed85, bottom85, top85 = lowess_with_confidence_bounds(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], eval_x, lowess_kw={"frac": 0.5}, N=48, conf_interval=0.95)

# Plotar os resultados
plt.figure(figsize=(14, 7))
plt.plot(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], 'o', label='Data Points')
plt.plot(eval_x, smoothed, 'r-', label='LOWESS')

# Preencher a área do intervalo de confiança
plt.fill_between(eval_x, bottom, top, color='r', alpha=0.2, label='95% IC')
plt.fill_between(eval_x, bottom90, top90, color='g', alpha=0.2, label='90% IC')
plt.fill_between(eval_x, bottom85, top85, color='b', alpha=0.2, label='85% IC')

plt.xlabel('Genótipos')
plt.ylabel('dx/dy')
plt.legend()
plt.grid(True)
plt.xticks(np.arange(dGroups_dGenotypes['NumberGenotypes'].min(), dGroups_dGenotypes['NumberGenotypes'].max() + 1, 1))
plt.tight_layout()
plt.savefig("panRNAomeGroupsTrajectory_0.8_firstDerivative_ic.png", dpi=300) 
plt.show()

'''
# Plotar a primeira derivada (dx/dy) com a área do desvio padrão da LOWESS
plt.figure(figsize=(14, 7))
plt.plot(dGroups_dGenotypes['NumberGenotypes'], dGroups_dGenotypes['dGroups'], 'o')#, label='Data Points')
plt.plot(lowess_x, lowess_y, 'r-', label='LOWESS')

# Preencher a área do desvio padrão
plt.fill_between(lowess_x, lowess_y - std_dev, lowess_y + std_dev, color='r', alpha=0.2, label='desvio padrão')

# OBS: A primeira derivada representa o coeficiente angular da tangente em cada um dos pontos na curva.
plt.title('Primeira derivada: Grupos do pan-transcriptoma')
plt.xlabel('Genótipos')
plt.ylabel('dx/dy')
plt.legend()
plt.grid(True)
plt.xticks(np.arange(dGroups_dGenotypes['NumberGenotypes'].min(), dGroups_dGenotypes['NumberGenotypes'].max() + 1, 1))
plt.tight_layout()

plt.savefig("panRNAomeGroupsTrajectory_0.8_firstDerivative_stddev.png", dpi=300) 
#plt.savefig("panTranscriptomeGroupsTrajectory_I2.8_firstDerivative.png", dpi=300)
plt.show()
'''
