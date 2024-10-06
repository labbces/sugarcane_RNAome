import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

file_path = 'GOFrequency/total/total_lncRNAs_GOFrequencies.tsv'
df = pd.read_csv(file_path, delimiter='\t')
df_top = df.nlargest(70, 'Frequency')
df_sorted = df_top.sort_values('Frequency', ascending=True)

plt.figure(figsize=(8, 12))
plt.plot(df_sorted['Frequency'], range(len(df_sorted)), marker='', linestyle='-', color='blue', linewidth=2)

# configuração dos eixos e labels
plt.yticks(range(len(df_sorted)), df_sorted['Name'], fontsize=8)
plt.xticks(fontsize=7)
plt.gca().yaxis.tick_right()

plt.xlabel('Frequência')
#plt.title('Top 70 Termos GO nos módulos enriquecidos')
#plt.title('Top 70 Termos GO nos lncRNAs anotados')
plt.title('Top 70 Termos GO no total dos lncRNAs anotados')
#plt.title('Top 70 Termos GO na interseção dos lncRNAs anotados')

plt.grid(False) 

# mostrar apenas inteiros no eixo X
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig('GOFrequency/total/top70_total_lncRNAs_GOFrequencies.pdf', bbox_inches='tight')
plt.savefig('GOFrequency/total/top70_total_lncRNAs_GOFrequencies.png', bbox_inches='tight')
plt.show()