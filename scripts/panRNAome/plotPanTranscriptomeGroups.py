from turtle import color
from matplotlib import markers, style
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#data = pd.read_csv("teste_boxplot.tsv", delimiter="\t")
#data = pd.read_csv("pan_RNAome.tsv", delimiter="\t")
#data = pd.read_csv("pan_RNAome_hard_soft.tsv", delimiter="\t")
data = pd.read_csv("pan_RNAome_pan_hard_soft.tsv", delimiter="\t")

#print(data)

# x = genotypes
# y = pan_transcriptome and core_transcriptome
# palette="binary" ## black and white
# palette="tab10"  ## blue and orange

#sns.boxplot(y="Número de Grupos", x = "Número de Genótipos", data = data, hue="Classificação", dodge=False)#, 
            #showmeans=True, meanline=True, meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            #medianprops={'visible': False}, whiskerprops={'visible': False},
            #zorder=10,
            #showfliers = False, showbox=False, showcaps=False)

plt.figure(figsize=(20,13))

#sns.scatterplot(y="Numero de Grupos", x = "Numero de Genotipos", data = data, hue="Classificacao",
#              style="Classificacao", alpha=0.5, palette="binary")

#atualizado 05/07/2023 - pan_rnaome
sns.stripplot(y="Grupos", x = "Genotipos", data = data, hue="Classe",
              jitter=True, marker="o", alpha=0.3, palette="tab10")
plt.yscale('log')



#plt.savefig("strip_pan_grey.png", format="png")              
#plt.savefig("grupos_pan_colorido_v2.png", format="png") #svg

#sns.scatterplot(y="Transcritos", x = "Genotipos", data = data, 
#              hue="Classe", alpha=0.5, 
#              palette="binary", style="Classe")# marker="o"

##sns.stripplot(y="Transcritos", x = "Genotipos", data = data, 
##              hue="Classe", jitter=True, marker="o", alpha=0.5)#, 
#              palette="binary")
#plt.yscale('log')

#plt.savefig("strip_transcritos_grey.png", format="png")
#plt.savefig("strip_transcritos_grey.svg", format="svg")

#sns.violinplot(y="Orthogroups", x = "genotypes", data = data,inner=None, color=".8")

plt.show()


