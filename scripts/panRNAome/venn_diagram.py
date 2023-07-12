#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Ler os arquivos
with open('CPC2_ncRNAs_list.txt', 'r') as f:
	set1 = set(f.read().strip().split(','))

with open('CPC2_PLncPRO_ncRNAs_list.txt', 'r') as f:
	set2 = set(f.read().strip().split(','))

with open('RNAplonc_ncRNAs_list.txt', 'r') as f:
	set3 = set(f.read().strip().split(','))

# Gerar o diagrama de Venn
venn_diagram = venn3([set1, set2, set3], ('CPC2', 'PLncPRO', 'RNAplonc'))

# Configurar o texto das interseções
venn_diagram.get_label_by_id('100').set_text(f"CPC2: {len(set1 - set2 - set3)}")
if venn_diagram.get_label_by_id('010') is not None:
	venn_diagram.get_label_by_id('010').set_text(f"PLncPRO: {len(set2 - set1 - set3)}")
if venn_diagram.get_label_by_id('001') is not None:
	venn_diagram.get_label_by_id('001').set_text(f"RNAplonc: {len(set3 - set1 - set2)}")
if venn_diagram.get_label_by_id('110') is not None:
	venn_diagram.get_label_by_id('110').set_text(f"CPC2 e PLncPRO: {len(set1 & set2 - set3)}")
if venn_diagram.get_label_by_id('101') is not None:
	venn_diagram.get_label_by_id('101').set_text(f"CPC2 e RNAplonc: {len(set1 & set3 - set2)}")
if venn_diagram.get_label_by_id('011') is not None:
	venn_diagram.get_label_by_id('011').set_text(f"PLncPRO e RNAplonc: {len(set2 & set3 - set1)}")
if venn_diagram.get_label_by_id('111') is not None:
	venn_diagram.get_label_by_id('111').set_text(f"CPC2, PLncPRO e RNAplonc: {len(set1 & set2 & set3)}")

# Exibir o diagrama de Venn
plt.savefig("diagrama_venn2.png")

