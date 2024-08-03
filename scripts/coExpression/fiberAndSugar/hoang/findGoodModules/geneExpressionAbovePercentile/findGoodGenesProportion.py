#!/usr/bin/env python

import pandas as pd

file_path = 'Hoang2017_counts_filters_VST_CNC_CV_above1.2.txt'
data = pd.read_csv(file_path, sep='\t')

# Transpor o dataframe para calcular o percentil com base em todos os valores de expressão
all_values = data.iloc[:, 1:].values.flatten()

# Função para calcular a proporção de valores baixos
def calculate_low_value_proportion(row, threshold):
    return (row < threshold).mean()

# Função principal para calcular e salvar proporções de valores baixos
def process_data(percentile):
    # Calcular o percentil para definir o valor "baixo"
    low_threshold = pd.Series(all_values).quantile(percentile / 100)
    print(f"Valor de corte ({percentile}º percentil): {low_threshold}")

    # Nome da coluna com o valor do corte
    column_name = f'low_value_proportion (< {low_threshold:.6f})'

    # Calcular a proporção de valores baixos para cada gene
    data[column_name] = data.iloc[:, 1:].apply(lambda row: calculate_low_value_proportion(row, low_threshold), axis=1)

    # Ordenar os genes pela proporção de valores baixos (do menor para o maior)
    sorted_genes = data.sort_values(by=column_name, ascending=True)

    # Salvar o nome do gene e a proporção de valores baixos em um arquivo CSV
    output_file = f'proportionOfExpressedGenes_{percentile}percentile.csv'
    sorted_genes[['Name', column_name]].to_csv(output_file, index=False)

    print(f"Lista de genes ordenados com proporção de valores baixos salva em '{output_file}'")

# Exemplo de uso com diferentes percentis
percentiles_to_explore = [25, 50, 75, 90]  
for percentile in percentiles_to_explore:
    process_data(percentile)

