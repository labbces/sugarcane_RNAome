#!/usr/bin/env python3

# Caminhos dos arquivos
caminho_arquivo1 = 'deepgraph_500_network.txt'
caminho_arquivo2 = 'lstrap_500_network.txt'


def extrair_correlacoes(linha):
    partes = linha.split(':')
    gene = partes[0].strip()
    correlacoes = [(item.split('(')[0].strip(), float(item.split('(')[1].replace(')', '').strip())) for item in partes[1].split()]
    return gene, correlacoes

def comparar_correlacoes(arquivo1, arquivo2):
    correlacoes_arquivo1 = {gene: correlacoes for gene, correlacoes in [extrair_correlacoes(linha) for linha in arquivo1]}
    correlacoes_arquivo2 = {gene: correlacoes for gene, correlacoes in [extrair_correlacoes(linha) for linha in arquivo2]}

    # Verificar se os genes do arquivo 1 existem no arquivo 2 e se os valores de correlação são iguais
    for gene, correlacoes1 in correlacoes_arquivo1.items():
        correlacoes2 = correlacoes_arquivo2.get(gene, [])

        if not correlacoes2:
            print(f"Gene {gene}: Não encontrado no arquivo 2.")
            continue

        for (gene1, valor1) in correlacoes1:
            # Procurar a correspondência no arquivo 2
            correspondencia = next((valor2 for gene2, valor2 in correlacoes2 if gene1 == gene2), None)

            if correspondencia is None:
                print(f"Gene {gene}: Valor de correlacao para {gene1} nao encontrado no LSTrAP.")
            elif valor1 != correspondencia:
                print(f"Gene {gene}: Valor de correlacao diferente para {gene1} DeepGraph {valor1}, LSTrAP {correspondencia}")

# Ler os dados dos arquivos
with open(caminho_arquivo1, 'r') as arquivo1, open(caminho_arquivo2, 'r') as arquivo2:
    dados_arquivo1 = arquivo1.readlines()
    dados_arquivo2 = arquivo2.readlines()

# Comparar as correlações
comparar_correlacoes(dados_arquivo1, dados_arquivo2)

