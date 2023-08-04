#!/usr/bin/env python

import argparse

# Uso
#./findConsensus.py lista1.txt lista2.txt lista3.txt saida.txt

def ler_identificadores_de_arquivo(nome_arquivo):
    with open(nome_arquivo, 'r') as arquivo:
        # Le o conteudo do arquivo e retorna uma lista com os identificadores
        identificadores = arquivo.read().strip().split(',')
    return identificadores

def encontrar_consenso(lista1, lista2, lista3):
    # Cria um conjunto com os identificadores de cada lista para encontrar o consenso de forma eficiente
    set1 = set(lista1)
    set2 = set(lista2)
    set3 = set(lista3)

    # Encontra a intersecao dos tres conjuntos
    consenso_set = set1.intersection(set2, set3)

    # Cria a lista consenso mantendo a ordem original dos identificadores
    consenso = [identificador for identificador in lista1 if identificador in consenso_set]

    return consenso

def salvar_identificadores_em_arquivo(identificadores, nome_arquivo):
    with open(nome_arquivo, 'w') as arquivo:
        for identificador in identificadores:
            arquivo.write(identificador + '\n')

def main():
    parser = argparse.ArgumentParser(description='Encontra o consenso entre as tres listas de ncRNAs identificados - output do CPC2, PLncPRO e RNAplonc')
    parser.add_argument('lista1', help='Arquivo contendo os identificadores da lista 1 (formato de lista separada por virgula)')
    parser.add_argument('lista2', help='Arquivo contendo os identificadores da lista 2 (formato de lista separada por virgula)')
    parser.add_argument('lista3', help='Arquivo contendo os identificadores da lista 3 (formato de lista separada por virgula)')
    parser.add_argument('saida', help='Nome do arquivo de saida para salvar os identificadores em consenso')

    args = parser.parse_args()

    lista1 = ler_identificadores_de_arquivo(args.lista1)
    lista2 = ler_identificadores_de_arquivo(args.lista2)
    lista3 = ler_identificadores_de_arquivo(args.lista3)

    identificadores_consenso = encontrar_consenso(lista1, lista2, lista3)

    salvar_identificadores_em_arquivo(identificadores_consenso, args.saida)

if __name__ == '__main__':
    main()
