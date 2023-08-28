import sys
import os

if len(sys.argv) <= 1:
    print('argumentos insuficientes')
    print('insira a saida do Rfam formatada')
    sys.exit()

def contar_classes_rna(arquivo, prefixos_rna):
    count = {prefixo: 0 for prefixo in prefixos_rna}

    # lê o arquivo e armazena as contagens de classes de RNA
    with open(arquivo, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            rnas = fields[1].split(',')
            for rna in rnas:
                rna = rna.lower()
                for sufixo in prefixos_rna:
                    if rna.startswith(sufixo.lower()):
                        if sufixo not in count:
                            count[sufixo] = 0
                        count[sufixo] += 1

    return count

arquivo = sys.argv[1]

classes_rna = ['sno','trna','mir','intron','u', 'tpp', 'enod40']

count = contar_classes_rna(arquivo, classes_rna)

# imprime as contagens de classes de RNA
print([i for i in count.values()])
