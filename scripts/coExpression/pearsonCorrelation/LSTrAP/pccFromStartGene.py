#!/usr/bin/env python3

import argparse
import numpy as np
import sys

def pcc(filename, mcl_output, start_gene):
    """
    Reads an htseq-count matrix, calculates the PCC (Pearson Correlation) for pairs of genes starting from a specific gene.
    It will return a text file with all correlated genes and a mcl but also genemania/cytoscape compatible file
    containing all gene pairs with a correlation of 0.7 or better.

    :param filename: path to input, an htseq-count matrix file
    :param mcl_output: Mcl compatible output
    :param start_gene: gene name to start calculations from
    """
    # Read Matrix and store nominators and denominators
    with open(filename, 'r') as fin:
        nominators, denominators, genes = [], [], []

        header = fin.readline()
        size = len(header.strip().split('\t'))

        for line in fin:
            parts = line.rstrip().split("\t")

            if size != len(parts):
                print("Warning! Unequal number of columns found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                quit()

            if len(parts) == size:
                temp = []
                for j in range(1, len(parts)):
                    try:
                        temp.append(float(parts[j]))
                    except ValueError:
                        print("Warning! Non-number character found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line, file=sys.stderr)
                        quit()

                row_values = np.array(temp)
                nomi = row_values-(sum(row_values)/len(row_values))
                denomi = np.sqrt(sum(nomi**2))

                if denomi != 0.0:
                    nominators.append(nomi)
                    denominators.append(denomi)
                    genes.append(parts[0])

    nominators = np.array(nominators)
    denominators = np.array(denominators)

    # Find index of specific gene
    try:
        start_index = genes.index(start_gene)
    except ValueError:
        print("Error: The specified gene '{}' is not found in the matrix.".format(start_gene), file=sys.stderr)
        quit()

    # Calculate PCC and write mcl_output
    with open(mcl_output, 'w') as mcl_out:
        print("Database OK.\nCalculating Pearson Correlation Coefficient and ranks.\n")
        for i, (nom, denom, gene) in enumerate(zip(nominators[start_index:], denominators[start_index:], genes[start_index:]), start=1):
            print("Calculated PCC values for sequence:%s, %d out of %d." % (gene, i, len(nominators)))

            nominator = np.dot(nominators, nom)
            denominator = np.dot(denominators, denom)
            pcc_values = nominator/denominator

            data = [{'score': p,
                     'gene': g,
                     'string': g + '(' + str(p) + ')'} for g, p in zip(genes, pcc_values) if g != gene]

            # sort by absolute pcc value
            data.sort(key=lambda x: x['score'], reverse=True)

            for s in data:
                # Keep scores > 0.7 subtract 0.7 from result to remap values to [0,0.3] as this is important for mcl
                if s['score'] > 0.7:
                    print(gene, s['gene'], s['score'] - 0.7, sep='\t', file=mcl_out)

    print("PCCs calculated and saved as %s." % (mcl_output), file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="./pcc.py")

    parser.add_argument('input', help='path to input')
    parser.add_argument('mcl_output', help='path to mcl compatible output')
    parser.add_argument('start_gene', help='gene to start calculations from')

    args = parser.parse_args()

    pcc(args.input, args.mcl_output, args.start_gene)

