#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(prog = "calculateDegree", description="calculate gene degree from network", add_help=True)
parser.add_argument("-i", dest="network_file", metavar="<network_mcl_out.txt>", help="PCC output file", required=True)

args = parser.parse_args()
network_file = args.network_file
output_file = network_file[:-4] + "_degree.tsv"

# calculate degree from each gene
def calculate_degree(network, output):
    degree_count = {}

    # read file line by line
    with open(network, 'r') as file:
        for line in file:
            source, target, weight = line.strip().split()

            if source in degree_count:
                degree_count[source] += 1
            else:
                degree_count[source] = 1

            if target in degree_count:
                degree_count[target] += 1
            else:
                degree_count[target] = 1

    # write degree from each gene
    with open(output, 'w') as out_file:
        for gene, degree in degree_count.items(): # sum degree_count[source] and degree_count[target]
            out_file.write(f"{gene}\t{degree}\n")
        #for source in degree_count:
            #out_file.write(f"{source}\t{degree_count[source]}\n")

calculate_degree(network_file, output_file)
