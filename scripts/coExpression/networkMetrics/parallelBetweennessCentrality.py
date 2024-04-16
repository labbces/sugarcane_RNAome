#!/usr/bin/env python3

from multiprocessing import Pool
import time
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import csv
import argparse

parser = argparse.ArgumentParser(description="parallel betweenness centrality with networkX.")
parser.add_argument("--edge_list", dest="edge_list", metavar="Correr_clusteringSimilarities.tsv", required=True, help="Path to the edge list")
parser.add_argument("--output", dest="output_file", metavar="betweenness.csv", required=True, help="Path to the output file")
parser.add_argument("--threads", dest="threads", metavar="2", type=int, required=True, help="threads")
parser.add_argument("--columns", dest="columns", nargs='+', type=int, metavar="0 1 2", default=[0, 1, 2], required=True, help="columns indices (source, target, weight)")
args = parser.parse_args()

def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x

def betweenness_centrality_parallel(G, output_file, processes=None):
    """Calculate and write parallel betweenness centrality"""
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['node', 'betweenness_centrality']
        #print(fieldnames)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        p = Pool(processes=processes)
        node_divisor = len(p._pool) * 4
        node_chunks = list(chunks(G.nodes(), G.order() // node_divisor))
        num_chunks = len(node_chunks)
        bt_sc = p.starmap(
            nx.betweenness_centrality_subset,
            zip(
                [G] * num_chunks,
                node_chunks,
                [list(G)] * num_chunks,
                [True] * num_chunks,
                [None] * num_chunks,
            ),
        )
        
        # Reduce the partial solutions
        bt_c = bt_sc[0]
        for bt in bt_sc[1:]:
            for n in bt:
                bt_c[n] += bt[n]
                
        for node in bt_c:
            writer.writerow({'node': node, 'betweenness_centrality': f"{bt_c[node]:.15f}"})
            
        return bt_c

df = pd.read_csv(args.edge_list, sep='\t', usecols=args.columns)
df.columns = ['source', 'target', 'weight']
G_my = nx.from_pandas_edgelist(df)

for G in [G_my]:

    print("")
    print("Computing betweenness centrality for:")
    print(G)

    # parallel version
    start = time.time()
    bt = betweenness_centrality_parallel(G, args.output_file ,processes=args.threads)
    print("\tParallel version")
    print(f"\t\tTime: {(time.time() - start):.4F} seconds")
    #print(f"\t\tBetweenness centrality for node out.5.8_1: {bt['out.5.8_1']:.5f}")
print("")