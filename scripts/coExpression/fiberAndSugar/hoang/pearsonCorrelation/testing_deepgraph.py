import os
from multiprocessing import Pool
import numpy as np
import pandas as pd
import deepgraph as dg

# whiten variables for fast parallel computation later on
counterrr=0

X_df = pd.read_csv("/home/renato/tmp/Svi_matrix_tpm_subsample.txt", sep="\t")
X_df.set_index("Name", inplace=True)
X_numpy_array = X_df.to_numpy()
n_samples = len(X_df.columns.to_list())
n_features = len(X_df.index.to_list())
print(f'We have {n_samples} samples and {n_features} features')
X = (X_numpy_array - X_numpy_array.mean(axis=1, keepdims=True)) / X_numpy_array.std(axis=1, keepdims=True)
np.save('samples', X)

# parameters (change these to control RAM usage)
step_size = 1e2
n_processes = 6
# load samples as memory-map
X = np.load('samples.npy', mmap_mode='r')
# create node table that stores references to the mem-mapped samples
v = pd.DataFrame({'index': range(X.shape[0])})
# connector function to compute pairwise pearson correlations
def corr(index_s, index_t):
    features_s = X[index_s]
    features_t = X[index_t]
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    global counterrr
    counterrr += 1
    if (counterrr % 100) == 0:
        print(f'Calculated correlation for a pair of nodes in deepgraph edge. {counterrr} Corrs done.')
    return corr
# index array for parallelization
pos_array = np.array(np.linspace(0, n_features*(n_features-1)//2, n_processes), dtype=int)
# parallel computation
def create_ei(i):
    from_pos = pos_array[i]
    to_pos = pos_array[i+1]
    # initiate DeepGraph
    g = dg.DeepGraph(v)
    print('Created graph. Will calculate correlations')
    # create edges
    g.create_edges(connectors=corr, step_size=step_size,
                   from_pos=from_pos, to_pos=to_pos)
    # store edge table
    g.e.to_pickle('tmp/correlations/{}.pickle'.format(str(i).zfill(3)))
# computation
if __name__ == '__main__':
    os.makedirs("tmp/correlations", exist_ok=True)
    indices = np.arange(0, n_processes - 1)
    p = Pool()
    for _ in p.imap_unordered(create_ei, indices):
        pass
