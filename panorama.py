from bin.process import load_names
from time import time
import numpy as np
import sys
from scipy.io import mmwrite
from scipy.sparse import vstack
from scanorama.scanorama import *

if __name__ == '__main__':
    from bin.config import data_names, names, namespace, path, output, metadata

    datasets, genes_list, cells_list, n_cells = load_names(data_names)

    t0 = time()
    datasets_dimred, datasets_norm, datasets, genes = correct(
        datasets, genes_list, ds_names=names,
        sigma=150, return_dimred=True
    )
    if VERBOSE:
        print('Integrated and batch corrected panoramas in {:.3f}s'
              .format(time() - t0))

    mmwrite(output + '%s_datasets_norm.mtx' % namespace, vstack(datasets_norm))
    mmwrite(output + '%s_datasets_dimred.mtx' % namespace, vstack(datasets_dimred))

    labels = []
    curr_label = 0
    for i, a in enumerate(datasets):
        labels += list(np.zeros(a.shape[0]) + curr_label)
        curr_label += 1
    labels = np.array(labels, dtype=int)
    
    embedding = visualize(datasets_dimred,
                          labels, path + namespace, names,
                          multicore_tsne=False, viz_cluster=True)

    metadata_into_file(embedding, labels, names, output, cells_list, namespace, metadata)

    # Uncorrected.
    datasets, genes_list, cells_list, n_cells = load_names(data_names)
    datasets, genes = merge_datasets(datasets, genes_list)
    datasets_dimred = dimensionality_reduce(datasets)
    
    labels = []
    names = []
    curr_label = 0
    for i, a in enumerate(datasets):
        labels += list(np.zeros(a.shape[0]) + curr_label)
        names.append(data_names[i])
        curr_label += 1
    labels = np.array(labels, dtype=int)

    embedding = visualize(datasets_dimred, labels,
                          path + namespace + '_uncorrected', names)