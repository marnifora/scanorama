from bin.process import load_names
from scanorama.scanorama import *
from scipy.io import mmwrite
from time import time
from bin.config import data_names, names, namespace, output, write, tsne, dimred

t0 = time()

print('Loading matrices from files')
datasets, genes_list, cells_list, n_cells, counts = load_names(data_names)

print('Selecting common genes across all datasets')
datasets, genes = merge_datasets(datasets, genes_list, ds_names=names)

if write and counts:
    t1 = time()
    mmwrite('{}{}_matrix_counts.mtx'.format(output, namespace), vstack(datasets))
    timew = time() - t1
elif write and not counts:
    timew = 0.0
    print('Matrix of counts has not been written as not all given matrices are matrix of counts!')

print('Normalization and dimension reduction')
datasets_dimred, datasets_norm = process_data(datasets, genes, dimred=dimred)

if write:
    t1 = time()
    mmwrite('{}{}_matrix_norm.mtx'.format(output, namespace), vstack(datasets_norm))
    mmwrite('{}{}_matrix_dimred.mtx'.format(output, namespace), vstack(datasets_dimred))
    timew += time() - t1

print('Scanorama adjusting process')
datasets_adjusted = assemble(datasets_dimred, ds_names=names)

if write or tsne:
    t1 = time()
    cells = []
    for c, name in zip(cells_list, names):
        for cell in c:
            cells.append('%s:%s' % (cell, name))
    timew += time() - t1

if write:
    t1 = time()
    mmwrite('{}{}_matrix_adjusted.mtx'.format(output, namespace), vstack(datasets_adjusted))
    with open(output + '{}_genes_list.txt'.format(namespace), 'w') as o:
        o.write('\n'.join(genes))
    with open(output + '{}_cells_list.txt'.format(namespace), 'w') as o:
        o.write('\n'.join(cells))
    timew += time() - t1

if write:
    print('Integrated and batch corrected panoramas in {:.3f} minutes (including {:.3f}'.format((time()-t0)/60, timew/60)
          + ' min for writing matrices into files)')
else:
    print('Integrated and batch corrected panoramas in {:.3f} minutes'.format((time() - t0)/60))

if tsne:
    print('Caculating t-SNE')
    calculate_tsne(vstack(datasets_adjusted), cells, namespace, output, 10)
