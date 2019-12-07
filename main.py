from bin.process import load_names
from scanorama.scanorama import *
from scipy.io import mmwrite
from time import time
from sklearn.preprocessing import normalize
import os
from bin.config import data_names, names, namespace, output, write, tsne, pc, metadata, scanpy, adjust, norm, dimred, datasets, path

t0 = time()

if not datasets and data_names:
    print('Loading matrices from files')
    datasets, genes_list, cells_list, n_cells, counts = load_names(data_names)

    print('Selecting common genes across all datasets')
    datasets, genes = merge_datasets(datasets, genes_list, ds_names=names)

    if write or tsne:
        t1 = time()
        cells = []
        for c, name in zip(cells_list, names):
            for cell in c:
                cells.append('%s:%s' % (cell, name))
        timew = time() - t1
        print("Cells' names has been established in {} s".format(timew))

    if write:
        t1 = time()
        with open(output + '{}_genes.txt'.format(namespace), 'w') as o:
            o.write('\n'.join(genes))
        with open(output + '{}_cells.txt'.format(namespace), 'w') as o:
            o.write('\n'.join(cells))
        if metadata is None:
            with open(output + '{}_metadata.tsv'.format(namespace), 'w') as o:
                o.write('Cell ID\tDataset\n')
                for cell in cells:
                    o.write('{}\t{}\n'.format(cell, cell.split(':')[-1]))
        print('Cells, genes and metadata have been written into files!')
        if all(counts):
            mmwrite('{}{}_matrix_counts.mtx'.format(output, namespace), vstack(datasets))
            timew = time() - t1
            print('Matrix of counts has been written in {:.3f} min'.format(timew/60))
        elif write and not all(counts):
            print('Matrix of counts has not been written as not all given matrices are matrix of counts!')

if len(datasets) == 1:
    print('WARNING: only one dataset given - Scanorama will not do anything!')

print('Normalization and dimension reduction')
if norm:
    datasets_norm = []
    for ds, c, n in zip(datasets, counts, names):
        if c:
            print('Normalizing {}'.format(n))
            datasets_norm.append(np.log1p(normalize(ds.copy(), axis=1)))
        else:
            print('{} already normalized'.format(n))
            datasets_norm.append(ds)
if dimred and pc > 0:
    print('Reducing dimension...')
    if 'datasets_norm' not in globals():
        datasets_norm = datasets
    datasets_dimred = dimensionality_reduce([ds.copy() for ds in datasets_norm], dimred=pc, scanpy=scanpy)
# datasets_dimred, datasets_norm = process_data(datasets, genes, dimred=pc)

if (norm or dimred) and write:
    t1 = time()
    mmwrite('{}{}_matrix_norm.mtx'.format(output, namespace), vstack(datasets_norm))
    if scanpy:
        mmwrite('{}{}_matrix_dimred-scanpy.mtx'.format(output, namespace), vstack(datasets_dimred))
    else:
        mmwrite('{}{}_matrix_dimred.mtx'.format(output, namespace), vstack(datasets_dimred))
    timew = time() - t1
    print('Norm and dimred matrices have been written in {:.3f} min'.format(timew/60))

if adjust:
    print('Scanorama adjusting process')
    if 'datasets_dimred' not in globals():
        datasets_dimred = [el.toarray() for el in datasets]
    datasets_adjusted = assemble(datasets_dimred, ds_names=names)

    if write:
        t1 = time()
        if scanpy:
            mmwrite('{}{}_matrix_adjusted-scanpy.mtx'.format(output, namespace), vstack(datasets_adjusted))
        else:
            mmwrite('{}{}_matrix_adjusted.mtx'.format(output, namespace), vstack(datasets_adjusted))
        timew = time() - t1
        print('Adjusted matrix has been written in {} min'.format(timew/60))

    if write:
        print('Integrated and batch corrected panoramas in {:.3f} minutes'.format((time()-t0)/60))
    else:
        print('Integrated and batch corrected panoramas in {:.3f} minutes'.format((time() - t0)/60))

if tsne:
    print('Caculating t-SNE')
    if scanpy:
        calculate_tsne(vstack(datasets_adjusted), cells, namespace + '_tsne_adjusted-scanpy.tsv', output, 10)
    else:
        calculate_tsne(vstack(datasets_adjusted), cells, namespace + '_tsne_adjusted.tsv', output, 10)
