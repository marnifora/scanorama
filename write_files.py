from bin.config import namespace, output, data_names, names
from bin.process import load_names
from scipy.io import mmwrite
from scipy.sparse import vstack
from scanorama.scanorama import merge_datasets

datasets, genes_list, cells_list, n_cells = load_names(data_names, norm=False)

datasets, genes = merge_datasets(datasets, genes_list, ds_names=names, union=False)

mmwrite(output + '%s_datasets.mtx' % namespace, vstack(datasets), field='integer')
