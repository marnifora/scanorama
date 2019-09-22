from bin.config import namespace, output, data_names, names
from bin.process import load_names
from scipy.io import mmwrite
from scipy.sparse import vstack
from scanorama.scanorama import merge_datasets, process_data

datasets, genes_list, cells_list, n_cells = load_names(data_names, norm=False)

datasets, genes = merge_datasets(datasets, genes_list, ds_names=names, union=False)

datasets_dimred, datasets_norm, genes = process_data(datasets, genes, dimred=0)

mmwrite(output + '%s_datasets_norm-new.mtx' % namespace, vstack(datasets_norm))

mmwrite(output + '%s_datasets-new.mtx' % namespace, vstack(datasets), field='integer')
