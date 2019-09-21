from bin.config import namespace, output, data_names
from bin.process import load_names
from scipy.io import mmwrite
from scipy.sparse import vstack

datasets, _, _, _ = load_names(data_names, norm=False)

mmwrite(output + '%s_datasets.mtx' % namespace, vstack(datasets))
