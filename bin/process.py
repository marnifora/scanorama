import gzip
import numpy as np
import os.path
import scipy.sparse
from scipy.io import mmread
from scipy.sparse import csr_matrix
from sklearn.preprocessing import normalize
import sys

MIN_TRANSCRIPTS = 600

def load_tab(fname, max_genes=40000):
    if fname.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
        
    with opener(fname, 'r') as f:
        if fname.endswith('.gz'):
            header = f.readline().decode('utf-8').rstrip().split('\t')
        else:
            header = f.readline().rstrip().split('\t')
            
        cells = [c.strip() for c in header[1:] if c]
        X = np.zeros((len(cells), max_genes))
        genes = []
        for i, line in enumerate(f):
            if i > max_genes:
                break
            if fname.endswith('.gz'):
                line = line.decode('utf-8')
            fields = line.strip().split('\t')
            genes.append(fields[0])
            X[:, i] = [ float(f) for f in fields[1:] ]
    return X[:, range(len(genes))], np.array(cells), np.array(genes)

def load_mtx(dname):
    with open(dname + '/matrix.mtx', 'r') as f:
        while True:
            header = f.readline()
            if not header.startswith('%'):
                break
        header = header.rstrip().split()
        n_genes, n_cells = int(header[0]), int(header[1])

        data, i, j = [], [], []
        for line in f:
            fields = line.rstrip().split()
            data.append(float(fields[2]))
            i.append(int(fields[1])-1)
            j.append(int(fields[0])-1)
        X = csr_matrix((data, (i, j)), shape=(n_cells, n_genes))

    genes = []
    with open(dname + '/genes.tsv', 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            genes.append(fields[1])
    assert(len(genes) == n_genes)

    with open(dname + '/barcodes.tsv', 'r') as f:
        cells = f.read().split('\n')

    return X, np.array(cells), np.array(genes)

def load_h5(fname, genome='mm10'):
    try:
        import tables
    except ImportError:
        sys.stderr.write('Please install PyTables to read .h5 files: '
                         'https://www.pytables.org/usersguide/installation.html\n')
        exit(1)
    
    # Adapted from scanpy's read_10x_h5() method.
    with tables.open_file(str(fname), 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()

            n_genes, n_cells = dsets['shape']
            data = dsets['data']
            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']

            X = csr_matrix((data, dsets['indices'], dsets['indptr']),
                           shape=(n_cells, n_genes))
            genes = [ gene for gene in dsets['genes'].astype(str) ]
            assert(len(genes) == n_genes)
            assert(len(genes) == X.shape[1])
            
        except tables.NoSuchNodeError:
            raise Exception('Genome %s does not exist in this file.' % genome)
        except KeyError:
            raise Exception('File is missing one or more required datasets.')    

    return X, np.array(genes)


def load_h5ad(fname):
    import anndata
    adata = anndata.read_h5ad(fname)
    X = adata.X
    cells = np.array(adata.obs.index)
    genes = np.array([el.upper() for el in adata.var.index])
    return X, genes, cells


def process_tab(fname, min_trans=MIN_TRANSCRIPTS):
    X, cells, genes = load_tab(fname)

    gt_idx = [ i for i, s in enumerate(np.sum(X != 0, axis=1))
               if s >= min_trans ]
    X = X[gt_idx, :]
    cells = cells[gt_idx]
    if len(gt_idx) == 0:
        print('Warning: 0 cells passed QC in {}'.format(fname))
    if fname.endswith('.txt'):
        cache_prefix = '.'.join(fname.split('.')[:-1])
    elif fname.endswith('.txt.gz'):
        cache_prefix = '.'.join(fname.split('.')[:-2])
    elif fname.endswith('.tsv'):
        cache_prefix = '.'.join(fname.split('.')[:-1])
    elif fname.endswith('.tsv.gz'):
        cache_prefix = '.'.join(fname.split('.')[:-2])
    else:
        sys.stderr.write('Tab files should end with ".txt" or ".tsv"\n')
        exit(1)
        
    cache_fname = cache_prefix + '.npz'
    np.savez(cache_fname, X=X, genes=genes, cells=cells)

    return X, cells, genes


def process_mtx(dname, min_trans=MIN_TRANSCRIPTS):
    X, cells, genes = load_mtx(dname)

    gt_idx = [ i for i, s in enumerate(np.sum(X != 0, axis=1))
               if s >= min_trans ]
    X = X[gt_idx, :]
    cells = cells[gt_idx]
    if len(gt_idx) == 0:
        print('Warning: 0 cells passed QC in {}'.format(dname))
    
    cache_fname = dname + '/tab.npz'
    scipy.sparse.save_npz(cache_fname, X, compressed=False)

    with open(dname + '/tab.genes.txt', 'w') as of:
        of.write('\n'.join(genes) + '\n')

    with open(dname + '/tab.cells.txt', 'w') as of:
        of.write('\n'.join(cells) + '\n')

    return X, cells, genes


def process_h5(fname, min_trans=MIN_TRANSCRIPTS):
    X, genes = load_h5(fname)

    gt_idx = [ i for i, s in enumerate(np.sum(X != 0, axis=1))
               if s >= min_trans ]
    X = X[gt_idx, :]
    if len(gt_idx) == 0:
        print('Warning: 0 cells passed QC in {}'.format(fname))
    
    if fname.endswith('.h5'):
        cache_prefix = '.'.join(fname.split('.')[:-1])
        
    cache_fname = cache_prefix + '.h5.npz'
    scipy.sparse.save_npz(cache_fname, X, compressed=False)

    with open(cache_prefix + '.h5.genes.txt', 'w') as of:
        of.write('\n'.join(genes) + '\n')

    return X, genes


def process_h5ad(fname, min_trans=MIN_TRANSCRIPTS):
    X, genes, cells = load_h5ad(fname)

    gt_idx = [i for i, s in enumerate(np.sum(X != 0, axis=1))
              if s >= min_trans]
    X = X[gt_idx, :]
    cells = cells[gt_idx]
    if len(gt_idx) == 0:
        print('Warning: 0 cells passed QC in {}'.format(fname))

    if fname.endswith('.h5ad'):
        cache_prefix = '.'.join(fname.split('.')[:-1])

    cache_fname = cache_prefix + '.h5ad.npz'
    scipy.sparse.save_npz(cache_fname, X, compressed=False)

    with open(cache_prefix + '.h5ad.genes.txt', 'w') as of:
        of.write('\n'.join(genes) + '\n')

    with open(cache_prefix + '.h5ad.cells.txt', 'w') as of:
        of.write('\n'.join(cells) + '\n')

    return X, genes, cells


def load_data(dname):
    read = False
    if os.path.isfile(dname + '.h5.npz') or dname.endswith('.h5.npz'):
        if not dname.endswith('.h5.npz'):
            dname = dname + '.h5.npz'
        X = scipy.sparse.load_npz(dname)
        counts = check_sparse(X, dname)
        with open(dname.replace('.h5.npz', '.h5.genes.txt')) as f:
            genes = np.array(f.read().rstrip().split())
        cells = None
    elif os.path.isfile(dname + '.npz') and dname.endswith('.npz'):
        if not dname.endswith('.npz'):
            dname = dname + '.npz'
        data = np.load(dname)
        X = data['X']
        counts = check_ndarray(X, dname)
        genes = data['genes']
        cells = data['cells']
        data.close()
    elif os.path.isfile(dname + '/tab.npz'):
        X = scipy.sparse.load_npz(dname + '/tab.npz')
        counts = check_sparse(X, dname)
        with open(dname + '/tab.genes.txt') as f:
            genes = np.array(f.read().rstrip().split())
        with open(dname + '/tab.cells.txt', 'r') as f:
            cells = np.array(f.read().rstrip().split())
    elif os.path.isfile(dname + '.raw.dge.txt') or dname.endswith('.raw.dge.txt'):
        if not dname.endswith('.raw.dge.txt'):
            dname = dname + '.raw.dge.txt'
        X = mmread(dname)
        read = True
    elif os.path.isfile(dname + '.h5ad') or dname.endswith('.h5ad'):
        import anndata
        if not dname.endswith('.h5ad'):
            dname = dname + '.h5ad'
        adata = anndata.read_h5ad(dname)
        X = adata.X
        counts = check_ndarray(X, dname)
        cells = np.array(adata.obs.index)
        genes = adata.var.index
    elif os.path.isfile(dname + '.mtx') or dname.endswith('.mtx'):
        if not dname.endswith('.mtx'):
            dname = dname + '.mtx'
        X = mmread(dname)
        read = True
    else:
        sys.stderr.write('Could not find: {}\n'.format(dname))
        exit(1)

    if read:
        counts = check_sparse(X, dname)
        with open('/'.join(dname.split('/')[:-1]) + '/genes.txt', 'r') as f:
            genes = np.array(f.read().rstrip().split())
        with open('/'.join(dname.split('/')[:-1]) + '/cells.txt', 'r') as f:
            cells = np.array(f.read().rstrip().split())
    if len(cells) != X.shape[0]:
        X = X.transpose()
    assert (len(cells), len(genes)) == X.shape

    genes = np.array([ gene.upper() for gene in genes ])
    return X, cells, genes, counts


def check_sparse(X, name):
    cx = X.tocoo()
    dt = str(cx.dtype)
    if not dt.startswith('int'):
        for el in cx.data:
            if not el.is_integer():
                print("WARNING: input matrix for dataset {} is not a matrix of counts".format(name))
                return False
    return True


def check_ndarray(X, name):
    for el in X.flat:
        if el % 1 != 0:
            print("WARNING: input matrix for dataset {} is not a matrix of counts".format(name))
            return False
    return True


def load_names(data_names, norm=False, log1p=False, verbose=True):
    # Load datasets.
    datasets = []
    genes_list = []
    cells_list = []
    n_cells = 0
    counts = True
    for dname in data_names:
        X_i, cells_i, genes_i, c = load_data(dname)
        if not c:
            counts = False
        if norm:
            X_i = normalize(X_i, axis=1)
        if log1p:
            X_i = np.log1p(X_i)
        X_i = csr_matrix(X_i)
            
        datasets.append(X_i)
        genes_list.append(genes_i)
        cells_list.append(cells_i)
        n_cells += X_i.shape[0]
        if verbose:
            print('Loaded {} with {} genes and {} cells'.
                  format(dname, X_i.shape[1], X_i.shape[0]))
    if verbose:
        print('Found {} cells among all datasets'
              .format(n_cells))

    return datasets, genes_list, cells_list, n_cells, counts


def save_datasets(datasets, genes, data_names, verbose=True,
                  truncate_neg=False):
    for i in range(len(datasets)):
        dataset = datasets[i].toarray()
        name = data_names[i]

        if truncate_neg:
            dataset[dataset < 0] = 0

        with open(name + '.scanorama_corrected.txt', 'w') as of:
            # Save header.
            of.write('Genes\t')
            of.write('\t'.join(
                [ 'cell' + str(cell) for cell in range(dataset.shape[0]) ]
            ) + '\n')

            for g in range(dataset.shape[1]):
                of.write(genes[g] + '\t')
                of.write('\t'.join(
                    [ str(expr) for expr in dataset[:, g] ]
                ) + '\n')


def process(data_names, min_trans=MIN_TRANSCRIPTS):
    for name in data_names:
        if os.path.isdir(name):
            process_mtx(name, min_trans=min_trans)
        elif os.path.isfile(name) and name.endswith('.h5'):
            process_h5(name, min_trans=min_trans)
        elif os.path.isfile(name + '.h5'):
            process_h5(name + '.h5', min_trans=min_trans)
        elif os.path.isfile(name + '.h5ad'):
            process_h5ad(name, min_trans=min_trans)
        elif os.path.isfile(name) and name.endswith('.h5ad'):
            process_h5ad(name, min_trans=min_trans)
        elif os.path.isfile(name) and name.endswith('.txt'):
            process_tab(name, min_trans=min_trans)
        elif os.path.isfile(name + '.txt'):
            process_tab(name + '.txt', min_trans=min_trans)
        elif os.path.isfile(name + '.txt.gz'):
            process_tab(name + '.txt.gz', min_trans=min_trans)
        elif os.path.isfile(name + '.tsv'):
            process_tab(name + '.tsv', min_trans=min_trans)
        elif os.path.isfile(name + '.tsv.gz'):
            process_tab(name + '.tsv.gz', min_trans=min_trans)
        else:
            sys.stderr.write('Warning: Could not find {}\n'.format(name))
            continue
        print('Successfully processed {}'.format(name))

if __name__ == '__main__':
    from .config import data_names

    process(data_names)
