import fileinput
import argparse
import os, sys

parser = argparse.ArgumentParser(description='Run Scanorama based on given datasets.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-mtx', action='store', metavar='FILE', type=str, default=None,
                    help='Mtx file - matrix of counts')
group1.add_argument('-conf', action='store', metavar='FILE', type=str, default=None,
                    help='Name of the configuration file (containing names and paths to the datasets)')
parser.add_argument('-p', '--path', action='store', metavar='PATH', type=str, required=False, default='./',
                    help='Path to the directory with the given datasets to process')
parser.add_argument('--batches', action='store', metavar='FILE', type=str, default=None,
                    help='File with batches - if mtx file containing merged datasets is provided')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    default=None, help='Directory for the results, if not given - it is [PATH]/results/')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default=None, help='Namespace for the result files, if not given - the same as input file')
parser.add_argument('-m', '--metadata', action='store', metavar='METADATA', type=str, required=False,
                    default=None, help='Optional file with metadata')
parser.add_argument('-pc', action='store', metavar='NUMBER', type=int, required=False,
                    default=100, help='Number of dimensions after reduction, default = 100')
parser.add_argument('-w', '--write', action='store_true')
parser.add_argument('--tsne', action='store_true')
parser.add_argument('--uncorrected', action='store_true')
parser.add_argument('--scanpy_pca', action='store_true')
parser.add_argument('--no_norm', action='store_true')
parser.add_argument('--no_dimred', action='store_true')
parser.add_argument('--no_adjust', action='store_true')

args = parser.parse_args()

path, write, tsne, uncorrected, pc, metadata, scanpy = \
    args.path, args.write, args.tsne, args.uncorrected, args.pc, args.metadata, \
    args.scanpy_pca

if args.namespace is None:
    if args.mtx is not None:
        string = args.mtx.strip('.mtx')
    elif args.conf is not None:
        string = args.conf.strip('.txt')
    namespace = string.split('/')[-1].split('_')[0]
else:
    namespace = args.namespace

if args.no_adjust:
    adjust = False
else:
    adjust = True

if args.no_dimred:
    dimred = False
else:
    dimred = True

if args.no_norm:
    norm = False
else:
    norm = True

if args.output is None:
    output = os.path.join(path, '/results')
else:
    output = args.output

if not os.path.isdir(output):
    os.mkdir(output)

if args.metadata is not None:
    metafile = args.metadata
else:
    metafile = '{}_metadata.tsv'.format(namespace)

datasets = []
names, data_names = [], []
if args.mtx is not None:
    from scipy.io import mmread
    import numpy as np

    batches = {}
    with open(os.path.join(path, metafile), 'r') as f:
        col = f.readline().strip().split('\t').index('Dataset')
        for i, line in enumerate(f):
            line = line.strip().split('\t')
            batches[line[col]] = batches.setdefault(line[col], []) + [i]

    genes = open(os.path.join(path, '{}_genes.txt'.format(namespace)), 'r').read().strip().split('\n')
    cells = open(os.path.join(path, '{}_cells.txt'.format(namespace)), 'r').read().strip().split('\n')
    matrix = mmread(os.path.join(path, args.mtx))
    genes_list = []
    cells_list = []
    for b in batches.values():
        b = np.array(b)
        m = matrix[b, :]
        datasets.append(m)
        genes_list.append(np.array(genes[b]))
        cells_list.append(np.array(cells[b]))
    print('{} datasets based on {} have been loaded'.format(len(batches), args.mtx))

else:
    with open(args.conf, 'r') as file:
        for line in file:
            l = line.rstrip().split('\t')
            if len(l) == 2:
                names.append(l[0])
                dn = os.path.join(path, l[1])
                data_names.append(dn)
            else:
                fields = line.rstrip().split(',')
                for f in fields:
                    if f.strip() == '':
                        continue
                    dn = os.path.join(path, f)
                    data_names.append(dn)
        if not names and data_names:
            names = []
            for el in data_names:
                n = el.strip().split('/')[-1].split('.')
                if len(n) > 1:
                    n = '.'.join(n[:-1])
                else:
                    n = n[0]
                names.append(n)
            if len(set(names)) < len(names):
                names = [el.strip().split('/')[-2] for el in data_names]
                if len(set(names)) < len(names):
                    print('ERROR: names of datasets are not unique')
                    sys.exit()
        print('Data names loaded')

if os.path.isfile(metafile):
    with open(metafile, 'r') as f:
        f.readline()
        metadata = {}
        for line in f:
            line = line.strip().split('\t')
            if line[0] in names:
                metadata[line[0]] = line[1:]
else:
    metadata = None
