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
                    default='./results/', help='Directory for the results')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default='panorama', help='Namespace for the result files')
parser.add_argument('-m', '--metadata', action='store', metavar='METADATA', type=str, required=False,
                    default=None, help='Optional file with metadata')
parser.add_argument('-pc', '--dimred', action='store', metavar='NUMBER', type=int, required=False,
                    default=100, help='Number of dimensions after reduction, default = 100')
parser.add_argument('-w', '--write', action='store_true')
parser.add_argument('--tsne', action='store_true')
parser.add_argument('--uncorrected', action='store_true')
parser.add_argument('--scanpy_pca', action='store_true')
parser.add_argument('--no_adjust', action='store_true')
args = parser.parse_args()

path, output, namespace, write, tsne, uncorrected, dimred, metadata, scanpy = \
    args.path, args.output, args.namespace, args.write, args.tsne, args.uncorrected, args.dimred, args.metadata, \
    args.scanpy_pca

if args.no_adjust:
    adjust = False
else:
    adjust = True

if not os.path.isdir(output):
    os.mkdir(output)

datasets = []
if args.mtx is not None:
    from scipy.io import mmread
    import numpy as np

    clusters = {}
    with open(os.path.join(path, args.clusters), 'r') as f:
        for i, line in enumerate(f):
            clusters[line.strip()] = clusters.setdefault(line.strip(), []) + [i]

    matrix = mmread(os.path.join(path, args.mtx))
    for cl in clusters.values():
        m = matrix[np.array(cl), :]
        datasets.append(m)
else:
    names, data_names = [], []
    with open(args.file, 'r') as file:
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

if args.metadata is not None:
    metafile = args.metadata
else:
    metafile = '.'.join(args.file.split('.')[:-1]) + '_metadata.txt'
if os.path.isfile(metafile):
    with open(metafile, 'r') as f:
        metadata = {}
        for line in f:
            line = line.strip().split('\t')
            if line[0] in names:
                metadata[line[0]] = line[1:]
else:
    metadata = None
