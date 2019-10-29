import argparse
import os
from scipy.io import mmread
import umap
from time import time

parser = argparse.ArgumentParser(description='Calculating UMAP embeddings based on given sparse matrix.')
parser.add_argument('-m', '--matrix', action='store', metavar='MTX_FILE', type=str, required=True,
                    help='Matrix file')
parser.add_argument('--cells', action='store', metavar='FILE', type=str,
                    help='Names of cells')
parser.add_argument('-p', '--path', action='store', metavar='PATH', type=str, help='Directory of the given data')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, help='Directory for the results')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default=None, help='Namespace for the analysis')
args = parser.parse_args()
path = args.path

if args.output is None:
    output = args.path
else:
    output = args.output

matrix = mmread(path + args.matrix).toarray()
namespace = args.namespace
if namespace is None:
    namespace = args.matrix.split('_')[0]
    outfile = args.matrix.replace('matrix', 'umap').replace('.mtx', '.tsv')
else:
    outfile = namespace + '_umap.tsv'

if args.cells is None:
    cells = '{}_cells.txt'.format(namespace)
else:
    cells = args.cells

t0 = time()
print('Starting calculating UMAP embeddings')
embeddings = umap.UMAP().fit_transform(matrix)
print('UMAP calculated in {:.3f} minutes'.format((time() - t0)/60))

cells_file = open(path + cells, 'r').read().strip().split('\n')
with open(output + outfile, 'w') as o:
    o.write('Cell ID\tUMAP-x\tUMAP-y\n')
    for cell, (x, y) in zip(cells_file, embeddings):
        o.write('%s\t%.5f\t%.5f\n' % (cell, x, y))
print('Embeddings written into {} file'.format(outfile))
