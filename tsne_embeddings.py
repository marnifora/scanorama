from scanorama.scanorama import calculate_tsne
import argparse
import os
from scipy.io import mmread

parser = argparse.ArgumentParser(description='Counting t-SNE embeddings based on given sparse matrix.')
parser.add_argument('matrix', action='store', metavar='MTX_FILE', type=str,
                    help='Matrix file')
parser.add_argument('-c', '--cells', action='store', metavar='FILE', type=str,
                    help='File with names of cells, default [NAMESPACE]_cells.txt')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    help='Directory for the results, if not given the same as PATH')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    help='Namespace for the analysis')
parser.add_argument('-p', '--path', action='store', metavar='PATH', type=str, help='Directory of the given data')
parser.add_argument('--n_jobs', action='store', metavar='NUMBER', type=int,
                    default=5, help='Number of jobs for calculating tSNE, default is 5')
parser.add_argument('--nomultitsne', action='store_true', help='If MulticoreTSNE package should NOT be used')
args = parser.parse_args()

if args.nomultitsne:
    multicore_tsne = False
else:
    multicore_tsne = True

path = args.path

if args.output is None:
    output = path
else:
    output = args.output

if args.namespace is None:
    namespace = args.matrix.split('_')[0]
    outfile = args.matrix.replace('matrix', 'tsne').replace('.mtx', '.tsv')
else:
    namespace = args.namespace
    outfile = '{}_tsne.tsv'.format(namespace)

if args.cells is None:
    cells_file = '{}_cells.txt'.format(namespace)
else:
    cells_file = args.cells

matrix = mmread(path + args.matrix)
cells = open(path + cells_file, 'r').read().strip().split('\n')

calculate_tsne(matrix, cells, outfile, output, n_jobs=args.n_jobs, multicore_tsne=multicore_tsne)
