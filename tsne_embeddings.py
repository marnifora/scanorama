from scanorama.scanorama import calculate_tsne
import argparse
import os
from scipy.io import mmread

parser = argparse.ArgumentParser(description='Counting t-SNE embeddings based on given sparse matrix.')
parser.add_argument('matrix', action='store', metavar='MTX_FILE', type=str,
                    help='Matrix file')
parser.add_argument('cells', action='store', metavar='FILE', type=str,
                    help='Names of cells')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    default='./results/', help='Directory for the results')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default=None, help='Namespace for the analysis')
args = parser.parse_args()

output = args.output
if not os.path.isdir(output):
    os.mkdir(output)
matrix = mmread(args.matrix)
cells = open(args.cells, 'r').read().strip().split('\n')
namespace = args.namespace
if namespace is None:
    namespace = args.matrix.split('/')[-1].split('_')[0]

calculate_tsne(matrix, cells, namespace, output)
