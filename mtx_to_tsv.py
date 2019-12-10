import argparse
from scipy.io import mmread
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Rewrite mtx matrix into tsv file')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-mtx', action='store', metavar='FILE', type=str, default=None,
                    help='Mtx file - matrix of counts')
parser.add_argument('-p', '--path', action='store', metavar='PATH', type=str, required=False, default='./',
                    help='Path to the directory with the given datasets to process')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    default=None, help='Directory for the results, if not given - it is [PATH]/results/')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default=None, help='Namespace for the result files, if not given - the same as input file')
parser.add_argument('-m', '--metadata', action='store', metavar='METADATA', type=str, required=False,
                    default=None, help='Optional file with metadata')
parser.add_argument('-g', '--genes', action='store', metavar='FILE', type=str, default=None,
                    help='Name of the genes file, deafult: [NAMESPACE]_genes.txt')
parser.add_argument('-c', '--cells', action='store', metavar='FILE', type=str, default=None,
                    help='Name of the cells file, deafult: [NAMESPACE]_cells.txt')
args = parser.parse_args()

if args.path is not None:
    path = args.path
    file = args.mtx
else:
    s = args.mtx.strip('/').split('/')
    path = '/' + '/'.join(s[:-1])
    file = s[-1]

if args.output is not None:
    output = args.output
else:
    output = path

if args.namespace is not None:
    namespace = args.namespace
else:
    namespace = file.split('_')[0]

if args.cells is None:
    cells = os.path.join(path, '{}_cells.txt'.format(namespace))
else:
    cells = args.cells

if args.genes is None:
    genes = os.path.join(path, '{}_genes.txt'.format(namespace))
else:
    genes = args.genes

m = mmread(os.path.join(path, file)).toarray()
index = open(cells, 'r').read().strip().split('\n')
header = open(genes, 'r').read().strip().split('\n')

df = pd.DataFrame(m, index=index, columns=header)
df.to_csv(os.path.join(output, file.replace('.mtx', '.tsv')), sep='\t', float_format='%.3f')
