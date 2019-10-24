from scanorama.scanorama import calculate_tsne
import argparse
import os
from scipy.io import mmread
import pandas as pd
from sklearn.metrics import silhouette_samples as silh_samples
from statistics import mean, median
from time import time

parser = argparse.ArgumentParser(description='Counting t-SNE embeddings based on given sparse matrix.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-m', '--matrix', action='store', metavar='MTX_FILE', type=str,
                    help='Matrix file')
group1.add_argument('-tsne', '--tsne', action='store', metavar='MTX_FILE', type=str,
                    help='File with tSNE coords')
parser.add_argument('-r', '--real', action='store', metavar='FILE', type=str,
                    help='Assigned clusters')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    default='./results/', help='Directory for the results')
args = parser.parse_args()

if args.matrix is not None:
    input = mmread(args.matrix).toarray()
    outfile = args.matrix.split('/')[-1].replace('matrix', 'silh').replace('.mtx', '.txt')
elif args.tsne is not None:
    input = pd.read_csv(args.tsne, sep='\t', header=0, index_col=0).values
    outfile = args.tsne.split('/')[-1].replace('tsne_', 'silh_tsne-').replace('.tsv', '.txt')
print('Starting calculating silh scores for {} file'.format(outfile))

with open(args.real, 'r') as f:
    real = f.read().split('\n')
print('Real labels loaded')

t0 = time()
scores = silh_samples(input, real)
print('Silh scores calculated in {:.2f} minutes'.format((time()-t0)/60))

with open(args.output + outfile, 'w') as f:
    for s in scores:
        f.write('{:.4f}\n'.format(s))

print('Mean of silh scores: {:.3f}'.format(mean(scores)))
print('Median of silh scores: {:.3f}'.format(median(scores)))
