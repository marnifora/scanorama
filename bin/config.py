import fileinput
import argparse
import os

parser = argparse.ArgumentParser(description='Run Scanorama based on given datasets.')
parser.add_argument('file', action='store', metavar='FILE', type=str,
                    help='Name of the configuration file (containing names and paths to the datasets)')
parser.add_argument('-p', '--path', action='store', metavar='PATH', type=str, required=False, default='./',
                    help='Path to the directory with the given datasets to process')
parser.add_argument('-o', '--output', action='store', metavar='OUT', type=str, required=False,
                    default='./results/', help='Directory for the results')
parser.add_argument('-n', '--namespace', action='store', metavar='NAMESPACE', type=str, required=False,
                    default='panorama', help='Namespace for the result files')
parser.add_argument('-m', '--metadata', action='store', metavar='METADATA', type=str, required=False,
                    default=None, help='Optional file with metadata')
parser.add_argument('-pc', '--dimred', action='store', metavar='NUMBER', type=int, required=False,
                    default=100, help='Number of dimensions after reduction')
parser.add_argument('-w', '--write', action='store_true')
parser.add_argument('--tsne', action='store_true')
parser.add_argument('--uncorrected', action='store_true')
args = parser.parse_args()

path, output, namespace, write, tsne, uncorrected, dimred = args.path, args.output, args.namespace, args.write, \
                                                            args.tsne, args.uncorrected, args.dimred

if not os.path.isdir(output):
    os.mkdir(output)

names, data_names = [], []
with open(args.file, 'r') as file:
    for line in file:
        l = line.rstrip().split('\t')
        if len(l) == 2:
            names.append(l[0])
            dn = (path + '/' + l[1]).replace('///', '/').replace('//', '/')
            data_names.append(dn)
        else:
            fields = line.rstrip().split(',')
            for f in fields:
                if f.strip() == '':
                    continue
                dn = (path + '/' + f).replace('///', '/').replace('//', '/')
                data_names.append(dn)
    if not names and data_names:
        names = [el.strip().split('/')[-1].split('.')[0] for el in data_names]
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
