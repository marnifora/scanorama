from bin.config import namespace, output

embeddings = []
o = open(output + namespace + '_metadata.tsv', 'r')
w = open(output + namespace + '_truemeta.tsv', 'w')
with open(output + namespace + '_cells_list.txt', 'r') as oo:
    cells = oo.read().strip().split()
print('number of cells: %d' % len(cells))

w.write(o.readline())
for cell, meta in zip(cells, o):
    w.write('%s\t%s\n' % (cell, '\t'.join(meta.strip().split()[1:])))
o.close()
w.close()

