from bin.config import data_names, names, namespace, output
from scanorama.scanorama import check_datasets, merge_datasets
from bin.process import load_names

datasets, genes_list, cells_list, n_cells = load_names(data_names)
datasets_full = check_datasets(datasets)
datasets, genes = merge_datasets(datasets_full, genes_list,
                                 ds_names=names, union=False)

with open(output + '%s_genes_list.txt' % namespace, 'w') as o:
    o.write('\n'.join(genes))
with open(output + '%s_cells_list.txt' % namespace, 'w') as o:
    o.write('\n'.join([ll for l in cells_list for ll in l]))

