import numpy as np
from scanorama.scanorama import *
from scipy.stats import ttest_ind
from sklearn.metrics import silhouette_samples as sil
from scipy.io import mmwrite

from bin.process import load_names

if __name__ == '__main__':
    from bin.config import data_names, namespace, path, output

    labels = np.array(
        open(path + 'cell_labels/all.txt').read().rstrip().split()
    )
    idx = range(labels.shape[0])
    #idx = np.random.choice(len(labels), size=int(len(labels)/2), replace=False)
    
    datasets, genes_list, cells_list, n_cells = load_names(data_names)
    datasets, genes = merge_datasets(datasets, genes_list)
    datasets_dimred, datasets_norm, genes = process_data(datasets, genes)

    mmwrite(output + 'panorama_silh_matrix.mtx', vstack(datasets_dimred))

    # Baseline without correction.
    X = np.concatenate(datasets_dimred)
    sil_non = sil(X[idx, :], labels[idx])
    print(np.median(sil_non))

    # scran MNN.
    X = np.loadtxt(path + 'corrected_mnn.txt')
    sil_mnn = sil(X[idx, :], labels[idx])
    print(np.median(sil_mnn))

    # Seurat CCA.
    X = np.loadtxt(path + 'corrected_seurat.txt')
    sil_cca = sil(X[idx, :], labels[idx])
    print(np.median(sil_cca))

    # Scanorama.
    X = np.concatenate(assemble(datasets_dimred, sigma=150))
    sil_pan = sil(X[idx, :], labels[idx])
    print(np.median(sil_pan))

    print(ttest_ind(sil_pan, sil_non))
    print(ttest_ind(sil_pan, sil_cca))
    print(ttest_ind(sil_pan, sil_mnn))
    
    plt.figure()
    plt.boxplot([ sil_non, sil_mnn, sil_cca, sil_pan ], showmeans=True, whis='range')
    plt.ylim([ -1, 1 ])
    plt.title('Distributions of Silhouette Coefficients')
    plt.xticks(range(1, 5), [ 'No correction', 'scran MNN', 'Seurat CCA', 'Scanorama' ])
    plt.ylabel('Silhouette Coefficient')
    plt.savefig(output + '{}_silhouette.svg'.format(namespace))
