#!/usr/bin/env python3

import sys

import numpy
import matplotlib.pyplot as plt
from matplotlib import cm, colors, markers
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages


def main():
    in_fname, out_fname = sys.argv[1:3]
    dtypes = [("Genome", "<U4"), ("CT", "<U30"), ("R2", numpy.float32),
              ('Experiment', "<U8" ), ('Condition', "<U10"), ('Features', "<U8"), ("Rep", numpy.int32)]
    alldata = numpy.loadtxt(in_fname, dtype=numpy.dtype(dtypes), skiprows=1, usecols=(1,2,3,4,5,6,7))
    alldata = alldata[numpy.where(alldata['Rep'] == 0)]
    fig, all_ax = plt.subplots(2, 1, figsize=(8,6),
                               gridspec_kw={'height_ratios': (7, 1)})
    ax = all_ax[0]
    all_CTs = numpy.unique(alldata['CT'])
    cmap = cm.get_cmap('tab20')
    color_dict = {}
    for i, ct in enumerate(all_CTs):
        color_dict[ct] = cmap((i + 0.5) / all_CTs.shape[0])
    shape_dict = {'mm10': "o", 'hg38': "D"}
    for g in ['mm10', 'hg38']:
        data = alldata[numpy.where(alldata['Genome'] == g)]
        pcolor = []
        pshape = []
        X = []
        Y = []
        offset0 = {'standard:' 0.5, 'nearest':6.5}
        offset1 = {'both': 0.5, 'promoter': 2.5, 'cre': 4.5}
        offset2 = {'shufflecre': 0.1, 'shuffletss': 0.45, 'control': 0.8}
        for i in range(data.shape[0]):
            Y.append(data['R2'][i])
            X.append(offset1[data['Features'][i]] + offset2[data['Condition'][i]] + offset0[data['Experiment'][i]])
            X[-1] += (numpy.random.random() - 0.5) * 0.5
            pcolor.append(color_dict[data['CT'][i]])
        ax.scatter(X, Y, color=pcolor, marker=shape_dict[g])
    ax.set_ylabel(r'Adjusted R^2')
    ax.set_xticks([1.1, 1.45, 1.8, 3.1, 3.45, 3.8, 5.1, 5.45, 5.8, 7.1, 7.45, 7.8, 9.1, 9.45, 9.8, 11.1, 11.45, 11.8])
    ax.set_xticklabels(numpy.rep(['ShC', 'ShP', 'C'], 6))
    ax.text(2.5, -1, 'Standard Refinement', horizontalalignment='center')
    ax.text(7.5, -1, 'Nearest Gene', horizontalalignment='center')
    ax.text(1, -0.1, "Promoters+CREs", horizontalalignment='center')
    ax.text(7, -0.1, "Promoters+CREs", horizontalalignment='center')
    ax.text(9, -0.1, "Promoters", horizontalalignment='center')
    ax.text(11, -.1, "CREs", horizontalalignment='center')
    ax = all_ax[1]
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 3)
    ax.scatter([0.1], [2], marker='o')
    ax.annotate('mm10', (0.17, 1.75))
    ax.scatter([0.1], [1], marker='D')
    ax.annotate('hg38', (0.17, 0.75))
    X = numpy.arange(all_CTs.shape[0]) // 3 + 1
    Y = numpy.arange(all_CTs.shape[0]) % 3 + 0.5
    pcolors = []
    for ct in all_CTs:
        pcolors.append(color_dict[ct])
    ax.scatter(X, Y, marker='s', color=pcolors)
    for i, txt in enumerate(all_CTs):
        ax.annotate(txt, (X[i] + 0.07, Y[i] - 0.25), color=pcolors[i])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.tight_layout()
    fig.savefig(out_fname)


main()
