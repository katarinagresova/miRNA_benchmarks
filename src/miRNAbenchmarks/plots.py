import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
import argparse

from miRNAbenchmarks.utils import TOOLS, DATASETS_CONFIG

def plot_precision_recall_single(data, tools, title, save_dir):

    plt.figure(figsize=(5, 5))

    for tool in tools:

        if tool not in data.columns:
            print(f'{tool} not in data')
            continue

        if tool == 'seed':
            p8, r8, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer8'].values, average='binary')
            p7, r7, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer7'].values, average='binary')
            p6, r6, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer6'].values, average='binary')
            p6b, r6b, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer6_bulge_or_mismatch'].values, average='binary')

            plt.plot(r8, p8, 'o', label='8mer')
            plt.plot(r7, p7, 'o', label='7mer')
            plt.plot(r6, p6, 'o', label='6mer')
            plt.plot(r6b, p6b, 'o', label='6mer_bulge_mismatch')

        precision, recall, _ = precision_recall_curve(data['label'], data[tool])
        plt.plot(recall, precision, label=TOOLS[tool])

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # save publication quality figures
    plt.savefig(save_dir + '/' + title + '.png', dpi=300, bbox_inches='tight')

def plot_precision_recall_row(dataset, tools, title, save_dir):

    ratios = DATASETS_CONFIG[dataset]['ratios']
    columns = len(ratios)
    dir_path = DATASETS_CONFIG[dataset]['path']

    fig, ax = plt.subplots(1, columns, figsize=(4*columns, 3.5), dpi=300)
    fig.suptitle(title)

    for i, ratio in enumerate(ratios):
        path = dir_path + f"miRNA_test_set_{ratio}_predictions.tsv"
        data = pd.read_csv(path, sep='\t')

        for tool in tools:

            if tool not in data.columns and tool != 'seed':
                print(f'{tool} not in data')
                continue

            if tool == 'seed':
                p8, r8, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer8'].values, average='binary')
                p7, r7, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer7'].values, average='binary')
                p6, r6, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer6'].values, average='binary')
                p6b, r6b, _, _ = precision_recall_fscore_support(data['label'].values, data['kmer6_bulge_or_mismatch'].values, average='binary')

                ax[i].plot(r8, p8, 'o', label='8mer')
                ax[i].plot(r7, p7, 'o', label='7mer')
                ax[i].plot(r6, p6, 'o', label='6mer')
                ax[i].plot(r6b, p6b, 'o', label='6mer_bulge_mismatch')

            else:
                precision, recall, _ = precision_recall_curve(data['label'], data[tool])
                ax[i].plot(recall, precision, label=TOOLS[tool])

            ax[i].set_xlabel('Recall')
            ax[i].set_ylabel('Precision')

            ax[i].set_xlim(0, 1)
            ax[i].set_ylim(0, 1)

    ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Tight layout often produces nice results
    # but requires the title to be spaced accordingly
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)

    # save publication quality figures
    plt.savefig(save_dir + '/' + title + '.png', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Plot precision-recall curves for different tools')
    parser.add_argument('-t', '--tools', 
                        action='store', 
                        type=str, 
                        help='Names of tools separated by a space. Use "all" to run all tools.', 
                        nargs="+",
                        default=['all'],
                        choices=list(TOOLS.keys()) + ['all'])
    parser.add_argument('-d', '--dataset', 
                        action='store', 
                        type=str, 
                        help='Names of dataset.', 
                        required=True,
                        choices=list(DATASETS_CONFIG.keys()))
    parser.add_argument('--title', help='Title of the plot', required=True)
    parser.add_argument('--save_dir', help='Directory to save the plot', required=True)

    args = parser.parse_args()

    if args.tools == ['all']:
        args.tools = list(TOOLS.keys())

    plot_precision_recall_row(args.dataset, args.tools, args.title, args.save_dir)