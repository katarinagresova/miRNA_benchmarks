import argparse
import numpy as np

def parse_args(tool_name):
    parser = argparse.ArgumentParser(description=tool_name + ' prediction.')
    parser.add_argument('--input', type=str, help='Path to the input file - miRNA and a gene sequence in a tab-separated format.')
    parser.add_argument('--miRNA_column', type=str, help='Name of the column containing miRNA sequences')
    parser.add_argument('--gene_column', type=str, help='Name of the column containing gene sequences')
    parser.add_argument('--output', type=str, help='Path to the output file')
    return parser.parse_args()

def one_hot_encoding(df, miRNA_col, gene_col, tensor_dim=(50, 20, 1)):
    """
    fun encodes miRNAs and mRNAs in df into binding matrices
    :param df: dataframe containing 'gene' and 'miRNA' columns
    :param tensor_dim: output shape of the matrix
    :return: numpy array of predictions
    """

    miRNA_length = 20
    gene_length = 50

    # alphabet for watson-crick interactions.
    alphabet = {"AT": 1., "TA": 1., "GC": 1., "CG": 1., "AU": 1., "UA": 1.}
    # create empty main 2d matrix array
    N = df.shape[0]  # number of samples in df
    shape_matrix_2d = (N, *tensor_dim)  # 2d matrix shape
    # initialize dot matrix with zeros
    ohe_matrix_2d = np.zeros(shape_matrix_2d, dtype="float32")

    # compile matrix with watson-crick interactions.
    for index, row in df.iterrows():
        for bind_index, bind_nt in enumerate(row[gene_col][:gene_length].upper()):
            for mirna_index, mirna_nt in enumerate(row[miRNA_col][:miRNA_length].upper()):
                base_pairs = bind_nt + mirna_nt
                ohe_matrix_2d[index, bind_index, mirna_index, 0] = alphabet.get(base_pairs, 0)

    return ohe_matrix_2d