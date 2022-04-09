#!/usr/bin/env python

import io
import os
import pandas as pd
from pickle import LIST, NONE
from posixpath import split
import allel
from sklearn.decomposition import PCA
import numpy as np

file_name = './for_aley.vcf'

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

dictionary = {}
callset = read_vcf(file_name)
important_names = []

for index, row in callset.iterrows():
    for col in callset.columns:

        dictionary[row["POS"]] = []

        for col in row.keys():
            # print(col)
            # print(row[col])
            if col.isdigit():
                a, b = row[col].split('|')
                
                dictionary[row["POS"]].append(int(a) + int(b))

SUPERLIST = []
for pos in dictionary:
    SUPERLIST.append(dictionary[pos])
NUMPYMATRIX = np.matrix(SUPERLIST)
print(NUMPYMATRIX)


k = 1 # target dimension(s)
pca = PCA(k) # Create a new PCA instance

data = np.array(SUPERLIST) # 2x2 data matrix
print("Data: ", data)
print("Reduced: ", pca.fit_transform(data)) # fit and transform

data = data - data.mean(axis=0) # Center data points
print("Centered Matrix: ", data)
cov = np.cov(data.T) / data.shape[0] # Get covariance matrix
print("Covariance matrix: ", cov)
v, w = np.linalg.eig(cov)

idx = v.argsort()[::-1] # Sort descending and get sorted indices
v = v[idx] # Use indices on eigv vector
w = w[:,idx] # 

print("Eigenvalue vektoru: ", v)
print("Eigenvektorler: ", w)


print("Sonuc: ", data.dot(w[:, :k]))


