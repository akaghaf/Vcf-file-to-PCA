from pickle import LIST, NONE
from posixpath import split
import allel
import numpy as np

callset = allel.read_vcf('for_aley.vcf')

gt = allel.GenotypeArray(callset['calldata/GT'])
STRLIST = str(gt)
LISTS = STRLIST.split("\n")
split_list = [i.replace("\'", "").split(' ', len(LISTS)) for i in LISTS]
# for pair in  split_list:
#     pairsum = int(pair[0]) + int(pair[2])

#for ploidyset in range(0, len(LISTS):
#NUMPYMATRIX = np.matrix(split_list)
#print(NUMPYMATRIX)
print(split_list)
#print(pairsum)