import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from Bio import SeqIO
from collections import Counter
from scipy import interpolate
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import pydot
import pylab
%matplotlib inline

def power_spectrum(sequences):
    pspec = []
    for seq in sequences:
        #indicator sequences
        ind = pd.get_dummies(pd.Series(list(seq))) 
        del ind['$']
        ind = np.fft.fft(ind,axis=0)
        ind = np.absolute(ind)**2
        ind = np.sum(ind,axis=1)
        pspec.append(ind)
    return pspec

def linear_scaling(pspec):
    scaled_ps = []
    M = len(max(pspec,key=len))
    for ps in pspec:
        N = len(ps)
        bound = np.arange(N)
        ip = interpolate.interp1d(bound,ps,fill_value='extrapolate')
        r = np.arange(0,N,N/M)
        y = ip(r)
        scaled_ps.append(y)
    return scaled_ps

def cubic_scaling(pspec):
    scaled_ps = []
    M = len(max(pspec,key=len))
    for ps in pspec:
        N = len(ps)
        bound = np.arange(N)
        ip = interpolate.CubicSpline(bound,ps,extrapolate=True)
        r = np.arange(0,N,N/M)
        y = ip(r)
        scaled_ps.append(y)
    return scaled_ps

def get_euclidean_distance(scaled_ps):
    dist_mat = []
    for i in range(1,len(scaled_ps)+1):
        row = []
        for j in range(i):
            if((i-1)==j):
                row.append(0)
            else:
                row.append(np.linalg.norm(scaled_ps[j]-scaled_ps[i-1]))
        dist_mat.append(row)
    return dist_mat

if __name__ == "__main__":
    names = []
    sequences = []
    with open("kmercoded.fasta", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq
            N = len(seq)
            c = Counter(seq)
            if(c['$']==N):
                continue
            names.append(record.description)
            sequences.append(record.seq)
        pspec = power_spectrum(sequences)
        scaled_ps = linear_scaling(pspec) #cubic_scaling(pspec)
        dist_mat = get_euclidean_distance(scaled_ps)
        constructor = DistanceTreeConstructor()
        distance_matrix_10 = _DistanceMatrix(names=names[0:10], matrix=dist_arr[0:10])
        tree_upgma_10 = constructor.upgma(distance_matrix_10)
        Phylo.draw(tree_upgma_10)
        tree_nj_10 = constructor.nj(distance_matrix_10)
        Phylo.draw(tree_nj_10)
    

