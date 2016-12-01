import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from Bio import SeqIO
from collections import Counter
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
        
    

