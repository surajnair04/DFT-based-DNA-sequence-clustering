import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.cluster import KMeans
#%matplotlib inline

# naive way of getting sequences
def fasta_read(file_name):
    f = open(file_name,"r")
    line = f.read()
    line2 = line.split('\n')
    sequences = []
    for idx,i in enumerate(line2):
        if(idx%4==1):
            sequences.append(i)
    f.close()
    return sequences

def power_spectrum(sequences):
    pspec = []
    for seq in sequences:
        #indicator sequences
        ind = pd.get_dummies(pd.Series(list(seq))) 

        #fft on each nucleotide
        U_A = np.fft.fft(ind['A'])
        U_C = np.fft.fft(ind['C'])
        U_G = np.fft.fft(ind['G'])
        U_T = np.fft.fft(ind['T'])

        #power spectrum of each nucleotide
        PS_A = np.square(np.absolute(U_A))
        PS_C = np.square(np.absolute(U_C))
        PS_G = np.square(np.absolute(U_G))
        PS_T = np.square(np.absolute(U_T))

        #combined power spectrum
        PS = PS_A+PS_C+PS_G+PS_T
        pspec.append(PS[1:len(PS)])         #Exclude the entry corresponding to index 0
    return pspec

def even_scaling(pspec):
    scaled_ps = []
    M = len(max(pspec))
    for ps in pspec:
        N = len(ps)
        temp = []
        temp[0] = ps[0]
        for k in range(1,M,1):
            p = round(k*N/M)
            q = round((k-1)*N/M)
            sc = 0
            for j in range(q,p+1,1):
                sc += ps[j]
            sc /= (p-q+1)
            temp[k] = sc
        scaled_ps.append(temp)
        
    return scaled_ps

if __name__ == "__main__":
    kmeans = KMeans(n_clusters=3, random_state=1)       #No. of clusters=?
    pow_spect = []
    for file_name in os.listdir('readsdata'):
        #print(file_name)
        sequences = fasta_read('readsdata/' + file_name)
        pow_spect.append(power_spectrum(sequences))
    
    pow_spect = [item for sublist in pow_spect for item in sublist]
    #print(pow_spect)
    
    kmeans.fit(pow_spect)
    #print(kmeans.labels_)
    #print(kmeans.cluster_centers_)                     #Cluster Centers Identified

    #Test on given sequences
    new_sequence = fasta_read('readsdata/reads_02_3.fq')
    print(kmeans.predict(power_spectrum(new_sequence)))

    for i in range(len(pow_spect)):
        plt.plot(pow_spect[i][1:200])
    plt.show()

    

