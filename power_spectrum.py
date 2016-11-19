import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

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
        pspec.append(PS)
    return pspec

if __name__ == "__main__":
    sequences = fasta_read('readsdata/reads_02_3.fq')
    pspec = power_spectrum(sequences)

    # plt.plot(pspec[500][1:200])
    # plt.show()

    
