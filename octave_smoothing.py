import numpy as np 
import matplotlib.pyplot as plt


def octave_smooth(F, freq, width=float(1.0/3.0), fs = 44100, W=None): 
    ##F is the Fourier transform of a signal
    ##freq is the array of frequency bins corresponding to the FFT signal
    #octave_width self explanatory, defaults to 1/3 octave 
    #W is the window or weights to be used in the averaging 
    
    S = np.array(np.shape(F)) #make the smoothed array the same shape as the input
    nbins = len(freq)

    prev_lo = 0 
    prev_hi = 0

    for i in xrange(len(S)): 
        # s[i] = sum(abs(S[i-octave_width:i+octave_width]))   
         
        fcentre = freq[i]

        flower = fcentre*2.0**(-width/2.0)
        lowerbin = int(np.floor(flower*nbins/fs))

        fupper = fcentre*2.0**(width/2.0)
        upperbin = int(np.ceil(fupper*nbins/fs))

        pwr_sum = abs(F[i]*F[i])

        for j in xrange(prev_lo,lowerbin): 
            pwr_sum -= abs(F[j]*F[j])
        
        for k in xrange(prev_hi, upperbin): 
            pwr_sum += abs(F[k]*F[k])

        S[i] = np.sqrt(pwr_sum/(upperbin-lowerbin+1))

        prev_lo = lowerbin
        prev_hi = upperbin
x
    ##need to still ensure S is conjugate even 

    #and that the nyquist value is real 

    return S

