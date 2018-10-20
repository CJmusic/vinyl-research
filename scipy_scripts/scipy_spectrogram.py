from scipy import signal
from scipy.io import wavfile


import matplotlib.pyplot as plt
import numpy as np



audio_dir = '/Users/cz/OneDrive\ -\ University\ of\ Waterloo/Vinyl_Project/audio_files/'
file_dir = '015_18_LiteToneTest/'

# def pwroctsmooth(Spec_a, octave_width = 0.33):
    # N = len(Spec_a)
    # an = 2.0**(octave_width/2.0)

    # old_lo = 0.0
    # old_hi = 0.0 

    # for i in xrange(N/2 + 1):
    #     lo_bin = np.round((i-1)/an) + 1
    #     hi_bin = np.round((i-1)*an)) + 1
    #     if hi_bin == lo_bin: 
    #         pwr_sum = abs(Spec_a[i]))**2.0
    #     else: 
    #         for J in xrange(old_lo, lo_bin-1): 
    #             pwr_sum = pwr_sum + abs(Spec_a[J])**2.0
    #         for K in xrange(old_hi+1, hi_bin)
    #             pwr_sum = pwr_sum + abs(Spec_a[K])**2.0
    #     smoothed_spec(range(np.floor(N/2)+2,N)) = np.conj(smoothed_spec[(np.ceil(N/2)]) 
    #     if N/2==floor(N/2):
    #         smoothed_spec(N/2+1)=abs(smoothed_tf(N/2+1))

def octave_smooth(Spec, octave_wdith = 0.33):
    N = len(Spec)
    for i in xrange(2,N-2,1): 
        Spec[i] = (Spec[i-2] + Spec[i-1] + Spec[i] + Spec[i+1] + Spec[i+2])/5.0
    return Spec
# fs, x = wavfile.read('110918UW_stuckontop_r1B.wav')
# fs, x = wavfile.read('110918UW_stuckontop_r1A.wav')


# fs, x = wavfile.read('sweep16kHz.wav')
fs, x = wavfile.read('13.5.wav')

# fs, x = wavfile.read('2.5.wav')
# fs, x = wavfile.read('tone1000.wav')

L = x.T[0]
# L = L.tolist()
# L = np.array(L[len(x)/2  + len(x)/3:len(x)/2 +len(x)/3+ 2**16])
print 'len(L): ', len(L)
f, t, Sxx = signal.spectrogram(L, fs, nperseg = 2**12, nfft = 2**12, noverlap=60)
##Default nperseg = 256


# f, t, Sxx = signal.spectrogram(snd_block, RATE)   
dBS = 10 * np.log10(Sxx)  # convert to dB
# plt.pcolormesh(t, f, dBS)

time = np.arange(0.0,len(L),1)
time = time/fs

# Sxx = octave_smooth(Sxx)

plt.subplot(211)
plt.pcolormesh(t, f, dBS, cmap='inferno')
plt.ylabel('Frequency [Hz]')
plt.colorbar()
# plt.yscale('symlog')
plt.xlabel('Time [sec]')
plt.subplot(212)
plt.plot(time,L)
plt.show()