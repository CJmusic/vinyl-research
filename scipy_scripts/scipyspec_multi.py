from scipy import signal
from scipy.io import wavfile
import matplotlib.pyplot as plt
import numpy as np
import os

audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
file_dir = '1015_18_LiteToneTest/'

for filename in os.listdir(audio_dir+file_dir):
    if filename.endswith(".wav"):
        fs, x = wavfile.read(audio_dir+file_dir+filename)
        L = x.T[0]
        f, t, Sxx = signal.spectrogram(L, fs)#, nfft = 2**16)

        time = np.arange(0.0,len(L),1)
        time = time*fs

        plt.figure()
        plt.subplot(211)
        plt.pcolormesh(t, f, (Sxx))
        plt.ylabel('Frequency [Hz]')
        # plt.yscale('symlog')
        plt.xlabel('Time [sec]')
        plt.colorbar()
        plt.subplot(212)
        plt.plot(time,L)
        plt.xlabel('Time [sec]')
        plt.ylabel('Amplitude')
        plt.savefig('plots/'+filename+'.png')
        f = None 
        t = None 
        x = None 
        fs = None 
        L = None 
        time = None 
        Sxx = None