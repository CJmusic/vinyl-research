from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

import soundfile as sf



def fft_audio(audio_file, phase = False):
    data = audio_file['data'] 
    fs = audio_file['fs'] 
    bit_depth = audio_file['bit_depth']
    # if stereo == False: 
        # data = data[0] ##by default take the first channel

    fft_a = fft(data)

    # if phase == True: 
    #     fft_a = fft_a[0:(len(fft_a)/2 + 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
    #     freq = np.linspace(0.0,fs/2.0,len(fft_a)/2 + 1)
    #     phi = np.real(fft_a[0])/(np.real(fft_a[0]) - np.imag(fft_a[0]))*np.pi ##haven't tested out the phase part of the code yet

    #     return fft_a, freq, phi

    # else:
    fft_a = fft_a[0:(len(fft_a)/2 - 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
    freq = np.linspace(0.0,fs/2.0,len(fft_a)/2 - 1)
    return np.array(freq_a), np.array(fft_a)


def fft_array(data, fs = 44100, phase = False):
    print 'MAX: ', max(data)
    fft_o = fft(data)

    fft_a = fft_o[0:(len(fft_o)/2 - 1)]
    freq_a = np.linspace(0.0,fs/2.0,len(fft_o)/2 - 1)
    # print 'FFT LENGTH: ', np.shape(fft_a)
    # print 'FREQ LENGTH: ', np.shape(freq_a)
    return np.array(freq_a), np.array(fft_a)