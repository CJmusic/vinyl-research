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

def stereo_normalize(block):
    #normalizes the left and right channels of a stereo file to max
    return block/max(block)


def octave_smooth(Spec, octave_wdith = 0.33):
    N = len(Spec)
    for i in xrange(2,N-2,1): 
        Spec[i] = (Spec[i-2] + Spec[i-1] + Spec[i] + Spec[i+1] + Spec[i+2])/5.0
    return Spec

def mono_spectrogram(file):
    fs, x = wavfile.read(file)
    L = x.T[0]
    time = np.arange(0.0,len(L),1) ##calculates a time array in order to plot the waveform of the audio file
    time = time/float(fs)

    f, t, Sxx = signal.spectrogram(L, fs, nperseg = 256, nfft = 256, noverlap=60) ##Default nperseg = 256, nfft =, noverlap = 
    dBS = 10 * np.log10(Sxx) ##converts the Spectrogram to deciBels 

    plt.figure(1)
    # Create room on the right
    plt.gcf().subplots_adjust(right=0.8)


    plt.subplot(211)
    plt.pcolormesh(t, f, dBS, cmap='inferno')
    plt.ylabel('Frequency [Hz]')
    # plt.colorbar()
    # plt.yscale('symlog') ##this logscales the frequency axis, however as of now it does not look good 
    plt.subplot(212)
    plt.plot(time,L)
    plt.xlabel('Time [sec]')
    plt.ylabel ('Amplitude')
    plt.xlim(xmin=0,xmax=max(time))
    cbar_ax = plt.gcf().add_axes([0.85, 0.15, 0.05, 0.7])
    # plt.xlim(min(t),max(t))
    cbar = plt.colorbar(cax=cbar_ax)
    cbar.set_label('Amplitude [dB]')
    plt.show()

def stereo_spectrogram(file, normalization = True):
    fs, x = wavfile.read(file)
    L = x.T[0]
    R = x.T[1]
    if normalization == True:
        L = stereo_normalize(L)
        R = stereo_normalize(R)

    time = np.arange(0.0,len(x),1) ##calculates a time array in order to plot the waveform of the audio file
    time = time/float(fs)   

    f, t, Sxx_L = signal.spectrogram(L, fs, nperseg = 256, nfft = 256, noverlap=0) ##Default nperseg = 256, nfft =, noverlap = 
    f, t, Sxx_R = signal.spectrogram(R, fs, nperseg = 256, nfft = 256, noverlap=0) ##Default nperseg = 256, nfft =, noverlap = 
    
    
    dBS_L = 10 * np.log10(Sxx_L) ##converts the Spectrogram to deciBels 
    dBS_R = 10 * np.log10(Sxx_R) ##converts the Spectrogram to deciBels 

    plt.figure(1)
    # Create room on the right
    plt.gcf().subplots_adjust(right=0.8)


    plt.subplot(311)
    plt.pcolormesh(t, f, dBS_L, cmap='summer')#, vmin = -80.0 , vmax = -40.0)
    plt.title('Left Channel')
    plt.ylabel('Frequency [Hz]')
    cbar_axL = plt.gcf().add_axes([0.85, 0.15, 0.01, 0.7])
    cbarL = plt.colorbar(cax=cbar_axL)

    plt.subplot(312)
    plt.pcolormesh(t, f, dBS_R, cmap='summer')#, vmin = -80.0 , vmax = -40.0)
    plt.title('Right Channel')
    plt.ylabel('Frequency [Hz]')
    # cbar_axR = plt.gcf().add_axes([0.95, 0.15, 0.005, 0.7])
    # cbarR = plt.colorbar(cax=cbar_axR)


    # plt.yscale('symlog') ##this logscales the frequency axis, however as of now it does not look good 
    plt.subplot(313)
    plt.plot(time, L, 'orange')
    plt.plot(time, R, 'blue')
    plt.title('Waveform')
    plt.xlabel('Time [sec]')
    plt.ylabel ('Amplitude')
    plt.xlim(xmin=0,xmax=max(time))
    # cbar_ax = plt.gcf().add_axes([0.9, 0.15, 0.025, 0.7])
    # plt.xlim(min(t),max(t))
    # cbar = plt.colorbar(cax=cbar_ax)
    cbarL.set_label('Amplitude [dB]')
    plt.show()



# mono_spectrogram('13.5.wav')
stereo_spectrogram('13.5.wav')
