from scipy import signal
from scipy.io import wavfile
import matplotlib.pyplot as plt
import numpy as np
import os

audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
file_dir = '1015_18_LiteToneTest/'


def stereo_normalize(block):
    #normalizes the left and right channels of a stereo file to max
    return block/max(block)

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

    plt.figure(figsize=(16,9),dpi=80)
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
    plt.savefig('plots/'+filename+'.png')


for filename in os.listdir(audio_dir+file_dir):
    if filename.endswith(".wav"):
        stereo_spectrogram(audio_dir + file_dir + filename)