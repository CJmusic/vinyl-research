from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

import soundfile as sf

from octave_smoothing import octave_smooth


def groove_map(file_path, filter = False): 
    SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
    soundfile = sf.SoundFile(file_path, 'r')
    # data_sf, fs_sf = sf.read(input, always_2d=True)
    bit_depth = SOUNDFILE_types[soundfile.subtype]
    print 'The bit depth of the audio is: ', bit_depth

    fig = plt.figure(1)
    fs, data = wavfile.read(file_path)

    L = data.T[0]/2**bit_depth
    R = data.T[1]/2**bit_depth

    time = np.arange(0.0,len(data),1) ##calculates a time array in order to plot the waveform of the audio file
    time = time/float(fs)  

    T = 1.8 #seconds, the period of rotation of a record
    T_s = T*fs #period as a number of samples
    # fs = 44100


    num_grooves = int(time[len(time)-1]/T)
    print 'length of files: ', len(time), len(data)
    print 'max_time: ',  time[len(time)-1]
    print 'num_grooves: ', num_grooves
    # num_grooves = 2
    num_grooves = 16

    for i in xrange(1,num_grooves): ##skip the first groove as this typically contains the needle drop
        start = int(i*T*fs)
        end = int(i*T*fs + T*fs)
        block_L = L[start:end]
        block_R = R[start:end]
        time = np.linspace(0.0,float(T*fs),T*fs)/float(fs) ##calculates a time array in order to plot the waveform of the audio file


        if filter == True:
            nyq = 0.5 * fs
            b, a = signal.butter(1, 20.0/nyq, btype='low', analog=False)
            block_L = signal.lfilter(b, a, block_L)
            block_R = signal.lfilter(b, a, block_R)


        fft_block_L = fft(block_L)
        fft_block_L = fft_block_L[0:(len(block_L)/2 + 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
        freq = np.linspace(0.0,fs/2.0,len(block_L)/2 + 1)

        octave_smooth(fft_block_L, freq)

        plt.subplot(221)    
        # plt.plot(time, block_L, label = 'groove %i' % i)
        plt.plot(block_L, label = 'groove %i' % i)

        # plt.grid(which = 'both')

        plt.subplot(223)    
        plt.plot(freq, 20*np.log10(fft_block_L/(2**bit_depth)), label = 'groove %i' % i)
        # plt.xscale('log')
        # plt.grid(which = 'both')

        if num_grooves < 6:
            plt.legend()

        
        fft_block_R = fft(block_R)
        fft_block_R = fft_block_R[0:(len(block_R)/2 + 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
        freq = np.linspace(0.0,fs/2.0,len(block_R)/2 + 1)
        
        plt.subplot(222)    
        # plt.plot(time, block_R, label = 'groove %i' % i)
        plt.plot(block_R, label = 'groove %i' % i)

        # plt.grid(which = 'both')

        plt.subplot(224)    
        plt.plot(freq, 20*np.log10(fft_block_R/(2**bit_depth)), label = 'groove %i' % i)
        plt.xscale('log')
        # plt.grid(which = 'both')
        if num_grooves < 6:
            plt.legend()




        
    plt.subplot(221)
    plt.grid(which = 'both')
    plt.title('Left Channel')
    # plt.ylim(-1,1)
    # plt.xlim(0,1.8)
    plt.ylabel('Amplitude')
    plt.xlabel('Time [s]')

    plt.subplot(222)
    plt.grid(which = 'both')
    plt.title('Right Channel')
    # plt.ylim(-1,1)
    # plt.xlim(0,1.8)
    plt.ylabel('Amplitude')
    plt.xlabel('Time [s]')

    plt.subplot(223)
    plt.grid(which = 'both')
    plt.xscale('log')
    plt.ylabel('Power [dB]')
    plt.xlabel('Freq [Hz]')

    plt.subplot(224)
    plt.grid(which = 'both')
    plt.xscale('log')
    plt.ylabel('Power [dB]')
    plt.xlabel('Freq [Hz]')


    plt.show()


audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
file_dir = '1015_18_LiteToneTest/'
file_name = '5.3declicked.wav'
# file_dir = '1101_18_LiteTone45rpm/'
# file_name = '45-5.2.wav'
# file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
# file_name = 'Silence 0_05.wav'
# file_name = 'Dec14-TestA1n-tone1000.wav'

# file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
# file_name = 'Silence 7_45.wav'
# file_dir = '00_digital_files/'
# file_name = 'tone100.wav'
file_dir = '1108_18_LiteToneMusicSample/'
file_name = '1-intro.wav'

file_path = audio_dir + file_dir + file_name
groove_map(file_path, filter = True)
groove_map(file_path, filter = False)


