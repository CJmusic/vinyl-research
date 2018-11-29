from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

import soundfile as sf

from octave_smoothing import octave_smooth
from fft_audio import fft_audio, fft_array


def groove_map(audio_file, T=1.8): 
    data = audio_file['data'] 
    fs = audio_file['fs'] 
    bit_depth = audio_file['bit_depth']

    L = data.T[0]#/2**bit_depth
    R = data.T[1]#/2**bit_depth

    time = np.arange(0.0,len(data),1)/fs ##calculates a time array in order to plot the waveform of the audio file
    num_grooves = int(len(data)/fs/T)

    print 'length of files: ', len(time), len(data)
    print 'max_time: ',  time[len(time)-1]
    print 'num_grooves: ', num_grooves

    groove_array = np.zeros([num_grooves,2,int(T*fs)]) #
    time_array = np.zeros([num_grooves,int(T*fs)])

    for i in xrange(0,num_grooves): ##skip the first groove as this typically contains the needle drop
        start = int(i*T*fs)
        end = int(i*T*fs + T*fs)
        print 'start: ', start
        print 'end: ', end
        block_L = L[start:end]
        block_R = R[start:end]
        time_array[i] = np.linspace(start/fs,end/fs,T*fs)##calculates a time array in order to plot the waveform of the audio file  
        try:      
            groove_array[i] = np.array([block_L,block_R])
        except: 
            continue

    return groove_array, time_array


if __name__ == '__main__':

    audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'

    # file_dir = '1015_18_LiteToneTest/'
    # file_name = '5.3declicked.wav'
    # file_name = '5.1.wav'

    # # file_dir = '1101_18_LiteTone45rpm/'
    # # file_name = '45-5.2.wav'
    # # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    # # file_name = 'Silence 0_05.wav'
    # # file_name = 'Dec14-TestA1n-tone1000.wav'

    # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    # file_name = 'Silence 7_45.wav'
    # # file_dir = '00_digital_files/'
    # # file_name = 'tone100.wav'
    # # file_dir = '1108_18_LiteToneMusicSample/'
    # # file_name = '1-intro.wav'
    file_dir = '1129-18_KingGizzard/'
    file_name = 'A33.wav'

    file_path = audio_dir + file_dir + file_name


    SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
    soundfile = sf.SoundFile(file_path, 'r')
    # data_sf, fs_sf = sf.read(input, always_2d=True)
    bit_depth = SOUNDFILE_types[soundfile.subtype]
    print 'The bit depth of the audio is: ', bit_depth

    fs, data = wavfile.read(file_path)
    audio_file = {}
    audio_file['data'] = data 
    audio_file['fs'] = fs 
    audio_file['bit_depth'] = bit_depth

    groove_array, time_array = groove_map(audio_file)
    print 'np.shape(groove_array)', np.shape(groove_array)
    fig = plt.figure(1)

    for i in xrange(6,len(groove_array)-2): 
        # print time_array[i]
        plt.subplot(211)
        plt.plot(time_array[0],groove_array[i][0],label = 'groove %i' % i)
        plt.subplot(212)
        freq, spec = fft_array(groove_array[i][0])
        plt.plot(freq,10.0*np.log10(spec/2**float(16)))

    plt.subplot(211)
    # plt.grid('on')
    plt.ylabel('Amplitude')
    plt.xlabel('Time (s)')
    plt.grid(which='both')
    plt.legend()
    plt.subplot(212)
    plt.ylabel('Power (dB)')
    plt.xlabel('Frequency (Hz)')
    plt.xscale('log')
    plt.grid(which='both')


    plt.legend()
    plt.show()

