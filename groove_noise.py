# from scipy import signal
# from scipy.io import wavfile
# from scipy.fftpack import fft
# # import scipy.signal as sig

# import matplotlib.pyplot as plt
# import numpy as np

# audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
# file_dir = '1015_18_LiteToneTest/'
# file_name = '5.3declicked.wav'
# file_dir = '1101_18_LiteTone45rpm/'
# file_name = '45-5.2.wav'

# # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
# # file_name = 'Silence 7_45.wav'


# plt.figure(1)
# plt.subplot(211)
# plt.subplot(212)

# audio_file = audio_dir + file_dir + file_name
# fs, data = wavfile.read(audio_file)
# L = data.T[0]
# R = data.T[1]

# time = np.arange(0.0,len(data),1) ##calculates a time array in order to plot the waveform of the audio file
# time = time/float(fs)  

# T = 1.8 #seconds, the period of rotation of a record
# T = 1.0 + 1.0/3.0
# T_s = T*fs #period as a number of samples
# fs = 44100


# num_grooves = int(time[len(time)-1]/T)
# print 'length of files: ', len(time), len(data)
# print 'max_time: ',  time[len(time)-1]
# print 'num_grooves: ', num_grooves
# # num_grooves = 3

# for i in xrange(num_grooves): 
#     # if i == 0: 
#         # continue
#     start = int(i*T*fs)
#     end = int(i*T*fs + T*fs)
#     block = L[start:end]

#     nyq = 0.5 * fs
#     normal_cutoff = 100 / nyq
#     b, a = signal.butter(1, 1000.0/nyq, btype='low', analog=False)
#     block = signal.lfilter(b, a, block)


#     fft_block = fft(block)
#     fft_block = abs(fft_block[0:(len(fft_block)/2)]) ##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
#     freq = np.arange(len(fft_block))*float(len(block))/float(fs)

#     plt.subplot(211)    
#     plt.plot(block, label = 'groove %i' % i)

#     plt.subplot(212)    
#     plt.plot(freq, 20*np.log10(fft_block/(2**16)), label = 'groove %i' % i)
    


# file_dir = '1101_18_LiteToneMiddle/'
# file_name = '5.3m.wav'

# audio_file = audio_dir + file_dir + file_name
# fs, data = wavfile.read(audio_file)
# L = data.T[0]
# R = data.T[1]

# time = np.arange(0.0,len(data),1) ##calculates a time array in order to plot the waveform of the audio file
# time = time/float(fs)  

# T = 1.8 #seconds, the period of rotation of a record
# T_s = T*fs #period as a number of samples

# num_grooves = int(time[len(time)-1]/T)
# print 'length of files: ', len(time), len(data)
# print 'max_time: ',  time[len(time)-1]
# print 'num_grooves: ', num_grooves


# # for i in xrange(num_grooves): 
# #     start = int(i*T*fs)
# #     end = int(i*T*fs + T*fs)
# #     block = L[start:end]

# #     fft_block = fft(block)
# #     fft_block = abs(fft_block[0:(len(fft_block)/2)]) ##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
# #     freq = np.arange(len(fft_block))*float(len(block))/float(fs)

# #     plt.subplot(211)    
# #     plt.plot(block,'b--', label = 'groove %im' % i)

# #     plt.subplot(212)    
# #     plt.plot(freq, 20*np.log10(fft_block/(2**16)),'b--', label = 'groove %im' % i)
    

# plt.subplot(211)
# plt.legend()
# plt.subplot(212)
# plt.legend()
# plt.xscale('log')


# plt.show()

"~~~~~~~~"
from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

import soundfile as sf

from octave_smoothing import octave_smooth


def groove_map(audio_file, start_groove,T=1.8 ,filter = False): 
    data = audio_file['data'] 
    fs = audio_file['fs'] 
    bit_depth = audio_file['bit_depth']

    L = data.T[0]#/2**bit_depth
    R = data.T[1]#/2**bit_depth

    time = np.arange(0.0,len(data),1)/fs ##calculates a time array in order to plot the waveform of the audio file
    # time = time/float(fs)  

    T_s = int(T*fs) #period as a number of samples

    num_grooves = int(time[len(time)-1]/T)
    print 'length of files: ', len(time), len(data)
    print 'max_time: ',  time[len(time)-1]
    print 'num_grooves: ', num_grooves
    # num_grooves = 2
    # num_grooves = 16
    groove_array = np.zeros([num_grooves,2,int(T*fs)])#,2,T*fs)
    for i in xrange(0,num_grooves): ##skip the first groove as this typically contains the needle drop
        start = int(i*T*fs)
        end = int(i*T*fs + T*fs)
        print 'start: ', start
        print 'end: ', end
        block_L = L[start:end]
        block_R = R[start:end]
        time = np.linspace(0.0,float(T*fs),T*fs)/float(fs) ##calculates a time array in order to plot the waveform of the audio file


        if filter == True:
            nyq = 0.5 * fs
            b, a = signal.butter(1, 20.0/nyq, btype='high', analog=False)
            block_L = signal.lfilter(b, a, block_L)
            block_R = signal.lfilter(b, a, block_R)


        # fft_block_L = fft(block_L)
        # fft_block_L = fft_block_L[0:(len(block_L)/2 + 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
        # # freq = np.linspace(0.0,fs/2.0,len(block_L)/2 + 1)
        # fft_block_R = fft(block_R)
        # fft_block_R = fft_block_R[0:(len(block_R)/2 + 1)]##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
        # freq = np.linspace(0.0,fs/2.0,len(block_R)/2 + 1)
        groove_array[i] = np.array([block_L,block_R])
    return groove_array, time



if __name__ == '__main__':

    audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
    file_dir = '1015_18_LiteToneTest/'
    file_name = '5.3declicked.wav'
    file_name = '5.1.wav'

    # file_dir = '1101_18_LiteTone45rpm/'
    # file_name = '45-5.2.wav'
    # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    # file_name = 'Silence 0_05.wav'
    # file_name = 'Dec14-TestA1n-tone1000.wav'

    # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    # file_name = 'Silence 7_45.wav'
    # file_dir = '00_digital_files/'
    # file_name = 'tone100.wav'
    # file_dir = '1108_18_LiteToneMusicSample/'
    # file_name = '1-intro.wav'

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

    groove_array, time = groove_map(audio_file, 1, filter = False)
    print 'np.shape(groove_array)', np.shape(groove_array)
    plt.figure(1)

    for i in xrange(1,len(groove_array)): 
        RMS_level = 20.0*np.log10(np.sqrt(sum(groove_array[i][0])**2)/len(groove_array[i][0]))
        print 'average: ', RMS_level
        if RMS_level < -30.0 : 
            plt.plot(groove_array[i][0])
    plt.show()

