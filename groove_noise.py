from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig

import matplotlib.pyplot as plt
import numpy as np

audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
file_dir = '1015_18_LiteToneTest/'
file_name = '5.3declicked.wav'
file_dir = '1101_18_LiteTone45rpm/'
file_name = '45-5.2.wav'

# file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
# file_name = 'Silence 7_45.wav'


plt.figure(1)
plt.subplot(211)
plt.subplot(212)

audio_file = audio_dir + file_dir + file_name
fs, data = wavfile.read(audio_file)
L = data.T[0]
R = data.T[1]

time = np.arange(0.0,len(data),1) ##calculates a time array in order to plot the waveform of the audio file
time = time/float(fs)  

T = 1.8 #seconds, the period of rotation of a record
T = 1.0 + 1.0/3.0
T_s = T*fs #period as a number of samples
fs = 44100


num_grooves = int(time[len(time)-1]/T)
print 'length of files: ', len(time), len(data)
print 'max_time: ',  time[len(time)-1]
print 'num_grooves: ', num_grooves
# num_grooves = 3

for i in xrange(num_grooves): 
    # if i == 0: 
        # continue
    start = int(i*T*fs)
    end = int(i*T*fs + T*fs)
    block = L[start:end]

    nyq = 0.5 * fs
    normal_cutoff = 100 / nyq
    b, a = signal.butter(1, 1000.0/nyq, btype='low', analog=False)
    block = signal.lfilter(b, a, block)


    fft_block = fft(block)
    fft_block = abs(fft_block[0:(len(fft_block)/2)]) ##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
    freq = np.arange(len(fft_block))*float(len(block))/float(fs)

    plt.subplot(211)    
    plt.plot(block, label = 'groove %i' % i)

    plt.subplot(212)    
    plt.plot(freq, 20*np.log10(fft_block/(2**16)), label = 'groove %i' % i)
    


file_dir = '1101_18_LiteToneMiddle/'
file_name = '5.3m.wav'

audio_file = audio_dir + file_dir + file_name
fs, data = wavfile.read(audio_file)
L = data.T[0]
R = data.T[1]

time = np.arange(0.0,len(data),1) ##calculates a time array in order to plot the waveform of the audio file
time = time/float(fs)  

T = 1.8 #seconds, the period of rotation of a record
T_s = T*fs #period as a number of samples

num_grooves = int(time[len(time)-1]/T)
print 'length of files: ', len(time), len(data)
print 'max_time: ',  time[len(time)-1]
print 'num_grooves: ', num_grooves


# for i in xrange(num_grooves): 
#     start = int(i*T*fs)
#     end = int(i*T*fs + T*fs)
#     block = L[start:end]

#     fft_block = fft(block)
#     fft_block = abs(fft_block[0:(len(fft_block)/2)]) ##originally was 0:(len(data_fft_a)/2) - 1, no idea why I had the minus one
#     freq = np.arange(len(fft_block))*float(len(block))/float(fs)

#     plt.subplot(211)    
#     plt.plot(block,'b--', label = 'groove %im' % i)

#     plt.subplot(212)    
#     plt.plot(freq, 20*np.log10(fft_block/(2**16)),'b--', label = 'groove %im' % i)
    

plt.subplot(211)
plt.legend()
plt.subplot(212)
plt.legend()
plt.xscale('log')


plt.show()