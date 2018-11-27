from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
# import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

import soundfile as sf

from octave_smoothing import octave_smooth
from groove_noise_lacquer import groove_map



audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
file_dir = '1015_18_LiteToneTest/'
file_name = '5.1.wav'
plt.figure(1)

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

groove_map(audio_file, 1, filter = False)

file_dir = '1101_18_LiteTone45rpm/'
file_name = '45-5.1.wav'
plt.figure(2)

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
groove_map(audio_file, 1, T=0.45, filter = False)

file_dir = '1101_18_LiteTone78rpm/'
file_name = '78-5.1.wav'
plt.figure(3)

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
groove_map(audio_file, 1, T=1.3 , filter = False)

plt.show()