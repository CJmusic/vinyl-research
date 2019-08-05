import matplotlib.pyplot as plt
import numpy as np
from octave_smoothing import octave_smooth
from groove_noise import groove_map
from audio_load import load_audio
from fft_audio import fft

print 'hello'
fs = 44100

a_33rpm_loc = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.1.wav'
a_45rpm_loc = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone45rpm/45-5.1.wav'
a_78rpm_loc = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone78rpm/78-5.1.wav'

a_33rpm = load_audio(a_33rpm_loc)
a_45rpm = load_audio(a_45rpm_loc)
a_78rpm = load_audio(a_78rpm_loc)

a_33_groove, a_33_time = groove_map(a_33rpm, 0, T = 1.8)
a_45_groove, a_45_time = groove_map(a_45rpm, 0, T = 1.3333)
a_78_groove, a_78_time = groove_map(a_78rpm, 0, T = 0.769)

# freq_45, fft_45 = fft(a_45_groove[i][0])
# freq_78, fft_78 = fft(a_78_groove[i][0])


plt.figure()
for i in xrange(len(a_33_groove)):
    freq_33, fft_33 = fft(a_33_groove[i][0])
    fft_33 = fft_block_L[0:(len(a_33_groove[i][0])/2 + 1)]
    freq_33 = np.linspace(0.0,fs/2.0,len(a_33_groove[i][0])/2 + 1)
    plt.plot(fft(a_33_groove[i][0]),'blue')

for i in xrange(len(a_45_groove)): 
    freq_33, fft_33 = fft(a_33_groove[i][0])
    fft_33 = fft_block_L[0:(len(a_33_groove[i][0])/2 + 1)]
    freq_33 = np.linspace(0.0,fs/2.0,len(a_33_groove[i][0])/2 + 1)
    plt.plot(fft(a_45_groove[i][0]),'orange')
for i in xrange(len(a_78_groove)): 
    freq_33, fft_33 = fft(a_33_groove[i][0])
    fft_33 = fft_block_L[0:(len(a_33_groove[i][0])/2 + 1)]
    freq_33 = np.linspace(0.0,fs/2.0,len(a_33_groove[i][0])/2 + 1)
    plt.plot(fft(a_33_groove[i][0]),'blue')
    plt.plot(fft(a_78_groove[i][0]),'green')

plt.xlim(0,22500)
plt.yscale('log')
plt.show()