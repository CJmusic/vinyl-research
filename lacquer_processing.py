import matplotlib.pyplot as plt 
import numpy as np 
import soundfile as sf
from scipy.io import wavfile
from scipy.fftpack import fft
import soundfile as sf

from octave_smoothing import octave_smooth


from fft_audio import fft_audio
###In this file, you feed in a stereo wave file to a digital_signal() and it will determine the signal used 
###and produce a wave to cancel out the signal and isolate just the noise 

def rm_signal(audio_file): 
    data = audio_file['data']
    fs = audio_file['fs'] 
    bit_depth = audio_file['bit_depth']

    data_fft, freq, phi = fft_audio(audio_file, phase = True)

    print 'np.shape(data): ', np.shape(data)
    print 'phi: ', phi
    max_freq = data_fft.argmax(axis=0)
    max_bin = freq[max_freq]
    max_amp = 0.5#data_fft[max_freq]

    time = np.linspace(0.0,float(len(data))/fs,len(data))
    wt = 2.0*np.pi*time*fs*max_bin
    print 'shape(wt): ', np.shape(wt)
    wt += phi
    wave = max_amp*np.sin( wt )
    plt.subplot(211)
    plt.plot(time, data)
    plt.subplot(212)
    plt.plot(time, wave)
    plt.show()

    return 


if __name__ == '__main__':

    audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'
    # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    # file_name = 'Dec14-TestA1n-tone1000.wav'

    # file_dir = '4A-RecordedSoundFiles/Dec14-TestAnormalized/Dec14-A1n_files/'
    file_dir = '00_digital_files/'
    file_name = 'tone1000.wav'
    file_path = audio_dir + file_dir + file_name


    SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
    soundfile = sf.SoundFile(file_path, 'r')
    # data_sf, fs_sf = sf.read(input, always_2d=True)
    bit_depth = SOUNDFILE_types[soundfile.subtype]
    print 'The bit depth of the audio is: ', bit_depth
    fs, data = wavfile.read(file_path)
    audio_file = {}
    audio_file['data'] = data.T[0][0:2**16]
    audio_file['fs'] = fs 
    audio_file['bit_depth'] = bit_depth

    rm_signal(audio_file)
    


