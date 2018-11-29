
import soundfile as sf 
from scipy.io import wavfile


# file_path = audio_dir + file_dir + file_name

def load_audio(file_path): 
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
    return audio_file
