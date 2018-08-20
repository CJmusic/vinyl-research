import numpy as np 

import soundfile as sf

class SOUNDFILE: 
    def __init__(self, Input, stereo = True): 
        self.SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}

        self.data, self.fs = sf.read(INPUT)
        self.bit_depth = self.SOUNDFILE_types[self.Input.subtype]

    def file_read():
        self.Input = sf.SoundFile(INPUT, 'r')

    def file_info():
        print 'Number of Channels: ', self.Input.channels
        # print 'The sampling rate is: ', self.fs
        # print 'The type of audio is: ', self.Input.subtype       
        # print 'The bit depth is : ', self.bit_depth
        # print 'There are ', len(self.data), ' samples.'
