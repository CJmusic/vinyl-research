import numpy as np 

import soundfile as sf
import wave

class SOUNDFILE: ##the class that represents audio
    def __init__(self, input): 
        self.SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
        
        self.data = []
        self.data_L = []
        self.data_R = []

        self.time = []
        self.freq = []
        self.data_windowed = []
        self.data_DBFS = []
        self.data_normalized = []
        self.data_fft = []

        self.file_read(input)
        self.file_info(input)


    def file_read(self,input):
        self.Input = sf.SoundFile(input, 'r')
        self.data, self.fs = sf.read(input)
        
        self.data_L = self.data.T[0] 
        self.data_R = self.data.T[1]

        self.bit_depth = self.SOUNDFILE_types[self.Input.subtype]

        
    def file_info(self,input):
        print 'Number of Channels: ', self.Input.channels
        print 'The sampling rate is: ', self.fs
        print 'The type of audio is: ', self.Input.subtype       
        print 'The bit depth is : ', self.bit_depth
        print 'There are ', len(self.data), ' samples.'


def loadfile(input):
    file = SOUNDFILE(input)