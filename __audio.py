import numpy as np 

import soundfile as sf

class SOUNDFILE: ##the class that represents audio
    def __init__(self, input): 
        self.SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
        self.file_read(input)
        self.file_info()


    def file_read(self,input):
        self.data, self.fs = sf.read(input)
        self.bit_depth = self.SOUNDFILE_types[self.Input.subtype]



    def file_info():
        print 'Number of Channels: ', self.Input.channels
        print 'The sampling rate is: ', self.fs
        # print 'The type of audio is: ', self.Input.subtype       
        # print 'The bit depth is : ', self.bit_depth
        # print 'There are ', len(self.data), ' samples.'


def loadfile(input):
    file = SOUNDFILE(input)