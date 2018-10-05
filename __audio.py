from scipy.fftpack import fft

import numpy as np 

import soundfile as sf
import wave
from scipy import signal

class SOUNDFILE: ##the class that represents audio
    def __init__(self, input, dft_npoints = 2.0**16): 
        self.SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
        self.dft_npoints = dft_npoints
        self.fs = 0.0

        self.data_a = np.empty(0)
        self.dft_a = np.empty(0)
        self.time_a = np.empty(0)
        self.freq_a = np.empty(0)

        self.bit_depth_i = 0 

        self.__file_read(input)
        self.__calc_time()
        self.__calc_freqbins()
        self.__dft_init()

        self.leadin_start = False
        self.signal_start = False

        self.scipy_freq = False 
        self.scipy_time = False 
        self.scipy_spec = False

        self.file_name = input
        self.file_info(input)


    def __file_read(self,input):
        ##The following processes the wav file and imports it into python using the PySoundFile package
        self.Input = sf.SoundFile(input, 'r')
        self.data_a, self.fs = sf.read(input, always_2d=True)
        self.bit_depth_i = self.SOUNDFILE_types[self.Input.subtype]
        return
        
        ##This fills in all the values of a SOUNDFILE class for use in the analysis 
        'Keeping the below 2 lines for reference'
        # self.data_L = self.data.T[0] 
        # self.data_R = self.data.T[1]

    def __calc_time(self):
        self.time_a = np.linspace(1,float(len(self.data_a))/self.fs, len(self.data_a))
        return


    def __calc_freqbins(self):
        k = np.arange(self.dft_npoints/2-1)
        T = len(self.data_a)/self.fs
        self.freq_a = k/T
        return

    def __dft_init(self):
        ##initializes an empty dft_a array to be filled with data as its needed

        # dft_length = int((start_t - end_t)/dft_npoints)
        # dft_cutoff = ((start_t - end_t)/dft_npoints) - dft_length)*dft_npoints 

        # self.dft_a = np.zeros(dft_length) 
        # self.dft_time_a = np.zeros(dft_length) ##this creates an array of time stamps for each index of the dft_a array 
        # self.dft_freq_a = np.zeros(dft_length)

        if len(self.data_a)%self.dft_npoints > 0: 
            padding = self.dft_npoints - (len(self.data_a)%self.dft_npoints)
            self.dft_a = np.zeros(int((len(self.data_a) + padding)/self.dft_npoints))
            self.dft_time_a = np.zeros(len(self.dft_a)) ##this creates an array of time stamps for each index of the dft_a array 
            self.dft_freq_a = np.zeros(len(self.dft_a))

            return
        else:
            self.dft_a = np.zeros(len(self.data_a)/self.dft_npoints)
            self.dft_time_a = np.zeros(len(self.dft_a)) ##this creates an array of time stamps for each index of the dft_a array 
            self.dft_freq_a = np.zeros(len(self.dft_a))
            return

    def detect_signal(self): 
        return

    def detect_leadin(self, tstart = 10.0):
        leadin_detect = False
        signal_detect = False
        # print 'DATA ARRAY: ', self.data_a
        index = 0
        for sample in self.data_a[:int(tstart)*44100]:
            if (np.average(sample)) > 10.0**(-30.0/20) and leadin_detect == False: 
                self.leadin_start = index
                leadin_detect = True
            if (np.average(sample)) > 10.0**(-10.0/20.0) and signal_detect == False: 
                self.signal_start = index
                signal_detect = True
                return self.leadin_start, self.signal_start
            index += 1
        
        return self.leadin_start, self.signal_start

    def leadin_noise(self):
        if self.leadin_start == False or self.signal_start == False: 
            self.detect_leadin()
        self.leadin_noise = (RMS_level(self, self.leadin_start, self.signal_start-self.leadin_start))
        return self.leadin_noise
      
        

    def dft_audio(self, start, end):
        ''' the dft_a array has the following structure: 
        dft_a[ time[LEFT[freq[amplitude]], RIGHT[freq[amplitude]]],
               time[LEFT[freq[amplitude]], RIGHT[freq[amplitude]]],
                                        ...
        ]
        -the first index is sorted by time, matching this index with the dft_time_a array will give the timestamp that this 
        dft point occurs at
        -the next index is for stereo, 0 for L 1 for R
        -next is a list of amplitude values sorted by frequency, each index corresponds to a frequency values found by matching 
        the index with the dft_freq_a array. These amplitudes are not in dB
        '''
        ##I want to make the array from the start, but only fill in values that I want to, so I have the full array, but some spots that aren't being analyzed are empty.
        
        ###HANDLED BY DFT INIT####
        # dft_length = int((start_t - end_t)/dft_npoints)
        # dft_cutoff = ((start_t - end_t)/dft_npoints) - dft_length)*dft_npoints 

        # self.dft_a = np.zeros(dft_length) 
        # self.dft_time_a = np.zeros(dft_length) ##this creates an array of time stamps for each index of the dft_a array 
        # self.dft_freq_a = np.zeros(dft_length)
        #####
        if start%self.dft_npoints < self.dft_npoints/2: 
            start_dft = int(start/self.dft_npoints)
        else: 
            start_dft = int(start/self.dft_npoints) + 1
        
        if end%self.dft_npoints < self.dft_npoints/2: 
            end_dft = int(end/self.dft_npoints)
        else: 
            end_dft = int(end/self.dft_npoints) + 1
        

        for dft_step in xrange(start_dft,end_dft):
            self.dft_a[dft_step], self.dft_freq_a[dft_step] = fft_audio(dft_step*self.dft_npoints,self.dft_npoints)
            self.dft_time_a[dft_step] = dft_step*self.dft_npoints/self.fs

        # if self.dft_a[int(start/self.dft_npoints)] != 0:
        #     ##This function will take the Fast Fourier Transform of a self.data and time array 
        #     data_fft = fft(self.data[start:start+ Npoints]) # calculate fourier transform (complex numbers list)
        #     data_fft_Re = len(self.data_fft)/2  # you only need half of the fft list (real signal symmetry)
        #     self.dft_a[int(start/self.dft_npoints)] = abs(self.data_fft[0:(self.data_fft_Re-1)])

        #     k = np.arange(data_fft_Re - 1)   
        #     T = float(len(data))/float((self.fs)) 
        #     freq = k/T

        #     self.dft_time_a[] = index*self.dft_npoints/self.fs
        #     self.dft_freq_a[] = freq

        #     return self.dft_a[int(start/self.dft_npoints)]
        # else:
        #     return self.dft_a[int(start/self.dft_npoints)]

    def fft_audio(self, start, Npoints):
        ##This function will take the Fast Fourier Transform of a self.data and time array 
        data = self.data[start:start+ Npoints]
        data_fft = fft(data) # calculate fourier transform (complex numbers list)
        data_fft_Re = len(data_fft)/2  # you only need half of the fft list (real signal symmetry)
        data_fft = abs(data_fft[0:(data_fft_Re-1)])
        # data_fft_DBFS = 20.0*np.log10(4.5933528*abs(data_fft)/len(data_fft))
        # data_fft_DBFS = 20.0*np.log10(abs(data_fft)/2**16)#/len(data_fft))
        
        k = np.arange(data_fft_Re - 1)   
        T = float(len(data))/float((self.fs)) 
        freq = k/T
        # freq = freq
        return data_fft, freq

    def scipy_spectograph(self,start, Npoints):
        self.scipy_freq, self.scipy_time, self.scipy_spec = signal.spectogram(self.data_a[start:start+Npoints] , self.fs)  
        return self.scipy_freq, self.scipy_time, self.scipy_spec
        
    def file_info(self,input):
        print 'Number of Channels: ', self.Input.channels
        print 'The sampling rate is: ', self.fs
        print 'The type of audio is: ', self.Input.subtype       
        print 'The bit depth is : ', self.bit_depth_i
        print 'The file is ', self.time_a[len(self.time_a)-1], ' seconds long'
        print 'There are ', len(self.data_a), ' samples.'



def loadfile(input):
    file = SOUNDFILE(input)
    return file


'''
Below are the various audio tests that can be run on a SOUNDFILE object.


'''


def detect_signal(SF, start = 0):
    return sample_start_i

def signal_to_noise(SF, start, Npoints):
    return snr_f

def harmonic_distortion(SF,start,Npoints):
    return snr_f

def logsweep_spectrum(SF,start,Npoints):
    return sweep_spectrum_A

def noise_spectrum(SF, start, Npoints):
    return noise_spectrum_A

def stereo_seperation(SF, start, Npoints):
    return stereo_seperation_r #r as in ratio

def count_clicks(SF, start, Npoints):
    return num_clicks

def RMS_level(SF, start, Npoints):
    RMS_level_a = np.sqrt(np.mean(SF.data_a[start:start+Npoints]**2))
    return RMS_level_a

### LEGACY CODE####
#     def chop_audio(self, start, Npoints):
#     data = self.data[start : start + Npoints]
#     data_DBFS = self.data_DBFS[start:start + Npoints]
#     time = self.time[start:start + Npoints]
#     return data, data_DBFS, time

# def flattop_window(self, data):#, self.data):
#     self.window = np.array([sig.flattop(len(data))])
#     self.data_windowed = (self.window[0])*(data)  
#     return self.data_windowed

# def A_weighting(self,data):
#     # Fs = 48000;
#     #Analog A-weighting filter according to IEC/CD 1672.
#     f1 = 20.598997
#     f2 = 107.65265
#     f3 = 737.86223
#     f4 = 12194.217
#     A1000 = 1.9997
            
#     NUMs = [(2*np.pi * f4)**2 * (10**(A1000/20)), 0, 0, 0, 0]
#     DENs = np.convolve([1, 4*np.pi * f4, (2*np.pi * f4)**2],
#                     [1, 4*np.pi * f1, (2*np.pi * f1)**2], mode='full')
#     DENs = np.convolve(np.convolve(DENs, [1, 2*np.pi * f3], mode='full'),[1, 2*np.pi * f2], mode='full')

#     #Bilinear transformation of analog design to get the digital filter.
#     [b,a] = sig.bilinear(NUMs,DENs,self.fs)
#     return sig.lfilter(b,a,data)

# def unweighted_rumble(self,data):
#     b, a = sig.butter(2, 400.0 / (0.5 * self.fs), btype='low')
#     y = sig.lfilter(b, a, data)        
#     b, a = sig.butter(1, 9.0 / (0.5 * self.fs), btype='high')        
#     y = sig.lfilter(b, a, y)
#     return y

# def weighted_rumble(self,data):
#     b, a = sig.butter(2, 500.0 / (0.5 * self.fs), btype='low')
#     y = sig.lfilter(b, a, data)        
#     b, a = sig.butter(2, 200.0 / (0.5 * self.fs), btype='high')        
#     y = sig.lfilter(b, a, y)
#     return y

# def inv_riaa(self,data): 
#     b, a = sig.butter(1, 1000.0 / (0.5 * self.fs), btype='high')
#     y = 0.01*sig.lfilter(b, a, data)    
#     # y = data
#     b, a = sig.butter(1, [20000.0 / (0.5 * self.fs)], btype='low')        
#     y = 1000.0*sig.lfilter(b, a, y)
#     return y


# def RMS(self,data):
#     data_RMS = np.sqrt(np.mean(data*data))
#     return data_RMS

# def rm_pops(self,data):
#     print 'RM POPS Function: '
#     data_RMS = self.RMS(data)
#     print 'data_RMS: ', data_RMS
#     print 'max(data): ', max(data)
#     i = 0
#     while max(data) > 2.0*data_RMS: 
#         # print 'max: ', max(data)
#         # print 'index: ', data[np.where(data == max(data))]
#         data[np.where(data == max(data))] = data_RMS
#         i += 1
#     print 'Pops and clicks removed: ', i
#     return data

# def locate_click(self, region, start, Npoints): 
#     data_p = region[start:start+Npoints]
#     n = np.where(region == max(data_p))
    
#     print('Locate click')
#     print(n[0][0])
#     print('Value of Click')
#     print(region[n[0][0]])
#     print ('START: ', n[0][0] - 2**16/2)
#     print ('END: ', n[0][0] + 2**16/2)
#     print ('length: ', len(self.data))
#     # self.chop_audio(n[0][0] - 2**16/2,2**16)#,2**16)
#     self.chop_audio(n[0][0] ,2**16)#,2**16)
### END ###