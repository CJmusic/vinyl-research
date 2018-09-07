from scipy.fftpack import fft

import numpy as np 

import soundfile as sf
import wave

class SOUNDFILE: ##the class that represents audio
    def __init__(self, input, dft_npoints = 2.0**16): 
        self.SOUNDFILE_types = {'PCM_16': 16, 'PCM_24': 24, 'PCM_32': 32, 'FLOAT': 32}
        self.dft_npoints = dft_npoints

        self.data_a = np.empty(0)
        # self.data_DBFS_a = np.empty(0)
        self.dft_a = np.empty(0)
        self.time_a = np.empty(0)
        self.freq_a = np.empty(0)

        self.bit_depth_i = 0 

        self.file_read(input)
        self.file_info(input)
        self.calc_time()
        self.calc_freqbins()
        self.__dft_init()

    def file_read(self,input):
        ##The following processes the wav file and imports it into python using the PySoundFile package
        self.Input = sf.SoundFile(input, 'r')
        self.data, self.fs = sf.read(input)
        
        self.data_DBFS_a = 20*np.log10(abs(np.asarray(self.data_normalized)))

        ##This fills in all the values of a SOUNDFILE class for use in the analysis 
        'Keeping the below 2 lines for reference'
        # self.data_L = self.data.T[0] 
        # self.data_R = self.data.T[1]
        self.bit_depth_i = self.SOUNDFILE_types[self.Input.subtype]

    def calc_time(self):
        self.time_a = np.linspace(1,float(len(self.data))/self.fs, len(self.data))
        return


    def calc_freqbins(self):
        k = np.arange(self.dft_npoints/2-1)
        T = len(self.data_a)/self.fs
        self.freq_a = k/T
        return

    def __dft_init(self):
        self.dft


    def dft_audio(self, start, Npoints):
        ##This function will take the Fast Fourier Transform of a self.data and time array 
        self.data_fft = fft(self.data[start:start+ Npoints]) # calculate fourier transform (complex numbers list)
        self.data_fft_Re = len(self.data_fft)/2  # you only need half of the fft list (real signal symmetry)
        self.data_fft = abs(self.data_fft[0:(self.data_fft_Re-1)])
        # data_fft_DBFS = 20.0*np.log10(4.5933528*abs(data_fft)/len(data_fft))
        self.data_fft_DBFS = 20.0*np.log10(abs(self.data_fft)/2**16)#/len(data_fft))
        
        k = np.arange(self.data_fft_Re - 1)   
        T = float(len(self.data))/float((self.fs)) 
        freq = k/T
        self.freq = freq
        return self.data_fft, self.data_fft_DBFS, self.freq

    def fft_audio(self, start, Npoints):
        ##This function will take the Fast Fourier Transform of a self.data and time array 
        self.data_fft = fft(self.data[start:start+ Npoints]) # calculate fourier transform (complex numbers list)
        self.data_fft_Re = len(self.data_fft)/2  # you only need half of the fft list (real signal symmetry)
        self.data_fft = abs(self.data_fft[0:(self.data_fft_Re-1)])
        # data_fft_DBFS = 20.0*np.log10(4.5933528*abs(data_fft)/len(data_fft))
        self.data_fft_DBFS = 20.0*np.log10(abs(self.data_fft)/2**16)#/len(data_fft))
        
        k = np.arange(self.data_fft_Re - 1)   
        T = float(len(self.data))/float((self.fs)) 
        freq = k/T
        self.freq = freq
        return self.data_fft, self.data_fft_DBFS, self.freq


        
    def file_info(self,input):
        print 'Number of Channels: ', self.Input.channels
        print 'The sampling rate is: ', self.fs
        print 'The type of audio is: ', self.Input.subtype       
        print 'The bit depth is : ', self.bit_depth
        print 'There are ', len(self.data), ' samples.'


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

def count_clicks(SF, start, Npoints)
    return num_clicks


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