from scipy import signal
from scipy.io import wavfile
from scipy.fftpack import fft
import matplotlib.pyplot as plt
import numpy as np

audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/'

file_dir = '1015_18_LiteToneTest/'
file = '5.1.wav'

file_dir = '00_digital_files/'
file = 'sweep16khz.wav'

fs, data_a = wavfile.read(audio_dir+file_dir+file)

nfft = 2048 ##this is the number of points in the fft 

def spectrum(data_a, start = 0 , end = False, fs = 44100, nfft = 2048): 
    nbins = nfft/2 + 1 ##number of frequency bins 

    if end == False: ##if no ending point was specified do the whole length of data_a
        end = len(data_a)

    time_a = np.linspace(start,end,nfft)*float(nfft/fs)

    FFT_l = []
    print 'start: ', start
    print 'end: ', end
    print 'steps: ', (end-start)/nfft
    print 'np.shape(data_a): ', np.shape(data_a)
    # FFT_a = np.ze<os(len(time_a),len(freq_a),2)
    # fft_ai = None
    for i in xrange(start,end,nfft): 
        # print 'i: ', i
        if i-nfft < 0:
            # try: ##if start is somewhere in the middle of the audio file, this will just take the points preceding the start in the fft 
            #     print 'start try'
            #     print i - nfft/2
            #     fft_ai = (fft(data_a.take(indices=xrange(i - nfft/2, i + nfft/2 + 1), axis = 0)))
            #     fft_ai = fft_ai.take(indices=xrange(0,nfft/2 + 1),axis=0)
            # except: ##if there aren't enough preceding points, duplicate the array OR: should I pad with zeros? 
            print 'start except'
            x_i = data_a.take(indices=xrange(i, i + nfft/2), axis = 0)
            data_ai = np.concatenate((x_i,x_i),axis = 0)
            print np.shape(data_ai)
            fft_ai = (fft(data_ai))
            fft_ai = fft_ai.take(indices=xrange(0,nfft/2 + 1),axis=0)

        if i-nfft > start and i+nfft < end: 
            # print 'middle'
            fft_ai = (fft(data_a.take(indices=xrange(i - nfft/2, i + nfft/2), axis = 0)))
            print 'np.shape(fft_ai) : ', np.shape(fft_ai)
            print 'np.shape(data_a.take(indices=xrange(i - nfft/2, i + nfft/2), axis = 0))): ', np.shape(data_a.take(indices=xrange(i - nfft/2, i + nfft/2), axis = 0))
            fft_ai = fft_ai.take(indices=xrange(0,nfft/2+1),axis=0)
            print 'np.shape(fft_ai) after take : ', np.shape(fft_ai)

        if i+nfft > end: 
            try: #if start is somewhere in the middle of the audio file, this will just take the points preceding the start in the fft 
                print 'end try'
                fft_ai = (fft(data_ai.take(indices=xrange(i - n/2, i + nfft/2), axis = 0)))
                fft_ai = fft_ai.take(indices=xrange(0,nfft/2+1),axis=0)
            except: ##if there aren't enough preceding points, duplicate the array OR: should I pad with zeros? 
                print 'end except'
                x_i = data_a.take(indices=xrange(i-nfft/2, i), axis = 0)
                data_ai = np.concatenate((x_i,x_i),axis = 0)
                fft_ai = (fft(x_i))
                fft_ai = fft_ai.take(indices=xrange(0,nfft/2),axis=0)

        # print 'IN LOOP: ', np.shape(fft_ai)
        # fft_li = fft_ai.tolist()
        # print 'len(fft_li): ', len(fft_li)
        # print 'len(fft_li[i`]): ', len(fft_li[0])
        if i == start: 
            continue
        # FFT_l.append(fft_li)
        FFT_l.append(fft_ai.T)
        # print len(FFT_l)
    # print 'len(FFT_l): ', len(FFT_l)
    # print np.shape(fft_ai)
    freq_a = np.arange(0,len(FFT_l[0][0]))*fs/nfft
    time_a = np.arange(0,len(FFT_l))*nfft/fs
    return FFT_l, freq_a, time_a


FFT_l, freq_a, time_a = spectrum(data_a)


for i in xrange(len(FFT_l)): 
    if np.shape(FFT_l[i]) != (2, 1025):
        print 'FFT_l shape: ', np.shape(FFT_l[i])
        print 'i: ', i
FFT_a = np.dstack(FFT_l)


print np.shape(data_a)
print np.shape(FFT_a)
print np.shape(FFT_a.T)

print np.shape(time_a)
print np.shape(freq_a)


plt.pcolormesh(time_a, freq_a, np.real(FFT_a[0]))
plt.colorbar()
# plt.ylim(0,22050)
plt.show()