import matplotlib.pyplot as plt
import numpy as np
import __audio as audio
    
    
def plot_wave(file, start, Npoints):#, time, self.data):
    # print len(file.time_a)
    # print len(file.data_a)
    # print 'start: ', start 
    # print 'Npoints: ', Npoints
    x = file.time_a[start:start+Npoints]
    y = (file.data_a[start:start+Npoints])

    # y = 20*np.log10(file.data_a[start:start+Npoints])
    # print len(x), len(y)


    plt.plot(x,y, alpha = 0.7)#, linewidth = 0.7)
    plt.title(file.file_name)
    plt.ylabel('Volume (dB)')
    plt.xlabel('Time (s)')
    plt.grid(which='both')
    plt.show() 
    

def plot_fft(file, start, Npoints):
    file.data_fft, file.data_fft_DBFS, file.freq = file.fft_audio(start, Npoints)
    plt.plot(file.freq, file.data_fft_DBFS, alpha = 0.7)#, linewidth = 0.5)         
    plt.ylabel('Intensity (dB FS)')
    plt.xlabel('Freq (Hz)')
    plt.xscale('log')
    plt.grid(which='both')
    plt.show()


def scipy_plotspectograph(file):
    plt.pcolormesh(file.scipy_time, file.scipy_freq, file.scipy_spec)
    plt.imshow(file.scipy_spec)
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.show()

def plot_dft(file, start, Npoints):
    if start%file.dft_npoints < file.dft_npoints/2: 
        start_dft = int(start/file.dft_npoints)
    else: 
        start_dft = int(start/file.dft_npoints) + 1
    
    if Npoints%file.dft_npoints < file.dft_npoints/2: 
        npoints_dft = int(Npoints/file.dft_npoints)
    else: 
        npoints_dft = int(Npoints/file.dft_npoints) + 1

    print 'np.shape(file.dft_a): ', np.shape(file.dft_a)
    print 'np.shape(file.dft_a.T): ', np.shape(file.dft_a.T)
    print 'np.shape(file.dft_a.T[0])', np.shape(file.dft_a.T[0])


    dft_L = file.dft_a.T[0]
    dft_R = file.dft_a.T[1].T

    print 'dft_time shape: ', np.shape(file.dft_time_a)
    print 'dft_freq shape: ', np.shape(file.dft_freq_a)
    print 'dft_L shape: ', np.shape(dft_L)
    freq_a = file.dft_freq_a.T[0]
    print 'freq_a shape: ', np.shape(freq_a)
    plt.pcolormesh(file.dft_time_a, file.dft_freq_a.T[0], dft_L)

    plt.yscale('log')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.show()


if __name__ == '__main__': 

    file = audio.loadfile('/Users/christopherzaworski/Documents/University/Vinyl_Project/audio_files/Dec14-TestA1-normalized.wav')
    plot_wave(file,1000,int(1000))
    plot_fft(file, 1000, int(2e4))