import matplotlib.pyplot as plt
    
    
def plot_wave(file, start, Npoints):#, time, self.data):
    print len(file.time)
    print len(file.data)

    x = file.time[start:start+Npoints]
    y = file.data[start:start+Npoints]
    print len(x), len(y)


    plt.plot(x,y, alpha = 0.7)#, linewidth = 0.7)
    plt.ylabel('Intensity')
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


if __name__ == '__main__': 
    import __audio as audio

    file = audio.loadfile('/Users/christopherzaworski/Documents/University/Vinyl_Project/audio_files/Dec14-TestA1-normalized.wav')
    plot_wave(file,1000,int(1000))
    plot_fft(file, 1000, int(2e4))