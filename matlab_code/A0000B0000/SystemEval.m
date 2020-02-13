%This file takes the noise recordings of the system and plots it and analyzes it. 
%
%
%Last edit: Wed Mar 6, 2019
%


addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/');





audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/FocusriteGen1/'

AUDIO_FILES = {'FocusriteGen1lineXLRshort.wav', 'FocusriteGen1lineTSshort.wav'};

%AUDIO_FILES = {[fsinst_shorted, fsinst_shortedgained, gain10reference, gain10recordnoise, recordnoise, reference, system_noise, system_gained] };


for i = (1:length(AUDIO_FILES)); 
    strcat(audio_dir,AUDIO_FILES{i});
    [data, fs] = audioread(strcat(audio_dir,AUDIO_FILES{i}));
    time = (0:length(data)-1)/fs;
    disp(AUDIO_FILES{i})
    disp(20.0*log10(rms(data))) %the default matlab rms function works fine
    %[data_fft, freq] = audio_spectrum(data, fs, 2^16, 2^16);
    %figure(1);
    %audio_plotspectrum(freq,data_fft);
    
    figure(2);
    grid on; hold on;
    %[data_psd, freq_psd] = audio_psd(data, 2^16, fs); 
    %audio_plotspectrum(freq_psd, data_psd)
    nfft = 2^16;
    %pwelch(data_L, 2^16, fs)
    [Pxx,freq]=pwelch(data,hanning(nfft,'periodic'),nfft/2,nfft,fs,'');
    plot(freq, 20.0*log10(Pxx(:,1)))  
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  


    saveas(gcf,strcat(AUDIO_FILES{i},'.png'))
end 
title('TS vs XLR input spectrum Focusrite')
legend(AUDIO_FILES)
    
    
    
    

