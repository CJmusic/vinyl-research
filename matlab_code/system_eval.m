%This file takes the noise recordings of the system and plots it and analyzes it. 
%
%
%Last edit: Wed Mar 6, 2019
%


addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');





audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/'

AUDIO_FILES = {'fsinst-shorted.wav';'fsint-shortedgained.wav';'gain10recordnoise.wav';'gain10reference.wav';'recordnoise.wav';'reference.wav';'system.wav';'systemgained.wav'};

%AUDIO_FILES = {[fsinst_shorted, fsinst_shortedgained, gain10reference, gain10recordnoise, recordnoise, reference, system_noise, system_gained] };


for i = (1:length(AUDIO_FILES)); 
    strcat(audio_dir,AUDIO_FILES{i});
    [data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{i}));
    
    disp(AUDIO_FILES{i})
    disp(20.0*log10(rms(data))) %the default matlab rms function works fine
    %[data_fft, freq] = audio_spectrum(data, fs, 2^16, 2^16);
    %figure(1);
    %audio_plotspectrum(freq,data_fft);
    
    figure(2);
    grid on;
    %[data_psd, freq_psd] = audio_psd(data, 2^16, fs); 
    %audio_plotspectrum(freq_psd, data_psd)
    nfft = 2^16;
    %pwelch(data_L, 2^16, fs)
    [Pxx,f]=pwelch(data,hanning(nfft,'periodic'),nfft/2,nfft,fs,'');
    audio_plotspectrum(f, Pxx, strcat(AUDIO_FILES{i},' PSD'));
    saveas(gcf,strcat(AUDIO_FILES{i},'.png'))
end 
    
    
    
    
    

