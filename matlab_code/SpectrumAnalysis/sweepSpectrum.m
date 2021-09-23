clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')

    tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/300b1600.957.wav');
end 

[dig_sweep, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/reference_files/digital_signals/sweep16kHz.wav');


figure(1)
time = (1:length(dig_sweep))/fs;
plot(time, dig_sweep)
title_string = '';
start_sam = 1; 
n_sam = 35*fs;

[data_fft, freq_fft] = audio_spectrum(dig_sweep, fs, start_sam, n_sam);

figure(2)
audio_plotspectrum(freq_fft, data_fft, title_string);
saveas(figure(2), 'spectrumdigitalsweep.png')

dig_sweepA = audio_Aweighting(data_fft);

[data_fft, freq_fft] = audio_spectrum(dig_sweep, fs, start_sam, n_sam);

figure(3)
audio_plotspectrum(freq_fft, data_fft, title_string);
saveas(figure(3), 'spectrumdigitalsweep.png')

seg = tracks('sweep');
[data_fft, freq_fft] = audio_spectrum(seg(:,1), fs, 1, length(seg));
figure(3)
audio_plotspectrum(freq_fft, data_fft, title_string);
saveas(figure(3), 'spectrumrecordsweep.png')


