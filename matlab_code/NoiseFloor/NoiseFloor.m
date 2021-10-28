clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')

[data, fs] = audioread('/Volumes/AUDIOBANK/digitalnoisefloor.wav');

[data_fft, freq] = audio_spectrum(data, 96000, 1, 2^16);
data_fft = pwroctsmooth_singlesided(data_fft,0.33);


plot(freq, 20.0*log10(abs(data_fft)),'k') 
grid on 
set(gca, 'XScale', 'log');
title('Noise floor of Recording System')
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

