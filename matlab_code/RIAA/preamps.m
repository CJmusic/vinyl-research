clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Volumes/AUDIOBANK/audio_files/RIAA/')

[technics, fs] = audioread('/Volumes/AUDIOBANK/audio_files/RIAA/technicssweep.wav');

[tec_data_fft, tec_freq_fft] = audio_spectrum(technics(:,1), fs, 1, length(technics));

fig = figure(1)
audio_plotspectrum(tec_freq_fft, tec_data_fft, 'Frequency response of the Technics SU-9070 preamp') 
xlim([1,25000])
saveas(fig, 'plots/technics.png')

[optonica, fs] = audioread('/Volumes/AUDIOBANK/audio_files/RIAA/optonicasweep.wav');

[opt_data_fft, opt_freq_fft] = audio_spectrum(optonica(:,1), fs, 1, length(optonica));

fig = figure(2)
audio_plotspectrum(opt_freq_fft, opt_data_fft, 'Frequency response of the Optonica SM-1400 preamp') 
xlim([1,25000])
saveas(fig, 'plots/optonica.png')



fig = figure(3)
plot(tec_freq_fft, 20.0*log10(abs(tec_data_fft)), 'k') 
hold on;
plot(opt_freq_fft, 20.0*log10(abs(opt_data_fft)), 'b') 
legend('Technics SU 9070', 'Optonica SM-1400')
grid on 
set(gca, 'XScale', 'log');
title('Frequency response of the Technics and Optonica preamps')
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
xlim([1,25000])
saveas(fig, 'plots/together.png')




