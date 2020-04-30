clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
addpath('D:\Code\vinyl-research\matlab_code\')
addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')

record1 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');


csig = record1('transition');
plot(csig);

Aw = audio_Aweighting(csig(:,1));
CCIRw = audio_CCIRweighting(csig(:,1));


[data_fft, freq] = audio_spectrum(csig, 96000, 5*96000, 2^16);
[data_ffta, freq] = audio_spectrum(Aw, 96000, 5*96000, 2^16);
[data_fftccir, freq] = audio_spectrum(CCIRw, 96000, 5*96000, 2^16);
size(freq)
size(data_fft)

plot(freq, 20.0*log10(data_fft),'k--')
hold on; grid on;
plot(freq, 20.0*log10(data_ffta),'k-.')  
plot(freq, 20.0*log10(data_fftccir),'k-')  
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
legend('non-weighted','A-weighted','CCIR-weighted')


A_L = 20.0*log10(rms_response(Aw(1,:)));
A_R = 20.0*log10(rms_response(Aw(2,:)));

CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));