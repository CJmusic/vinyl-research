clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Volumes/AUDIOBANK/audio_files/')

N = 2^16;
fs = 96000;
time = (0:N-1)/fs;

% seg = [0.5*sin(1000*2*pi*time); 0.5*sin(1000*2*pi*time)];
% seg = seg.';

tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')
seg = tracks('1kHz');
seg = seg(0.33*length(seg):0.33*length(seg) + N - 1,:);



% seg = seg(:,1);
win = flattopwin(N);
windowfactor = 0.2155774;
winseg = seg.*win;

[seg_fft, freq_fft] = audio_spectrum(winseg, fs, 1, N);

normalization = max(abs(seg_fft))/(sqrt(2)*windowfactor);
segnorm = seg./normalization;
winsegnorm = segnorm.*win/windowfactor;
[seg_fftnorm, freq_fft] = audio_spectrum(winsegnorm, fs, 1, N);


disp(strcat('max fft before norm... ', num2str(max(abs(seg_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fft))))))
disp(strcat('max fft after norm... ', num2str(max(abs(seg_fftnorm)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fftnorm))))))

normalization_L = normalization(1);
normalization_R = normalization(2);
