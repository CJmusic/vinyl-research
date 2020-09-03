clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('D:\Code\vinyl-research\matlab_code\audio_functions')

%generate test signal
L = 2^16;
fs = 96000;
t = (0:L-1)/fs;
A = 0.5;% this amplitude should come out of the analysis
seg = [A*sin(1000*(L/fs)*2*pi*t); A*sin(1000*2*pi*t)];
seg = seg';%%%%%%%%%%%%%%%%%%%%

lr=1;%%%%%%%%%%%%%%%%% left or right chosen first
seg = seg(:,lr);%%%%%%%%%%%%%%%%

%plot test signal
figure(1)
plot(t,seg); grid on;
xlabel('Time [s]')
ylabel('Signal')
title('test signal')

%take and plot spectrum
figure(2)
[data_fft, freq_fft] = audio_spectrum(seg, fs, 1, L);%%%%%%%%%%%%%
audio_plotspectrum(freq_fft, data_fft, 'spectrum of test signal');

figure(5)
win = flattopwin(L);
%win = window(@flattopwin,L,'periodic');
plot(win); grid on;
title('flattop window')
seg2 =seg.*win;%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, L);
audio_plotspectrum(freq_fft, data_fft, 'spectrum after window');
disp(strcat('after window... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))

% fft normalization

% normalization = max(abs(data_fft)), when file is unity amplitude, flattop window
% this value turns out to be 0.215577429446209 if freq falls exactly on a bin

normalization = 0.2156;% for flattop window.  0.5 for hann window
seg3 = seg2/normalization;%%%%%%%%%%%%%

figure(7)
[data_fft, freq_fft] = audio_spectrum(seg3, fs, 1, L);
audio_plotspectrum(freq_fft, data_fft, 'normalized after window');
disp(strcat('after window norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))

amplitude = max(abs(data_fft))% this is normalization when A=1

