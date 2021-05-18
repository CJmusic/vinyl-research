close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Volumes/AUDIOBANK/audio_files/')


% fs = 96000;
% N = 3*fs;


fs = 96000;
N = 9*fs;

time = (1:N)/fs;
data = zeros(length(time),1);
data(1) = 1; 

data = 0.3333*sin(1000*2*pi*time);
data = data.';
% size(y)

seg = data;

% tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')
% data = tracks('1kHz');
time= (0:length(data)-1)/fs;

figure(1)
plot(time, data)
title('data')
% seg = data(0.33*length(data):0.33*length(data) + N - 1,:);

win = flattopwin(N);
windowfactor = 0.2155774;
% winseg = seg.*win/windowfactor;
winseg = seg;
[seg_fft, freq_fft] = audio_spectrum(winseg, fs, 1, N); 
figure(2)
audio_plotspectrum(freq_fft, seg_fft, '');
title('data spectrum')

% seg_fft(1,:) = [0,0];
% seg_fft(length(seg_fft)-1,:) = [0,0];
% normalization = max(abs(seg_fft))/(sqrt(2));
normalization = rms_response(seg);
segnorm = seg./normalization;



timenorm= (0:length(segnorm)-1)/fs;
figure(3)
plot(timenorm, segnorm)
title('normalized data')

winsegnorm = segnorm.*win/windowfactor;
[seg_fftnorm, freq_fft] = audio_spectrum(winsegnorm, fs, 1, N);
figure(4)
audio_plotspectrum(freq_fft, seg_fftnorm, '');
% title('normalized data spectrum')
xlim([970,1030])
ylim([-80,0])
saveas(figure(4),'normalizedspectrumflattop.png')

[seg_fftnorm, freq_fft] = audio_spectrum(segnorm, fs, 1, N);
figure(5)
audio_plotspectrum(freq_fft, seg_fftnorm, '');
% title('no window normalized data spectrum')
xlim([970,1030])
ylim([-80,0])
saveas(figure(5),'normalizedspectrumnowindow.png')

RMSNORM = rms_response(segnorm);
disp(strcat('rms normalized data....', num2str(RMSNORM)))

% plot(y)
% [data_fft, freq_fft] = audio_spectrum(y, fs, 1, N);
% audio_plotspectrum(freq_fft, data_fft, '');
