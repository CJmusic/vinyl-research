clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Volumes/AUDIOBANK/audio_files/')

%generate test signal
L = 2^16;
fs = 96000;
t = (0:L-1)/fs;
seg = [0.5*sin(1000*2*pi*t); 0.5*sin(1000*2*pi*t)];
seg = seg.';
%generate test signal ends

%use real test signal
tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')

seg = tracks('1kHz');
seg = seg(0.33*length(seg):0.33*length(seg) + L - 1,:);

size(seg)


%use real test signal ends




lr=1;%%%%%%%%%%%%%%%%% left or right chosen first
seg = seg(:,lr);%%%%%%%%%%%%%%%%

%plot test signal
figure(1)
plot(seg); grid on;
title('test signal')
xlabel('Time [s]')
ylabel('Signal')
title('test signal')

%take and plot spectrum
figure(2)
[data_fft, freq_fft] = audio_spectrum(seg, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, 'spectrum of test signal');

% fft normalization

normalization = max(abs(data_fft));
seg2 = seg/normalization;

disp(strcat('max fft before norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))

% rms normalization

% segRMS = rms(seg);
% normalization = segRMS*sqrt(2);
% disp(strcat('segRMS * root 2... ', num2str(segRMS*sqrt(2))))

figure(3)
plot(seg2); grid on; 
title('test signal after normalization')

figure(4)
[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, 'spectrum after normalization');
disp(strcat('max fft after norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))

figure(10)
win = flattopwin(L);
% win = window(@flattopwin,L,'periodic');
plot(win); grid on;
title('flattop window')
% seg3 = seg2.*win*(length(win)/sum(win));
size(seg)
size(win)

seg3 = seg.*win;
seg3 = seg3/0.2156;

figure(5)
[data_fft, freq_fft] = audio_spectrum(seg3, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, 'spectrum after window');
disp(strcat('after window... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))


