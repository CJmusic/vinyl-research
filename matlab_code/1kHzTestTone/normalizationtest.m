clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)



L = 2^16;
fs = 96000;
t = (1:L)/fs;
seg = [0.5*sin(1000*2*pi*t); 0.5*sin(1000*2*pi*t)];
seg = seg.';

figure(1)
plot(seg)

figure(2)
[data_fft, freq_fft] = audio_spectrum(seg, fs, 1, 2^12);
% disp(strcat('max(abs(data_fft))... ', num2str(20*log10(max(abs(data_fft))))))
audio_plotspectrum(freq_fft, data_fft, 'audio_spectrum');

disp(strcat('before norm... ', num2str(max(abs(data_fft)))))

segRMS = rms(seg);
disp(strcat('segRMS root 2... ', num2str(segRMS*sqrt(2))))

normalization = segRMS*sqrt(2);
seg2 = seg/normalization;

figure(3)
plot(seg2)

figure(4)
[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, 'audio_spectrum');
disp(strcat('after norm... ', num2str(max(abs(data_fft)))))

win = flattopwin(L);
figure(10)
plot(win)
seg3 = seg2.*win;

figure(5)
[data_fft, freq_fft] = audio_spectrum(seg3, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, 'audio_spectrum');
disp(strcat('after window... ', num2str(max(abs(data_fft)))))

