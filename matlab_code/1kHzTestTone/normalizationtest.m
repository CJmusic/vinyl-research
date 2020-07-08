clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)



L = 2^16;
fs = 96000;
t = (1:L)/fs;
seg = [sin(1000*2*pi*t); sin(1000*2*pi*t)];
seg = seg.';

figure(1)
plot(seg)

size(seg)

figure(3)
% win = flattopwin(L);
% seg2 = seg.*win;
seg2 = seg;
N = 2^12;
xdft = fft(seg2, N)/N;
xdft = xdft(1:N/2+1);
xdft(2:end-1) = xdft(2:end-1);
fdft = fs*(0:(N/2))/N;

size(xdft)
size(fdft)

plot(fdft,20*log10(abs(xdft)))
grid on 
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  


% max(seg)
% max(abs(xdft))
% 20*log10(max(abs(xdft)))

[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, 2^12);
figure(2)
audio_plotspectrum(freq_fft, data_fft, 'audio_spectrum');

max(abs(data_fft))
max(abs(xdft))

segRMS = rms(seg)
segRMS*sqrt(2)

% win = flattopwin(L);
% seg = seg.*win;

% size(seg)
% size(win)
% seg = seg.*win;
% fftsigL = fft(seg(:,1));
% % fftsigL = fft(seg(:,1));
% fftsigL = fftsigL(1:L/2+1)/L;
% fftsigR = fft(seg(:,2));
% % fftsigR = fft(seg(:,2))
% fftsigR = fftsigR(1:L/2+1)/L;
% fftfreq = fs*(0:(L/2))/L;

% figure(4)
% audio_plotspectrum(fftfreq, [fftsigL, fftsigR], 'norm')

% figure(5)
% plot(fftfreq, abs([fftsigL, fftsigR]))

% peak_L = max(abs(fftsigL));
% peak_R = max(abs(fftsigR));

% % sigRMS= [peak_L, peak_R]
% % normalization=sqrt(2)*sigRMS*40/7; %digital value of peak level
% normalization = [peak_L, peak_R]; %digital value of peak level


