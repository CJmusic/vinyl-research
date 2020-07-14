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

size(seg)

figure(2)
% win = flattopwin(L);
% seg2 = seg.*win;
seg2 = seg;
N = 2^12;
xdft = fft(seg2, N)/N;
xdft = xdft(1:N/2+1);
psdx = (1/(9600*N)) * abs(xdft).^2;
xdft(2:end-1) = 2*xdft(2:end-1);
fdft = fs*(0:(N/2))/N;

size(xdft)
size(fdft)

plot(fdft,20*log10(abs(xdft)))
% plot(fdft,10*log10(abs(psdx)))
grid on 
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

disp(strcat('max xdft... ', num2str(20*log10(max(abs(xdft))))))


% max(seg)
% max(abs(xdft))
% 20*log10(max(abs(xdft)))

[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, 2^12);
% disp(strcat('max(abs(data_fft))... ', num2str(20*log10(max(abs(data_fft))))))
figure(3)
audio_plotspectrum(freq_fft, data_fft, 'audio_spectrum');

disp(strcat('max(abs(data_fft))', num2str(max(abs(data_fft)))))
disp(strcat('max(abs(xdft))', num2str(max(abs(xdft)))))

segRMS = rms(seg);
disp(strcat('segRMS root 2... ', num2str(segRMS*sqrt(2))))




