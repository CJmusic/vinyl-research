clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Volumes/AUDIOBANK/audio_files/')

L = 2^16;
fs = 96000;
t = (0:L-1)/fs;


%generate test signal
% seg = [0.1*sin(1000*2*pi*t); 0.1*sin(1000*2*pi*t)];
% seg = seg.';
%generate test signal ends

% %use real test signal
tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')
seg = tracks('1kHz');
seg = seg(0.33*length(seg):0.33*length(seg) + L - 1,:);

%use real test signal ends



%%%%%%%%%%%%%%%%% left or right chosen first %%%%%%%%%%%%%%%%
lr=1;
seg = seg(:,lr);
seg = seg;
%plot test signal
figure(1)
time = (0:length(seg)-1)/fs;
plot(seg,'k'); grid on;
xlabel('Time [s]')
ylabel('Signal')
% title('test signal')
saveas(figure(1),'plots/testsignal.png')

%take and plot spectrum
L = 2^16;
win = flattopwin(L);
winseg = seg.*win;
windowfactor = 0.2155774;%4.6433 matlab answer, 4.858 my correction factor, 0.2155774 for flattop window, 0.5 for hann window
winseg = winseg/windowfactor;


figure(2)
[data_fft, freq_fft] = audio_spectrum(winseg, fs, 1, L);
audio_plotspectrum(freq_fft, data_fft, '');
saveas(figure(2),'plots/testsignalspectrum.png')



%~~~~~~~~~~~ fft normalization ~~~~~~~~~~~%
 
normalization = max(abs(data_fft))/sqrt(2);

%~~~~~~~~~~~ fft normalization ~~~~~~~~~~~%

segRMS = rms_response(seg);
normalization2 = segRMS/sqrt(2);
disp(strcat('segRMS * root 2... ', num2str(segRMS*sqrt(2))))

[data_fft3, freq_fft3] = audio_spectrum(seg, fs, 1, L);
audio_plotspectrum(freq_fft3, data_fft3, '');
normalization3 = max(abs(data_fft3))/sqrt(2);


segnorm = seg/normalization;
segnorm2 = seg/normalization2;
segnorm3 = seg/normalization3;


disp(strcat('max fft before norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))


figure(3)
timenorm = (0:length(segnorm)-1)/fs;
plot(timenorm,segnorm,'k'); grid on; 
xlabel('Time [s]')
ylabel('Signal')
% title('test signal after normalization')
saveas(figure(3),'plots/testsignalnorm.png')

figure(4)
[data_fft, freq_fft] = audio_spectrum(segnorm, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, '');
disp(strcat('norm... ', num2str(max(abs(normalization)))))
disp(strcat('max fft after norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))
saveas(figure(4),'plots/spectrumnorm.png')

figure(5)
timenorm2 = (0:length(segnorm2)-1)/fs;
plot(timenorm2,segnorm2,'k'); grid on; 
xlabel('Time [s]')
ylabel('Signal')
% title('test signal after normalization')
saveas(figure(3),'plots/testsignalnorm2.png')

figure(6)
[data_fft2, freq_fft2] = audio_spectrum(segnorm2, fs, 1, 2^12);
audio_plotspectrum(freq_fft2, data_fft2, '');
disp(strcat('norm2... ', num2str(max(abs(normalization2)))))
disp(strcat('max fft after norm2... ', num2str(max(abs(data_fft2)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft2))))))
saveas(figure(4),'plots/spectrumnorm2.png')


figure(7)
timenorm3 = (0:length(segnorm3)-1)/fs;
plot(timenorm3,segnorm3,'k'); grid on; 
xlabel('Time [s]')
ylabel('Signal')
% title('test signal after normalization')
saveas(figure(3),'plots/testsignalnorm3.png')

figure(8)
[data_fft3, freq_fft3] = audio_spectrum(segnorm3, fs, 1, 2^12);
audio_plotspectrum(freq_fft3, data_fft3, '');
disp(strcat('norm2... ', num2str(max(abs(normalization3)))))
disp(strcat('max fft after norm2... ', num2str(max(abs(data_fft3)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft3))))))
saveas(figure(4),'plots/spectrumnorm3.png')
