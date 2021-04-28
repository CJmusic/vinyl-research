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
seg = [0.1*sin(1000*2*pi*t); 0.1*sin(1000*2*pi*t)];
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
seg = 0.5*seg;
%plot test signal
figure(1)
time = (0:length(seg)-1)/fs;
plot(seg,'k'); grid on;
% title('test signal')
xlabel('Time [s]')
ylabel('Signal')
% title('test signal')
saveas(figure(1),'plots/testsignal.png')

%take and plot spectrum
figure(2)
[data_fft, freq_fft] = audio_spectrum(seg, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, '');
saveas(figure(2),'plots/testsignalspectrum.png')



% fft normalization

normalization = max(abs(data_fft));
normalization2 = rms_response(seg);

% normalization = normalization2

seg2 = seg/normalization;

disp(strcat('max fft before norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))

% rms normalization

% segRMS = rms(seg);
% normalization = segRMS*sqrt(2);
% disp(strcat('segRMS * root 2... ', num2str(segRMS*sqrt(2))))

figure(3)
time2 = (0:length(seg2)-1)/fs;
plot(time2,seg2,'k'); grid on; 
xlabel('Time [s]')
ylabel('Signal')
% title('test signal after normalization')
saveas(figure(3),'plots/testsignalnorm.png')

figure(4)
[data_fft, freq_fft] = audio_spectrum(seg2, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, '');
disp(strcat('max fft after norm... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))
saveas(figure(4),'plots/spectrumnorm.png')


figure(10)
win = flattopwin(L);
% win = window(@flattopwin,L,'periodic');
plot(win,'k'); grid on;
xlabel('Sample number')
ylabel('Level')
% title('flattop window')
% seg3 = seg2.*win*(length(win)/sum(win));
size(seg)
size(win)

seg3 = seg.*win;
seg3 = seg3/0.2156;
saveas(figure(10),'plots/window.png')


figure(5)
[data_fft, freq_fft] = audio_spectrum(seg3, fs, 1, 2^12);
audio_plotspectrum(freq_fft, data_fft, '');
disp(strcat('after window... ', num2str(max(abs(data_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(data_fft))))))
saveas(figure(5),'plots/spectrumwindowed.png')


normalization
normalization2