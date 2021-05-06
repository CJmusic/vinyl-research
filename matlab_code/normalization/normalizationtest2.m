clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Volumes/AUDIOBANK/audio_files/')

fs = 96000;
N = 0.05*fs;
time = (0:N-1)/fs;

% seg = [0.5*sin(1000*2*pi*time); 0.5*sin(1000*2*pi*time)];
% seg = seg.';

tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-03b1558.270.wav')
seg = tracks('1kHz');
seg = seg(0.33*length(seg):0.33*length(seg) + N - 1,:);

figure(2)
plot(seg)

% seg = seg(:,1);
win = flattopwin(N);
windowfactor = 0.2155774;
winseg = seg.*win/windowfactor;

[seg_fft, freq_fft] = audio_spectrum(winseg, fs, 1, N);

normalization = max(abs(seg_fft))/(sqrt(2));
segnorm = seg./normalization;
winsegnorm = segnorm.*win/windowfactor;
[seg_fftnorm, freq_fft] = audio_spectrum(winsegnorm, fs, 1, N);


disp(strcat('max fft before norm... ', num2str(max(abs(seg_fft)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fft))))))
disp(strcat('max fft after norm... ', num2str(max(abs(seg_fftnorm)))))
disp(strcat('dB... ', num2str(20*log10(max(abs(seg_fftnorm))))))

normalization_L = normalization(1);
normalization_R = normalization(2);

%plot test signal
figure(1)
time = (0:length(seg)-1)/fs;
plot(time, seg); grid on;
xlabel('Time [s]')
ylabel('Signal')
% title('test signal')
legend(['left channel', 'right channel'])
saveas(figure(1),'plots/testsignal.png')

figure(2)
audio_plotspectrum(freq_fft, seg_fft, '');
legend(['left channel', 'right channel'])
saveas(figure(2),'plots/testsignalspectrum.png')

figure(3)
timenorm = (0:length(segnorm)-1)/fs;
plot(timenorm,segnorm); grid on; 
xlabel('Time [s]')
ylabel('Signal')
legend(['left channel', 'right channel'])
saveas(figure(3),'plots/testsignalnorm.png')

figure(4)
audio_plotspectrum(freq_fft, seg_fftnorm, '');
legend(['left channel', 'right channel'])
saveas(figure(4),'plots/spectrumnorm.png')


figure(5)
plot(win,'k')
grid on;
ylabel('Amplitude')
xlabel('Sample')
saveas(figure(5),'plots/window.png')
