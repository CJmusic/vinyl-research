clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

% the 1kHz test signal on the test pressings is 7.07cm/s peak or 5.00 cm/s rms
% we want this to plot to 0 dB on our windowed spectra
%generate test signal
L = 2^16;
fs = 96000;
t = (0:L-1)'/fs;
A = 0.3456;% this amplitude should come out of the analysis
data = A*sin(1000*(L/fs)*2*pi*t);

%plot raw test signal
figure(1)
plot(t,data); grid on;
xlabel('Time [s]')
ylabel('amplitude')
title('1kHz raw test signal')

%take and plot spectrum
figure(2)
[wincdata_fft, freq_fft] = audio_spectrum(data, fs, 1, L);
audio_plotspectrum(freq_fft, wincdata_fft, 'unwindowed spectrum of raw test signal');

figure(5)
win = flattopwin(L);
%win = window(@flattopwin,L,'periodic');
plot(win); grid on;
title('flattop window')

windata =data.*win;% this windowed data will have a lower rms value

% When file is unity pk amplitude, flattop window gives 0.215577429446209
% for max(abs(data_fft)), if freq falls exactly on a bin

figure(6)
windowfactor = 0.2156;% 0.2155774 for flattop window, 0.5 for hann window
%inverse of this factor will correct spectral peak amplitude for windowed data
%window-corrected-data is windata/windowfactor%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wincdata_fft, freq_fft] = audio_spectrum(windata/windowfactor, fs, 1, L);
% wincdata_fft is now corrected for the window
audio_plotspectrum(freq_fft, wincdata_fft, 'unnormalized window-corrected spectrum');
disp(strcat('pk after window... ', num2str(max(abs(wincdata_fft)))))
disp(strcat('pk dB... ', num2str(20*log10(max(abs(wincdata_fft))))))

% the data is now windowed but still has the original rms amplitude
% we want the rms amplitude to be 1.0 (i.e. 0 dB)

pk_freq_fft=max(abs(wincdata_fft))%normalizing factor %%%%%%%%%%%%%%%%%%%%%%%

%now we use this peak to normalize the corrected windowed data

% normdata=windata/(windowfactor*pk_freq_fft);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normdata=data/(windowfactor*pk_freq_fft);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
[normdata_fft, freq_fft] = audio_spectrum(normdata, fs, 1, L);
% data_fft is single sided
audio_plotspectrum(freq_fft, normdata_fft, 'normalized windowed spectrum');
disp(strcat('pk after window norm... ', num2str(max(abs(normdata_fft)))))
disp(strcat('pk dB... ', num2str(20*log10(max(abs(normdata_fft))))))

amplitude = max(abs(normdata_fft))% this is final normalized value

figure(1)
plot(t,normdata); grid on;
xlabel('Time [s]')
ylabel('amplitude')
title('1kHz normalized test signal')
