clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')

% the 1kHz test signal on the test pressings is 7.071 cm/s peak or 5.00 cm/s rms
% we want this to be 0 dB for our analysis

%generate test signal
L = 2^16;
fs = 96000;
t=(0:L-1)'/fs;
A = 0.3456;
data = A*sin(1000*(L/fs)*2*pi*t);

rmsrawdata=rms_response(data) %initial rms

%window, find peak, correct for window and rms
windata =data.*flattopwin(L);
[windata_fft, freq_fft] = audio_spectrum(windata, fs, 1, L);
windowfactor = 0.2155774;
pk_freq_fft=max(abs(windata_fft))/windowfactor;
rms_freq_fft=pk_freq_fft/sqrt(2);

normdata=data/rms_freq_fft;%this is now normalized to rms of unity

rmsnormdata=rms_response(normdata) %check rms


%%CZ NORMALIZATION BELOW 


%~~~~~~~~~~~~~~~~~~~~ NORMALIZATION ~~~~~~~~~~~~~~~~~~~~%

%~~~~ separate out the 1 kHz track
% t = 1;
% sigtime = timedata(floor(timestamps(t,1)*fs):floor(timestamps(t,2)*fs));
% sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);

% sigtime = timedata(floor(timestamps(t,1)*fs) - lagdiff :floor(timestamps(t,2)*fs) - lagdiff);
% sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
data(:,1) = data;
data(:,2) = data;
size(data)
% sig = data.';
% size(data)

L = 2^16;
seg = data;
% seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);


figure(100)
plot(seg)
title('unnormalized data')

figure(101)
[fftrawseg, fftfreq]=audio_spectrum(seg,fs,1,L);
audio_plotspectrum(fftfreq, fftrawseg, 'before normalization')

%~~~~ window the data
win = flattopwin(L);
disp(strcat('size seg...', num2str(size(seg))))
disp(strcat('size win...', num2str(size(win))))

winseg = seg(:,1).*win; %% NEED TO NORMALIZE EACH CHANNEL -cz
disp(strcat('size winseg...', num2str(size(winseg))))

%------correct for window to get amplitude
[fftseg, fftfreq] = audio_spectrum(winseg, fs, 1, L);
windowfactor = 0.2155774;% 0.2155774 for flattop window, 0.5 for hann window
fftseg=fftseg/windowfactor;

figure(102)
audio_plotspectrum(fftfreq, fftseg, 'after windowing, correction')

%~~~ take peak of fft as pk amplitude, convert to rms to normalize
amplitude = max(abs(fftseg)); % I've already implemented, now rms
normalization = [amplitude, amplitude]/sqrt(2); %CHANGE need to normalize each channel individually

disp(strcat('max data before...', num2str(max(data))))
disp(strcat('max data before dB...', num2str(20.0*log10(max(abs(fftseg))))))
disp(strcat('rms data before...', num2str(rms(data))))
disp(strcat('rms data before dB...', num2str(20.0*log10(rms(data)))))

normalization

normalization_L = normalization(1);
normalization_R = normalization(2);
disp(strcat('normalization_L...', num2str(normalization_L)))
disp(strcat('normalization_R...', num2str(normalization_R)))

%~~ normalize the data 
data(:,1)=data(:,1)./normalization_L;% now normalized to 5cm/s rms    
data(:,2)=data(:,2)./normalization_R;% now normalized to 5cm/s rms 

sig = data;

% This is for seperating out the normalized 1 kHz track, not needed here

% t = 1;
% sig = data(floor(timestamps(t,1)*fs) - lagdiff : floor(timestamps(t,2)*fs) - lagdiff,:);
L = 2^16;        
% seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);
seg = sig;

[fftseg, fftfreq] = audio_spectrum(seg, fs, 1, L);
disp(strcat('max data after...', num2str(max(data))))
disp(strcat('max data after dB...', num2str(20.0*log10(max(abs(fftseg))))))
disp(strcat('rms data after...', num2str(rms(data))))
disp(strcat('rms data after dB...', num2str(20.0*log10(rms(data)))))

% sig = data(floor(timestamps(t,1)*fs): floor(timestamps(t,2)*fs),:);
[spec, fftfreq] = audio_spectrum(seg, fs, 1, 2^16);

figure(103)
audio_plotspectrum(fftfreq, spec, 'spectrum after normalization')

figure(104)
plot(data);
title('normalized raw data')
