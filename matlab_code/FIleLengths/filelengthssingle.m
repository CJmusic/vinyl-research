clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    % addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/')

    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/3.wav');
    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/samerecordtest/2.wav');

    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/042820_A0000B0000r1a.wav');
    % record1 = SeperateTracks('/AUDIOBANK/audio_files/A0137B0137/003a.wav')
    disp('Calling seperatetracks')
    record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav')
    record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3b.wav')

    % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/040318_A0000B0000r001a.wav');

end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')
    addpath('')

    record1 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/003a.wav');
    record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0137B0137/000a.wav');    
end

data1 = record1('1kHz');
data1 = data1(1:length(data1)-1,1);
time = (1:length(data1))/96000;

data1 = data1/max(data1);
plot(time, data1);

% plot(data1)
hold on; grid on;
data2 = record2('1kHz');
time = (1:length(data1))/96000;
data2 = data2(1:length(data1),1);
data2 = data2/max(data2);
size(data1)
size(data2)
plot(time, data2);
% plot(data2)


figure(2)
data3 = record1('leadout');
data3 = data3(:,1);
time = (1:length(data3))/96000;
plot(time,data3)
hold on; grid on;
data4 = record2('leadout');
data4 = data4(:,1);
time = (1:length(data4))/96000;
plot(time,data4)

% size(data1)
% size(data2)
% figure(3)
% plot((data1 - data2))
% sum((data1 - data2).^ 2)
% correlation = sqrt(sum((data1 - data2).^ 2))/length(data1)  % >= R2016b: auto-expand
% correlation = correlation / max(correlation)  % Normalize to [0, 1]

figure(3)
data3 = record1('sweepV2');
data3 = data3(:,1);
time = (1:length(data3))/96000;
plot(time,data3)
hold on; grid on;
data4 = record2('sweepV2');
data4 = data4(:,1);
time = (1:length(data4))/96000;
plot(time,data4)

fs = 96000;
datastart = 58*fs;
dataend = 60*fs;
data1 = data1(datastart:dataend);

figure(4)
[b,a]=butter(2,[700*2/fs 1400*2/fs]);
data1=filtfilt(b,a,data1);
signal = hilbert(data1);
time = (1:length(data1))/96000;
plot(time,abs(signal))
N = length(data1)
startsum=round(0.01*N);stopsum=round(0.25*N);
amplitde=sum(abs(signal(startsum:stopsum)))/(stopsum-startsum+1)

figure(5)
[acor_L,lags_L] = xcorr(data3,data4);
[M_L,I_L] = max(abs(acor_L));
lagdiff_L = lags_L(I_L);
lagdiff = lagdiff_L;
plot(lags_L, acor_L)