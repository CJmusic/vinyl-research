clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/')

tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/110b.wav')

signal_names = tracks.keys;
signals = tracks.values;
fs = 96000;
for i = 1:length(signal_names);
    trackname = signal_names{i};
    sig = signals{i};
    time = (1:length(sig))/fs;

    figure(i)
    title(trackname)
    subplot(2,1,1)
    plot(time, sig(:,1))
    title('left channel')
    grid on;
    subplot(2,1,2)
    plot(time, sig(:,2))
    title('right channel')
    grid on;
    
    

    [spec, fftfreq] = audio_spectrum(sig, fs, 1, 2^16);
    % fig = figure('Visible', 'off')
    figure(i+100)
    audio_plotspectrum(fftfreq, spec, strcat(trackname,' spectrum'));
    
end


% fs = 96000;
% seg = tracks('1kHz');
% % seg = seg(1:2^16,:);
% L = 2^16;
% seg = seg(floor(length(seg)/2) - L/2:floor(length(seg)/2) + L/2 - 1,:);
% size(seg)

% [spec, fftfreq] = audio_spectrum(seg, fs, 1, 2^16);

% figure(2)
% audio_plotspectrum(fftfreq, spec, 'after seperate tracks')
