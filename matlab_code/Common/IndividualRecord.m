clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/Code/vinyl-research/matlab_code/')

tracks = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/039b.wav')

signal_names = tracks.keys;
signals = tracks.values;
fs = 96000;
for i = 1:length(signal_names);
    trackname = signal_names{i};
    sig = signals{i};
    time = (1:length(sig))/fs;

    plot_track(time,sig,trackname,trackname)
    

    % [spec, fftfreq] = audio_spectrum(sig, fs, 1, 2^16);
    % % fig = figure('Visible', 'off')
    % figure(i+100)
    % audio_plotspectrum(fftfreq, spec, strcat(trackname,' spectrum'));
    
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
function plot_track(x, y, trackname, filename)
  
    fig = figure('Visible', 'off')
    subplot(2,1,1)
    title(strcat(trackname, ' left channel'))
    plot(x,y(:,1),'k')
    subplot(2,1,2)
    plot(x,y(:,2),'k')
    title(strcat(trackname, ' right channel'))
    grid on;
    ylim([-1,1])
    xlabel('time [s]')
    ylabel('signal level')
    plotname = strcat('IndividualRecord/', filename, '.png')
    saveas(fig, plotname);




end