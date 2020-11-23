 clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

file_ref = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav';
file1 = '/Volumes/AUDIOBANK/audio_files/A0000B0000/042119_A0000B0000r070a.wav';
file2 = '/Volumes/AUDIOBANK/audio_files/A0000B0000/042119_A0000B0000r071a.wav';
file3 = '/Volumes/AUDIOBANK/audio_files/A0000B0000/042419_A0000B0000r074a.wav';

filename_ref = 'transition A0000B0000 28a'
filename1 = 'transition A0000B0000 70a';
filename2 = 'transition A0000B0000 71a';
filename3 = 'transition A0000B0000 74a';

[reference, ~] = SeperateTracks(file_ref)
[record1, ~] = SeperateTracks(file1); 
[record2, ~] = SeperateTracks(file2);
[record3, ~] = SeperateTracks(file3);

ref = reference('transition');
sig = record1('transition');
sig2 = record2('transition');
sig3 = record3('transition');

[~, clicksref_L] = ClickDetect(ref(:,1))
[~, clicksref_R] = ClickDetect(ref(:,2))
[~, clicks_L] = ClickDetect(sig(:,1))
[~, clicks_R] = ClickDetect(sig(:,2))
[~, clicks2_L] = ClickDetect(sig2(:,1))
[~, clicks2_R] = ClickDetect(sig2(:,2))
[~, clicks3_L] = ClickDetect(sig3(:,1))
[~, clicks3_R] = ClickDetect(sig3(:,2))

plot_waveclicks(ref, clicksref_L, clicksref_R, filename_ref)
plot_waveclicks(sig, clicks_L, clicks_R, filename1)
plot_waveclicks(sig2, clicks2_L, clicks2_R, filename2)
plot_waveclicks(sig3, clicks3_L, clicks3_R, filename3)

CommonClicks(clicksref_L, clicksref_L)
CommonClicks(clicksref_R, clicksref_R)
CommonClicks(clicks_L, clicksref_L)
CommonClicks(clicks_R, clicksref_R)
CommonClicks(clicks2_L, clicksref_L)
CommonClicks(clicks2_R, clicksref_R)
CommonClicks(clicks3_L, clicksref_L)
CommonClicks(clicks3_R, clicksref_R)

function plot_waveclicks(sig, clicks_L, clicks_R, filename)
    fig = figure('Visible', 'off')
    subplot(2,1,1)
    plot(sig(:,1),'k')
    ylim([-0.2, 0.2])
    hold on; grid on;
    % plot(sig2)
    title(strcat(filename, ' left channel'))
    
    for xi = 1:length(clicks_L)
        x1 = (clicks_L(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end

    subplot(2,1,2)
    plot(sig(:,1),'k')
    ylim([-0.2, 0.2])
    hold on; grid on;
    % plot(sig2)
    title(strcat(filename, ' right channel'))

    for xi = 1:length(clicks_R)
        x1 = (clicks_R(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end

    plotname = strcat('plots/ClickCompare/', filename, '.png')
    saveas(fig, plotname);
end