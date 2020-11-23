clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

file_ref = '/Volumes/AUDIOBANK/audio_files/A0137B0137/074a.wav';

[reference, ~] = SeperateTracks(file_ref);

sig = reference('transition');
fs = 96000;
%plot the click and the spectrum, as well as the corrected audio

[csig_L, clicks_L] = ClickDetect(sig(:,1));
[csig_R, clicks_R] = ClickDetect(sig(:,2));

size(csig_L)
size(csig_R)

csig = [csig_L, csig_R];

size(csig)


for i = (1:length(clicks_L))

    start_sam = clicks_L(i) - 0.0005*fs;
    if start_sam < 0
        start_sam = 1;
    end

    end_sam = clicks_L(i) + 0.002*fs;
    if end_sam > length(sig)
        end_sam = length(sig);
    end


    click = sig(start_sam:end_sam,1);
    plot_click(click, strcat('Click ',num2str(i) ,' detected in left channel on record A0137B0137 74 side a'))

    cclick = csig(start_sam:end_sam,1);
    plot_click(cclick, strcat('Corrected audio of Click ',num2str(i) ,' detected in left channel on record A0137B0137 74 side a'))

    [click_spec, freq] = audio_spectrum(click,96000, 1, length(click));
    


end

for i = (1:length(clicks_R))
    start_sam = clicks_L(i) - 0.0005*fs;
    if start_sam < 0
        start_sam = 1;
    end

    end_sam = clicks_L(i) + 0.002*fs;
    if end_sam > length(sig)
        end_sam = length(sig);
    end

    click = sig(start_sam:end_sam,2);
    plot_click(click, strcat('Click ',num2str(i) ,' detected in right channel on record A0137B0137 74 side a'))

    cclick = csig(start_sam:end_sam,2);
    plot_click(cclick, strcat('Corrected audio of Click ',num2str(i) ,' detected in right channel on record A0137B0137 74 side a'))
end



function plot_click(click, titlestring)
    fs = 96000;
    time = (1:length(click))/fs;
    time = time*1000;

    fig = figure('Visible', 'off')
    plot(time,click,'k')
    grid on;
    ylim([-0.5,0.5])
    title(titlestring)
    xlabel('time [ms]')
    ylabel('signal level')
    plotname = strcat('plots/PlotClicks/', titlestring, '.png')
    saveas(fig, plotname);

end