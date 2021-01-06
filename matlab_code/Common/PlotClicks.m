clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

file_ref = '/Volumes/AUDIOBANK/audio_files/A0137B0137/-01a1558.884.wav';
file_ref = '/Volumes/AUDIOBANK/audio_files/A0137B0137/-02a1558.981.wav';
[reference, ~] = SeperateTracks(file_ref);

sig = reference('10kHz');
figure(1)
plot(sig)
fs = 96000;
%plot the click and the spectrum, as well as the corrected audio

[csig_L, clicks_L] = ClickDetect(sig(:,1));
[csig_R, clicks_R] = ClickDetect(sig(:,2));

disp('NUMBER OF CLICKS')
size(clicks_L)
size(clicks_R)

csig = [csig_L, csig_R];

size(csig)


for i = (1:length(clicks_L))

    start_sam = clicks_L(i) - 0.02*fs;
    if start_sam < 0
        start_sam = 1;
    end

    end_sam = clicks_L(i) + 0.02*fs;
    if end_sam > length(sig)
        end_sam = length(sig);
    end


    click = sig(start_sam:end_sam,1);
    plot_click(click, strcat('click ', num2str(i), ' test1'))
    [click_spec, freq] = audio_spectrum(click,96000, 1, length(click));
    


end

for i = (1:length(clicks_R))
    % start_sam = clicks_L(i) - 0.0005*fs;
    start_sam = clicks_L(i) - 0.02*fs;
    if start_sam < 0
        start_sam = 1;
    end

    end_sam = clicks_L(i) + 0.02*fs;
    if end_sam > length(sig)
        end_sam = length(sig);
    end

    click = sig(start_sam:end_sam,2);
    plot_click(click, strcat('Click ',num2str(i) ,' test2'))

end

plot_waveclicks(sig, clicks_L, clicks_R, 'TESTINGWAVES')
length(clicks_L)
length(clicks_R)

function plot_click(click, titlestring)
    disp('PLOT CLICKS CALLED')
    fs = 96000;
    time = (1:length(click))/fs;
    time = time*1000;

    fig = figure('Visible', 'off')
    plot(time,click,'k')
    grid on;
    % ylim([-0.5,0.5])
    title(titlestring)
    xlabel('time [ms]')
    ylabel('signal level')
    plotname = strcat('plots/PlotClicks/1kHz/', titlestring, '.png')
    saveas(fig, plotname);

end


function plot_waveclicks(sig, clicks_L, clicks_R, filename)
    fig = figure('Visible', 'off')
    subplot(2,1,1)
    plot(sig(:,1),'k')
    ylim([-2, 2])
    hold on; grid on;
    % plot(sig2)
    title(strcat(filename, ' left channel'))
    
    for xi = 1:length(clicks_L)
        x1 = (clicks_L(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
    end

    subplot(2,1,2)
    plot(sig(:,2),'k')
    hold on; grid on;
    ylim([-2, 2])
    % plot(sig2)
    title(strcat(filename, ' right channel'))

    for xi = 1:length(clicks_R)
        x1 = (clicks_R(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
    end


    plotname = strcat('plots/PlotClicks/1kHz/', filename, '.png')
    saveas(fig, plotname);
end

