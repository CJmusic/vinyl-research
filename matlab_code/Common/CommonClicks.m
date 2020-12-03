clear all;close all;clc
set(0,'DefaultLineLineWidth',0.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - U   niversity of Waterloo/School/Vinyl_Project/audio_files/')

[record1, ~] = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav');
[record2, ~] = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/042119_A0000B0000r070a.wav');

sig = record1('transition');
sig2 = record2('transition');

[~, clicks_L] = ClickDetect(sig(:,1));
[~, clicks_R] = ClickDetect(sig(:,2));
[~, clicks2_L] = ClickDetect(sig2(:,1))
[~, clicks2_R] = ClickDetect(sig2(:,2))

plot_waveclicks(sig, clicks_L, clicks_R, 'Clicks detected in transition track 28a')
plot_waveclicks(sig2, clicks2_L, clicks2_R, 'Clicks detected in transition track 70a')


RELAX = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000];
COM_L = [];
COM_R = [];
for i = 1:length(RELAX); 
    common_L = num_comclicksTest(clicks_L, clicks2_L, RELAX(i))
    common_R = num_comclicksTest(clicks_R, clicks2_R, RELAX(i))
    COM_L = [COM_L, common_L];
    COM_R = [COM_R, common_R];
end

COM_L
COM_R

fig = figure('Visible', 'off')
scatter(RELAX, COM_L,'ko')
hold on;
scatter(RELAX, COM_R,'kx')
grid on;
legend({'left channel', 'right channel'}) 
title('Common clicks vs relaxation')
xlabel('relaxation')
ylabel('# of common clicks')
plotname = strcat('plots/CommonClicks/RelaxationPlotScatter.png')
saveas(fig, plotname);

fig = figure('Visible', 'off')
plot(RELAX, COM_L,'k')
hold on;
plot(RELAX, COM_R,'k--')
grid on;
legend({'left channel', 'right channel'}) 
title('Common clicks vs relaxation')
xlabel('relaxation')
ylabel('# of common clicks')
plotname = strcat('plots/CommonClicks/RelaxationPlotPlot.png')
saveas(fig, plotname);
% disp(strcat('common .....' , num2str
clicks_L
clicks_R
clicks2_L
clicks2_R


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
        line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
    end

    subplot(2,1,2)
    plot(sig(:,2),'k')
    ylim([-0.2, 0.2])
    hold on; grid on;
    % plot(sig2)
    title(strcat(filename, ' right channel'))

    for xi = 1:length(clicks_R)
        x1 = (clicks_R(xi));
        line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
    end

    plotname = strcat('plots/CommonClicks/', filename, '.png')
    saveas(fig, plotname);
end

% function num_comclicks = CommonClicks(clicks, clicks_ref);
function num_comclicks = num_comclicksTest(clicks, clicks_ref, relaxation);
    % relaxation = 750;
    % relaxation = 1024;
    comclicks = [];
    for i = (1:length(clicks_ref))
        upper_bound = clicks < clicks_ref(i) + relaxation;
        lower_bound  = clicks(upper_bound) > clicks_ref(i) - relaxation;
        if sum(lower_bound) > 0;
            comclicks = [comclicks, clicks_ref(i)];
            continue
        end
    end

    num_comclicks = length(comclicks);

end
