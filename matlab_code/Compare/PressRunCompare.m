clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

c = gray(20);


% Load in A0000B0000 table/audio table
% load in A0137B0137 table/audio table

% do histograms for RMS of each

A0137B0137 = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137-AudioTableMay12.csv');
A0000B0000 = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTableMay26.csv');

head(A0137B0137)
head(A0000B0000)

tracks = unique(A0137B0137.track);
tracks

plotnum = 0;

for i = (1:length(tracks))
    track = tracks{i};
    
    

    if strcmp(track, 'quiet') || strcmp(track, 'transition') || strcmp(track, 'quiet2')


        statsLa = datastats(A0000B0000.RMS_L(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'a'),:))
        statsRa = datastats(A0000B0000.RMS_R(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'a'),:))
        
        
        lower_binLa = statsLa.mean - 10;
        lower_binRa = statsRa.mean - 10;
        upper_binLa = statsLa.mean + 10;
        upper_binRa = statsRa.mean + 10;


        statsLb = datastats(A0000B0000.RMS_L(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'b'),:))
        statsRb = datastats(A0000B0000.RMS_R(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'b'),:))

        lower_binLb = statsLb.mean - 10;
        lower_binRb = statsRb.mean - 10;
        upper_binLb = statsLb.mean + 10;
        upper_binRb = statsRb.mean + 10;


        plotnum = plotnum + 1;
        figure(plotnum)
        histogram(A0000B0000.RMS_L(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'a'),:),50,'BinLimits',[lower_binLa,upper_binLa], 'facecolor',[0.3 0.3 0.3])
        hold on; grid on;
        histogram(A0000B0000.RMS_R(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'a'),:),50,'BinLimits',[lower_binLa,upper_binLa], 'facecolor',[0.3 0.3 0.3])
        title(strcat('A0000B0000 ', track, ' RMS noise side a'))
        legend('left channel', 'right channel')
        dim = [0.2 0.5 0.3 0.3];
        str = {strcat('left mean :',num2str(statsLa.mean)),strcat('left std :',num2str(statsLa.std)),strcat('right mean :',num2str(statsRa.mean)),strcat('right std :',num2str(statsRa.std))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        saveas(figure(plotnum),strcat(track, 'A0000B0000aRMS.png'))
        

        title(strcat('A0000B0000 ', track, ' RMS noise side a'))
        
        plotnum = plotnum + 1;
        figure(plotnum)
        histogram(A0000B0000.RMS_L(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'b'),:),50,'BinLimits',[lower_binLb,upper_binLb], 'facecolor',[0.3 0.3 0.3])
        hold on; grid on;
        histogram(A0000B0000.RMS_R(strcmp(A0000B0000.track, track) & strcmp(A0000B0000.side, 'b'),:),50,'BinLimits',[lower_binLb,upper_binLb], 'facecolor',[0.3 0.3 0.3])
        title(strcat('A0000B0000 ', track, ' RMS noise side b'))
        legend('left channel', 'right channel')
        dim = [0.2 0.5 0.3 0.3];
        str = {strcat('left mean :',num2str(statsLa.mean)),strcat('left std :',num2str(statsLa.std)),strcat('right mean :',num2str(statsRa.mean)),strcat('right std :',num2str(statsRa.std))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        saveas(figure(plotnum),strcat(track, 'A0000B0000bRMS.png'))

        % plotnum += 1;
        % figure(plotnum)

        statsLa = datastats(A0137B0137.RMS_L(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'a'),:))
        statsRa = datastats(A0137B0137.RMS_R(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'a'),:))

        lower_binLa = statsLa.mean - 10;
        lower_binRa = statsRa.mean - 10;

        upper_binLa = statsLa.mean + 10;
        upper_binRa = statsRa.mean + 10;

        plotnum = plotnum + 1;
        figure(plotnum)
        histogram(A0137B0137.RMS_L(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'a'),:),50,'BinLimits',[lower_binLa, upper_binLa], 'facecolor',[0.3 0.3 0.3])
        hold on; grid on;
        histogram(A0137B0137.RMS_R(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'a'),:),50,'BinLimits',[lower_binLa, upper_binLa], 'facecolor',[0.3 0.3 0.3])
        legend('left channel', 'right channel')
        dim = [0.2 0.5 0.3 0.3];
        str = {strcat('left mean :',num2str(statsLa.mean)),strcat('left std :',num2str(statsLa.std)),strcat('right mean :',num2str(statsRa.mean)),strcat('right std :',num2str(statsRa.std))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title(strcat('A0137B0137', track, ' RMS noise side a'))
        saveas(figure(plotnum),strcat(track, 'A0137B0137aRMS.png'))

        statsLb = datastats(A0137B0137.RMS_L(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'b'),:))
        statsRb = datastats(A0137B0137.RMS_R(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'b'),:))

        lower_binLb = statsLb.mean - 10;
        lower_binRb = statsRb.mean - 10;

        upper_binLb = statsLb.mean + 10;
        upper_binRb = statsRb.mean + 10;


        plotnum = plotnum + 1;
        figure(plotnum)
        histogram(A0137B0137.RMS_L(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'b'),:),50,'BinLimits',[lower_binLb,upper_binLb], 'facecolor',[0.3 0.3 0.3])
        hold on; grid on;
        histogram(A0137B0137.RMS_R(strcmp(A0137B0137.track, track) & strcmp(A0137B0137.side, 'b'),:),50,'BinLimits',[lower_binLb,upper_binLb], 'facecolor',[0.3 0.3 0.3])
        legend('left channel', 'right channel')
        dim = [0.2 0.5 0.3 0.3];
        str = {strcat('left mean :',num2str(statsLa.mean)),strcat('left std :',num2str(statsLa.std)),strcat('right mean :',num2str(statsRa.mean)),strcat('right std :',num2str(statsRa.std))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title(strcat('A0137B0137 ', track, ' RMS noise side b'))
        saveas(figure(plotnum),strcat(track, 'A0137B0137bRMS.png'))


    end

end

% figure(1)
% histogram(A0000B0000.RMS_L(strcmp(A0000B0000.track,'transition'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.3 0.3 0.3])


% figure(2)
% histogram(A0137B0137.RMS_L(strcmp(A0137B0137.track,'transition'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.3 0.3 0.3])
