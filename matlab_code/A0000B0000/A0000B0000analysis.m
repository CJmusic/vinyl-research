clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

if ismac()
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/')
    data_folder = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/A0000B0000/'
    AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTable.csv') 
end
if ispc()
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\')
    data_folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'
    AudioTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\A0000B0000-AudioTable.csv');
end


plotnum = 0;


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[-50,-30])
histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[-50,-30])
legend('RMS Left Channel', 'RMS Right Channel')
ylabel('number of records')
xlabel('RMS level [dB]')
title('RMS noise in quiet track')
xlim([-50,-30])
saveas(figure(plotnum),'RMSquiet.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500])
histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500])
legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
ylabel('number of records')
xlabel('number of clicks')
title('number of clicks in quiet track')
xlim([0,500])
saveas(figure(plotnum),'clicksquiet.png')


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30])
histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30])
legend('RMS Left Channel', 'RMS Right Channel')
ylabel('number of records')
xlabel('RMS level [dB]')
title('RMS noise in quiet2 track')
xlim([-50,-30])
saveas(figure(plotnum),'RMSquiet2.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500])
histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500])
legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
ylabel('number of records')
xlabel('RMS level [dB]')
title('RMS noise in quiet2 track')
% xlim([0,500])
saveas(figure(plotnum),'clicksquiet2.png')



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30])
histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30])
legend('RMS Left Channel', 'RMS Right Channel')
ylabel('number of records')
xlabel('RMS level [dB]')
title('RMS noise in transition track')
xlim([-50,-30])
saveas(figure(plotnum),'rmstransition.png')


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000])
histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000])
legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
ylabel('number of Records')
xlabel('number of clicks')
title('num clicks in transition track')
xlim([0,3000])
saveas(figure(plotnum),'clickstransition.png')


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50])
histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50])
legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
ylabel('number of Records')
xlabel('wow and flutter')
title('Wow and flutter in 3150Hz track')
saveas(figure(plotnum),'wow1.png')



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50])
histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50])
legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
ylabel('number of Records')
xlabel('wow and flutter')
title('Wow and flutter in 3150Hz2 track')
saveas(figure(plotnum),'wow2.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30])
histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30])
legend('thd Left Channel', 'thd Right Channel')
ylabel('number of Records')
xlabel('thd [dB]')
title('THD in the 1kHz track')
saveas(figure(plotnum),'thd1.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30])
histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30])
legend('thd Left Channel', 'thd Right Channel')
ylabel('number of Records')
xlabel('thd [dB]')
title('THD in the 1kHz2 track')
saveas(figure(plotnum),'thd2.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL'),:),50,'BinLimits',[-50,0])
ylabel('number of Records')
xlabel('stereo bleed')
title('stereo bleed in the 1kHzL track')
saveas(figure(plotnum),'stereoL1.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR'),:),50,'BinLimits',[-50,0])
ylabel('number of Records')
xlabel('stereo bleed')
title('stereo bleed in the 1kHzR track')
saveas(figure(plotnum),'stereoR1.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL2'),:),50,'BinLimits',[-50,0])
ylabel('number of Records')
xlabel('stereo bleed')
title('stereo bleed in the 1kHzL2 track')
saveas(figure(plotnum),'stereoL2.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR2'),:),50,'BinLimits',[-50,0])
ylabel('number of Records')
xlabel('stereo bleed')
title('stereo bleed in the 1kHzR2 track')
saveas(figure(plotnum),'stereoR2.png')
