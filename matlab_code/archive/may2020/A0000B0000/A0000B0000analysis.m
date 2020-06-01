clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
% set(0,'DefaultColorOrder','gray')
% set(0,'DefaultFigureColormap',feval('gray'));
% colormap('gray');


if ismac()
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/')
    data_folder = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/A0000B0000/'
    AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTableMay26.csv') 
    AudioStats = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_AudioStats.csv')
    SensorTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SensorTable.csv')
end
if ispc()
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\')
    data_folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'
    AudioTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\A0000B0000-AudioTable.csv');
end

%remove row 76 from SensorValues, press closed but no steam was pumped through 
SensorTable([76],:) = [];
AudioTable([76],:) = [];
AudioTable = sortrows(AudioTable,7);



AudioStats
AudioStats(strcmp(AudioStats.track,'transition') & strcmp(AudioStats.measurement,'RMS_L'),:)


plotnum = 0;

% %~~~~~~~~~~~~~HISTOGRAM PLOTS~~~~~~~~~~~~~~~~~%

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.5 0.5 0.5])
% legend('RMS Left Channel', 'RMS Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in quiet track')
% xlim([-50,-30])
% saveas(figure(plotnum),'RMSquiet.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500], 'facecolor',[0.5 0.5 0.5])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of records')
% xlabel('number of clicks')
% title('number of clicks in quiet track')
% xlim([0,500])
% saveas(figure(plotnum),'clicksquiet.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.5 0.5 0.5])
% legend('RMS Left Channel', 'RMS Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in quiet2 track')
% xlim([-50,-30])
% saveas(figure(plotnum),'RMSquiet2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500], 'facecolor',[0.5 0.5 0.5])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in quiet2 track')
% % xlim([0,500])
% saveas(figure(plotnum),'clicksquiet2.png')



% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30], 'facecolor',[0.5 0.5 0.5])
% legend('RMS Left Channel', 'RMS Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in transition track')
% xlim([-50,-30])
% saveas(figure(plotnum),'rmstransition.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000], 'facecolor',[0.5 0.5 0.5])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of Records')
% xlabel('number of clicks')
% title('num clicks in transition track')
% xlim([0,3000])
% saveas(figure(plotnum),'clickstransition.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50], 'facecolor',[0.5 0.5 0.5])
% legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
% ylabel('number of Records')
% xlabel('wow and flutter')
% title('Wow and flutter in 3150Hz track')
% saveas(figure(plotnum),'wow1.png')



% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50], 'facecolor',[0.5 0.5 0.5])
% legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
% ylabel('number of Records')
% xlabel('wow and flutter')
% title('Wow and flutter in 3150Hz2 track')
% saveas(figure(plotnum),'wow2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30], 'facecolor',[0.5 0.5 0.5])
% legend('thd Left Channel', 'thd Right Channel')
% ylabel('number of Records')
% xlabel('thd [dB]')
% title('THD in the 1kHz track')
% saveas(figure(plotnum),'thd1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30], 'facecolor',[0.3 0.3 0.3])
% histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30], 'facecolor',[0.5 0.5 0.5])
% legend('thd Left Channel', 'thd Right Channel')
% ylabel('number of Records')
% xlabel('thd [dB]')
% title('THD in the 1kHz2 track')
% saveas(figure(plotnum),'thd2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL'),:),50,'BinLimits',[-50,0], 'facecolor',[0.3 0.3 0.3])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzL track')
% saveas(figure(plotnum),'stereoL1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR'),:),50,'BinLimits',[-50,0], 'facecolor',[0.3 0.3 0.3])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzR track')
% saveas(figure(plotnum),'stereoR1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL2'),:),50,'BinLimits',[-50,0], 'facecolor',[0.3 0.3 0.3])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzL2 track')
% saveas(figure(plotnum),'stereoL2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR2'),:),50,'BinLimits',[-50,0], 'facecolor',[0.3 0.3 0.3])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzR2 track')
% saveas(figure(plotnum),'stereoR2.png')

% %~~~~~~~~~~~~~HISTOGRAM PLOTS END~~~~~~~~~~~~~~~~~%


% %~~~~~~~~~~~~~SENSOR TABLE PLOTS~~~~~~~~~~~~~~~~~%

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(SensorTable.recordNumber,SensorTable.maxMouldSteamIn_F,'k-o')
% ylabel('maxMouldSteamIn [F]')
% xlabel('record number')
% title('record number vs maxMouldSteamIn_F')
% saveas(figure(plotnum),'maxMouldSteamIn_F.png')

% % SensorTable

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(SensorTable.recordNumber,SensorTable.maxMouldSteamOutBottom_F,'k-o')
% ylabel('MouldSteamOutBottom_F [F]')
% xlabel('record number')
% title('record number vs MouldSteamOutBottom_F')
% saveas(figure(plotnum),'MouldSteamOutBottom_F.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(SensorTable.recordNumber,SensorTable.maxMouldSteamOutTop_F,'k-o')
% ylabel('MouldSteamOutTop_F [F]')
% xlabel('record number')
% title('record number vs MouldSteamOutTop_F')
% saveas(figure(plotnum),'MouldSteamOutTop_F.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(SensorTable.recordNumber, SensorTable.maxPressForce_Ton,'k-o')
% ylabel('maxPressForce_Ton')
% xlabel('record number')
% title('record number vs maxPressForce_Ton')
% saveas(figure(plotnum),'maxPressForce_Ton.png')


% %~~~~~~~~~~~~~SENSOR TABLE PLOTS END~~~~~~~~~~~~~%

% %~~~~~~~~~~~~~AUDIO TABLE PLOTS~~~~~~~~~~~~~%


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), AudioTable.RMS_L(strcmp(AudioTable.track,'quiet'),:),'ko')
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), AudioTable.RMS_R(strcmp(AudioTable.track,'quiet'),:),'kx')
ylabel('RMS [dB]')
xlabel('record number')
title('record number vs RMS in quiet track')
saveas(figure(plotnum),'RMSquiet.png')




plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), AudioTable.clicks_L(strcmp(AudioTable.track,'quiet'),:),'ko')
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), AudioTable.clicks_R(strcmp(AudioTable.track,'quiet'),:),'kx')
ylabel('RMS [dB]')
xlabel('record number')
title('record number vs RMS in quiet track')
saveas(figure(plotnum),'RMSquiet.png')






%~~~~~~~~~~~~~AUDIO TABLE PLOTS END~~~~~~~~~~~~~%
Tbl = SensorTable(ismember(SensorTable.recordNumber,AudioTable.record),:)
Tbl.Properties.VariableNames([1]) = {'record'};
Tbl = outerjoin(Tbl, AudioTable)
%~~~~~~~~~~~~~MIXED TABLE PLOTS~~~~~~~~~~~~~%
plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'quiet'),:),'ko')
plot(AudioTable.record(strcmp(AudioTable.track,'quiet'),:), Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'quiet'),:),'kx')
ylabel('RMS [dB]')
xlabel('record number')
title('record number vs Mould Steam Out in quiet track')
saveas(figure(plotnum),'mouldsteamout.png')


plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(Tbl.maxMouldSteamIn_F(strcmp(Tbl.track,'quiet'),:), Tbl.RMS_L(strcmp(Tbl.track,'quiet'),:),'ko')
ylabel('RMS [dB]')
xlabel('rmaxMouldSteamIn_F')
title('maxMouldSteamIn_F vs RMS in quiet track')
saveas(figure(plotnum),'maxMouldSteamIn_F.png')



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'quiet'),:), Tbl.RMS_L(strcmp(Tbl.track,'quiet'),:),'ko')
ylabel('RMS [dB]')
xlabel('maxPressForce_Ton')
title('maxPressForce_Ton vs RMS in quiet track')
saveas(figure(plotnum),'maxPressForce_Ton.png')



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'quiet'),:), Tbl.RMS_L(strcmp(Tbl.track,'quiet'),:),'ko')
ylabel('RMS [dB]')
xlabel('maxMouldSteamOutBottom_F')
title('maxMouldSteamOutBottom_F vs RMS in quiet track')
saveas(figure(plotnum),'maxMouldSteamOutBottom_F.png')



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;
plot(Tbl.maxExtruderBarrelZone2Temp_F(strcmp(Tbl.track,'quiet'),:), Tbl.RMS_L(strcmp(Tbl.track,'quiet'),:),'ko')
ylabel('RMS [dB]')
xlabel('maxExtruderBarrelZone2Temp_F')
title('maxExtruderBarrelZone2Temp_F vs RMS in quiet track')
saveas(figure(plotnum),'maxExtruderBarrelZone2Temp_F.png')



%~~~~~~~~~~~~~MIXED TABLE PLOTS END~~~~~~~~~~~~~%

