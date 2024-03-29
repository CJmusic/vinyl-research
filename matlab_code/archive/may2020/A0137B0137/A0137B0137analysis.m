clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

c = gray(20);


if ismac()
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/')
    ('/Volumes/AUDIOBANK/audio_files/A0137B0137/A0137B0137-AudioTable.csv');
    data_folder = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/A0137B0137/';
    % AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTable.csv') 
    % AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137-AudioTableApr14.csv');
    % AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137-AudioTableApr28.csv');
    AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137-AudioTableMay12.csv');


    SensorTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_SensorTable.csv');
    RecordTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_RecordNumbers.csv');
    AudioError = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_AudioStats.csv')
    SensorError = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SensorStats.csv')
end
if ispc()
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\')
    data_folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\';
    % AudioTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\A0137B0137-AudioTableApr28.csv');
    AudioTable = readtable('E:\audio_files\A0137B0137\A0137B0137-AudioTable.csv');


    SensorTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\A0137B0137_SensorTable.csv');
    RecordTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0137B0137\A0137B0137_RecordTable.csv');
    AudioError = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\A0000B0000_AudioStats.csv')
    SensorError = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\A0000B0000_SensorStats.csv')
end

pressruns = unique(AudioTable.pressing);
tracks = unique(AudioTable.track);

% head(AudioTable)
AudioTable.RecordID = erase(AudioTable.record,'.wav');
AudioTable.RecordID = erase(AudioTable.RecordID,'a');
AudioTable.RecordID = erase(AudioTable.RecordID,'b');
AudioTable.RecordID =  cellfun(@str2num, AudioTable.RecordID);

% AudioTable.RecordID = cell(AudioTable.RecordID);

% AudioStats = struct(pressruns)

StatNames2 = {'pressing', 'track', 'AvgRMS_L','StdRMS_L', 'AvgRMS_R', 'StdRMS_R','AvgClicks_L', 'StdevClicks_L', 'AvgClicks_R', 'StdevClicks_R', 'AvgWow_L', 'StdWow_L', 'AvgStereobleed', 'StdevStereobleed'};
SensorStats = cell(0,14);

StatNames = {'pressing', 'track', 'AvgRMS_L','StdRMS_L', 'AvgRMS_R', 'StdRMS_R', 'AvgA_L', 'StdA_L', 'AvgA_R', 'StdA_R','AvgClicks_L', 'StdevClicks_L', 'AvgClicks_R', 'StdevClicks_R', 'AvgWow_L', 'StdWow_L', 'AvgWow_R', 'StdWow_R', 'AvgStereobleed', 'StdevStereobleed'};
AudioStats = cell(0,20);


for j = (1:length(tracks))
    for i = (1:length(pressruns))
        RMS_L = struct2table(datastats(AudioTable.RMS_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
        RMS_R = struct2table(datastats(AudioTable.RMS_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
        A_L = struct2table(datastats(AudioTable.A_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
        A_R = struct2table(datastats(AudioTable.A_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
        clicks_L = struct2table(datastats(AudioTable.clicks_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
        clicks_R = struct2table(datastats(AudioTable.clicks_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
        wow_L= struct2table(datastats(AudioTable.wow_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
        wow_R= struct2table(datastats(AudioTable.wow_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
        stereo_bleed = struct2table(datastats(AudioTable.stereo_bleed(strcmp(AudioTable.pressing,pressruns{i}),:)));
        

        

        AudioStats = [AudioStats ; cell2table({pressruns{i}, tracks{j}, RMS_L.mean, RMS_L.std, RMS_R.mean, RMS_R.std,A_L.mean, A_L.std, A_R.mean, A_R.std,clicks_L.mean, clicks_L.std, clicks_R.mean, clicks_R.std, wow_L.mean, wow_L.std, wow_R.mean, wow_R.std, stereo_bleed.mean, stereo_bleed.std})];
      
    end
end

AudioStats.Properties.VariableNames{'Var1'}='pressrun';
AudioStats.Properties.VariableNames{'Var2'}='track';
AudioStats.Properties.VariableNames{'Var3'}='AvgRMS_L';
AudioStats.Properties.VariableNames{'Var4'}='StdRMS_L';
AudioStats.Properties.VariableNames{'Var5'}='AvgRMS_R';
AudioStats.Properties.VariableNames{'Var6'}='StdRMS_R';
AudioStats.Properties.VariableNames{'Var7'}='AvgA_L';
AudioStats.Properties.VariableNames{'Var8'}='StdA_L';
AudioStats.Properties.VariableNames{'Var9'}='AvgA_R';
AudioStats.Properties.VariableNames{'Var10'}='StdA_R';
AudioStats.Properties.VariableNames{'Var11'}='AvgClicks_L';
AudioStats.Properties.VariableNames{'Var12'}='StdevClicks_L';
AudioStats.Properties.VariableNames{'Var13'}='AvgClicks_R';
AudioStats.Properties.VariableNames{'Var14'}='StdevClicks_R';
AudioStats.Properties.VariableNames{'Var15'}='AvgWow_L';
AudioStats.Properties.VariableNames{'Var16'}='StdWow_L';
AudioStats.Properties.VariableNames{'Var17'}='AvgWow_R';
AudioStats.Properties.VariableNames{'Var18'}='StdWow_R';
AudioStats.Properties.VariableNames{'Var19'}='AvgStereobleed';
AudioStats.Properties.VariableNames{'Var20'}='StdevStereobleed';



% AudioStats
Tbl = outerjoin(RecordTable, SensorTable);%,'VariableNames', 'RecordNumber')%;, SensorTable)
% Tbl = innerjoin(RecordTable, SensorTable);%,'VariableNames', 'RecordNumber')%;, SensorTable)
% Tbl
Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
writetable(Tbl,'Tbl1.csv')

head(Tbl)
head(AudioTable)

Tbl = outerjoin(Tbl, AudioTable);%, 'VariableNames', 'RecordNumber')%;, SensorTable)
% Tbl = innerjoin(Tbl, AudioTable);%, 'VariableNames', 'RecordNumber')%;, SensorTable)
% Tbl
Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
writetable(Tbl,'Tbl2.csv')

Tbl = outerjoin(Tbl, AudioStats);
% Tbl = innerjoin(Tbl, AudioStats);
% writetable(Tbl,'Tbl.csv')
writetable(Tbl,'Tbl3.csv')

Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
Tbl.Properties.VariableNames([33]) = {'track'};


% Tbl = Tbl(strcmp(Tbl.side,'a'),:);
writetable(Tbl,'Tbl4.csv')

Tbl = Tbl(strcmp(Tbl.side,'a'),:);

for i = (1:length(Tbl.Properties.VariableNames))
    disp(string(Tbl.Properties.VariableNames(i)));
end

AudioError.ste = AudioError.std/sqrt(5);
AudioError(strcmp(AudioError.track,'transition') & strcmp(AudioError.measurement,'RMS_L'),:);
SensorError.ste = SensorError.std/sqrt(5);
% SensorError(strcmp(SensorError.track,'') & strcmp(SensorError.measurement,'RMS_L'),:)


plotnum = 0;



plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;

%% PULL ERROR VALUE
tb1 = AudioError(strcmp(AudioError.track,'quiet2'),:)
tb1 = AudioError(strcmp(tb1.measurement,'RMS_L'),:)
% tb1.ste
figure(plotnum); grid on; hold on;
H = bar([Tbl.AvgRMS_L(strcmp(AudioStats.track,'quiet2'),:), Tbl.AvgRMS_R(strcmp(AudioStats.track,'quiet2'),:)], 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet2'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgRMS_L(strcmp(AudioStats.track,'quiet2'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';
numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet2'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgRMS_R(strcmp(AudioStats.track,'quiet2'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';

set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)


legend('RMS Left Channel', 'RMS Right Channel')
xlabel('number of records')
ylabel('RMS level [dB]')
title('RMS noise in quiet track')
saveas(figure(plotnum),'RMSquiet.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;

%% PULL ERROR VALUE
% tb1 = AudioError(strcmp(AudioError.track,'quiet2'),:)
% tb1 = AudioError(strcmp(tb1.measurement,'A_L'),:)
% tb1.ste
figure(plotnum); grid on; hold on;
H = bar([Tbl.AvgA_L(strcmp(AudioStats.track,'quiet2'),:), Tbl.AvgA_R(strcmp(AudioStats.track,'quiet2'),:)], 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
% numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet2'),:)))
% er = errorbar(double(categorical(pressruns)),AudioStats.AvgRMS_L(strcmp(AudioStats.track,'quiet2'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
% er.LineStyle = 'none';
% numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet2'),:)))
% er = errorbar(double(categorical(pressruns)),AudioStats.AvgRMS_R(strcmp(AudioStats.track,'quiet2'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
% er.LineStyle = 'none';

set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)


legend('RMS Left Channel', 'RMS Right Channel')
xlabel('number of records')
ylabel('RMS level [dB]')
title('A weighted RMS noise in quiet track')
saveas(figure(plotnum),'ARMSquiet.png')


plotnum = plotnum + 1;

figure(plotnum); grid on; hold on;
H = bar([Tbl.AvgWow_L(strcmp(AudioStats.track,'3150Hz'),:), Tbl.AvgWow_R(strcmp(AudioStats.track,'3150Hz'),:)], 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
numcats = (height(AudioStats(strcmp(AudioStats.track,'3150Hz'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgWow_L(strcmp(AudioStats.track,'3150Hz'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';
numcats = (height(AudioStats(strcmp(AudioStats.track,'3150Hz'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgWow_R(strcmp(AudioStats.track,'3150Hz'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';

set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)


legend('Wow Left Channel', 'Wow Right Channel')
xlabel('number of records')
ylabel('Wow [Hz]')
title('Avg Wow in 3150 Hz track')
saveas(figure(plotnum),'WowAvg.png')

plotnum = plotnum + 1;

figure(plotnum); grid on; hold on;
H = bar([Tbl.StdWow_L(strcmp(AudioStats.track,'3150Hz'),:), Tbl.StdWow_R(strcmp(AudioStats.track,'3150Hz'),:)], 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [.9 .9 .9];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
numcats = (height(AudioStats(strcmp(AudioStats.track,'3150Hz'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.StdWow_L(strcmp(AudioStats.track,'3150Hz'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';
numcats = (height(AudioStats(strcmp(AudioStats.track,'3150Hz'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.StdWow_R(strcmp(AudioStats.track,'3150Hz'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';

%% SET X VALUES
% xticklabels(pressruns)
% xtickangle(45)
% er.LineStyle = 'none';

set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)


legend('Wow Left Channel', 'Wow Right Channel')
xlabel('number of records')
ylabel('Wow [Hz]')
title('Std Wow in 3150 Hz track')
saveas(figure(plotnum),'WowStd.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;

% H = bar(Tbl.AvgStereobleed(strcmp(AudioStats.track,'1kHz2'),:), 'LineWidth', 2)
% H(1).FaceColor = [0.6 0.6 0.6];
% % H(2).FaceColor = [.9 .9 .9];
% numcats = (height(AudioStats(strcmp(AudioStats.track,'1kHz2'),:)))
% er = errorbar(double(categorical(pressruns)),AudioStats.AvgStereobleed(strcmp(AudioStats.track,'1kHz2'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);    
% er.LineStyle = 'none';
% grid on; hold on;

% xticklabels(pressruns)
% xtickangle(45)

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;

H = bar(Tbl.AvgStereobleed(strcmp(AudioStats.track,'1kHzL'),:), 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
numcats = (height(AudioStats(strcmp(AudioStats.track,'1kHzL'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgStereobleed(strcmp(AudioStats.track,'1kHzL'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';

%% SET X VALUES
% xticklabels(pressruns)
% xtickangle(45)
% er.LineStyle = 'none';

set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)

xlabel('number of records')
ylabel('stereo bleed [dB]')
title('Stereo bleed in  1HzL track')
saveas(figure(plotnum),'Stereobleedin1HzLtrack.png')

plotnum = plotnum + 1;
figure(plotnum); grid on; hold on;

H = bar([Tbl.AvgClicks_L(strcmp(AudioStats.track,'quiet'),:), Tbl.AvgClicks_R(strcmp(AudioStats.track,'quiet'),:)], 'LineWidth', 2)
H(1).FaceColor = [0.6 0.6 0.6];
H(2).FaceColor = [0.9 0.9 0.9];

%% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgClicks_L(strcmp(AudioStats.track,'quiet'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';
numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet'),:)))
er = errorbar(double(categorical(pressruns)),AudioStats.AvgClicks_R(strcmp(AudioStats.track,'quiet'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
er.LineStyle = 'none';

%% SET X VALUES
xticklabels(pressruns)
xtickangle(45)
er.LineStyle = 'none';

% set(gca,'xticklabel',pressruns)
set(gca,'xticklabel',pressruns)
ax=gca;
ax.FontSize=8;
ax.XTick = (1:length(pressruns))   %THIS WAY, YOU SET HOW MANY XTICKS YOU WANT FOR YOUR XTICKLABELS
xtickangle(45)

xlabel('number of records')
ylabel('number of clicks')
title('Number of clicks in quiet tracks')
saveas(figure(plotnum),'clicksquiet.png')


% Tbl.A_L
% Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet'))
plotnum = plotnum + 1;
figure(plotnum);  
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),'ko')
grid on; hold on;
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),'kx')
legend('left channel', 'right channel')
title('PressingNumber vs minMouldSteamIn')
saveas(figure(plotnum),'PressingNumber vs minMouldSteamIn.png')

plotnum = plotnum + 1;
figure(plotnum);  
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.RMS_L(strcmp(Tbl.track,'quiet')),'ko')
grid on; hold on;
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.RMS_R(strcmp(Tbl.track,'quiet')),'kx')
legend('left channel', 'right channel')
title('PressingNumber vs RMS')
saveas(figure(plotnum),'PressingNumber vs RMS.png')


plotnum = plotnum + 1;
figure(plotnum);  
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.A_L(strcmp(Tbl.track,'quiet')),'ko')
grid on; hold on;
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.A_R(strcmp(Tbl.track,'quiet')),'kx')
title('PressingNumber vs ARMS')
legend('left channel', 'right channel')
saveas(figure(plotnum),'PressingNumber vs ARMS.png')

plotnum = plotnum + 1;
figure(plotnum);  
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.clicks_L(strcmp(Tbl.track,'quiet')),'ko')
grid on; hold on;
plot(Tbl.PressingNumber(strcmp(Tbl.track,'quiet')), Tbl.clicks_R(strcmp(Tbl.track,'quiet')),'kx')
title('PressingNumber vs Clicks')
legend('left channel', 'right channel')
saveas(figure(plotnum),'PressingNumber vs Clicks.png')


plotnum = plotnum + 1;
figure(plotnum);  
plot(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.RMS_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
plot(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.RMS_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('RMS vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs RMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title('Wow vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs Wow.png')


plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'ko')
grid on; hold on;
title('stereo bleed vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs stereo_bleed.png')

plotnum = plotnum + 1;
figure(plotnum);  
scatter(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
scatter(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('ARMS vs maxPressForce')
saveas(figure(plotnum),'maxPressForce vs ARMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title(' Wow vs maxPressForce')
saveas(figure(plotnum),'maxPressForce vs Wow.png')



plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.maxPressForce_Ton(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'ko')
grid on; hold on;
title('stereo bleed vs maxPressForce')
saveas(figure(plotnum),'maxPressForce vs stereo_bleed.png')

plotnum = plotnum + 1;
figure(plotnum);  
scatter(Tbl.minExtruderMeltTemp_F(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
scatter(Tbl.minExtruderMeltTemp_F(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('RMS vs minExtruderMeltTemp')
saveas(figure(plotnum),'minExtruderMeltTemp vs RMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minExtruderMeltTemp_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.minExtruderMeltTemp_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title('Wow vs minExtruderMeltTemp')
saveas(figure(plotnum),'minExtruderMeltTemp vs Wow.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minExtruderMeltTemp_F(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'ko')
grid on; hold on;
title('stereo bleed vs minExtruderMeltTemp')
saveas(figure(plotnum),'minExtruderMeltTemp vs stereo_bleed.png')

plotnum = plotnum + 1;
figure(plotnum);  
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('RMS vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs RMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title('Wow vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs Wow.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'ko')
grid on; hold on;
legend('left channel', 'right channel')
title('stereo bleed vs minMouldSteamIn')
saveas(figure(plotnum),'minMouldSteamIn vs stereo_bleed.png')

plotnum = plotnum + 1;
figure(plotnum);  
scatter(Tbl.maxMouldSteamOutTop_F(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
scatter(Tbl.maxMouldSteamOutTop_F(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('RMS vs maxMouldSteamOutTop')
saveas(figure(plotnum),'maxMouldSteamOutTop vs RMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title('Wow vs maxMouldSteamOutBottom')
saveas(figure(plotnum),'maxMouldSteamOutBottom vs Wow.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.maxMouldSteamOutBottom_F(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'ko')
grid on; hold on;
title('stereo bleed vs maxMouldSteamOutBottom')
saveas(figure(plotnum),'maxMouldSteamOutBottom vs stereo bleed.png')


plotnum = plotnum + 1;
figure(plotnum);  
scatter(Tbl.minMouldSteamOutBottom_F(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),'ko')
grid on; hold on;
scatter(Tbl.minMouldSteamOutBottom_F(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'kx')
legend('left channel', 'right channel')
title('RMS vs minMouldSteamOutBottom')
saveas(figure(plotnum),'minMouldSteamOutBottom vs RMS.png')

plotnum = plotnum + 1;
figure(plotnum); 
scatter(Tbl.minMouldSteamOutBottom_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:),'ko')
grid on; hold on;
scatter(Tbl.minMouldSteamOutBottom_F(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'kx')
legend('left channel', 'right channel')
title(' Wow vs minMouldSteamOutBottom')
saveas(figure(plotnum),'minMouldSteamOutBottom vs Wow.png')

plotnum = plotnum + 1;
plot_scatter(plotnum, Tbl.minMouldSteamOutBottom_F(strcmp(Tbl.track,'1kHzL')),Tbl.stereo_bleed(strcmp(Tbl.track,'1kHzL'),:),'minMouldSteamOutBottom vs stereo bleed')

plotnum = plotnum + 1;
plot_scatter2(plotnum,Tbl.PressingNumber(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet'),:),Tbl.PressingNumber(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'RMS vs Pressing number')

plotnum = plotnum + 1;
plot_scatter2(plotnum,Tbl.PressingNumber(strcmp(Tbl.track,'quiet')),Tbl.clicks_L(strcmp(Tbl.track,'quiet')), Tbl.PressingNumber(strcmp(Tbl.track,'quiet')),Tbl.clicks_R(strcmp(Tbl.track,'quiet'),:),'Clicks vs Pressing number')

plotnum = plotnum + 1;
figure(plotnum);  
plot_scatter2(plotnum,Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.clicks_L(strcmp(Tbl.track,'quiet'),:),Tbl.minMouldSteamIn_F(strcmp(Tbl.track,'quiet')),Tbl.clicks_R(strcmp(Tbl.track,'quiet'),:),'clicks vs minMouldSteamIn')


plotnum = plotnum + 1;
plot_scatter2(plotnum,Tbl.maxExtruderBarrelZone3Temp_F(strcmp(Tbl.track,'quiet')),Tbl.A_L(strcmp(Tbl.track,'quiet')),Tbl.maxExtruderBarrelZone3Temp_F(strcmp(Tbl.track,'quiet')),Tbl.A_R(strcmp(Tbl.track,'quiet'),:),'maxExtruderBarrelZone3Temp_Fvs RMS.png')

plotnum = plotnum + 1;
plot_scatter2(plotnum,Tbl.PressingNumber(strcmp(Tbl.track,'3150Hz')),Tbl.wow_L(strcmp(Tbl.track,'3150Hz'),:), Tbl.PressingNumber(strcmp(Tbl.track,'3150Hz')),Tbl.wow_R(strcmp(Tbl.track,'3150Hz'),:),'Wow vs PressingNumber')

function plot_scatter2(plotnum, x1, y1, x2, y2, titlestring)
    figure(plotnum);  
    scatter(x1,y1,'ko')
    grid on; hold on;
    scatter(x2,y2,'kx')
    legend('left channel', 'right channel')
    title(titlestring)
    saveas(figure(plotnum), titlestring)
end

function plot_scatter(plotnum, x1, y1, titlestring)
    figure(plotnum);  
    scatter(x1,y1,'ko')
    grid on; hold on;
    legend('left channel', 'right channel')
    title(titlestring)
    saveas(figure(plotnum), titlestring)
end


% function plot_histogram(Tbl,plotnum,track,measurement)

%     % plotnum = plotnum + 1;
%     figure(plotnum); grid on; hold on;

%     %% PULL ERROR VALUE
%     tb1 = AudioError(strcmp(AudioError.track,track,:))
%     tb1 = AudioError(strcmp(tb1.measurement,measurement),:)
%     tb1.ste

%     %% PLOT BAR GRAPH WITH COLORS 
%     H = bar([Tbl.AvgRMS_L(strcmp(AudioStats.track,'quiet'),:), 'LineWidth', 2)
%     H(1).FaceColor = [0.6 0.6 0.6];

%     %% TWO SETS OF ERROR BARS FOR L AND R CHANNELS
%     numcats = (height(AudioStats(strcmp(AudioStats.track,'quiet'),:)));
%     er = errorbar(double(categorical(pressruns)),AudioStats.AvgRMS_L(strcmp(AudioStats.track,'quiet'),:), ones(numcats,1)*tb1.ste, ones(numcats,1)*tb1.ste,'k','LineWidth', 2);   
%     er.LineStyle = 'none';

%     %% SET X VALUES
%     xticklabels(pressruns)
%     xtickangle(45)

%     legend('RMS Left Channel', 'RMS Right Channel')
%     xlabel('pressing')
%     ylabel('RMS level [dB]')
%     title('RMS noise in quiet track')
%     saveas(figure(plotnum),'RMSquiet.png')

% end
% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(SensorTable.PressingNumber,SensorTable.minMouldSteamOutBottom_F,'ro')

% % plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% plot(Tbl.PressingNumber,Tbl.minMouldSteamOutBottom_F,'b.')


% grid on; hold on;


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet'),:),50,'BinLimits',[0,500])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of records')
% xlabel('number of clicks')
% title('number of clicks in quiet track')
% xlim([0,500])
% saveas(figure(plotnum),'clicksquiet.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30])
% histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[-50,-30])
% legend('RMS Left Channel', 'RMS Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in quiet2 track')
% xlim([-50,-30])
% saveas(figure(plotnum),'RMSquiet2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'quiet2'),:),50,'BinLimits',[0,500])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in quiet2 track')
% % xlim([0,500])
% saveas(figure(plotnum),'clicksquiet2.png')



% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.RMS_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30])
% histogram(AudioTable.RMS_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[-50,-30])
% legend('RMS Left Channel', 'RMS Right Channel')
% ylabel('number of records')
% xlabel('RMS level [dB]')
% title('RMS noise in transition track')
% xlim([-50,-30])
% saveas(figure(plotnum),'rmstransition.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.clicks_L(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000])
% histogram(AudioTable.clicks_R(strcmp(AudioTable.track,'transition'),:),50,'BinLimits',[0,3000])
% legend('Num Clicks Left Channel', 'Num Clicks Right Channel')
% ylabel('number of Records')
% xlabel('number of clicks')
% title('num clicks in transition track')
% xlim([0,3000])
% saveas(figure(plotnum),'clickstransition.png')


% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50])
% histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz'),:),50,'BinLimits',[10,50])
% legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
% ylabel('number of Records')
% xlabel('wow and flutter')
% title('Wow and flutter in 3150Hz track')
% saveas(figure(plotnum),'wow1.png')



% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.wow_R(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50])
% histogram(AudioTable.wow_L(strcmp(AudioTable.track,'3150Hz2'),:),50,'BinLimits',[10,50])
% legend('wow and flutter Left Channel', 'wow and flutter Right Channel')
% ylabel('number of Records')
% xlabel('wow and flutter')
% title('Wow and flutter in 3150Hz2 track')
% saveas(figure(plotnum),'wow2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30])
% histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz'),:),50,'BinLimits',[-60,-30])
% legend('thd Left Channel', 'thd Right Channel')
% ylabel('number of Records')
% xlabel('thd [dB]')
% title('THD in the 1kHz track')
% saveas(figure(plotnum),'thd1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.THD_L(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30])
% histogram(AudioTable.THD_R(strcmp(AudioTable.track,'1kHz2'),:),50,'BinLimits',[-60,-30])
% legend('thd Left Channel', 'thd Right Channel')
% ylabel('number of Records')
% xlabel('thd [dB]')
% title('THD in the 1kHz2 track')
% saveas(figure(plotnum),'thd2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL'),:),50,'BinLimits',[-50,0])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzL track')
% saveas(figure(plotnum),'stereoL1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR'),:),50,'BinLimits',[-50,0])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzR track')
% saveas(figure(plotnum),'stereoR1.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzL2'),:),50,'BinLimits',[-50,0])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzL2 track')
% saveas(figure(plotnum),'stereoL2.png')

% plotnum = plotnum + 1;
% figure(plotnum); grid on; hold on;
% histogram(AudioTable.stereo_bleed(strcmp(AudioTable.track,'1kHzR2'),:),50,'BinLimits',[-50,0])
% ylabel('number of Records')
% xlabel('stereo bleed')
% title('stereo bleed in the 1kHzR2 track')
% saveas(figure(plotnum),'stereoR2.png')
