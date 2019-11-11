%% DataProcess.m 
% inputs: AudioData, SensorValues 
% outputs: StatsTable

close all

addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')

dataFile = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/A0000B0000-data.csv')
dataFolder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'

% dataFile = ('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\A0000B0000-data.csv')


byN = readtable(dataFile);
for i = (1:length(byN.Properties.VariableNames))
    disp(byN.Properties.VariableNames(i))
end
byN

%make a table with columns : 
% track , measurement , mean , median .... , std

col_names = {'track','measurement','max', 'min', 'mean', 'median', 'range', 'std'};
% statsTable = table('VariableNames', col_names)
statsTable  = cell2table(cell(0,8), 'VariableNames', col_names);
statsTable
% each row is the table filtered by track and then written to a new table as measurements 

measurements = {'normalization_L', 'normalization_R', 'RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};

signal_names = {'leadin',    % 1
'1kHz',      % 2
'10kHz',     % 3
'100Hz',     % 4
'sweep',     % 5
'quiet',     % 6
'3150Hz',    % 7 
'1kHzL',     % 8
'sweepL',    % 9
'1kHzR',     % 10
'sweepR',    % 11
'1kHzV',     % 12
'sweepV',    % 13
'transition',% 14
'1kHz2',     % 15
'10kHz2',    % 16
'100Hz2',    % 17
'freqsweep2',% 18
'quiet2',    % 19
'3150Hz2',   % 20
'1kHzL2',    % 21
'sweepL2',   % 22
'1kHzR2',    % 23
'sweepR2',   % 24
'1kHzV2',    % 25
'sweepV2',   % 26
'leadout'    % 27
};

for i = (1:length(measurements))
    for j = (1:length(signal_names))
        % intTable = byN(strcmp(byN(:,9),signal_names{j}),:);
        int_stats = datastats(getData(byN, signal_names{j}, measurements{i}));

        intTable = cell2table({signal_names{j}, measurements{i},int_stats.max,int_stats.min,int_stats.mean,int_stats.median,int_stats.range,int_stats.std},'VariableNames', col_names);

        % intTable
        statsTable = [statsTable ; intTable];
    end
end
statsTable

writetable(statsTable,strcat(dataFolder,'statsTable.csv'));


%DO TOTAL CLICKS 
data_return = (table2array(byN(:,'clicks_L')));
data_return
% if strcmp(class(data_return),'cell')
    % data_return = str2double(data_return);
% end

% datastats(getData(byN, signal_names{j}, measurements{i}));



%~~~~~ OLD CODE ~~~~~%
% figure(1); hold on; grid on;
% % scatter(getData(byN,'transition', 'PressForce_Ton'), getData(byN,'transition', 'clicks_L'))
% scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_L'))
% scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_R'))
% title('RMS levels in transition track')
% xlabel('record #')
% ylabel('RMS level [dB]')
% legend(['Left', 'Right'])

% figure(2); hold on; grid on;
% hist(getData(byN,'transition','RMS_L'),20)
% hist(getData(byN,'transition','RMS_R'),20)
% title('RMS level Histogram')
% xlabel('num records')
% title('RMS level [dB]')

% figure(4); hold on; grid on;
% % scatter(getData(byN,'transition', 'PressForce_Ton'), getData(byN,'transition', 'clicks_L'))
% scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_L'))
% scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_R'))
% title('RMS levels in transition track')
% xlabel('record #')
% ylabel('RMS level [dB]')
% legend(['Left', 'Right'])

% %% find the 10 records with the min and max noise levels

% byRMS = sortrows(byN);
% disp('Highest RMS_L')
% mins = getData(byRMS,'transition','record');
% maxs = getData(byRMS,'transition','record');

% mins(1:10)
% maxs(end-10:end)

% stats_RMS_L = datastats(getData(byRMS,'transition','RMS_L'))
% stats_RMS_R = datastats(getData(byRMS,'transition','RMS_R'))
% stats_clicks_R = datastats(getData(byRMS,'transition','clicks_R'))
% stats_clicks_L = datastats(getData(byRMS,'transition','clicks_L'))

% disp('standard errors')
% stats_RMS_R.std/sqrt(height(byN))
% stats_RMS_L.std/sqrt(height(byN))
%~~~~~ OLD CODE ENDS ~~~~~% 






function data_return = getData(Tbl, track, param) 
    data_return = (table2array(Tbl(strcmp(Tbl.track, track), param)));
    if strcmp(class(data_return),'cell')
        data_return = str2double(data_return);
    end
end 