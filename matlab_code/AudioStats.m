%% DataProcess.m 
% inputs: AudioData, SensorValues 
% outputs: AudioStats

close all

addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')

AudioFile = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/A0000B0000-data.csv')
SensorFile = ('')

% dataFolder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'

% AudioFile = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_bin\A0000B0000\A0000B0000-data.csv')


audio_byN = readtable(AudioFile);
for i = (1:length(audio_byN.Properties.VariableNames))
    disp(audio_byN.Properties.VariableNames(i))
end
audio_byN

sensor_names = {'recordNumber',	
'id',	
'RecordTimeStamp',
'PressPosition_Inches',
'PressForce_Ton',
'MouldSteamIn_PSI',
'MouldSteamIn_F',
'MouldSteamOutTop_F',
'MouldSteamOutBottom_F',
'ExtruderFeedthroatTemp_F',
'ExtruderBarrelZone1Temp_F',	
'ExtruderBarrelZone2Temp_F',
'ExtruderBarrelZone3Temp_F',
'ExtruderDieZoneTemp_F',
'ExtruderPremouldTemp_F',
'ExtruderMeltTemp_F'};

sensor_variables = {'RecordNumber',
'maxPressForce_Ton',
'minPressForce_Ton',
'maxMouldSteamIn_PSI',
'minMouldSteamIn_PSI',
'maxMouldSteamIn_F',
'minMouldSteamIn_F',
'maxMouldSteamOutTop_F',
'minMouldSteamOutTop_F',
'maxMouldSteamOutBottom_F',
'minMouldSteamOutBottom_F',
'maxExtruderFeedthroatTemp_F',
'minExtruderFeedthroatTemp_F',
'maxExtruderBarrelZone1Temp_F',
'minExtruderBarrelZone1Temp_F',
'maxExtruderBarrelZone2Temp_F',
'minExtruderBarrelZone2Temp_F',
'maxExtruderBarrelZone3Temp_F',
'minExtruderBarrelZone3Temp_F',
'maxExtruderDieZoneTemp_F',
'minExtruderDieZoneTemp_F',
'maxExtruderPremouldTemp_F',
'minExtruderPremouldTemp_F',
'maxExtruderMeltTemp_F',
'minExtruderMeltTemp_F'};


%now gotta do the sensor values 
sensor_raw = readtable(SensorFile)
for i = (1:max(sensor_raw.RecordNumber))
    sensor_table.RecordNumber(i) = sensor_raw.RecordNumber;
    sensor_table.maxPressForce_Ton(i) = max(sensor_raw.PressForce_Ton);
    sensor_table.minPressForce_Ton(i) = min(sensor_raw.PressForce_Ton);
    sensor_table.maxMouldSteamIn_PSI(i) = max(sensor_raw.MouldSteamIn_PSI);
    sensor_table.minMouldSteamIn_PSI(i) = min(sensor_raw.MouldSteamIn_PSI);
    sensor_table.maxMouldSteamIn_F(i) = max(sensor_raw.MouldSteamIn_F);
    sensor_table.minMouldSteamIn_F(i) = min(sensor_raw.MouldSteamIn_F);
    sensor_table.maxMouldSteamOutTop_F(i) = max(sensor_raw.MouldSteamOutTop_F);
    sensor_table.minMouldSteamOutTop_F(i) = min(sensor_raw.MouldSteamOutTop_F);
    sensor_table.maxMouldSteamOutBottom_F(i) = max(sensor_raw.MouldSteamOutBottom_F);
    sensor_table.minMouldSteamOutBottom_F(i) = min(sensor_raw.MouldSteamOutBottom_F);
    sensor_table.maxExtruderFeedthroatTemp_F(i) = max(sensor_raw.ExtruderFeedthroatTemp_F);
    sensor_table.minExtruderFeedthroatTemp_F(i) = min(sensor_raw.ExtruderFeedthroatTemp_F);
    sensor_table.maxExtruderBarrelZone1Temp_F(i) = max(sensor_raw.ExtruderBarrelZone1Temp_F);
    sensor_table.minExtruderBarrelZone1Temp_F(i) = min(sensor_raw.ExtruderBarrelZone1Temp_F);
    sensor_table.maxExtruderBarrelZone2Temp_F(i) = max(sensor_raw.ExtruderBarrelZone2Temp_F);
    sensor_table.minExtruderBarrelZone2Temp_F(i) = min(sensor_raw.ExtruderBarrelZone2Temp_F);
    sensor_table.maxExtruderBarrelZone3Temp_F(i) = max(sensor_raw.ExtruderBarrelZone3Temp_F);
    sensor_table.minExtruderBarrelZone3Temp_F(i) = min(sensor_raw.ExtruderBarrelZone3Temp_F);
    sensor_table.maxExtruderDieZoneTemp_F(i) = max(sensor_raw.ExtruderDieZoneTemp_F);
    sensor_table.minExtruderDieZoneTemp_F(i) = min(sensor_raw.ExtruderDieZoneTemp_F);
    sensor_table.maxExtruderPremouldTemp_F(i) = max(sensor_raw.ExtruderPremouldTemp_F);
    sensor_table.minExtruderPremouldTemp_F(i) = min(sensor_raw.ExtruderPremouldTemp_F);
    sensor_table.maxExtruderMeltTemp_F(i) = max(sensor_raw.ExtruderMeltTemp_F);
    sensor_table.minExtruderMeltTemp_F(i) = min(sensor_raw.ExtruderMeltTemp_F);
    
end 




%make a table with columns : 
% track , measurement , mean , median .... , std

col_names = {'track','measurement','max', 'min', 'mean', 'median', 'range', 'std'};
% AudioStats = table('VariableNames', col_names)
AudioStats  = cell2table(cell(0,8), 'VariableNames', col_names);
AudioStats
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
% datastats([1,2,3,1,2,6])

for i = (1:length(measurements))
    for j = (1:length(signal_names))
        % intTable = byN(strcmp(byN(:,9),signal_names{j}),:);
        int_stats = datastats(getData(byN, signal_names{j}, measurements{i}));

        intTable = cell2table({signal_names{j}, measurements{i},int_stats.max,int_stats.min,int_stats.mean,int_stats.median,int_stats.range,int_stats.std},'VariableNames', col_names);

        % intTable
        AudioStats = [AudioStats ; intTable];
    end
end
AudioStats

writetable(AudioStats,strcat(dataFolder,'A0000B0000_AudioStats.csv'));


% %DO TOTAL CLICKS 
% data_return = (table2array(byN(:,'clicks_L')));
% data_return


function data_return = getData(Tbl, track, param) 
    data_return = (table2array(Tbl(strcmp(Tbl.track, track), param)));
    if strcmp(class(data_return),'cell')
        data_return = str2double(data_return);
    end
end 