%% DataProcess.m 
% inputs: AudioData, SensorValues 
% outputs: StatsTable

close all
addpath('E:\audio_files\A0000B0000\')
addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\');

% folder = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\rawpress_data\');
% dataFolder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
folder = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/rawpress_data/')
dataFolder = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/')

SensorTable = strcat(dataFolder, 'A0000B0000_PressTable.csv')

clc
disp('~~~~~~~~~~~~~~~~SENSOR PROCESS~~~~~~~~~~~~~~~~~~')


% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
% SensorTable = ('')

%make a table with columns : 
% track , measurement , mean , median .... , std
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

sensor_variables = {'recordNumber',
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
size(sensor_variables)
sensor_raw = readtable(SensorTable);
sensor_table = cell2table(cell(0,25), 'VariableNames', sensor_variables)


for i = (1:max(sensor_raw.recordNumber))
    int_table = sensor_raw(sensor_raw.recordNumber == i,:);

    % FinalTable = alldata(strcmp(alldata.HOMETOWNS, 'London') & strcmp(alldata.CLASSES', 'Class B2'), :)


    sensor_table.recordNumber(i) = i;
    sensor_table.maxPressForce_Ton(i) = max(int_table.PressForce_Ton);
    sensor_table.minPressForce_Ton(i) = min(int_table.PressForce_Ton);
    sensor_table.maxMouldSteamIn_PSI(i) = max(int_table.MouldSteamIn_PSI);
    sensor_table.minMouldSteamIn_PSI(i) = min(int_table.MouldSteamIn_PSI);
    sensor_table.maxMouldSteamIn_F(i) = max(int_table.MouldSteamIn_F);
    sensor_table.minMouldSteamIn_F(i) = min(int_table.MouldSteamIn_F);
    sensor_table.maxMouldSteamOutTop_F(i) = max(int_table.MouldSteamOutTop_F);
    sensor_table.minMouldSteamOutTop_F(i) = min(int_table.MouldSteamOutTop_F);
    sensor_table.maxMouldSteamOutBottom_F(i) = max(int_table.MouldSteamOutBottom_F);
    sensor_table.minMouldSteamOutBottom_F(i) = min(int_table.MouldSteamOutBottom_F);
    sensor_table.maxExtruderFeedthroatTemp_F(i) = max(int_table.ExtruderFeedthroatTemp_F);
    sensor_table.minExtruderFeedthroatTemp_F(i) = min(int_table.ExtruderFeedthroatTemp_F);
    sensor_table.maxExtruderBarrelZone1Temp_F(i) = max(int_table.ExtruderBarrelZone1Temp_F);
    sensor_table.minExtruderBarrelZone1Temp_F(i) = min(int_table.ExtruderBarrelZone1Temp_F);
    sensor_table.maxExtruderBarrelZone2Temp_F(i) = max(int_table.ExtruderBarrelZone2Temp_F);
    sensor_table.minExtruderBarrelZone2Temp_F(i) = min(int_table.ExtruderBarrelZone2Temp_F);
    sensor_table.maxExtruderBarrelZone3Temp_F(i) = max(int_table.ExtruderBarrelZone3Temp_F);
    sensor_table.minExtruderBarrelZone3Temp_F(i) = min(int_table.ExtruderBarrelZone3Temp_F);
    sensor_table.maxExtruderDieZoneTemp_F(i) = max(int_table.ExtruderDieZoneTemp_F);
    sensor_table.minExtruderDieZoneTemp_F(i) = min(int_table.ExtruderDieZoneTemp_F);
    sensor_table.maxExtruderPremouldTemp_F(i) = max(int_table.ExtruderPremouldTemp_F);
    sensor_table.minExtruderPremouldTemp_F(i) = min(int_table.ExtruderPremouldTemp_F);
    sensor_table.maxExtruderMeltTemp_F(i) = max(int_table.ExtruderMeltTemp_F);
    sensor_table.minExtruderMeltTemp_F(i) = min(int_table.ExtruderMeltTemp_F);
    
end 

sensor_table
writetable(sensor_table,strcat(dataFolder,'A0000B000_SensorTable.csv'));


col_names = {'sensor_measurement','max', 'min', 'mean', 'median', 'range', 'std'};

SensorStats  = cell2table(cell(0,7), 'VariableNames', col_names);
SensorStats

for i = (1:width(sensor_table))
        i
        % intTable = byN(strcmp(byN(:,9),sensor_variables{j}),:);
        sensor_table(:,i)
        int_stats = datastats(table2array(sensor_table(:,i)));

        intTable = cell2table({sensor_variables(i),int_stats.max,int_stats.min,int_stats.mean,int_stats.median,int_stats.range,int_stats.std},'VariableNames', col_names);

        % intTable
        SensorStats = [SensorStats ; intTable];
    end

SensorStats

writetable(SensorStats,strcat(dataFolder,'SensorStats.csv'));


%DO TOTAL CLICKS 
% data_return = (table2array(byN(:,'clicks_L')));
% data_return
% if strcmp(class(data_return),'cell')
    % data_return = str2double(data_return);
% end

% datastats(getData(byN, signal_names{j}, measurements{i}));





function data_return = getData(Tbl, rec_num, param) 
    %modified for sensor values
    disp('inside data_return')
    
    rec_num
    Tbl.recordNumber
    Tbl
    data_return = (table2array(Tbl(Tbl.recordNumber == rec_num, param)));
    if strcmp(class(data_return),'cell')
        data_return = str2double(data_return);
    end
end 