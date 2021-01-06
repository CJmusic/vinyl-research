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
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')
    data_folder = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/A0000B0000/'
    % AudioTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000-AudioTableOct26.csv') ;
    AudioTable = readtable('/Volumes/AUDIOBANK/audio_files/A0000B0000/A0000B0000-AudioTable.csv');
    SensorTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SensorTable.csv');
    % RecordTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SensorTable.csv');
    AudioStats = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_AudioStats.csv');

    SettingsTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_SettingsTable.csv');
    RecordTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/A0000B0000_RecordTable.csv');

    head(AudioTable)
    head(SensorTable)
    % head(RecordTable)
    head(AudioStats)


end
if ispc()
    addpath('D:\Code\vinyl-research\matlab_code\Common')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\')
    data_folder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'
    AudioTable = readtable('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\A0000B0000-AudioTable.csv');
end

%remove row 76 from SensorValues, press closed but no steam was pumped through 
SensorTable([76],:) = [];
AudioTable([76],:) = [];
AudioTable = sortrows(AudioTable,7);



plotnum = 0;

pressruns = unique(AudioTable.pressing);
tracks = unique(AudioTable.track);


% clean up AudioTable
% AudioTable.RecordID = erase(AudioTable.record,'.wav')
% AudioTable.RecordID = erase(AudioTable.RecordID,'a')
% AudioTable.RecordID = erase(AudioTable.RecordID,'b')
% AudioTable.RecordID =  cellfun(@str2num, AudioTable.RecordID)

% AudioTable.RecordID = AudioTable.record;
% AudioTable.RecordID = eraseBetween(AudioTable.record, 'a', 'v');
AudioTable.RecordID = extractBetween(AudioTable.record, 19, 21);
AudioTable.RecordID =  cellfun(@str2num, AudioTable.RecordID)



head(AudioTable)


%~~~~~~~~~~~~~~~ GENERATE AUDIOSTATS TABLE ~~~~~~~~~~~~~~~~~~%

    % StatNames2 = {'pressing', 'track', 'AvgRMS_L','StdRMS_L', 'AvgRMS_R', 'StdRMS_R','AvgClicks_L', 'StdevClicks_L', 'AvgClicks_R', 'StdevClicks_R', 'AvgWow_L', 'StdWow_L', 'AvgStereobleed', 'StdevStereobleed'};
    % SensorStats = cell(0,14);

    % StatNames = {'pressing', 'track', 'AvgRMS_L','StdRMS_L', 'AvgRMS_R', 'StdRMS_R', 'AvgA_L', 'StdA_L', 'AvgA_R', 'StdA_R','AvgClicks_L', 'StdevClicks_L', 'AvgClicks_R', 'StdevClicks_R', 'AvgWow_L', 'StdWow_L', 'AvgWow_R', 'StdWow_R', 'AvgStereobleed', 'StdevStereobleed'};
    % AudioStats = cell(0,20);


    % for j = (1:length(tracks))
    %     for i = (1:length(pressruns))
    %         RMS_L = struct2table(datastats(AudioTable.RMS_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         RMS_R = struct2table(datastats(AudioTable.RMS_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         A_L = struct2table(datastats(AudioTable.A_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         A_R = struct2table(datastats(AudioTable.A_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         clicks_L = struct2table(datastats(AudioTable.clicks_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         clicks_R = struct2table(datastats(AudioTable.clicks_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         wow_L= struct2table(datastats(AudioTable.wow_L(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         wow_R= struct2table(datastats(AudioTable.wow_R(strcmp(AudioTable.pressing,pressruns{i}),:)));
    %         stereo_bleed = struct2table(datastats(AudioTable.stereo_bleed(strcmp(AudioTable.pressing,pressruns{i}),:)));
            

            

    %         AudioStats = [AudioStats ; cell2table({pressruns{i}, tracks{j}, RMS_L.mean, RMS_L.std, RMS_R.mean, RMS_R.std,A_L.mean, A_L.std, A_R.mean, A_R.std,clicks_L.mean, clicks_L.std, clicks_R.mean, clicks_R.std, wow_L.mean, wow_L.std, wow_R.mean, wow_R.std, stereo_bleed.mean, stereo_bleed.std})];
        
    %     end
    % end

    % AudioStats.Properties.VariableNames{'Var1'}='pressrun';
    % AudioStats.Properties.VariableNames{'Var2'}='track';
    % AudioStats.Properties.VariableNames{'Var3'}='AvgRMS_L';
    % AudioStats.Properties.VariableNames{'Var4'}='StdRMS_L';
    % AudioStats.Properties.VariableNames{'Var5'}='AvgRMS_R';
    % AudioStats.Properties.VariableNames{'Var6'}='StdRMS_R';
    % AudioStats.Properties.VariableNames{'Var7'}='AvgA_L';
    % AudioStats.Properties.VariableNames{'Var8'}='StdA_L';
    % AudioStats.Properties.VariableNames{'Var9'}='AvgA_R';
    % AudioStats.Properties.VariableNames{'Var10'}='StdA_R';
    % AudioStats.Properties.VariableNames{'Var11'}='AvgClicks_L';
    % AudioStats.Properties.VariableNames{'Var12'}='StdevClicks_L';
    % AudioStats.Properties.VariableNames{'Var13'}='AvgClicks_R';
    % AudioStats.Properties.VariableNames{'Var14'}='StdevClicks_R';
    % AudioStats.Properties.VariableNames{'Var15'}='AvgWow_L';
    % AudioStats.Properties.VariableNames{'Var16'}='StdWow_L';
    % AudioStats.Properties.VariableNames{'Var17'}='AvgWow_R';
    % AudioStats.Properties.VariableNames{'Var18'}='StdWow_R';
    % AudioStats.Properties.VariableNames{'Var19'}='AvgStereobleed';
    % AudioStats.Properties.VariableNames{'Var20'}='StdevStereobleed';
%~~~~~~~~~~~~~~~ GENERATE AUDIOSTATS TABLE ENDS ~~~~~~~~~~~~~~~~~~%

%~~~~~~~~~~~~~~~ MERGING TABLES STARTS ~~~~~~~~~~~~~~%
    % disp('AudioStats')
    % head(AudioStats)
    disp('AudioTable')
    head(AudioTable)
    % disp('RecordTable')
    % head(RecordTable)
    disp('SensorTable')
    head(SensorTable)

    % join RecordTable and SensorTable
    % Tbl = outerjoin(RecordTable, SensorTable);%,'VariableNames', 'RecordNumber')%;, SensorTable)
    % Tbl.Properties.VariableNames([1]) = {'RecordID'};
    % writetable(Tbl,'Tbl1.csv')

    % head(Tbl)

    Tbl = SensorTable;
    Tbl.Properties.VariableNames([1]) = {'RecordID'};

    % join in AudioTable
    Tbl = outerjoin(Tbl, AudioTable, 'Keys', {'RecordID', 'RecordID'});%, 'VariableNames', 'RecordNumber')%;, SensorTable)
    Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
    writetable(Tbl,'Tbl.csv')

    head(Tbl)

    Tbl2 = outerjoin(SettingsTable, RecordTable, 'Keys', {'pressing', 'pressing'});%, 'VariableNames', 'RecordNumber')%;, SensorTable)
    writetable(Tbl,'Tbl2.csv')

    Tbl = outerjoin(Tbl,Tbl2,'Keys', {'PressingNumber', 'PressingNumber'});
    writetable(Tbl,'Tbl3.csv')
    % join in AudioStats
    % Tbl = outerjoin(Tbl, AudioStats);
    % writetable(Tbl,'Tbl3.csv')
    % Tbl.Properties.VariableNames([3]) = {'pressing'};


    Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
    % Tbl.Properties.VariableNames([33]) = {'track'};
    Tbl.pressing = [];
    Tbl.Properties.VariableNames([58]) = {'pressing'};
    

    Tbl.RecordID = Tbl.PressingNumber;
    Tbl.RecordNumber = Tbl.PressingNumber;
    Tbl.RecordID_AudioTable = [];
    Tbl.PressingNumber_Tbl2 = [];

    Tbl.pressing_RecordTable = [];
    Tbl.TimeStamp = [];
    Tbl.timestamp = [];
    writetable(Tbl,'Tbl4.csv')

    Tbl = movevars(Tbl, 'RecordID', 'After', 'PressingNumber');
    Tbl = movevars(Tbl, 'pressing', 'After', 'RecordID');
    Tbl = movevars(Tbl, 'RecordNumber', 'After', 'pressing');
    

    writetable(Tbl,'A0000B0000.csv')