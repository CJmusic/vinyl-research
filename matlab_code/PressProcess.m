
addpath('E:\audio_files\A0000B0000\')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
folder = ('D:\OneDrive - University of Waterloo\Vinyl_Project\data\121918_A0000B0000\')

% folder = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/')
% rawFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/A0000B0000.csv')
% timeFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/oct10A0000B0000.csv')


% read in the press file and the time file 
% join the two arrays 

date_tag = '2012-12-19'


% find the cycle of the press 

%#####_SensorValues
% tracks the data from the sensors in the press via timestamps
selectedColumns ={'id',
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
'ExtruderMeltTemp_F',}

opts = detectImportOptions(strcat(folder,date_tag,'_SensorValues.csv'));
getvaropts(opts,selectedColumns);
opts = setvartype(opts,selectedColumns,'string');
opts.SelectedVariableNames = selectedColumns;

SensorValues = readtable(strcat(folder,date_tag,'_SensorValues.csv'),opts);

%#####_JobDetailsCurrent
% tracks cycle times



JobDetails = readtable(strcat(folder,date_tag,'_JobDetailsCurrent.csv'));


%#####_ADAPT_DATA
% tracks the press options changes via timestamps
ADAPT = readtable(strcat(folder,'121918_CZA0000B0000_ADAPT_DATA.xlsx'));

%#####_TimeStamps
% Recorded by hand timestamps
TimeStamps = readtable(strcat(folder,'2012-12-20_TimeStamps.xlsx'));
SensorValues.RecordTimeStamp = datetime(SensorValues.RecordTimeStamp, 'InputFormat', 'yyyy-MM-dd HH:mm')
TimeStamps.TimeStamp = datetime(TimeStamps.TimeStamp, 'InputFormat', 'MMMM d, yyyy hh:mm:ss')

for i = (1:length(TimeStamps.TimeStamp))
    [closestTimeStamp,closestIndex] = min(abs(SensorValues.RecordTimeStamp-TimeStamps.TimeStamp(i)));
    if isnan(closestTimeStamp)
        continue
    end

    % TimeStamps.TimeStamp(i)
    % SensorValues.RecordTimeStamp(closestIndex)
    % SensorValues.PressPosition_Inches(closestIndex)
    SensorValues.PressForce_Ton(closestIndex)

    
    if 

    % find the max or min PressForce in a nearby area, 
    %write needed values to table  

end


SensorValues.Properties.VariableNames;

% get the hand recorded timestamp
% identify a pressing cycle by press position (one open and close)
% -> need to keep in mind times when parameter's changed 

% RECORDID is added by this code 
