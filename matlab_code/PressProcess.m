
addpath('E:\audio_files\A0000B0000\')
folder = ('D:\OneDrive - University of Waterloo\Vinyl_Project\data\121918_A0000B0000\')

% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
% folder = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/121918_A0000B0000/')



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
SensorValues = [SensorValues; readtable(strcat(folder,'2012-12-20_SensorValues.csv'),opts)];

%#####_JobDetailsCurrent
% tracks cycle times



JobDetails = readtable(strcat(folder,date_tag,'_JobDetailsCurrent.csv'));
JobDetails = [JobDetails; readtable(strcat(folder,'2012-12-20_JobDetailsCurrent.csv'))];


%#####_ADAPT_DATA
% tracks the press options changes via timestamps
ADAPT = readtable(strcat(folder,'121918_CZA0000B0000_ADAPT_DATA.xlsx'));
ADAPT = [ADAPT; readtable(strcat(folder,'122018_CZA0000B0000_ADAPT_DATA.xlsx'))];

%#####_TimeStamps
% Recorded by hand timestamps
TimeStamps = readtable(strcat(folder,'2012-12-20_TimeStamps.xlsx'));
SensorValues.RecordTimeStamp = datetime(SensorValues.RecordTimeStamp, 'InputFormat', 'yyyy-MM-dd HH:mm');6

TimeStamps.TimeStamp = datetime(TimeStamps.TimeStamp, 'InputFormat', 'MMMM d, yyyy hh:mm:ss');

for i = (1:length(TimeStamps.TimeStamp))
    % if TimeStamps.TimeStamp(i) > datetime(date_tag)
    %     date_tag
    %     continue
    % end
    % disp('FOUND MATCH')
    [closestTimeStamp,closestIndex] = min(abs(SensorValues.RecordTimeStamp-TimeStamps.TimeStamp(i)));
    if isnan(closestTimeStamp)
        continue
    end
    disp(strcat(string(SensorValues.RecordTimeStamp(closestIndex)),'.....',string(TimeStamps.TimeStamp(i))))

    
    % SensorValues.RecordTimeStamp(closestIndex)
    % SensorValues.PressPosition_Inches(closestIndex)
    % SensorValues.PressForce_Ton(closestIndex)

    % find the max or min PressForce in a nearby area, 
    %write needed values to table  
    if str2num(SensorValues.PressPosition_Inches(closestIndex)) > 1 
        % disp('Press closed')
        SensorValues.PressPosition_Inches(closestIndex);
    else
        % disp('Press open') 
        SensorValues.PressPosition_Inches(closestIndex);
        for k = (1:10)
            pl = str2double(SensorValues.PressPosition_Inches(closestIndex - k));
            pr = str2double(SensorValues.PressPosition_Inches(closestIndex + k));
            
            if pl > 1
                closestIndex = closestIndex - k;
                % SensorValues.PressPosition_Inches(closestIndex) 
                break
            end
            if pr > 1
                closestIndex = closestIndex + k;
                % SensorValues.PressPosition_Inches(closestIndex)
                break
            end
        end
    end 

    j = closestIndex
    while SensorValues.PressPosition_Inches(j) > 1
        j = j - 1
    end
    press_open = j; 
    while SensorValues.PressPosition_Inches(1) < 1
        j = j - 1
    end
    upper_bound = j;


    j = closestIndex
    while SensorValues.PressPosition_Inches(j) > 1
        j = j + 1
    end
    press_open = j; 
    while SensorValues.PressPosition_Inches(1) < 1
        j = j + 1
    end
    upper_bound = j;

end


SensorValues.Properties.VariableNames;

% get the hand recorded timestamp
% identify a pressing cycle by press position (one open and close)
% -> need to keep in mind times when parameter's changed 

% RECORDID is added by this code 
