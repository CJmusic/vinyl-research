
% addpath('E:\audio_files\A0000B0000\')
% folder = ('D:\OneDrive - University of Waterloo\Vinyl_Project\data\121918_A0000B0000\');

addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
folder = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/data/121918_A0000B0000/')



date_tags = {'121918', '122018'};

AggTable = [];

for i=(1:length(date_tags))
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
    'ExtruderMeltTemp_F'};

    opts = detectImportOptions(strcat(folder,date_tags{i},'_SensorValues.csv'));
    getvaropts(opts,selectedColumns);
    opts = setvartype(opts,selectedColumns,'string');
    opts.SelectedVariableNames = selectedColumns;
    
    if i == 1
        SensorValues = readtable(strcat(folder,date_tags{i},'_SensorValues.csv'),opts);
    else
        SensorValues = [SensorValues; readtable(strcat(folder,date_tags{1},'_SensorValues.csv'),opts)];
    end

%#####_JobDetailsCurrent
% tracks cycle times

    if i == 1
        JobDetails = readtable(strcat(folder,date_tags{i},'_JobDetailsCurrent.csv'));
    else
        JobDetails = [JobDetails; readtable(strcat(folder, date_tags{i},'_JobDetailsCurrent.csv'))];
    end

    %#####_ADAPT_DATA
    % tracks the press options changes via timestamps
    if i == 1
        ADAPT = readtable(strcat(folder, date_tags{i}, '_CZA0000B0000_ADAPT_DATA.xlsx'));
    else
        ADAPT = [ADAPT; readtable(strcat(folder,date_tags{i},'_CZA0000B0000_ADAPT_DATA.xlsx'))];
    end


end

SensorValues.PressPosition_Inches = str2double(SensorValues.PressPosition_Inches);
SensorValues.PressForce_Ton = str2double(SensorValues.PressForce_Ton);
SensorValues.MouldSteamIn_F = str2double(SensorValues.MouldSteamIn_F);
SensorValues.MouldSteamOutTop_F = str2double(SensorValues.MouldSteamOutTop_F);
SensorValues.MouldSteamOutBottom_F = str2double(SensorValues.MouldSteamOutBottom_F);
SensorValues.ExtruderFeedthroatTemp_F = str2double(SensorValues.ExtruderFeedthroatTemp_F);
SensorValues.ExtruderBarrelZone1Temp_F = str2double(SensorValues.ExtruderBarrelZone1Temp_F);
SensorValues.ExtruderBarrelZone2Temp_F = str2double(SensorValues.ExtruderBarrelZone2Temp_F);
SensorValues.ExtruderBarrelZone3Temp_F = str2double(SensorValues.ExtruderBarrelZone3Temp_F);
SensorValues.ExtruderDieZoneTemp_F = str2double(SensorValues.ExtruderDieZoneTemp_F);
% SensorValues.ExtruderFeedthroatTemp_F = str2double(SensorValues.ExtruderFeedthroatTemp_F);
SensorValues.ExtruderPremouldTemp_F = str2double(SensorValues.ExtruderPremouldTemp_F);
SensorValues.ExtruderMeltTemp_F = str2double(SensorValues.ExtruderMeltTemp_F);


%#####_TimeStamps
% Recorded by hand timestamps

TimeStamps = readtable(strcat(folder,'121918_TimeStamps.xlsx'));
SensorValues.RecordTimeStamp = datetime(SensorValues.RecordTimeStamp, 'InputFormat', 'yyyy-MM-dd HH:mm');
TimeStamps.TimeStamp = datetime(TimeStamps.TimeStamp, 'InputFormat', 'MMMM d, yyyy hh:mm:ss');

for i = (1:length(TimeStamps.TimeStamp))
    % if TimeStamps.TimeStamp(i) > datetime(date_tag)
    %     date_tag
    %     continue
    % end
    % disp('FOUND MATCH')

    %% convert all the necessary columns to numbers
    [closestTimeStamp,closestIndex] = min(abs(SensorValues.RecordTimeStamp-TimeStamps.TimeStamp(i)));
    if isnan(closestTimeStamp)
        continue
    end
    % disp(strcat(string(SensorValues.RecordTimeStamp(closestIndex)),'.....',string(TimeStamps.TimeStamp(i))))

    
    % SensorValues.RecordTimeStamp(closestIndex)
    % SensorValues.PressPosition_Inches(closestIndex)
    % SensorValues.PressForce_Ton(closestIndex)

    % find the max or min PressForce in a nearby area, 
    %write needed values to table  
    % SensorValues.PressPosition_Inches(closestIndex)
    if str2double(SensorValues.PressPosition_Inches(closestIndex)) > 1 
        % disp('Press closed')
        SensorValues.PressPosition_Inches(closestIndex);
    else
        % disp('Press open') 
        SensorValues.PressPosition_Inches(closestIndex);

        % if press is open then find the first value where the press is closed
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

    %% closestIndex points to a position where the press was closed 

    % look backwards in time to find when the press was last open
    j = closestIndex;
    while str2double(SensorValues.PressPosition_Inches(j)) > 1 && j > (closestIndex - 10) % && j < closestIndex 
        j = j - 1;
    end
    press_closed = j + 1; %this represents the table index when the press was closed

    %look backwards to find when the press 
    while str2double(SensorValues.PressPosition_Inches(j)) < 1 && j < closestIndex - 10
        j = j - 1;
    end
    lower_bound = j + 1; %

    % look forwards in time to find the next time the press closes
    j = closestIndex;
    
    while str2double(SensorValues.PressPosition_Inches(j)) > 1 && j < closestIndex + 10
        j = j + 1;
    end
    
    press_open = j - 1; 

    try  
        while str2double(SensorValues.PressPosition_Inches(j)) < 1 && j < closestIndex + 10
            j = j + 1;
        end
    catch 
        disp('END OF FILE')
    end
    
    upper_bound = j - 1;


    % pull the necessary adapt settings

    % perform all the needed measurements 
    % SensorValues(lower_bound:upper_bound)

    %% record number and record ID will be different because of Junk records 
    % determine from the time stamps when the nearest record 

    % TimeStamps.TimeStamp(i)
    % TimeStamps.RECORDID{i}
    lower_bound
    upper_bound
    SensorValues.PressPosition_Inches(lower_bound:upper_bound)


end
size(SensorValues)
% get the hand recorded timestamp
% identify a pressing cycle by press position (one open and close)
% -> need to keep in mind times when parameter's changed 

% RECORDID is added by this code 
