

%~~~~~~~~~ LOADING TABLES~~~~~~~~~~~%
    addpath('E:\audio_files\A0000B0000\')
    addpath('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\');

    % folder = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\rawpress_data\');
    % dataFolder = 'D:\OneDrive - University of Waterloo\School\Vinyl_Project\data\A0000B0000\'

    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/') 
    folder = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/rawpress_data/')
    dataFolder = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0000B0000/')

    clc
    disp('~~~~~~~~~~~~~~~~SENSOR PROCESS~~~~~~~~~~~~~~~~~~')
    date_tags = {'121918', '122018'};

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
            SensorValues = [SensorValues; readtable(strcat(folder,date_tags{i},'_SensorValues.csv'),opts)];
        end
        % disp(strcat('SensorValues',num2str(height(SensorValues))))

        %#####_JobDetailsCurrent
        % tracks cycle times

        if i == 1
            JobDetails = readtable(strcat(folder,date_tags{i},'_JobDetailsCurrent.csv'));
        else
            JobDetails = [JobDetails; readtable(strcat(folder, date_tags{i},'_JobDetailsCurrent.csv'))];
        end
        % disp(strcat('JobDetails',num2str(height(JobDetails))))

        %#####_ADAPT_DATA
        % tracks the press options changes via timestamps
        if i == 1
            ADAPT = readtable(strcat(folder, date_tags{i}, '_CZA0000B0000_ADAPT_DATA.xlsx'));
        else
            ADAPT = [ADAPT; readtable(strcat(folder,date_tags{i},'_CZA0000B0000_ADAPT_DATA.xlsx'))];
        end
        % disp(strcat('ADAPT',num2str(height(ADAPT))))


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



i = 1;
press_close = 1;
press_open = 1;
lower_bound = 1;
upper_bound = 1;
record_number = 0;
recordNumbers = [];
recordArray = [];
height(SensorValues)
while i < height(SensorValues)
        % disp(strcat('pre i:    ', num2str(i)))
        % disp(strcat('pre position:  ', num2str(SensorValues.PressPosition_Inches(i))))
        %~~~ if press is open keep looking
        if SensorValues.PressPosition_Inches(i) < 1
            i = i + 1; 
            continue

        else        

            press_close = i;
            press_open = i;
            while SensorValues.PressPosition_Inches(press_close) > 1
                press_close = press_close + 1;
            end 
            % disp('press_close set')

            lower_bound = press_open - 1;


            while SensorValues.PressPosition_Inches(lower_bound) < 1 
                lower_bound = lower_bound - 1;

                if lower_bound < 1
                    lower_bound = 1;
                    break
                end
            end
            lower_bound = lower_bound + 1;
            % disp('lower_bound set')


            while SensorValues.PressPosition_Inches(press_open) > 1
                press_open = press_open + 1;
            end 

            % disp('press_open set')


            upper_bound = press_open;
            while SensorValues.PressPosition_Inches(upper_bound) < 1 
                upper_bound = upper_bound + 1;
                if upper_bound > height(SensorValues)
                    recordArray = [recordArray; [record_number,lower_bound, press_close, press_open, height(SensorValues)-1]];
                    disp('BREAKING')
                    break
                end
            end
            % disp('upper_bound set')

        end 
        i = upper_bound;
        record_number = record_number + 1;
        recordArray = [recordArray; [record_number, lower_bound, press_close, press_open, upper_bound]];
        
        
    end    

% recordArray
% disp('going into 4 loop')
% recordTable = table('VariableNames', {'recordNumber'})
recordTable  = array2table(zeros(height(SensorValues),1), 'VariableNames', {'recordNumber'});
% optionsChange = array2table(zeros(0,9), 'VariableNames',{...});
% recordTable
disp('HEIGHT RECORD')
height(recordTable)
recordArray(length(recordArray),5) = recordArray(length(recordArray),5) - 1
length(recordArray)

for i = (1:length(recordArray))
    % recordArray(1,1)
    % recordArray(1,2)
    % recordArray(1,3)
    % recordArray(1,4)
    % recordArray(1,5)
    for j = (recordArray(i,2):recordArray(i,4))
        % recordTable(j,'recordNumber')
        % recordArray(i,1) 
        recordTable(j,'recordNumber') = num2cell(recordArray(i,1));
    end
end
% recordTable
% recordTable = table(record,'VariableNames', {'recordNumber'});
height(recordTable)
height(SensorValues)

recordTable = [recordTable SensorValues];
writetable(recordTable, strcat(dataFolder,'A0000B0000_PressTable.csv'));
% recordNumbers
