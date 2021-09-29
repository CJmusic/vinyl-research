% SensorProcess.m 
% Created by Chris Zaworski 
% Input is the ########_SensorValues.csv
% Output is the SensorTable.csv

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

opts = detectImportOptions('########_SensorValues.csv');
getvaropts(opts,selectedColumns);
opts = setvartype(opts,selectedColumns,'string');
opts.SelectedVariableNames = selectedColumns;
if i == 1
    disp('i EQUALED 1')
    SensorValues = readtable(strcat(folder,date_tags{i},'_SensorValues.csv'),opts)
else
    SensorValues = [SensorValues; readtable(strcat(folder,date_tags{i},'_SensorValues.csv'),opts)];
end

height(SensorValues)
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


i = 1;
press_close = 1;
press_open = 1;
lower_bound = 1;
upper_bound = 1;
record_number = 0;
recordNumbers = [];
recordArray = [];
while i < height(SensorValues)
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
        lower_bound = press_open - 1;
        while SensorValues.PressPosition_Inches(lower_bound) < 1 
            lower_bound = lower_bound - 1;
            if lower_bound < 1
                lower_bound = 1;
                break
            end
        end
        lower_bound = lower_bound + 1;
        while SensorValues.PressPosition_Inches(press_open) > 1
            press_open = press_open + 1;
        end 
        upper_bound = press_open;
        while SensorValues.PressPosition_Inches(upper_bound) < 1 
            upper_bound = upper_bound + 1;
            if upper_bound > height(SensorValues)
                recordArray = [recordArray; [record_number,lower_bound, press_close, press_open, height(SensorValues)-1]];
                break
            end
        end
    end 
    i = upper_bound;
    record_number = record_number + 1;
    recordArray = [recordArray; [record_number, lower_bound, press_close, press_open, upper_bound]];
end    

recordTable  = array2table(zeros(height(SensorValues),1), 'VariableNames', {'recordNumber'});
recordArray(length(recordArray),5) = recordArray(length(recordArray),5) - 1

for i = (1:length(recordArray))
    for j = (recordArray(i,2):recordArray(i,4))
         recordTable(j,'recordNumber') = num2cell(recordArray(i,1));
    end
end

recordTable = [recordTable SensorValues];
writetable(recordTable, strcat(dataFolder,'A0000B0000_SensorTable.csv'));
