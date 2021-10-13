
AudioTable = readtable('/Volumes/AUDIOBANK/audio_files/A0137B0137/A0137B0137-AudioTable.csv');
SensorTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_SensorTable.csv');
RecordTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_RecordNumbers.csv');

SettingsTable = readtable('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/data/A0137B0137/A0137B0137_SettingsTable.csv')

pressruns = unique(AudioTable.pressing);
tracks = unique(AudioTable.track);

% clean up AudioTable
AudioTable.RecordID = erase(AudioTable.record,'.wav');

for i = 1:height(AudioTable)
    AudioTable.RecordID{i} = str2num(AudioTable.RecordID{i}(1:3));
end

AudioTable.RecordID = cell2table(AudioTable.RecordID);
AudioTable.RecordID = AudioTable.RecordID.Var1;



%~~~~~~~~~~~~~~~ MERGING TABLES STARTS ~~~~~~~~~~~~~~%

    Tbl = outerjoin(RecordTable, SensorTable);
    Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
    writetable(Tbl,'Tbl1.csv')

    head(Tbl)
    head(AudioTable)
    
    Tbl(isnan(Tbl.RecordID), :) = [];
    % join in AudioTable
    Tbl = outerjoin(Tbl, AudioTable, 'Keys', {'RecordID', 'RecordID'});
    Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
    writetable(Tbl,'Tbl2.csv')

    Tbl.Properties.VariableNames([3]) = {'pressing'};

    % join in SettingsTable
    Tbl = outerjoin(Tbl, SettingsTable);
    writetable(Tbl,'Tbl4.csv')

    Tbl.Properties.VariableNames([1]) = {'PressingNumber'};
    Tbl.Properties.VariableNames([33]) = {'track'};
    Tbl.Properties.VariableNames([3]) = {'pressing'};

    Tbl.Properties.VariableNames([2]) = {'RecordID'};

    Tbl.PressingNumber_SensorTable = [];
    Tbl.pressing_AudioTable = [];
    Tbl.pressing_SettingsTable = [];
    Tbl.RecordID_AudioTable = [];
    writetable(Tbl,'PRESSING.csv')

