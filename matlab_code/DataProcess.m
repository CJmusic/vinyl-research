close all


dataFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/A0000B0000-data.csv')

Tbl = readtable(dataFile);
for i = (1:length(Tbl.Properties.VariableNames))
    disp(Tbl.Properties.VariableNames(i))
end
Tbl

figure(1); hold on; grid on;
scatter(getData(Tbl,'transition', 'PressForce_Ton'), getData(Tbl,'transition', 'clicks_L'))



function data_return = getData(Tbl, track, param) 
    data_return = (table2array(Tbl(strcmp(Tbl.track, track), param)));
    if strcmp(class(data_return),'cell')
        data_return = str2double(data_return);
    end
end 