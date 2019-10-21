close all


% dataFile = ('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/A0000B0000-data.csv')
dataFile = ('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\A0000B0000\A0000B0000-data.csv')


byN = readtable(dataFile);
for i = (1:length(byN.Properties.VariableNames))
    disp(byN.Properties.VariableNames(i))
end
byN

figure(1); hold on; grid on;
% scatter(getData(byN,'transition', 'PressForce_Ton'), getData(byN,'transition', 'clicks_L'))
scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_L'))
scatter(getData(byN,'transition', 'record'), getData(byN,'transition', 'RMS_R'))
title('RMS levels in transition track')
xlabel('record #')
ylabel('RMS level [dB]')
legend(['Left', 'Right'])

figure(2); hold on; grid on;
hist(getData(byN,'transition','RMS_L'))
hist(getData(byN,'transition','RMS_R'))


%% find the 10 records with the min and max noise levels

byRMS = sortrows(byN);
disp('Highest RMS_L')
mins = getData(byRMS,'transition','record');
maxs = getData(byRMS,'transition','record');

mins(1:10)
maxs(end-10:end)

function data_return = getData(Tbl, track, param) 
    data_return = (table2array(Tbl(strcmp(Tbl.track, track), param)));
    if strcmp(class(data_return),'cell')
        data_return = str2double(data_return);
    end
end 