close all; clear all; clc;

folder = '';
% try 
    % addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\')
    % files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

% addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% folder = ('/Volumes/AUDIOBANK/audio_files/A0000B0000/')
% files = dir(strcat(folder,'*.wav'))

addpath('E:\audio_files\A0000B0000\')
folder = ('E:\audio_files\A0000B0000\')
files = dir(strcat(folder,'*.wav'))

csvdata = {'date', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'THD_L', 'THD_R', 'wow_L', 'wow_R', 'stereo_bleed'};
size(csvdata)



pressingID = files.folder
pressingID = pressingID(end-9:end)

% files.folder(1)
try

    T = readtable(strcat(folder,pressingID,'.csv'))
    T = readtable('myDataFile.csv')
catch
    disp('CANT FIND CSV FILTER')
    T  = cell2table(cell(0,length(csvdata)), 'VariableNames', csvdata);

    % T = cell2table(empt,'VariableNames',csvdata(1,:))
    % writetable(T,strcat(folder,pressingID,'.csv'));
end


for i = (1:length(files))
    filename = files(i).name
    T.Properties.VariableNames


    T.record
    str2num(filename(19:21))
    T

    if isempty(T.record)
    elseif ismember(str2num(filename(19:21)), T.record)
        disp('ALREADY PROCESSED')
        continue
    end

    file = strcat(files(i).folder,'/',files(i).name)

    date = 0;
    pressing = 0;
    top_stamper = 0;
    top_hits = 0;
    bottom_stamper = 0;
    bottom_hits = 0;
    record = 0;
    side  = 0;
    track  = 0;
    
    %STRIP RELEVANT INFO FROM NAME 
    date = str2num(filename(1:6))
    record = str2num(filename(19:21))

    top_stamper = filename(8)
    pressing = filename(8:16)
    top_hits = str2num(filename(9:12)) + record
    bottom_stamper = filename(13)
    bottom_hits = str2num(filename(14:17)) + record
    side = filename(22)




    output = RecordProcess(file)
    numrec = size(output);
    for i=(1:numrec)
        track = output(i,1)
        lagdiff = output(i,2)
        normalization_L = output(i,3)
        normalization_R = output(i,4)
        RMS_L = output(i,5)
        RMS_R = output(i,6)
        clicks_L = output(i,7)
        clicks_R = output(i,8)
        THD_L = output(i,9)
        THD_R = output(i,10)
        wow_L = output(i,11)
        wow_R = output(i,12)
        stereo_bleed = output(i,13)
      

        date
        % date{1}
        track{1}
        lagdiff{1,1}
        normalization_L{1,1}
        normalization_R{1,1}
        RMS_L{1,1}
        RMS_R{1,1}
        clicks_L{1,1}
        clicks_R{1,1}
        THD_L{1,1}
        THD_R{1,1}
        wow_L{1,1}
        wow_R{1,1}
        stereo_bleed{1,1}
        

        csvdata(end+1,:) = {date, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed};
        % N = {date, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow_L, wow_R, stereo_bleed};
        % dlmwrite(strcat(folder,pressingID,'.csv'),N,'delimiter',',','-append');
        numbers = randi(9, 10, 1);
        num_str = num2str(numbers);
        num_cell = mat2cell(num_str, ones(10, 1), 1);

        T1 = cell2table(csvdata(end,:),'VariableNames',csvdata(1,:))
        T
        T = [T;T1]

        writetable(T,'myDataFile.csv');
    end

end

% dlmwrite('test.csv',M,'delimiter',',');
% N = randn(4,4);
% dlmwrite('test.csv',N,'delimiter',',','-append');

% T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:))
% writetable(T,'myDataFile.csv');