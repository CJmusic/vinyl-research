
folder = '';
% try 
    addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\')
    files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

csvdata = {'date', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization_L', 'normalization_R','RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'THD_L', 'THD_R', 'wow', 'stereo_bleed'};

for i = (1:length(files))
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
    filename = files(i).name;
    date = filename(1:6);
    top_stamper = filename(8);
    pressing = filename(8:16);
    top_hits = filename(9:12);
    bottom_stamper = filename(13);
    bottom_hits = filename(14:17);
    record = filename(19:20);
    side = filename(21);

    top_hits = str2num(top_hits) + str2num(record);
    bottom_hits = str2num(bottom_hits) + str2num(record);



    output = RecordProcess(file)
    disp('size output')
    numrec = size(output)
    for i=(1:numrec)
        track = output(i,1);
        lagdiff = output(i,2);
        normalization_L = output(i,3);
        normalization_R = output(i,4);
        RMS_L = output(i,5);
        RMS_R = output(i,6);
        clicks_L = output(i,7);
        clicks_R = output(i,8);
        THD_L = output(i,9);
        THD_R = output(i,10);
        wow = output(i,11)
        stereo_bleed = output(i,12);

        csvdata(end+1,:) = {date, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization_L, normalization_R, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed};
    end

end


T = cell2table(csvdata(2:end,:),'VariableNames',csvdata(1,:));
writetable(T,'myDataFile.csv');