
folder = '';
% try 
    addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\')
    files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')
% catch
    % addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    % files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

csvdata = {'date', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization', 'RMS_L', 'RMS_R', 'clicks_L', 'clicks_R', 'THD_L', 'THD_R', 'wow', 'stereo_bleed'};
% size(headers)
% csvdata = {};
% csvdata = {headers};
size(csvdata)
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
    filename = files(i).name
    date = filename(1:6)
    pressing = filename(8:16)
    record = filename(19:20)
    side = filename(21)    
    top_stamper = filename(8)
    top_hits = filename(9:12)
    bottom_stamper = filename(13)
    bottom_hits = filename(14:17)



    [lagdiff, normalization, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed] = RecordProcess(file);
    csvdata
    csvdata(end+1,:) = {date, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed}

    
end