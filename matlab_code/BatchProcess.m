
folder = '';
% try 
%     addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\')
%     files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')
% catch
    addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/')
    files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
% end

headers = {'date', 'pressing', 'top_stamper', 'top_hits',   'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'normalization', 'RMS_L, RMS_R', 'clicksL','clicksR', 'thd', 'wow', 'bleed'};


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
    [lagdiff, normalization, RMS_L, RMS_R, clicks_L, clicks_R, THD_L, THD_R, wow, stereo_bleed] = RecordProcess(file);
    % [lagdiff,normalization,RMS_L,RMS_R,clicksL,clicksR,thd,wow,bleed] = RecordProcess(file)
    

    dataTable(i) = [date, pressing, top_stamper, top_hits, bottom_stamper, bottom_hits, record, side, track, lagdiff, normalization, RMS_L, RMS_R, clicks_L, clicks_R, thd, wow, bleed ]
    
end