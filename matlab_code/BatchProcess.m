



folder = '';

files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')

headers = {'date', 'pressing', 'top_stamper', 'top_hits', 'bottom_stamper', 'bottom_hits', 'record', 'side', 'track', 'lagdiff', 'norm', 'rms', 'clicksL','clicksR', 'thd', 'wow', 'bleed'}


for i = (1:length(files))
    file = strcat(files(i).folder,'/',files(i).name)

    
    %STRIP RELEVANT INFO FROM NAME 
    

    output = RecordProcess(file);
    
end