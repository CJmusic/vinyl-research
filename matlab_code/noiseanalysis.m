

%filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00-calibration_signals/*.wav';
path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00-calibration_files';

files = dir(strcat(path,'*.wav'));
files
%filepath = dir(filename);

%FOR LOOP to loop through the noise recordings
for file = files
    strcat(path,'\*.wav',file.name)
    file.name
    [data, fs] = audioread(file.name)%strcat(path,'\*.wav',file.name));
    
    
    n_sam = size(data);
    freq = fs*(0:(n_sam/2))/n_sam;
    
    fft_data = fft(data)/n_sam;
    fft_data = fft_data(1:size(fft_data)/2+1);
    
    clf(figure(1))
    figure(1);
    grid on; hold on; legend;
    set(gca, 'XScale', 'log'); 
    
    plot(freq,20*log10(real(fft_data)), 'DisplayName', filename);
    
end

