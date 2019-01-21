

%filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00-calibration_signals/*.wav';
%path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00-calibration_files';

%path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00-calibration_files';

path = 'FocusriteGen1'

%files = dir(strcat(path,'/*.wav'));
%files.name
%filepath = dir(filename);

clf(figure(1))
figure(1);
%FOR LOOP to loop through the noise recordings
%for file = files'
%    file
%    [pathstr,name,ext] = fileparts(strcat(path,'/',file.name));
%    %strcat(path,'/*.wav',file.name)
%    %file.name
%    [data, fs] = audioread(strcat(path,'/',file.name));
%    
%    
%    n_sam = length(data);
%    freq = fs*(0:(n_sam/2))/n_sam;
%    
%    fft_data = fft(data)/n_sam;
%    fft_data = fft_data(1:size(fft_data)/2+1);
%    
%    grid on; hold on; legend;
%    set(gca, 'XScale', 'log'); 
%    
%    plot(freq,20*log10(real(fft_data)), 'DisplayName', name);
%    
%end


%file
%strcat(path,'/*.wav',file.name)
%file.name

%%%THIS CODE BELOW WORKS, THIS IS NEED TO PUT IN THE LOOP


path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00_calibration_signals/FocusriteStereo1KhzCh26.2V.wav'

[data, fs] = audioread(path);

size(data(1))
data_L = transpose(data)(2); 
data_R = transpose(data)(1);

n_sam = length(data);
freq = fs*(0:(n_sam/2));
freq = freq/n_sam;

fft_data_L = fft(data_L)/n_sam;
fft_data_L = fft_data_L(1:size(fft_data_L)/2+1);

fft_data_R = fft(data_R)/n_sam;
fft_data_R = fft_data_R(1:size(fft_data_R)/2+1);

clf(figure(1))
figure(1);
grid on; hold on; legend;
set(gca, 'XScale', 'log'); 

plot(freq,20*log10(real(fft_data_L)));%, 'DisplayName', filename);
plot(freq,20*log10(real(fft_data_R)));%, 'DisplayName', filename);
    
    
