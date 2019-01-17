filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.3.wav';

[data, fs] = audioread(filename);

time = linspace(0,(length(data)-1)/fs,length(data));

n_block = 2048; 
n_sam = length(data);
data_fft = fft(data)/n_block;
freq = fs*(0:(n_block)/2)/n_block;
signal_indices = [];

%for xi = 1:size(data)
    %if abs(data(xi) - data(xi + 1)) > 0.2%data(xi) > threshold
        %signal_indices = [signal_indices,xi];
    %end
%end

RMS_Value = rms(data(:,1));
threshold = RMS_Value/2

for xi = 1:size(data)
    if abs(data(xi) - data(xi + 1)) > 0.2%data(xi) > threshold
        signal_indices = [signal_indices,xi];
    end
end
signal_indices=transpose(signal_indices);
size(signal_indices)
fig=figure(1); 
hax=axes; 
hold on;
plot(time,data)
for xi = 1:size(signal_indices)
    x1 = time(signal_indices(xi));
    line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
end
