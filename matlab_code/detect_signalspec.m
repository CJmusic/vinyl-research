filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.3.wav';

[data, fs] = audioread(filename);

% time = linspace(0,(length(data)-1)/fs,length(data));

% n_block = 2048; 
% n_sam = length(data);
% data_fft = fft(data)/n_block;
% freq = fs*(0:(n_block)/2)/n_block;
X = data(:,1);
size(X)
% dur=0.5;
winSize= 2^8;%round(fs*dur);
overlap=0;%round(winSize/2);
fftsize=winSize;
% figure
[s,f,t] = spectrogram(X,winSize,overlap,fftsize,fs,'yaxis');

% spectrogram(X,winSize,overlap,fftsize,fs,'yaxis')
size(s)
x = [];
signal_indices = [];

signal = [];
clicks = [];
index = 0;

for k = 1:size(s,2);
    low_freq = sum(s(1:5,k));
    hi_freq = sum(s(6:end,k));
    if abs(hi_freq)> 0.5;
        signal_index = k;
        signal_indices = [signal_indices, k*winSize];
        if size(signal_indices) == 0 || signal_indices(k-1) <= (k-1)*winSize;
            start_signal = k;
        end
    end

    length(signal_indices)
    if (k-1)*winSize == signal_indices(int_16(length(signal_indices))); 
        end_signal = k;
        if end_signal - start_signal == 1;
            click = [click, start_signal];
            continue
        end
        signal = [signal, [start_signal,end_signal]];
    end
end
% [s,f,t] = spectrogram(data,fs)
% imagesc (t, f, log(s));
figure(1);
clf();
plot(x,'red')
figure(2);
clf();
plot(X)
hold on;
size(data)
size(signal_indices)
signal_indices = transpose(signal_indices);
for xi = 1:size(signal_indices)
    x1 = (signal_indices(xi));
    line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
end

% spectrogram(data,1024);
% spectrogram(data,'yaxis')
