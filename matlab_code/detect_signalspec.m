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

signal_buffer = [];


signal_indices_Test = []; %%This array is meant to test the detection algorithm 
                          %%independent of the indexing etc.

for k = 1:size(s,2);
    low_freq = sum(s(1:5,k));
    hi_freq = sum(s(6:end,k));
    if abs(hi_freq) > 0.5;
        %signal_index = k;
        signal_indices_Test = [signal_indices_Test, k*winSize]; %%%THIS ARRAY TO TEST DETECTION
        signal_buffer = [signal_buffer, k*winSize];
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
signal_indices_Test = transpose(signal_indices_Test);
for xi = 1:size(signal_indices_Test)
    x1 = (signal_indices_Test(xi));
    line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
end
