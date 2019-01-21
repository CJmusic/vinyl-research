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
    if abs(hi_freq)/abs(low_freq) > 2.0;
        %If there is sufficient high frequency content in the signal then keep track of
        %these indexes in a buffer
        signal_indices_Test = [signal_indices_Test, k*winSize]; %%%THIS ARRAY TO TEST DETECTION
        signal_buffer = [signal_buffer, k*winSize];
        continue;
    end
    %if length(signal_buffer) == 1; 
    %    %If a signal is detected in one window, but not the next it is most likely a click
    %    %add it to the clicks array and clear the current buffer
    %    disp('CLICK BUFFER');
    %    clicks = [clicks, signal_buffer];
    %    signal_buffer = [];
    %    continue;
    %end
    if length(signal_buffer) > 1 && signal_buffer(length(signal_buffer))/winSize - k < -50;
        %if k*winSize - signal_buffer(1) > 10;  
            %If a signal was detected and it stretches over more than a few windows, then add 
            %it to the signal array as start and ending indices and clear the buffer
            start_signal = signal_buffer(1); 
            end_signal   = signal_buffer(length(signal_buffer)); 
            signal = [signal;[start_signal, end_signal]];
            signal_buffer = []; %clear the signal buffer array
            continue;
     elseif length(signal_buffer) > 1;
        continue; 
        %end;
    end
    continue;
end


length(signal)

%this loops through the signal array and amalgamates adjacent entries 
%for xi = 1:length(signal);
%    if length(signal) - xi <=  0;
%        disp('YA GOOFED')
%        disp(xi)
%        continue;
%    end;
%
%    disp('index: ')
%    disp(xi)
%    disp('length of signal')
%    disp(length(signal))
%    disp('start: ');
%    disp(signal(xi,2)/winSize);
%    disp( 'end: ');
%    disp( signal(xi+1,1)/winSize);
%    
%    if signal(xi + 1, 1)/winSize -  signal(xi, 2)/winSize < 30; 
%%        start_buffer = signal(x1,1);
%%        end_buffer = signal(x1,2);
%        disp('Yall good')
%        signal(xi,2) = signal(xi+1,1);           
%        signal(xi+1,:) = [];
%    end;
%end;
        
signal/winSize
%length(clicks)
length(signal)
%clicks
%signal/winSize
% [s,f,t] = spectrogram(data,fs)
% imagesc (t, f, log(s));
figure(2);
clf();
plot(X)
hold on;
signal_indices_Test = transpose(signal_indices_Test);
for xi = 1:size(signal_indices_Test)
    x1 = (signal_indices_Test(xi));
    line([x1 x1], get(gca, 'ylim'),'Color', [0,0,0] + 0.8 ,'LineStyle', '--');
end


for i = 1:size(clicks);
    x1 = (clicks(i));
    line([x1 x1], get(gca, 'ylim'),'Color', 'blue','LineStyle', '-');
end  
  
for i = 1:size(signal);
    x1 = (signal(i, 1));
    x2 = (signal(i, 2));
    y1 =  1;
    y2 = -1;
    %area([x1,x2]);
    signal_area = fill([x1 x1 x2 x2],[y1 y2 y2 y1],'r');
    alpha(signal_area, '0.5');
    %line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '-'); 
    %line([x2 x2], get(gca, 'ylim'),'Color', 'green','LineStyle', '-');
end  




