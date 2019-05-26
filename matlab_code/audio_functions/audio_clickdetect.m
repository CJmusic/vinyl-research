%{
Lining up clicks: 

clicks = [ sample numbers ] 

% a click has a sustained high derivative 
% the energy of a click is between 2-4 kHz

%} 

% TESTING THE CLICK DETECT FUNCTION
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/260219_noisereferenceinst/');
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/'

AUDIO_FILES = {'one.wav','two.wav','three.wav', 'four.wav', 'five.wav'};

clear all; clc; close all;
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/';

audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz.wav';
audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz.wav';

AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

%[data, time, fs] = audio_load(strcat(audio_dir,AUDIO_FILES{1}));
[data, time, fs] = audio_load(audio_bin);
data = data(5.1*fs : 10.1*fs,:);
audio_clickdetecttest(data, fs);

function clicks = audio_clickdetecttest(data, fs);
    data = data(:,1);
%~~~~~~~~~~~~~~~~PRE-FILTERS~~~~~~~~~~~~~~~~~~~
    fc = 50;
    [b,a] = butter(6,300/(fs/2),'low');
    data = filter(b, a, data); % use filtfilt

    [b,a] = butter(6,3000/(fs/2),'high');
    data = filter(b, a, data); % use filtfilt

    % [b,a] = butter(6,fc/(fs/2),'low');
    % data = filter(b, a, data); % use filtfilt
    %fc = 1000.0 
    %  [A,B,C,D] = butter(10,[1000 10000]/fs/2);
    % [b, a] = butter(1, [1000.0 4000.0]/(fs/2),'bandpass'); % this should 
    % % data_f = filter(b, a, data); % use filtfilt
    % data = filter(b, a, data); % use filtfilt
%~~~~~~~~~~~~~~~~PRE-FILTERS END~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~PLOT DATA~~~~~~~~~~~~~~~~~~~~~~~~~
    time = (1:length(data))/fs;
    % figure(1);
    % freqs(b,a)
    % title('bandpass filter')
    % figure(2); grid on; hold 
    % plot(time, data)
    % title('audio data')
    % find the start of the click 
    figure(1); grid on; hold on 
    plot(time, data)
    title('audio data')

    p_data = data.^2;
    figure(2); grid on; hold on;
    plot(time, p_data)
    title('audio power')
%~~~~~~~~~~~~~~~~PLOT DATA END~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~COMBI PEAKS/ENV METHOD~~~~~~~~~~~~~~

    % p_data = data.^2;
    % clicks = [];
    % zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    % zero_indices = zci(data);
    % size(zero_indices)
    % size(p_data)
    % avg_peak = 0.0;
    % prev_peak = 1.0;
    % peak_value = 0.0;

    % peak_values = zeros(length(zero_indices)-1);   
    % peak_indices = zeros(length(zero_indices)-1);
    
    % [peak_values, peak_indices] = findpeaks(p_data(:,1));
    % adata = peak_values;
    % % adata = abs(data);
    % lenClick = 6*96; % 6 ms
    % lenClick = 3;
    % prevPeak = 0.0;
    % peakLen = 0;
    % clicks = [];

    % figure(3); grid on;
    % plot(peak_values);



    % for i=(1:length(adata))
    %     comparator = true;
    %     % buff = data(i+lenClick);
    %     % if comparator == false
    %     %     prevPeak = 0.0;
    %     %     peakLen = 0;
    %     % elseif comparator == true
    %         % prevPeak = data(i);
    %         % peakLen = peakLen + 1;
    %     % end % if comparator
    %     if peakLen > lenClick
    %         %%% CLICK DETECTED %%%
    %         clicks = [clicks, i-lenClick];
    %         comparator = false;
    %         prevPeak = 0.0;
    %         peakLen = 0;
    %     end 

    %     threshold = 100;

    %     if adata(i) > threshold*prevPeak
    %         prevPeak = adata(i);
    %         peakLen = 0; 
    %         % comparator = true;
    %     elseif adata(i) < prevPeak
    %         peakLen = peakLen + 1;
    %     end
    % end % for adata




%~~~~~~~~~~~~COMBI PEAKS/ENV METHOD END~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~DECAY ENV METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % adata = abs(data);
    % lenClick = 6*96; % 6 ms
    % prevPeak = 0.0;
    % peakLen = 0;
    % clicks = [];
    % for i=(1:length(adata))
    %     comparator = true;
    %     % buff = data(i+lenClick);
    %     % if comparator == false
    %     %     prevPeak = 0.0;
    %     %     peakLen = 0;
    %     % elseif comparator == true
    %         % prevPeak = data(i);
    %         % peakLen = peakLen + 1;
    %     % end % if comparator
    %     if peakLen > lenClick
    %         %%% CLICK DETECTED %%%
    %         clicks = [clicks, i-lenClick];
    %         comparator = false;
    %         prevPeak = 0.0;
    %         peakLen = 0;
    %     end 

    %     threshold = 0.9;

    %     if adata(i) > threshold*prevPeak
    %         prevPeak = adata(i);
    %         peakLen = 0; 
    %         % comparator = true;
    %     elseif adata(i) < prevPeak
    %         peakLen = peakLen + 1;
    %     end
    % end % for adata

%~~~~~~~~~~~~~~~~~~~~~~~DECAY ENV METHOD END~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~AVG PEAK METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % p_data = data.^2;
    p_data = abs(data);
    clicks = [];
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    zero_indices = zci(data);
    size(zero_indices)
    size(p_data)
    avg_peak = 0.0;
    prev_peak = 1.0;
    peak_value = 0.0;

    peak_values = zeros(length(zero_indices)-1);   
    peak_indices = zeros(length(zero_indices)-1);
    
    [peak_values, peak_indices] = findpeaks(p_data(:,1));

    mAvgWidth = 2;
    for i = (mAvgWidth+1:length(peak_values)-mAvgWidth-1)
        if i == mAvgWidth + 1
            peak_values(i-mAvgWidth : i+mAvgWidth);
            iPeaks = peak_values(i-mAvgWidth : i+mAvgWidth);
            mAvg = (1/(mAvgWidth*2 + 1))*sum(iPeaks);
        end
        plow = peak_values(i+mAvgWidth); 
        phigh =  peak_values(i-mAvgWidth);
        mAvg = mAvg - (phigh - plow)/(21);

        threshold = 5;
        if peak_values(i) > threshold*mAvg
            if length(peak_values(peak_values(i : i+10) > 0.5*peak_values(i))) < 2;
                continue    
            end
            clicks = [clicks, peak_indices(i)];
        end
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~AVG PEAK METHOD END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    %     [peak_value, peak_index] = max(p_data(zero_indices(i): zero_indices(i+1)));
    %     peak_values = [peak_values, peak_value];
    %     peak_indices = [peak_indices, peak_index];
    %     %if peak_value > 10.0*prev_peak;
    %    if peak_value/prev_peak > 10.0;
    %       click = zero_indices(i) + peak_index;
    %       clicks = [clicks, click];
    %    else; 
    %       avg_peak = (avg_peak + prev_peak - )/2;
    %    end
    %    prev_peak = peak_value;
    % end
    % size(peak_values)
    
%~~~~~~~~~~~~~~~~PLOTTING~~~~~~~~~~~~~~~~~~~
    % figure(2); hold on;
    % for xi = 1:length(zero_indices)-1
    %     % xi
    %     x1 = time(zero_indices(xi));
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    % figure(1); hold on;
    % for xi = 1:length(peak_indices)-1
    %     x1 = peak_indices(xi)/fs;
    %     % x1 = peak_indices(xi);
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    for xi = 1:length(clicks)-1
        x1 = time(clicks(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');

        figure(2); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end
%~~~~~~~~~~~~~~~~PLOTTING END~~~~~~~~~~~~~~~
    % determine the click's polarity 
    
    % length and number of Oscillations 
    
end   


% CODE BELOW is old click method of first differences  


% function clicks = audio_clickdetect(data, fs);
%     time = (0:length(data)-1)/fs; 
%     d_data = diff(data,1)*fs;% since delta_t = 1/fs:w;
%     clicks = [];
%     threshold = 10*rms(d_data(1:1025));
%     disp('starting click detect')
%     for i = (1:length(data)-1);
%         if i > 512 && i + 512 < length(d_data);
%             % threshold = threshold - abs(d_data(i-511)) + abs(d_data(i+511));
%             threshold = 10*rms(d_data(i-511:i+512));
%         end
%         if d_data(i) > threshold;
%             click = i;
%             %click_timestamp = i/fs;   
%             clicks = [clicks, click];   
%         end
%     end
%     disp('Number of clicks');
%     disp('new file')
%     size(clicks) 
    
%     %plotting below
    
%     %clf(figure(1));clf(figure(2)); 
%     %figure(1); grid on; hold on;
%     %plot(time, data); 
%     % for xi = 1:length(clicks);
%     %      x1 = time(clicks(xi));
%     %      line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%     % end
%     % title('Click detection: Amplitude vs. Time'); 
%     % xlabel('Time [s]');
%     % ylabel('Amplitude');
    
%     %figure(2); grid on; hold on;
%     %t_diff = length(time) - length(d_data)
%     %plot(time(1:end-t_diff),d_data); grid on; hold on;
%     %for xi = 1:length(clicks);
%     %     x1 = time(clicks(xi));
%     %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
%     % end
    
%     % title('Click detection: 1st Derivative vs. Time'); 
%     % xlabel('Time [s]');
%     % ylabel('dA/dt (s^-1)');
% end 




% %}