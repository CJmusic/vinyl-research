%{
Lining up clicks: 

clicks = [ sample numbers ] 

% a click has a sustained high derivative 
% the energy of a click is between 2-4 kHz

%} 

% TESTING THE CLICK DETECT FUNCTION
clc; close all;
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/');
audio_dir = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/from_John/';

audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz.wav';
audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz.wav';

AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

[data, time, fs] = audio_load(audio_bin);
data = data(5.1*fs : 10.1*fs,:);
audio_clickdetecttest(data, fs);

function clicks = audio_clickdetecttest(data, fs)
    % data = data(:,1);
    time = (1:length(data))/fs;
%~~~~~~~~~~~~~~~~PRE-FILTERS~~~~~~~~~~~~~~~~~~~
    fc = 50;
    [b,a] = butter(6,2000/(fs/2),'low');
    data = filter(b, a, data); % use filtfilt

    [b,a] = butter(6,10000/(fs/2),'high');
    data = filter(b, a, data); % use filtfilt

    % [b,a] = butter(6,fc/(fs/2),'low');
    % data = filter(b, a, data); % use filtfilt
    %fc = 1000.0 
    %  [A,B,C,D] = butter(10,[1000 10000]/fs/2);
    % [b, a] = butter(1, [1000.0 4000.0]/(fs/2),'bandpass'); % this should 
    % % data_f = filter(b, a, data); % use filtfilt
    % data = filter(b, a, data); % use filtfilt
%~~~~~~~~~~~~~~~~PRE-FILTERS END~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~COMBI PEAKS/ENV METHOD~~~~~~~~~~~~~~
    % p_data = data.^2;
    p_data = abs(data);
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    zerosL = zci(data(:,1));
    zerosR = zci(data(:,2));

    peak_values = zeros(length(zerosL)-1,2);
    peak_indices = zeros(length(zerosL)-1,2);
    % length(zerosL)
    % length(p_data)
    % disp('last zero')
    % length(p_data) - zerosL(length(zerosL))
    % length(p_data) - length(zerosL)
    % zerosL(length(zerosL)-10:length(zerosL))
    size(p_data)
    for i = (1:length(zerosL)-1)
        buffL = p_data(zerosL(i):zerosL(i+1),1);
        buffR = p_data(zerosL(i):zerosL(i+1),2);
        [peakL, indexL] = max(buffL);
        [peakR, indexR] = max(buffR);
        peak_values(i,1)  = peakL;
        peak_indices(i,1) = indexL;
        peak_indices(i,2) = indexR;

    end
    for i = (1:length(zerosR)-1)

        peak_indices(i,2) = indexR;
        peak_values(i,2)  = peakR;
        [peakR, indexR] = max(buffR);
    end


    % [peak_values, peak_indices] = findpeaks(p_data(:,1));
    disp('Sizes')
    size(peak_indices)
    size(peak_values)
    size(zerosL)

    clicks = [];
    clicks_peaks = [];

    threshold = 1.2;
    lenClick = 1412;
    mAvgWidth = 2;
    avgBuffersize = 20; %% how many peaks to take into account for the moving avg
    
    N = 0;
    mAvg = 0;
    for i = (1:length(peak_values))
        N = i;
        % takes the moving average from the left
        mAvg = mAvg + peak_values(i,:)/N;

        if peak_values(i) > threshold*mAvg
            % disp('CLICK DETECTED')
            mAvg = mAvg - peak_values(i)/N; %% DON'T include clicks in the moving avg
            % if max(peak_values(i,1)) > max(peak_values(i,2))
            %     clpolarity = 0; %% click is in the left channel
            % else
            %     clpolarity = 1; %% click is in the right channel 
            % end

            if length(clicks) == 0;
                clicks = [clicks, peak_indices(i)];  
            elseif peak_indices(i) - clicks(length(clicks)) > lenClick %% makes sure the same click isn't recorded multiple times
                clicks = [clicks, peak_indices(i)];  
                % clicks = [clicks, zerosL(i)];  
            end
            % clicks = [clicks, [peak_indices(i-3:i+8)]]; %% What's recorded here is the indices of the peaks, perhaps the zeros make more sense
        end

        if N > avgBuffersize
            mAvg = mAvg - peak_values(i-avgBuffersize)/N;
        end
    end


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
%~~~~~~~~~~~~~~~~PLOTTING~~~~~~~~~~~~~~~~~~~
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

    % p_data = data.^2;
    % figure(2); grid on; hold on;
    % plot(time, p_data)
    % title('audio power')
    % figure(2); hold on;
    % for xi = 1:length(zerosL)-1
    %     % xi
    %     x1 = time(zerosL(xi));
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    % figure(1); hold on;
    % for xi = 1:length(peak_indices)-1
    %     x1 = peak_indices(xi)/fs;
    %     % x1 = peak_indices(xi);
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    disp('NUM CLICKS')
    clicks
    length(clicks)
    for xi = 1:length(clicks)
        x1 = time(clicks(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');

%         % figure(2); hold on;
%         % line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
% %~~~~~~~~~~~~~~~~PLOTTING END~~~~~~~~~~~~~~~
%         %% plot each individual clicks
%         figure(10+i); 
%         plot(time(clicks(i)-lenClick/2:clicks(i)+lenClick/2),data(clicks(i)-lenClick/2:clicks(i)+lenClick/2,:));
%         grid on;
    end
end %% function declaration  

