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
    % pData = data.^2;
    pData = abs(data);
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    zerosL = zci(data(:,1));
    zerosR = zci(data(:,2));

    peakValues = zeros(length(zerosL)-1,2);
    peakIndices = zeros(length(zerosL)-1,2);
    % length(zerosL)
    % length(pData)
    % disp('last zero')
    % length(pData) - zerosL(length(zerosL))
    % length(pData) - length(zerosL)
    % zerosL(length(zerosL)-10:length(zerosL))
    size(pData)
    for i = (1:length(zerosL)-1)
        buffL = pData(zerosL(i):zerosL(i+1),1);
        buffR = pData(zerosL(i):zerosL(i+1),2);
        [peakL, indexL] = max(buffL);
        [peakR, indexR] = max(buffR);
        peakValues(i,1)  = peakL;
        peakIndices(i,1) = indexL;
        peakIndices(i,2) = indexR;

    end
    for i = (1:length(zerosR)-1)

        peakIndices(i,2) = indexR;
        peakValues(i,2)  = peakR;
        [peakR, indexR] = max(buffR);
    end


    % [peakValues, peakIndices] = findpeaks(pData(:,1));
    disp('Sizes')
    size(peakIndices)
    size(peakValues)
    size(zerosL)

    clicks = [];
    clicks_peaks = [];

    threshold = 1.2;
    lenClick = 1412;
    mAvgWidth = 2;
    avgBuffersize = 20; %% how many peaks to take into account for the moving avg
    
    N = 0;
    mAvg = 0;
    for i = (1:length(peakValues))
        N = i;
        % takes the moving average from the left
        mAvg = mAvg + peakValues(i,:)/N;

        if peakValues(i) > threshold*mAvg
            % disp('CLICK DETECTED')
            mAvg = mAvg - peakValues(i)/N; %% DON'T include clicks in the moving avg
            % if max(peakValues(i,1)) > max(peakValues(i,2))
            %     clpolarity = 0; %% click is in the left channel
            % else
            %     clpolarity = 1; %% click is in the right channel 
            % end

            if length(clicks) == 0;
                clicks = [clicks, peakIndices(i)];  
            elseif peakIndices(i) - clicks(length(clicks)) > lenClick %% makes sure the same click isn't recorded multiple times
                clicks = [clicks, peakIndices(i)];  
                % clicks = [clicks, zerosL(i)];  
            end
            % clicks = [clicks, [peakIndices(i-3:i+8)]]; %% What's recorded here is the indices of the peaks, perhaps the zeros make more sense
        end

        if N > avgBuffersize
            mAvg = mAvg - peakValues(i-avgBuffersize)/N;
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

    % pData = data.^2;
    % figure(2); grid on; hold on;
    % plot(time, pData)
    % title('audio power')
    % figure(2); hold on;
    % for xi = 1:length(zerosL)-1
    %     % xi
    %     x1 = time(zerosL(xi));
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end
    % figure(1); hold on;
    % for xi = 1:length(peakIndices)-1
    %     x1 = peakIndices(xi)/fs;
    %     % x1 = peakIndices(xi);
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

