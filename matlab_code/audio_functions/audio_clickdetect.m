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

AUDIO_FILES = {'Bcorrelation_test_1.wav','Bcorrelation_test_2.wav','Bcorrelation_test_3.wav'};

[data, time, fs] = audio_load(audio_bin);
tStart = 5.1;
tEnd = 10.1;
data = data(tStart*fs : tEnd*fs,:);
audio_clickdetecttest(data, fs);

function clicks = audio_clickdetecttest(data, fs)
    % data = data(:,1);
    time = (1:length(data))/fs;
%~~~~~~~~~~~~~~~~PRE-FILTERS~~~~~~~~~~~~~~~~~~~
    % freqLow = 2000
    % freqHigh = 10000
    % [b,a] = butter(6,freqLow/(fs/2),'low');
    % data = filter(b, a, data); % use filtfilt

    % [b,a] = butter(6,freqHigh/(fs/2),'high');
    % data = filter(b, a, data); % use filtfilt
%~~~~~~~~~~~~~~~~PRE-FILTERS END~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~PEAKS METHOD~~~~~~~~~~~~~~~~~~~~~
    aData = abs(data);
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    zero_indices = zci(data);
    size(zero_indices)
    size(aData)
    avg_peak = 0.0;
    prev_peak = 1.0;
    peak_value = 0.0;

    peak_values = zeros(length(zero_indices)-1,1);   
    peak_indices = zeros(length(zero_indices)-1,1);
    
    [peak_values, peak_indices] = findpeaks(aData(:,1));

    clicks = [];
    clicks_peaks = [];

    threshold = 5;
    lenClick = 1412;
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

        if peak_values(i) > threshold*mAvg
            % if length(peak_values(peak_values(i : i+10) > 0.5*peak_values(i))) < 2;
            %     continue    
            % end

            if length(clicks) == 0;
                clicks = [clicks, peak_indices(i)];  
            elseif peak_indices(i) - clicks(length(clicks)) > lenClick %% makes sure the same click isn't recorded
                clicks = [clicks, peak_indices(i)];  
            end

            % clicks = [clicks, [peak_indices(i-3:i+8)]]; %% What's recorded here is the   
                                                        %  indices of the peaks, perhaps
                                                        %  the zeros make more sense
        end
    end
%~~~~~~~~~~~~~~~~~~~PEAKS METHOD END~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~COMBI PEAKS/ENV METHOD~~~~~~~~~~~~~~
    % aData = abs(data);
    % zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    % zerosL = zci(data(:,1));
    % zerosR = zci(data(:,2));

    % peakValues = zeros(length(zerosL)-1,2);
    % peakIndices = zeros(length(zerosL)-1,2);
    % size(aData)
    % for i = (1:length(zerosL)-1)
    %     buffL = aData(zerosL(i):zerosL(i+1),1);
    %     buffR = aData(zerosL(i):zerosL(i+1),2);

    %     [peakL, indexL] = max(buffL);
    %     [peakR, indexR] = max(buffR);
        
    %     peakValues(i,1)  = peakL;
    %     peakValues(i,1)  = peakL;

    %     peakIndices(i,1) = indexL;
    %     peakIndices(i,2) = indexR;

    % end
    % for i = (1:length(zerosR)-1)

    %     peakIndices(i,2) = indexR;
    %     peakValues(i,2)  = peakR;
    %     [peakR, indexR] = max(buffR);
    % end

    % clicks = [];

    % threshold = 1.2;
    % lenClick = 1412;
    % mAvgWidth = 2;
    % avgBuffersize = 200; %% how many peaks to take into account for the moving avg
    
    % N = 0;
    % mAvg = 0;
    % for i = (1:length(peakValues))
    %     N = i;
    %     if i > avgBuffersize;
    %         N = avgBuffersize;
    %     end
    %     % takes the moving average from the left
    %     mAvg = mAvg + peakValues(i,:)/N;

    %     if peakValues(i) > threshold*mAvg
    %         % disp('CLICK DETECTED')
    %         mAvg = mAvg - peakValues(i)/N; %% DON'T include clicks in the moving avg
    %         if length(clicks) == 0;
    %             clicks = [clicks, peakIndices(i)];  
    %         elseif peakIndices(i) - clicks(length(clicks)) > lenClick %% makes sure the same click isn't recorded multiple times
    %             clicks = [clicks, peakIndices(i)];  
    %             % clicks = [clicks, zerosL(i)];  
    %         end
    %         % clicks = [clicks, [peakIndices(i-3:i+8)]]; %% What's recorded here is the indices of the peaks, perhaps the zeros make more sense
    %     end

    %     if N > avgBuffersize
    %         mAvg = mAvg - peakValues(i-avgBuffersize)/N;
    %     end
    % end



%~~~~~~~~~~~~~~~~PLOTTING~~~~~~~~~~~~~~~~~~~
    

    figure(1); grid on; hold on 
    plot(time, data)
    title('audio data')
    ylim([-0.1 0.1])

    figure(4); grid on; hold on 
    plot(time, data)
    title('audio data')
    ylim([-0.1 0.1])

    figure(2); grid on; hold on;
    plot(time, aData)
    title('abs audio')

    audio_bin2 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/r26-96kHz-declicked.wav';
    [data2, time2, fs] = audio_load(audio_bin2);
    tStart = 5.1;
    tEnd = 10.1;
    data2 = data2(tStart*fs : tEnd*fs,:);
    time2 = (1:length(data2))/fs;

    figure(3); grid on; hold on 
    plot(time, data2)
    ylim([-0.1 0.1])
    title('audio data2')
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

    for xi = 1:length(clicks)
        x1 = time(clicks(xi));


        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');

        figure(2); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');

        % figure(10+xi); 
        % plot(time(clicks(xi)-lenClick/2:clicks(xi)+lenClick/2),data(clicks(xi)-lenClick/2:clicks(xi)+lenClick/2,:));
        % grid on;
    end
% %~~~~~~~~~~~~~~~~PLOTTING END~~~~~~~~~~~~~~~
end %% function declaration  


%~~~~~~~~~~~~~~~~~~~~~~~~~DECAY ENV METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{ 
an alternate method of detecting clicks based on the paper by 
Kinzie and Graveraux 1971

currently does not work

%}

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