%{
Lining up clicks: 

clicks = [ sample numbers ] 

% a click has a sustained high derivative 
% the energy of a click is between 2-4 kHz

%} 

% TESTING THE CLICK DETECT FUNCTION
clc; close all;

files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/*.wav')
files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')

%de_files = dir('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/declicked/')

for i = (1:length(files))
    file = strcat(files(i).folder,'/',files(i).name)
    de_file = strcat(files(i).folder,'/declicked/',files(i).name) 
    %audio_bin = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/4.wav';
    %audio_bin2 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/click_testing/declicked/4.wav';
    %[data, time, fs] = audio_load(audio_bin);
    %[data2, time2, fs] = audio_load(audio_bin2);

    [data, time, fs] = audio_load(file);
    [data2, time2, fs] = audio_load(de_file);
    % tStart = 5.5;
    % tEnd = 10.1;
    % dataL = data(tStart*fs : tEnd*fs,1);
    % dataR = data(tStart*fs : tEnd*fs,2);
    dataL = data(:,1);
    dataR = data(:,2);
    time = (1:length(dataL))/fs;
    [clicksL] = audio_clickdetecttest(dataL, fs);
    [clicksR] = audio_clickdetecttest(dataR, fs);


    %~~~~~~~~~~~~~~~~PLOTTING~~~~~~~~~~~~~~~~~~~

    figure(i); grid on; hold on 
    plot(time, dataL)
    %plot(time, dataR)
    title('audio data')
    %ylim([-0.1 0.1])

    size(data)
    size(data2)
    declick = data2(:,1) - data(:,1); 
    plot(time,declick,'r');


    %figure(2); grid on; hold on;
    %plot(time, abs(dataL))
    %title('abs audio')


    % tStart = 5.1;
    % tEnd = 6.1;
    % data2 = data2(tStart*fs : tEnd*fs,1);
    data2 = data2(:,1);
    time2 = (1:length(data2))/fs;

    %figure(2); grid on; hold on 
    %plot(time, data2)
    %%ylim([-0.1 0.1])
    %title('audio data2')


    %figure(3); grid on; hold on 
    %plot(time, data - data2)
    %%ylim([-0.1 0.1])
    %title('audio data2')


    % PLOT CLICKS
    for xi = 1:length(clicksL)
        x1 = time(clicksL(xi));
        figure(i); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end

    saveas(figure(i),strcat(files(i).folder,i,'.png'))

    for xi = 1:length(clicksR)
        x1 = time(clicksR(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '--');
    end
    clicksL
    clicksR
    % PLOT INDIVIDUAL CLICKS
    figure(10+xi); 
    plot(time(clicks(xi)-lenClick/2:clicks(xi)+lenClick/2),data(clicks(xi)-lenClick/2:clicks(xi)+lenClick/2,:));
    grid on;
% % %~~~~~~~~~~~~~~~~PLOTTING END~~~~~~~~~~~~~~~


end %files for loop


function [clicks] = audio_clickdetecttest(data, fs)
% function [clicks] = audio_clickdetect(data, fs)
    % data = data(:,1);
%~~~~~~~~~~~~~~~~PRE-FILTERS~~~~~~~~~~~~~~~~~~~
    % freqLow = 2000;
    % freqHigh = 10000;
    % [b,a] = butter(6, freqLow/(fs/2), 'low');
    % data = filter(b, a, data); % use filtfilt

    % [b,a] = butter(6,freqHigh/(fs/2),'high');
    % data = filter(b, a, data); % use filtfilt
%~~~~~~~~~~~~~~~~PRE-FILTERS END~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~PEAKS METHOD~~~~~~~~~~~~~~~~~~~~~
    aData = abs(data);
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
    zerosIndices = zci(data);
    % size(zerosIndices)
    % size(aData)
    avg_peak = 0.0;
    prev_peak = 1.0;
    peak_value = 0.0;

    peakValues = zeros(length(zerosIndices)-1,1);   
    peakIndices = zeros(length(zerosIndices)-1,1);
    
    for i = (1:length(zerosIndices)-1)
        [peak, index] = max(aData(zerosIndices(i):zerosIndices(i+1)));
        peakIndices(i) = zerosIndices(i) + index;
        peakValues(i)  = peak;
    end
    % [peakValues, peakIndices] = findpeaks(aData(:,1));

    clicks = [];
    clicks_peaks = [];

    threshold = 20;
    disp('threshold')
    % threshold = std(peakValues);
    lenClick = 1412;
    mAvgWidth = 10;
    for i = (mAvgWidth+1:length(peakValues)-mAvgWidth-1)

        if i == mAvgWidth + 1
            iPeaks = peakValues(i-mAvgWidth : i+mAvgWidth);
            mAvg = (1/(mAvgWidth*2 + 1))*sum(iPeaks);
        end

        if peakValues(i) > threshold*mAvg
            %%% I also want to take into account the width of the peak,
            %   preferebly this should be by sample and not just peaks
            %   as then it will catch clicks that appear within music 
            if length(clicks) == 0;
                clicks = [clicks, peakIndices(i)];  
            elseif peakIndices(i) - clicks(length(clicks)) > lenClick %% makes sure the same click isn't recorded
                clicks = [clicks, peakIndices(i)];  
            end
        end
    end
    disp('END OF CLICKDETECT')
    size(clicks)
%~~~~~~~~~~~~~~~~~~~PEAKS METHOD END~~~~~~~~~~~~~~~~~~~~

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

%~~~~~~~~~~~~COMBI PEAKS/ENV METHOD END~~~~~~~~~~~~