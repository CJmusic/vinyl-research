%% This file recreates Audacity's click removal in MATLAB. 
%
%
% christopher zaworski  
%
% last edit: may 19, 2019


%{

What I think is currently wrong:
   -  theres some messed up indexing, that's causing the clicks to be detected a 
      bit before their actual position

% bool EffectClickRemoval::ProcessOne(int count, WaveTrack * track, sampleCount start, sampleCount len)
% function testaudClickRemoval(data)


%}

% ~~~~~~~~~~~~~~~~~~~~~TESTING~~~~~~~~~~~~~~~~~~~~~ %{

clf(figure(1));clf(figure(2));
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/040319_r26fivetrials/');
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_functions/');
% addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/');
%record_dir = dir('');

% wave_files = dir(strcat(path_folder,'*.wav'));

file_path = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/030419_r26fivetrials/one.wav';

coh_start = 8.0;
coh_end = 18.0;
% wave_names{i} = sprintf(wave_files(i).name)
record = audio_recordclass(file_path);
% reference = record;
% clicks_ref = audio_clickdetect(reference.data, reference.fs);
%     clicks_ref = testaudClickRemoval(reference.data, reference.fs);
rec_cohere = record.data(coh_start*record.fs:coh_end*record.fs,:); 
recTime = (0:length(rec_cohere)-1)/record.fs; 

disp('BUFFER INFO')
[data_aud, clicks] = testaudClickRemoval(record.data(:,1));
size(data_aud)

fig1 = figure(1)
plot(record.time, record.data(:,1))
grid on;

fig2 = figure(2)
plot(record.time, data_aud)
grid on;

% ~~~~~~~~~~~~~~~~~~~TESTING END~~~~~~~~~~~~~~~~~~~ %

function [data, clicks] = testaudClickRemoval(data, fs)
    %%% Variables from other files

    % Param( Threshold, int,     wxT("Threshold"),  200,     0,       900,     1  );
    % Param( Width,     int,     wxT("Width"),      20,      0,       40,      1  );
    mThresholdLevel = 500;
    mClickWidth = 30;
    windowSize = 8192;
    sep = 2049;
    len = 2^16;
    % len = length(data); % len is the length of the buffer in audacity,  

    %% Audacity function RemoveClicks
    bResult = false; % says if the func did something,
                     % 
    % i = 0;
    % j = 0;
    left = 0;

    s2 = int16(sep/2); % originally assigned up above but might not need to do it, this is simply 
    % msw = 0.0;
    ww = 0;

    % my code to could clicks and list them in an array
    % num_clicks = 0;
    % clicks = [];

    data_power = zeros(len,2);  %%data_power is the signals energy

    % data_power is the power of the signal
    for i=(1:len) % calculate the rms level
        data_power(i) = data(i)*data(i);
    end

    % ms_seq = zeros(length(data),2);
    ms_seq = data_power;
    % for i = (1:len)
        % ms_seq(i]=b2[i];

    % for i=(1:i*i:sep) %% in c++ 
    i = 1;

    while i <= sep
        i = i*2;
        for j=(1:len-i)
            ms_seq(j) = ms_seq(j) + ms_seq(j+i);
        end
    end
    disp('after ms_seq set')
    % ms_seq

    sep = i;

    disp('GOING INTO 4 LOOP')
    % sep looks to be 2^8, probably the width of the samples being looked at
    disp('length data')
    size(data)
    disp('len-sep')
    len-sep

    clicks = [];


    for i=(1:len-sep)
        % size(ms_seq(i))
        % size(sep)
        % size(ms_seq)
        ms_seq(i) = ms_seq(i)/sep;
    end

    % truncate sep to the next lowest power of two
    % sep appears to be the len
    % sep = i; 
                % length of data/2  
    % wrc = 0;

    %%% this is a for loop in the C++ but probably needs to be converted to a while loop 
%    for(wrc=mClickWidth/4; wrc>=1; wrc /= 2) {
    wrc = mClickWidth/4;
    disp('Final While loop')
    while wrc >= 1
        %% wrc looks like a resolution, so it loops through the file multiple times at each 
        %  ClickWidth (ie: resolution) and tests it to see if there is a click there
        wrc = wrc/2
        
        %% ww is the variable ClickWidth without units
        ww = mClickWidth/wrc  

        for i=(1:len-sep)
            %% len-sep simply avoids the last iteration I believe
            msw = 0;
            for j=(1:ww+1)
                msw = msw + data_power(i+s2+j); %% msw seems to be a sum of the power at weird indices
                                        %  or using some weird c++ syntax 
            end
            msw = msw/ww;

            %% this is the point where clicks are detected and replaced in 
            %  data = the list of data 
            %  ms_seq = main test for clicks? 
            if msw >= mThresholdLevel * ms_seq(i)/10
                %%  if theres a click detected in the msw, then I believe this loop 
                %   goes through in more detail to locate the click more precisely
                disp('CLICK DETECTED?')
                data(i) = 0;
                if left == 0
                    % disp('LEFT 0')
                    % disp('CLICK DETECTED?')
                    left = i+s2; %% this is where left gets assigned its initial value
                end
                elseif (left ~= 0 && (i-left) <= ww*2); %% ww ClickWidth/4
                % disp('i-left')
                % i-left
                % disp('ww*2')
                % ww*2
                % if ((i-left) <= ww*2); %% ww ClickWidth/4
                    % disp('CLICK BEING REMOVED')
                    lv = data(left);
                    rv = data(i+ww+s2);
                    for j = (left:i+ww+s2)
                        bResult = true;
                        %% The line below replaces the sample in the data with I think
                        %  is the average of the remaining samples being analyzed
                        % data(j) = (rv*(j-left) + lv*(i+ww+s2-j))/(i+ww+s2-left);
                        data(j) = 1.0;
                        % data_power(j) = data(j)*data(j);
                    end
                left = 0;
                elseif(left ~= 0)
                    left = 0;
                % end % if left ~=
            end %if msw >=

            % data(i) = 0;
        end %while wrc loop
    end % end for loop
    % return bResult
    % Comment from Audacity writer
%    ww runs from about 4 to mClickWidth.  wrc is the reciprocal;
%    chosen so that integer roundoff doesn't clobber us.
    
















    % "below is the function ProcessOne translated, as far as I can tell it really just cleans up the wave data from Audacity and gets things ready to be processed. Not needed for our case" 
%     if(len <= windowSize/2);
%         disp("windowSize must be bigger than ", windowSize/2);
%     end

%     idealBlockLen = track.GetMaxBlockSize * 4; 

%     " as far as I can tell this is simply a maximum number of samples in the data "

% %    auto idealBlockLen = track->GetMaxBlockSize() * 4;
% %    if (idealBlockLen % windowSize ~= 0)
% %       idealBlockLen += (windowSize - (idealBlockLen % windowSize));

%     "the lines above confirm and simply set the block length to a factor of the max data size"

% %    bool bResult = true;
% %    decltype(len) s = 0;
% %    Floats data{ idealBlockLen };
% %    Floats datawindow{ windowSize };

%     "main iteration loop below: "

%     while (length(data) - s) > windowSize/2;
%         %% in the audacity code this is a call to a function, I've replaced it with the output of the actual function
%         block = [min( length( dataSize ), 
%         max( sampleCount(0), length(data - s))];
        
%         track()
%         data = 0;
%         floatSample;
%         start + s 
%         block

% %       track->Get((samplePtr) data.get(), floatSample, start + s, block);

%         "unless this Get command is doing something special, this is simply returning the values as listed"

%         %% I believe that the block is the main array of audio data
%         for i = (1:windowSize/2:block);
%         % for (1: i+windowSize/2 - 1 < block, i+= windowSize/2);

% %       {
% %          auto wcopy = std::min( windowSize, block - i );
%             wcopy = min(windowSize, block - i);
%             for (1:wcopy-1);
%                 datawindow[j] = data[i+j];
%             end
%             for (wcopy:windowSize);
%                datawindow[j] = 0; 
%             end

%         %    mbDidSomething = or()      
%         %    if mbDidSomething = false; 
%             if RemoveClicks(windowSize, datawindow) ~= false;
%                 mbDidSomething = true;
%             end

%             for j=(0:wcopy)
%                 data[i+j] = datawindow;
%             end

%             "this mbDidSomething variable is probably useless in our needs"
%             % if mbDidSomething = True;
%             %     track.set()
%             %     samplePtr) data.get(), floatSample, start + s, block)
%             % end

%             s += block;
%             if count == 0 | len == 0;
%                 bResult = false;
%             break;


end %function testaudClickRemoval