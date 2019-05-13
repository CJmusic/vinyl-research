% bool EffectClickRemoval::ProcessOne(int count, WaveTrack * track, sampleCount start, sampleCount len)
% function aud_clickremoval(buffer)
function buffer = aud_clickremoval(buffer, fs)
    %%% Variables from other files

    % Param( Threshold, int,     wxT("Threshold"),  200,     0,       900,     1  );
    % Param( Width,     int,     wxT("Width"),      20,      0,       40,      1  );
    mThresholdLevel = 500;
    mClickWidth = 30;
    windowSize = 8192;
    sep = 2049;


    %%%%

    % in the matlab code buffer is the array of buffer
    % buffer = buffer; 
    len = length(buffer); 
    % I believe the above are the two inputs in the audacity functions;


    %% Audacity function RemoveClicks
    bResult = false; %says if the func did something
    % i = 0;
    % j = 0;
    left = 0;

    % msw = 0.0;
    ww = 0;

    ms_seq = zeros(length(buffer),2);
    b2 = zeros(len,2);

    % b2 is the power of the signal
    for i=(1:len) % calculate the rms level
        b2(i) = buffer(i)*buffer(i);
    end

    for i=(1:i*i:sep) %% in c++ 
        for j=(1:len-i)
            ms_seq(j) = ms_seq(j) + ms_seq(j+i);
        end
    end

    disp('GOING INTO 4 LOOP')
    % sep looks to be 2^8, probably the width of the samples being looked at
    for i=(1:len-sep)
        % size(ms_seq(i))
        % size(sep)
        % size(ms_seq)
        ms_seq(i) = sep/ms_seq(i);
    end

    % truncate sep to the next lowest power of two
    % sep appears to be the len
    % sep = i; 
    s2 = int16(sep/2); % originally assigned up above but might not need to do it, this is simply 
                % length of buffer/2  
    % wrc = 0;

    %%% this is a for loop in the C++ but probably needs to be converted to a while loop 
%    for(wrc=mClickWidth/4; wrc>=1; wrc /= 2) {
    wrc = mClickWidth/4;
    while wrc >= 1
        wrc = wrc/2;
        ww = mClickWidth/wrc; % 
    end %while wrc loop

    for i=(1:len-sep);
        msw = 0;
        for j=(1:ww+1)
            msw = msw + b2(i+s2+j); %% msw seems to be a sum of the power at weird indeces
                                    %  or using some weird c++ syntax 
        end
        msw = msw/ww;


        %% this is the point where clicks are detected and replaced in 
        %  the buffer = the list of buffer 
        %  ms_seq = main test for clicks? 
        if msw >= mThresholdLevel * ms_seq(i)/10;
            if left == 0;
                left = i+s2; %% this is where left gets assigned its initial value
            end
        else;
            if(left ~= 0 && (i-left) <= ww*2); %% ww ClickWidth/4
                lv = buffer(left);
                rv = buffer(i+ww+s2);
                for j = (left:i+ww+s2)
                    bResult = true;
                    %% The line below replaces the sample in the buffer with I think
                    %  is the average of the remaining samples being analyzed
                    buffer(j) = (rv*(j-left) + lv*(i+ww+s2-j))/(i+ww+s2-left);
                    b2(j) = buffer(j)*buffer(j)
                end
            elseif(left ~= 0)
                left = 0;
            end % if left ~=
        end %if msw >=
    end % end for loop
    % return bResult
    % Comment from Audacity writer
%    ww runs from about 4 to mClickWidth.  wrc is the reciprocal;
%    chosen so that integer roundoff doesn't clobber us.
    
















    % "below is the function ProcessOne translated, as far as I can tell it really just cleans up the wave buffer from Audacity and gets things ready to be processed. Not needed for our case" 
%     if(len <= windowSize/2);
%         disp("windowSize must be bigger than ", windowSize/2);
%     end

%     idealBlockLen = track.GetMaxBlockSize * 4; 

%     " as far as I can tell this is simply a maximum number of samples in the buffer "

% %    auto idealBlockLen = track->GetMaxBlockSize() * 4;
% %    if (idealBlockLen % windowSize ~= 0)
% %       idealBlockLen += (windowSize - (idealBlockLen % windowSize));

%     "the lines above confirm and simply set the block length to a factor of the max buffer size"

% %    bool bResult = true;
% %    decltype(len) s = 0;
% %    Floats buffer{ idealBlockLen };
% %    Floats bufferwindow{ windowSize };

%     "main iteration loop below: "

%     while (length(buffer) - s) > windowSize/2;
%         %% in the audacity code this is a call to a function, I've replaced it with the output of the actual function
%         block = [min( length( bufferSize ), 
%         max( sampleCount(0), length(buffer - s))];
        
%         track()
%         buffer = 0;
%         floatSample;
%         start + s 
%         block

% %       track->Get((samplePtr) buffer.get(), floatSample, start + s, block);

%         "unless this Get command is doing something special, this is simply returning the values as listed"

%         %% I believe that the block is the main array of audio buffer
%         for i = (1:windowSize/2:block);
%         % for (1: i+windowSize/2 - 1 < block, i+= windowSize/2);

% %       {
% %          auto wcopy = std::min( windowSize, block - i );
%             wcopy = min(windowSize, block - i);
%             for (1:wcopy-1);
%                 bufferwindow[j] = buffer[i+j];
%             end
%             for (wcopy:windowSize);
%                bufferwindow[j] = 0; 
%             end

%         %    mbDidSomething = or()      
%         %    if mbDidSomething = false; 
%             if RemoveClicks(windowSize, bufferwindow) ~= false;
%                 mbDidSomething = true;
%             end

%             for j=(0:wcopy)
%                 buffer[i+j] = bufferwindow;
%             end

%             "this mbDidSomething variable is probably useless in our needs"
%             % if mbDidSomething = True;
%             %     track.set()
%             %     samplePtr) buffer.get(), floatSample, start + s, block)
%             % end

%             s += block;
%             if count == 0 | len == 0;
%                 bResult = false;
%             break;


end %function aud_clickremoval