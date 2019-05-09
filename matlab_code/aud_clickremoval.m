% bool EffectClickRemoval::ProcessOne(int count, WaveTrack * track, sampleCount start, sampleCount len)
function aud_clickremoval(data);
    buffer = data; 
    len = length(data); 
    % I believe the above are the two inputs in the audacity functions;


    %% Audacity function RemoveClicks
    bResult = false; %says if the func did something
    i = 0;
    j = 0;
    left = 0;

    msw = 0.0;
    ww = 0;
    s2 = sep/2; 

    ms_seq = zeros(data);
    b2 = zeros(len);

    for i=(0,len); % calculate the rms level
        b2[i] = buffer[i]*buffer[i];
    end


    for i=(1:i*i:sep); %% in c++ 
        for j=(0:len-i);
            ms_seq[j] += ms_seq[j+i];
        end

    % truncate sep to the next lowest power of two
    sep = i;

    for i=(0:len-sep);
        ms_seq[i] = sep/ms_seq;
    end

    wrc = 0;

    %%% this is a for loop in the C++ but probably needs to be converted to a while loop 
%    for(wrc=mClickWidth/4; wrc>=1; wrc /= 2) {
    while wrc >=1; 
        wrc /= 2;
    end %while wrc loop


    % Comment from Audacity writer
%    ww runs from about 4 to mClickWidth.  wrc is the reciprocal;
%    chosen so that integer roundoff doesn't clobber us.
    
















    % "below is the function ProcessOne translated, as far as I can tell it really just cleans up the wave data from Audacity and gets things ready to be processed. Not needed for our case" 
%     if(len <= windowSize/2);
%         disp("windowSize must be bigger than ", windowSize/2);
%     end

%     idealBlockLen = track.GetMaxBlockSize * 4; 

%     " as far as I can tell this is simply a maximum number of samples in the buffer "

% %    auto idealBlockLen = track->GetMaxBlockSize() * 4;
% %    if (idealBlockLen % windowSize != 0)
% %       idealBlockLen += (windowSize - (idealBlockLen % windowSize));

%     "the lines above confirm and simply set the block length to a factor of the max buffer size"

% %    bool bResult = true;
% %    decltype(len) s = 0;
% %    Floats buffer{ idealBlockLen };
% %    Floats datawindow{ windowSize };

%     "main iteration loop below: "

%     while (length(data) - s) > windowSize/2;
%         %% in the audacity code this is a call to a function, I've replaced it with the output of the actual function
%         block = [min( length( bufferSize ), 
%         max( sampleCount(0), length(data - s))];
        
%         track()
%         buffer = 0;
%         floatSample;
%         start + s 
%         block

% %       track->Get((samplePtr) buffer.get(), floatSample, start + s, block);

%         "unless this Get command is doing something special, this is simply returning the values as listed"

%         %% I believe that the block is the main array of audio data
%         for i = (1:windowSize/2:block);
%         % for (1: i+windowSize/2 - 1 < block, i+= windowSize/2);

% %       {
% %          auto wcopy = std::min( windowSize, block - i );
%             wcopy = min(windowSize, block - i);
%             for (1:wcopy-1);
%                 datawindow[j] = buffer[i+j];
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
%                 buffer[i+j] = datawindow;
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