% This file will find local peaks and average them to
% get the amplitude of the sine wave.
%
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function amplitude = audio_findamplitude(data,freq,fs);
    % frequency is the frequency of the sine wave 
    % you're looking for
    peak = [];
    for i = (1, length(data),freq*fs);
       peak(i) = max(data); 
    end % for loop
    amplitude = mean(peak);

end % function amplitude