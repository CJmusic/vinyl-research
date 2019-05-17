% This file will find local peaks and average them to
% get the amplitude of the sine wave.
%
%
% christopher zaworski
% last edit : march 31, 2019
%
%

function [amplitude] = audio_findamplitude(data,freq,fs)
    % frequency is the frequency of the sine wave 
    % you're looking for
    wav_peaks = [];
    for i = (1, length(data)-1,freq*fs);
       wav_peaks(i) = max(data(i:i+freq*fs)); 
    end % for loop
    amplitude = mean(peak);

end % function amplitude