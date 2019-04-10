
%~ Loop through all the files and line them up with the reference

%{ LINEUP %}
% This function needs to take two arrays, and return two arrays that are lined up with one another. 
% It should: 
%       -  
%
%This function doesn't like stereo files as they come so I'll need to clean them up, then realign and redo it
%
%
%
%
%function [cdata, ctime] = audio_lineup(data_file, fs_file, time_file, data_ref)
function lagDiff = audio_lineup(data_file, data_ref, fs_file)
    disp('audio_lineup function')

    [acor_L,lags_L] = xcorr(data_file(:,1),data_ref(:,1));
    [acor_R,lags_R] = xcorr(data_file(:,2),data_ref(:,2));

    [~,I_L] = max(abs(acor_L));
    lagDiff_L = lags_L(I_L);
    [~,I_R] = max(abs(acor_R));
    lagDiff_R = lags_R(I_R);
    
    disp('lagDiff_L')
    lagDiff_L
    disp('lagDiff_R')
    lagDiff_R
    lagDiff = lagDiff_L
    disp ('ending audio_lineup function')
end

