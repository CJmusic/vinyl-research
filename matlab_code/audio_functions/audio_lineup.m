
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

   % [acor,lag] = xcorr(data_file_L,data_ref_L);
    [acor_L,lags_L] = xcorr(data_file(:,1),data_ref(:,1));
    [acor_R,lags_R] = xcorr(data_file(:,2),data_ref(:,2));
    figure(50); hold on; grid on; 
    plot(acor_L)
    plot(acor_R,'--')

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
  %  cdata_file_L = data_file_L(lagDiff+1:end);
  %  ctime_file = (0:length(cdata_file_L)-1)/fs_file;

  %  %[acor,lag] = xcorr(data_file_R,data_ref_R);
  %  %[~,I] = max(abs(acor));
  %  %lagDiff = lag(I)
  %  %timeDiff = lagDiff/fs_file
  %  %cdata_file_R = data_file_R(lagDiff+1:end);
  %  %ctime_file = (0:length(cdata_file_R)-1)/fs_file;

  %  cdata = data_file(lagDiff+1:end,:);

  %  %cdata = [cdata_file_L, cdata_file_R];
  %  ctime = ctime_file;


%
%    size(data_file_L)
%    size(data_ref_R)
%    [acor_L,lag_L] = xcorr(data_file_L,data_ref_L);
%    [~,I_L] = max(abs(acor_L));
%    lagDiff_L = lag(I_L)
%    timeDiff_L = lagDiff_L/fs_file
%    cdata_file_L = data_file_L(lagDiff_L+1:end);
%    ctime_file_L = (0:length(cdata_file_L)-1)/fs_file;
%    
%    [acor_R,lag_R] = xcorr(data_file_R,data_ref_R);
%    [~,I_R] = max(abs(acor_R));
%    lagDiff_R = lag(I_R)
%    timeDiff_R = lagDiff_R/fs_file
%    cdata_file_R = data_file(lagDiff_R+1:end);
%    ctime_file_R = (0:length(cdata_file_R)-1)/fs_file;
%
%    cdata_file = [cdata_file_L, cdata_file_R];
%    ctime_file = ctime_file_L

end

