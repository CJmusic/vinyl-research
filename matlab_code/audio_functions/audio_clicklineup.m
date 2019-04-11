

% this function looks through the clicks in two audio files and generates the click matrix 
% last edit : march 17 2018

%function (cdata, cdata_ref) = audio_clicklineup(data, fs, data_ref, fs_ref);
%function lagdiff = audio_clicklineup(clicks, clicks_ref);
%
%     diff_array = []; % this array contains the distances between every click, each row 
%                      % represents a click in the referenc each column represents a click in the file being looked at 
%     for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
%                                     % files
%         diff_array = [diff_array; clicks - clicks_ref(xi)];
%     end
%     lagdiff = mode(diff_array(:)); % the time difference between the two signals is the most common distance between clicks
%end    


%    Code below left for reference
function [lagdiff] = audio_clicklineup(data, fs, clicks_ref, clicks);
   % diff_array = []; % this array contains the distances between every click, each row 
   %                  % represents a click in the referenc each column represents a click in the file being looked at 
   % for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
   %                                 % files
   %     diff_array = [diff_array; clicks - clicks_ref(xi)];
   % end
   % lagdiff = mode(diff_array(:)); % the time difference between the two signals is the most common distance between clicks

    %time_ref = (1:length(data_ref))/fs_ref;
    if nargin < 4; 
        [clicks] = audio_clickdetect(data, fs);
    end
        
        
    time = (1:length(data))/fs;
        %[clicks_ref] = audio_clickdetect(data_ref, fs_ref);

     diff_array = [] % this array contains the distances between every click, each row 
                      % represents a click in the referenc each column represents a click in the file being looked at 
     for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
                                     % files
         diff_array = [diff_array; clicks - clicks_ref(xi)];
     end
     lagdiff = mode(diff_array(:)) % the time difference between the two signals is the most common distance between clicks

    % cdata = circshift(data, lagdiff, 1);
    % ctime = (0:length(cdata)-1)/fs;
    % size_diff = length(ctime) - length(cdata)    


%    if lagdiff > 0;
%        % positive lagdiff means that the data array is delayed compared to the reference
%        cdata = data(lagdiff + 1:end,:); 
%%        ctime = time - lagdiff/fs; 
%    elseif lagdiff < 0; 
%        % negative lagdiff means that the data array is ahead of the reference  
%        cdata = data(1: end + lagdiff,:);   
%%        ctime = time - lagdiff/fs;
%    elseif lagdiff == 0; %then nothing needs to be corrected  
%        cdata = data; %no lag 
%%        ctime = time - lagdiff/fs; 
%    end


    % ctime = (0:length(cdata)-1)/fs;
    % size_diff = length(ctime) - length(cdata)    

    %size_diff = length(data_ref) - length(cdata)
%    cdata_ref = data_ref(abs(size_diff)+1:end,:); %need to      
%    ctime_ref = (1:length(data_ref))/fs_ref;

    %if size_diff > 0; %% the following is to ensure the data arrays are the same length for calculating the coherence
    %                  %% probably this isn't necessary after the lag diff is properly handled 
    %    disp('>0'); 
    %    cdata_ref = data_ref(size_diff+1:end,:); 
    %elseif size_diff < 0; 
    %    cdata = cdata(size_diff,:);
    %    cdata_ref = data_reff;
    %else;
    %    cdata_ref = data_ref;
    %end

   % disp('sizes corrected')
   % size(cdata)
   % size(cdata_ref)

   % [amp_coh, freq_coh] = audio_mscohere(cdata_ref, cdata, fs);

   % % plotting below : 
   % figure(1);
   % grid on; hold on;
   % plot(time,data(:,1));%, 'Color', [i/5,i/5,i/5] );
   % legend(AUDIO_FILES);
   % title('Click detection: Amplitude vs. Time'); 
   % xlabel('Time [s]');
   % ylabel('Amplitude');

   % figure(2); grid on; hold on;
   % %plot(freq_coh, amp_coh); 
   % %set(gca, 'XScale', 'log');
   % semilogx(freq_coh, amp_coh);
   % legend(AUDIO_FILES);
   % xlabel('frequency [Hz]')
   % title('Coherence compared to first recording')
end
