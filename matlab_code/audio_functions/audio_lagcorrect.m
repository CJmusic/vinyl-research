%{



%} 


function [cdata, ctime, cdata_ref, ctime_ref] = audio_lagcorrect(data, fs, data_ref, fs_ref, lagdiff);
     time = (0:length(data)); 
     time_ref = (0:length(data_ref));
%
%     diff_array = []; % this array contains the distances between every click, each row 
%                      % represents a click in the referenc each column represents a click in the file being looked at 
%     for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
%                                     % files
%         diff_array = [diff_array; clicks - clicks_ref(xi)];
%     end
%     lagdiff = mode(diff_array(:)); % the time difference between the two signals is the most common distance between clicks

    if lagdiff > 0;
        % positive lagdiff means that the data array is delayed compared to the reference
        cdata = data(lagdiff + 1:end,:); 
        ctime = time - lagdiff/fs; 
    elseif lagdiff < 0; 
        % negative lagdiff means that the data array is ahead of the reference  
        cdata = data(1: end + lagdiff,:);   
        ctime = time - lagdiff/fs;
    elseif lagdiff == 0; %then nothing needs to be corrected  
        cdata = data; %no lag 
        ctime = time - lagdiff/fs; 
    end

    size_diff = length(data_ref) - length(cdata)

    cdata_ref = data_ref(abs(size_diff)+1:end,:); %need to      
    %ctime = (0:length(cdata)); 
    ctime_ref = (0:length(cdata_ref));
end 
