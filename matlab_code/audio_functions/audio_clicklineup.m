

% this function looks through the clicks in two audio files and attempts to calculate the 
% time delay between them  

function lagdiff = audio_clicklineup(clicks_ref, clicks);
    diff_array = []; % this array contains the distances between every click, each row 
                     % represents a click in the referenc each column represents a click in the file being looked at 
    for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
                                    % files
        diff_array = [diff_array; clicks - clicks_ref(xi)];
    end
    lagdiff = mode(diff_array(:)); % the time difference between the two signals is the most common distance between clicks
end
