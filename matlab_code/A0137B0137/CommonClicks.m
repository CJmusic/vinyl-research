
function num_comclicks = CommonClicks(clicks, clicks_ref, lag_diff);
    
    if length(clicks) == 0; 
        num_comclicks = 0;
        return 
    end
    diff_array = []; % this array contains the distances between every click, each row 
                     % represents a click in the referenc each column represents a click in the file being looked at 
    for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
                                    % files
        diff_array(xi,:) = [clicks - clicks_ref(xi)];
    end
    lagDiff = mode(diff_array(:)); % the time difference between the two signals is the most common distance between clicks
    relaxation = floor(0.05*96000);

    diff_array
    lagDiff 
    for i = length(diff_array)
        if ismember(lag_diff, diff_array(i,:))
            CommonClicks = CommonClicks + 1;
        end
    end
    CommonClicks
end
