
function num_comclicks = num_comclicks(clicks, clicks_ref, lag_diff);
    %loop through clicks 
    %look for other clicks in
    disp('INSIDE COMMON CLICKS')    
    clicks
    clicks_ref
    
    size(clicks)
    size(clicks_ref)



    % find(clicks_ref(i) - relaxation < clicks < clicks_ref(i) + relaxation)
    num_comclicks = 0;
    relaxation = 100;
    for i = (1:length(clicks_ref))

        comclicks = (clicks_ref(i) - relaxation < clicks < clicks_ref(i) + relaxation);
        % if find(clicks_ref(i) - relaxation < clicks < clicks_ref(i) + relaxation) ~= 0

        % num_comclicks = num_comclicks + 1;
        num_comclicks = length(comclicks);
        % end 
        % num_comclicks
    end
   
    % disp('INSIDE COMMON CLICKS')    
    % num_comclicks = 0;
    % if length(clicks) == 0; 
    %     return 
    % end
    % diff_array = []; % this array contains the distances between every click, each row 
    %                  % represents a click in the referenc each column represents a click in the file being looked at 
    % for xi = (1:length(clicks_ref));%this makes an array with the distance between each click in the two
    %                                 % files
    %     diff_array(xi,:) = [clicks - clicks_ref(xi)];
    % end
    % lag_diff
    % lagDiff = mode(diff_array(:)) % the time difference 
    % diff_array
    % between the two signals is the most common distance between clicks
    % relaxation = floor(0.05*96000);
    % lag_diff

    % for i = size(diff_array,1)
    %     for j = size(diff_array,2)
    %          if (lag_diff) - relaxation < (diff_array(i,j)) < (lag_diff) + relaxation
    %         % if relaxation < diff_array(i,j) < relaxation
    %             disp('FOUND COMMON CLICK')
    %             diff_array(i,j)
    %             lag_diff
    %             num_comclicks = num_comclicks + 1;
    %             % continue
    %         end
    %     end
    % end
    % num_comclicks
    % disp('COMMON CLICKS ENDS')    
end
