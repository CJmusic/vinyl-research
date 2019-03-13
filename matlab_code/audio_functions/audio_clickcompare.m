% audio_clickcompare
% christopher zaworski
% last edit : march 12 2018
%
% This file counts and compares the number of clicks in two click matrices



function info = audio_clickcompare(clicks_ref, clicks);
    click_mat = audio_clickmatrix(clicks_ref, clicks);
    lagdiff = audio_clicklineup(clicks_ref, clicks);    
    ismember(click_mat, lagdiff); % find all occurances of the lag difference, hopefully sorted 
                                  % column row    
    clicks_com = []; % clicks common to both file
    info = ''    
end
