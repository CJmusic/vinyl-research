% function [csig, clicks] = ClickDetectTest(sig)
% function [csig, clicks] = ClickDetect(sig)
function [csig, clicks] = ClickDetect(sig)
    threshold = 200;
    clickwidth = 20;
    csig = sig;

    sep = 2048;
    s2 = sep/2;
    
    % replaces first two for loops in c++ code
    b2 = sig.^2; 
    ms_seq = b2; 

    % for i = (1:sep) %possibly start at 2? 
    for ii = (1:floor(log2(sep))) %possibly start at 2? 
        i = 2^ii;
        for j =(1:length(sig)-sep)
            % ms_seq(j,1) = ms_seq(j,1) + ms_seq(j+i,1);
            % ms_seq(j,2) = ms_seq(j,2) + ms_seq(j+i,2);

            ms_seq(j) = ms_seq(j) + ms_seq(j+i);

        end
    end
    % ms_seq(:,1) = ms_seq(:,1)./sep;
    % ms_seq(:,2) = ms_seq(:,2)./sep;

    ms_seq = ms_seq./sep;


    
    % threshold = 200; %in audacity runs from 200-900
    % clickwidth = 20; %in audacity runs from 20-40
    
    clicks = [];
    left = 0;
    % while len - s > windowSize/2
    % for wrc = (clickwidth/4:1)
    % for ww = (4:clickwidth) %% in audacity this runs from 4 
    %     wrc = clickwidth/ww;
    wrc = clickwidth/4;
    while wrc >= 1
        wrc = wrc/2;
        ww = clickwidth/wrc;
        for i = (1:length(sig)-2*sep); %% NEED TO TEST IF THIS FIXES CLICK AT END BUG 
        % for i = (1:length(sig)-sep);
            msw = 0;
            for j = (1:ww) 
                msw = msw + b2(i + s2 + j);
            end
            msw = msw/ww;
            if  msw >= threshold*ms_seq(i)/10
                clickdetected = 0;
                if left == 0
                    left = i + s2;
                    % clicks = [clicks, i + s2]; %%??
                end
            else 
                if(left ~= 0 && floor(i-left+s2) <= ww*2)
                    lv = sig(left);
                    rv = sig(i+ww+s2);
                    for j = (left:i+ww+s2)
                        % disp('REMOVING CLICK')
                        %%detected click?
                        if clickdetected == 0;
                            clicks = [clicks, j];
                            clickdetected = 1;
                        end
                        csig(j) = (rv*(j-left) + lv*(i+ww+s2-j))/(i+ww+s2-left);
                        %% perhaps I should have a cb2??
                        b2(j) = csig(j).^2;
                    end
                    left = 0;
                elseif left ~= 0
                    left = 0;
                end
            end
        end
    end

    % for i = (1:length(clicks)-1)
    %     abs(clicks(i)-clicks(i+1))
    %     abs(clicks(i)-length(sig))
    %     if abs(clicks(i) - clicks(i+1)) < 50000;
    %         clicks(i) = [];
    %     end 
    %     if abs(clicks(i) - length(sig)) < 50000;
    %         clicks(i) = [];
    %     end    
    % end

    % if length(clicks) > 0
    %     if abs(clicks(end) - length(sig)) < 50000;
    %         clicks(end) = [];
    %      end    
    % end



end
