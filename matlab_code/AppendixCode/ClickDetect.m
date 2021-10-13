function [csig, clicks] = ClickDetect(sig)
    threshold = 200;
    clickwidth = 20;
    csig = sig;

    sep = 2048;
    s2 = sep/2;
    
    % replaces first two for loops in c++ code
    b2 = sig.^2; 
    ms_seq = b2; 

    for ii = (1:floor(log2(sep))) %possibly start at 2? 
        i = 2^ii;
        for j =(1:length(sig)-sep)

            ms_seq(j) = ms_seq(j) + ms_seq(j+i);

        end
    end

    ms_seq = ms_seq./sep;

    
    % threshold = 200; %in audacity runs from 200-900
    % clickwidth = 20; %in audacity runs from 20-40
    
    clicks = [];
    left = 0;
    wrc = clickwidth/4;
    while wrc >= 1
        wrc = wrc/2;
        ww = clickwidth/wrc;
        for i = (1:length(sig)-2*sep); 
            msw = 0;
            for j = (1:ww) 
                msw = msw + b2(i + s2 + j);
            end
            msw = msw/ww;
            if  msw >= threshold*ms_seq(i)/10
                clickdetected = 0;
                if left == 0
                    left = i + s2;
                end
            else 
                if(left ~= 0 && floor(i-left+s2) <= ww*2)
                    lv = sig(left);
                    rv = sig(i+ww+s2);
                    for j = (left:i+ww+s2)
                        if clickdetected == 0;
                            clicks = [clicks, j];
                            clickdetected = 1;
                        end
                        csig(j) = (rv*(j-left) + lv*(i+ww+s2-j))/(i+ww+s2-left);
                        b2(j) = csig(j).^2;
                    end
                    left = 0;
                elseif left ~= 0
                    left = 0;
                end
            end
        end
    end

end
