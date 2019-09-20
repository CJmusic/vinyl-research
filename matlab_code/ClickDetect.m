

function csig = ClickDetect(sig)
    csig = sig;


    sep = 2048;
    s2 = sep/2;

    b2 = sig.^2;
    ms_seq = rms(sig);

    ms_seq = b2; 
    for i = (1:sep) %possibly start at 2? 
        for j =(1:length(sig)-1)
            ms_seq(j) = ms_seq(j) + ms_seq(j+i);
        end
    end
    ms_seq = ms_seq/sep;
    disp('sizes')
    size(sig)
    size(ms_seq)
    
    threshold = 500; %in audacity runs from 200-900
    clickwidth = 30; %in audacity runs from 20-40
    
    clicks = [];
    left = 0;
    %while len - s > windowSize/2
    for ww = (4:clickwidth) %% in audacity this runs from 4 
        for i = (1:length(sig)-sep);
            msw = 0;
            for j = (1:ww)
                msw = msw + b2(i + s2 + j);
            end
            msw = msw/ww;

            if  msw >= threshold*ms_seq(i)/10
                left = i + s2;
                % clicks = [clicks, i + s2]; %%??
            else 
                if(left ~= 0 && int(i-left+s2) <= ww*2)
                    lv = sig(left);
                    rv = sig(i+ww+s2);
                    for j = (left:i+ww+s2)
                        %%detected click?
                        csig = rv*(j-left) + lv*(i+ww+s2-j)/float(i+ww+s2-left)
                        %% perhaps I should have a cb2??
                        b2(j) = csig.^2;
                    end
                    left = 0;
                % elseif left ~= 0;
                %     left = 0;
                end
            end
        end
    end

end
