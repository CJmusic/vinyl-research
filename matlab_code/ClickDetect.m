

function csig = ClickDetect(sig)
    csig = sig;


    sep = 2048;
    s2 = sep/2;

    % replaces first two for loops in c++ code
    b2 = sig.^2; 
    
    figure(3); grid on;
    plot(b2)
    
    ms_seq = b2; 
    size(ms_seq)
    % for i = (1:sep) %possibly start at 2? 
    for ii = (1:floor(log2(sep))) %possibly start at 2? 
        i = 2^ii;
        for j =(1:length(sig)-sep)
            % ms_seq(j)
            % ms_seq(j+i)
            ms_seq(j,1) = ms_seq(j,1) + ms_seq(j+i,1);
            ms_seq(j,2) = ms_seq(j,2) + ms_seq(j+i,2);
        end
    end
    size(ms_seq)
    ms_seq(:,1) = ms_seq(:,1)./sep;
    ms_seq(:,2) = ms_seq(:,2)./sep;
    disp('sizes')
    size(sig)
    size(ms_seq)

    ms_seq(1000:1100)


    figure(4); grid on;
    plot(ms_seq)
    


    
    threshold = 500; %in audacity runs from 200-900
    clickwidth = 30; %in audacity runs from 20-40
    
    clicks = [];
    left = 0;
    % while len - s > windowSize/2
    % for ww = (4:clickwidth) %% in audacity this runs from 4 
        
    %     for i = (1:length(sig)-sep);
    %         msw = 0;
    %         for j = (1:ww) %% msw is an array that contains the RMS value for the next blocks
    %             msw = msw + b2(i + s2 + j);
    %         end
    %         % msw = msw/ww
    %         % threshold
    %         % ms_seq(i)
    %         % thresh = (threshold*ms_seq(i))/10
    %         % left
    %         if  msw >= threshold*ms_seq(i)/10
    %             disp('THRESHOLD TRIGGERED')
    %             left = i + s2;
    %             % clicks = [clicks, i + s2]; %%??
    %         else 
    %             if(left ~= 0 && floor(i-left+s2) <= ww*2)
    %                 lv = sig(left);
    %                 rv = sig(i+ww+s2);
    %                 for j = (left:i+ww+s2)
    %                     disp('REMOVING CLICK')
    %                     %%detected click?
    %                     csig = rv*(j-left) + lv*(i+ww+s2-j)/(i+ww+s2-left);
    %                     %% perhaps I should have a cb2??
    %                     b2(j) = csig.^2;
    %                 end
    %                 % left
    %                 left = 0;
    %             % elseif left ~= 0;
    %             %     left = 0;
    %             end
    %         end
    %     end
    % end

end
