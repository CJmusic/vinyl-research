% clear all;close all;clc
% set(0,'DefaultLineLineWidth',0.5);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)


% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
% addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')
% addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
% addpath('/Volumes/AUDIOBANK/audio_files/A0000B0000/')

% % file_ref = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav';
% % file_ref = '/Volumes/AUDIOBANK/audio_files/A0137B0137/074a.wav';
% file_ref = '/Volumes/AUDIOBANK/audio_files/A0137B0137/181b1558.466.wav';

% [reference, ~] = SeperateTracks(file_ref);

% sig = reference('sweep');

% % clickwidths = [5,10,15,20,25,30,35,40];
% clickwidths = 20;
% a = (100:10:200);
% b = (200:20:600);
% c = 1000;
% thresholds = [a,b,c];
% thresholds = (20:20:600)
% fig = figure('Visible','off')
% csigs_L = [];
% csigs_R = [];
% numclicks_L = [];
% numclicks_R = [];
% CLICKS_L = {};
% CLICKS_R = {};

% for c = 1:length(clickwidths)
%     numclicks_L = [];
%     disp(strcat('clickwidth... ', num2str(clickwidths(c))))
%     for t = (1:length(thresholds))
%         [csig_L, clicks_L] = ClickDetectTest(sig(:,1), thresholds(t), clickwidths(c));
%         [csig_R, clicks_R] = ClickDetectTest(sig(:,2), thresholds(t), clickwidths(c));
%         disp(strcat('clicksdetected...', num2str(length(clicks_L))))
%         numclicks_L = [numclicks_L, length(clicks_L)];
%         numclicks_R = [numclicks_R, length(clicks_R)];

%         % csig = [csig_L, csig_R];
%         % titlestring = strcat('Clicks detected at threshold=',num2str(thresholds(t)));
%         % plot_waveclicks(csig, clicks_L,clicks_R, titlestring)


%         % disp(strcat('numclicks... ', num2str(numclicks)))
%     end
%     % scatter(thresholds,numclicks_L);
%     % hold on; grid on;
%     % plot(thresholds,numclicks_R,'--')
%     % hold on; grid on;

%     scatter(thresholds,numclicks_L,'ko'); hold on; grid on;
%     scatter(thresholds,numclicks_R,'kx'); hold on; grid on;
% end

% % legend({'5','10','15','20','25','30','35','40'})
% legend({'left channel', 'right channel'})
% ylim([0,400])
% titlestring = 'Clicks vs. Threshold sweep track';
% title(titlestring)
% saveas(fig, strcat('plots/ClickDetect/',titlestring,'.png'))



% thresholds = 600:-50:100;

% prev_clicks = []
% for t = (1:length(thresholds))
%     [csig_L, clicks_L] = ClickDetectTest(sig(:,1), thresholds(t), 20);
%     disp(strcat('clicks detected.... ', num2str(length(clicks_L))))
%     new_clicks = setdiff(clicks_L, prev_clicks);
%     disp(strcat('new clicks.... ', num2str(length(new_clicks))))

%     prev_clicks = clicks_L;


%     for c = 1:length(new_clicks)
%         plot_click(new_clicks(c), sig, char(strcat('click',{' '}, num2str(c),' detected at threshold', {' '},num2str(thresholds(t)))))
%         if c > 5
%             break 
%         end
%     end

% end





% function plot_click(click, sig, titlestring)
%     fs = 96000;
%     start_sam = click - 0.0005*fs;
%     if start_sam < 0
%         start_sam = 1;
%     end

%     end_sam = click + 0.002*fs;
%     if end_sam > length(sig)
%         end_sam = length(sig);
%     end


%     click = sig(start_sam:end_sam,1);



%     fs = 96000;
%     time = (1:length(click))/fs;
%     time = time*1000;

%     fig = figure('Visible', 'off')
%     plot(time,click,'k')
%     grid on;
%     ylim([-0.5,0.5])
%     title(titlestring)
%     xlabel('time [ms]')
%     ylabel('signal level')
%     plotname = strcat('plots/ClickDetect/clickplots/', titlestring, '.png')
%     saveas(fig, plotname);

% end

% function plot_waveclicks(sig, clicks_L, clicks_R, filename)
%     fig = figure('Visible', 'off')
%     subplot(2,1,1)
%     plot(sig(:,1),'k')
%     ylim([-0.2, 0.2])
%     hold on; grid on;
%     % plot(sig2)
%     title(strcat(filename, ' left channel'))
    
%     for xi = 1:length(clicks_L)
%         x1 = (clicks_L(xi));
%         line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
%     end

%     subplot(2,1,2)
%     plot(sig(:,2),'k')
%     ylim([-0.2, 0.2])
%     hold on; grid on;
%     % plot(sig2)
%     title(strcat(filename, ' right channel'))

%     for xi = 1:length(clicks_R)
%         x1 = (clicks_R(xi));
%         line([x1 x1], get(gca, 'ylim'),'Color', 'green','LineStyle', '--');
%     end

%     plotname = strcat('plots/CommonClicks/', filename, '.png')
%     saveas(fig, plotname);
% end

filename = '/Volumes/AUDIOBANK/audio_files/A0137B0137/298b1600.386.wav';

filename = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a1558.066.wav'

[tracks, info_array] = SeperateTracks(filename);

signal_names = tracks.keys;
signals = tracks.values;

clicks_L = {};
clicks_R = {};

for t = 1:length(signal_names);
    sig = signals{t}; 

    [~, CLICKS_L] = ClickDetectTest(sig(:,1));
    [~, CLICKS_R] = ClickDetectTest(sig(:,2));

    clicks_L{t} = CLICKS_L;
    clicks_R{t} = CLICKS_R;
    %create map container of all the clicks

end

clicks_L = containers.Map(signal_names,clicks_R);
clicks_R = containers.Map(signal_names,clicks_R);

clicks_L

clicks_L('1kHz')


save('A0000B0000r028a1558.066clicks_L', 'clicks_L')
save('A0000B0000r028a1558.066clicks_R', 'clicks_R')

function [csig, clicks] = ClickDetectTest(sig)%, threshold, clickwidth)
% function [csig, clicks] = ClickDetect(sig)
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
