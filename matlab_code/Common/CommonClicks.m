% clear all;close all;clc
% set(0,'DefaultLineLineWidth',1.5);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)


% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')


% % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
% record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 

% % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3b.wav');  
% % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
% record2 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-01a.wav');
% % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/mincool5b.wav');
% % record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/minpucksize1a.wav');

% figure(1)
% plot(record1('100Hz'))
% hold on; grid on;
% plot(record2('100Hz'))

% [~, clicks1] = ClickDetect(record1('100Hz'));
% [~, clicks2] = ClickDetect(record2('100Hz'));
% common = num_comclicks(clicks1, clicks2);
  
% disp(strcat('clicks 1.....', num2str((length(clicks1)))))
% disp(strcat('clicks 2.....', num2str((length(clicks2)))))
% disp(strcat('common .....' , num2str(common)))


function num_comclicks = num_comclicks(clicks, clicks_ref);
    relaxation = 750;
    comclicks = [];
    for i = (1:length(clicks_ref))
        upper_bound = clicks < clicks_ref(i) + relaxation;
        lower_bound  = clicks(upper_bound) > clicks_ref(i) - relaxation;
        if sum(lower_bound) > 0;
            comclicks = [comclicks, clicks_ref(i)];
            continue
        end
    end

    num_comclicks = length(comclicks);


    % for xi = 1:length(clicks)
    %     x1 = (clicks(xi));
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '--');
    % end
    % for xi = 1:length(clicks_ref)
    %     x1 = (clicks_ref(xi));
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'blue','LineStyle', '--');
    % end

    % for xi = 1:length(comclicks)
    %     x1 = (comclicks(xi))
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black');
    % end

end
