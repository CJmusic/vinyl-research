% clear all;close all;clc
% set(0,'DefaultLineLineWidth',1.5);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultAxesFontWeight','bold')
% set(0,'DefaultAxesLineWidth',1.5)


% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')


% record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/mincool5b.wav');

% figure(1)
% plot(record1('transition'))
% hold on; grid on;
% plot(record2('transition'))

% [~, clicks1] = ClickDetect(record1('transition'),200,20);
% [~, clicks2] = ClickDetect(record2('transition'),200,20);
% num_comclicks(clicks1, clicks2)


function num_comclicks = num_comclicks(clicks, clicks_ref);
    % PLOT CLICKS
    % for xi = 1:length(clicks)
    %     x1 = (clicks(xi));
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end

    % % saveas(figure(1),strcat(files(i).folder,i,'.png'))

    % for xi = 1:length(clicks_ref)
    %     x1 = (clicks_ref(xi));
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '--');
    % end



    % below seems to work
    relaxation = 96000;
    comclicks = [];
    for i = (1:length(clicks_ref))
        upper_bound = clicks < clicks_ref(i) + relaxation;        
        lower_bound  = clicks(upper_bound) > clicks_ref(i) - relaxation;
        if sum(lower_bound) > 0;
            % disp('COMMON CLICK')
            comclicks = [comclicks, clicks_ref(i)];
            continue
        end
    end




    % num_comclicks = length(comclicks)

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
    %     x1 = (comclicks(xi));
    %     figure(1); hold on;
    %     line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    % end

   
end
