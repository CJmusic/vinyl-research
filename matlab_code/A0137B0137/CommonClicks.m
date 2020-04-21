clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')


% record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/mincool5b.wav');

figure(1)
plot(record1('transition'))
hold on; grid on;
plot(record2('transition'))

[~, clicks1] = ClickDetect(record1('transition'),200,20);
[~, clicks2] = ClickDetect(record2('transition'),200,20);
num_comclicks(clicks1, clicks2)


function num_comclicks = num_comclicks(clicks, clicks_ref);
    %loop through clicks 
    %look for other clicks in
    disp('INSIDE COMMON CLICKS')    
    % clicks;
    % clicks_ref;
    
    size(clicks)
    size(clicks_ref)


    
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

    % num_comclicks = 0;
    % relaxation = 240;
    % comclicks = [];
    % for i = (1:length(clicks_ref))
    %     index = find(clicks < clicks_ref(i) + relaxation & clicks > clicks_ref(i) - relaxation)
    %     if index > 0 
    %         comclicks = [comclicks, clicks_ref(i)];
    %     end
    % end


    % below seems to work
    clicks
    clicks_ref
    relaxation = 96000;
    comclicks = [];
    for i = (1:length(clicks_ref))
        % i
        upper_bound = clicks < clicks_ref(i) + relaxation;
        
        size(upper_bound)
        lower_bound  = clicks(upper_bound) > clicks_ref(i) - relaxation;
        % size(lower_bound)
        size(lower_bound)
        if sum(lower_bound) > 0;
            disp('COMMON CLICK')
            comclicks = [comclicks, clicks_ref(i)];
            continue
        end

        disp('Upper Bound')
        upper_bound
        size(upper_bound)
        disp('Lower Bound')
        lower_bound
        size(lower_bound)

        size(comclicks)
    end

        
        %%% below seems to work, but is very slow
    % comclicks = [];
    % relaxation = 9600;
    % for i = (1:length(clicks_ref))
    %     for k = (1:length(clicks))
    %         if clicks(k) >= clicks_ref(i) - relaxation & clicks(k) < clicks_ref(i) + relaxation;
    %             comclicks = [comclicks, clicks_ref(i)];
    %             disp('COMMON CLICK')
    %         end
    %     end
    % end


        % i
        % clicks_ref(i)
        % clicks(i)

        % clicks_ref(i) - relaxation
        % clicks_ref(i) + relaxation
        % clicks_ref(i) - relaxation
        % clicks_ref(i) + relaxation
        % try
        %     [ind,~] = find(clicks >= clicks_ref(i,1) - relaxation & clicks < clicks_ref(i,1) + relaxation)
        %     [ind,~] = find(clicks >= clicks_ref(i,2) - relaxation & clicks < clicks_ref(i,2) + relaxation)
        %     if ind != 
        % catch
        %     continue
        % end
        % if find(clicks_ref(i) - relaxation < clicks < clicks_ref(i) + relaxation) ~= 0

        % num_comclicks = num_comclicks + 1;
        % end 
        % num_comclicks
    % end

    num_comclicks = length(comclicks)

    for xi = 1:length(clicks)
        x1 = (clicks(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'red','LineStyle', '--');
    end
    for xi = 1:length(clicks_ref)
        x1 = (clicks_ref(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'blue','LineStyle', '--');
    end

    for xi = 1:length(comclicks)
        x1 = (comclicks(xi));
        figure(1); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
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
