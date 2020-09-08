clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)


addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')


% record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r028a.wav'); 

record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3b.wav');  
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
% record2 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/-01a.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/mincool5b.wav');
% record2 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/minpucksize1a.wav');

% figure(1)
% plot(record1('100Hz'))
% hold on; grid on;
% plot(record2('100Hz'))

[~, clicks1] = ClickDetect(record1('1kHz'))
[~, clicks2] = ClickDetect(record2('1kHz'))

figure(1)
plot(record1('1kHz'))
figure(2)
plot(record2('1kHz'))

disp(strcat('common clicks record 1... ', num2str(length(clicks1))))
disp(strcat('common clicks record 2... ', num2str(length(clicks2))))


RELAX = [0, 250, 500, 750, 1000, 1250, 1500];
COM = [];
for i = 1:length(RELAX); 
    common = num_comclicks(clicks1, clicks2, RELAX(i))
    COM = [COM, common];
end

figure(100)
plot(RELAX, COM)
grid on; 
title('Common clicks vs relaxation')
xlabel('relaxation')
ylabel('# of common clicks')

disp(strcat('clicks 1.....', num2str((length(clicks1)))))
disp(strcat('clicks 2.....', num2str((length(clicks2)))))
% disp(strcat('common .....' , num2str



function num_comclicks = num_comclicks(clicks, clicks_ref, relaxtion);
    % relaxation = 750;
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
