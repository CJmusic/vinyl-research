

% need to test various values of threshold and width and produce a plot of 
% number of clicks vs. those two values and look for a plateau


clc; close all;
addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\')
addpath('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\declicked')

files = dir('D:\OneDrive - University of Waterloo\Vinyl_Project\audio_bin\click_testing\*.wav')

for i = (1:length(files))
    file = strcat(files(i).folder,'/',files(i).name)
    de_file = strcat(files(i).folder,'/declicked/',files(i).name) 
    csig = [];
    [sig, fs] = audioread(file);
    [desig, fs] = audioread(de_file);

    sig = sig(500*fs:550*fs,:);
    desig = desig(500*fs:550*fs,:);


    time = (0:length(sig)-1)/fs;
    time2 = (0:length(sig)-1)/fs;

    [csig(:,1), clicksL] = ClickDetectTest(sig(:,1));
    [csig(:,2), clicksR] = ClickDetectTest(sig(:,2));


    %~~~~~~~~~~~~~~~~PLOTTING~~~~~~~~~~~~~~~~~~~

    figure(i); grid on; hold on 
    subplot(3,1,1)
    plot(time, sig)
    title('original')
    subplot(3,1,2)
    plot(time, csig)
    title('my click removal')
    subplot(3,1,3)
    plot(time, desig)
    title('audacitys click removal')
    for xi = 1:length(clicksL)
        x1 = time(clicksL(xi));
        figure(100 + i); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end
    for xi = 1:length(clicksR)
        x1 = time(clicksR(xi));
        figure(100 + i); hold on;
        line([x1 x1], get(gca, 'ylim'),'Color', 'black','LineStyle', '--');
    end

    figure(100 + i); grid on; 
    plot(time, desig-csig)
    title('audacity - mine')
end


