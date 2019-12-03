

% need to test various values of threshold and width and produce a plot of 
% number of clicks vs. those two values and look for a plateau

clc; close all;

addpath('D:\Code\vinyl-research\matlab_code\')
file = ('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\040318_A0000B0000r001b.wav')

%~~~~ MAC FILES ~~~~%

% addpath('/Users/cz/Code/vinyl-research/matlab_code')
% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/click_testing/')
% addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/click_testing/declicked')
% file = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031418_A0000B0000r27a.wav')

%~~~~ MAC FILES END ~~~~%

[sig, fs] = audioread(file);

sig = sig(418*fs:515*fs,:);
time = (0:length(sig))/fs;


%~~~~ TESTING DEPENDENCY ~~~~~%

% clickArr = [];
% thresholdArr = [];

% i = 1;
% for clickwidth = (20:5:40);
%     i
%     clickwidth
%     for threshold = (200:100:900);
%         threshold
%         thresholdArr = [thresholdArr, threshold];
%         [~, clicks] = ClickDetect(sig, threshold, clickwidth);
%         clickArr = [clickArr, length(clicks)];
%     end
%     figure(1); grid on; hold on;
%     plot(thresholdArr, clickArr)
%     grid on; hold on;
%     set(gca,'XTick',200:50:900);
%     title('num clicks vs threshold and width')
%     xlabel('threshold')
%     ylabel('num of clicks')
%     grid on; 
%     clickArr = [];
%     thresholdArr = [];
%     i = i + 1;
        
% end
   
% legend('20','25','30','35','40')

%~~~~~ TESTING DEPENDENCY ENDS ~~~~~~%

clickArr = [];
widthArr = [];
thresholdArr = [];

% ~~~~ testing each individually ~~~~ %



for threshold = (200:400:900);
    threshold
    thresholdArr = [thresholdArr, threshold];
    [~, clicks] = ClickDetect(sig, threshold, 30);
    clickArr = [clickArr, length(clicks)];

end
   

figure(1); grid on; 
plot(thresholdArr, clickArr)
set(gca,'XTick',1:150:1500);
title('Threshold varying between 200-900, width = 30')
xlabel('threshold')
ylabel('num of clicks')
grid on; 

clickArr = [];
widthArr = [];
for clickwidth = (1:100:1000);
    clickwidth
    [~, clicks] = ClickDetect(sig, 550, clickwidth);
    clickArr = [clickArr, length(clicks)];
    widthArr = [widthArr, clickwidth];


end

figure(2); 
plot(widthArr, clickArr)
set(gca,'XTick',20:1:40);
title('clickwidth varying between 20-40, threshold = 550')
xlabel('clickwidth')
ylabel('num of clicks')
grid on; 



%%~~~~ Looking at both parameters together in a 3d plot ~~~~~%

%% BROKEN %%

% thresholdArr = (200:100:900);
% widthArr = (20:2.5:37.5);


% for i = (1:length(thresholdArr))
%     threshold = thresholdArr(i)
%     for j = (1:length(widthArr))
%         clickwidth = widthArr(j)
%         [~, clicks] = ClickDetect(sig, threshold, clickwidth);
%         length(clicks)
%         clickArr(i,j) = length(clicks);
%     end
% end

% thresholdArr
% widthArr
% clickArr

% length()



% figure(3)
% surf(widthArr, thresholdArr, clickArr)
% grid on
% title('Surf plot')
% xlabel('clickwidth')
% ylabel('threshold')



% figure(4)
% contour(widthArr, thresholdArr, clickArr)
% grid on
% title('Contour plot')
% xlabel('clickwidth')
% ylabel('threshold')


% figure(5)
% scatter(widthArr,thresholdArr,clickArr)
% grid on
% title('Scatter plot')
% xlabel('clickwidth')
% ylabel('threshold')

