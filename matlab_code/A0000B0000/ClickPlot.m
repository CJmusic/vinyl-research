clc; clf(figure(1)); clf(figure(2));

[sig, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r28a.wav');

length(sig)
length(sig)/fs
disp('running click detect')
[csig(:,1), CLICKS_L] = ClickDetect(sig(:,1),200,20);
% [csig(:,2), CLICKS_R] = ClickDetect(sig(:,2),200,20);
size(CLICKS_L)



disp('calculating r')
r = 960/1.8 - floor(CLICKS_L/(fs*1.8)) + (length(sig)/3)/(fs*1.8);
disp('calculating theta')
theta = 2*pi*(CLICKS_L - r*1.8*fs)/(1.8*fs);

figure(1); 
polarplot(theta, r, 'o')
label = (length(sig)/3)/(fs*1.8);
hold on;
polarplot(linspace(0,2*pi,100),linspace(label,label,100),'k:')

% figure(2)
% r = floor(CLICKS_L/(fs*1.8)) + (length(sig)/3)/(fs*1.8);
% theta = 2*pi*(CLICKS_L - r*1.8*fs)/(1.8*fs);

% polarscatter(theta, r)
% hold on;
% polarplot(linspace(0,2*pi,100),linspace(label,label,100),'k:')
