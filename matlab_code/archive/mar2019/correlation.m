
% clr;


% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1129-18_KingGizzard/A33.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone45rpm/45-5.2.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.3.wav';
filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone78rpm/78-5.1.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone78rpm/1.4.wav';

%filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/misc/longleadout-trimmed.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/misc/longleadout-beginning.wav';

%filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/lacquernoise.wav';

filename = 'ChirpSilentLeadinCUT.wav'

[data, fs] = audioread(filename);


fs
time = linspace(0,(length(data)-1)/fs,length(data));
rotation_speed = 33.33333;%45;
T = 60/rotation_speed; %this is the length of one groove segment
num_segs = (floor(length(data)/fs/T))
n_sam = round(T*fs)
time_seg = time(1:n_sam);
seg_array = []; %need to 

num_segs = 20;

for ng = 1:num_segs
    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
end
size(seg_array)
corr_L = [];
corr_R = [];

%%This loop for correlation with the first groove in the series
for ns = 1:num_segs
    %corr_coeff = sum(seg_array(:,:,1).*seg_array(:,:,ns)) - sum(seg_array(:,:,1)).*sum(seg_array(:,:,ns));
    corr_coeff = sum(seg_array(:,:,20).*seg_array(:,:,ns)); 
    corr_L = [corr_L,corr_coeff(1)];
    corr_R = [corr_R,corr_coeff(2)];
end

clf(figure(1));
figure(1);
hold on;
plot(corr_R);
plot(corr_L);
title('Correlation with first groove')
xlabel('Groove number')
ylabel('Correlation Coeff')

corr_L = [];
corr_R = [];
%This loop for correlation with the next groove in the series
for ns = 1:num_segs-1
    %corr_coeff = sum(seg_array(:,:,ns).*seg_array(:,:,ns+1)) - sum(seg_array(:,:,ns)).*sum(seg_array(:,:,ns+1));
    corr_coeff = sum(seg_array(:,:,ns).*seg_array(:,:,ns+1));
    ns
    corr_L = [corr_L,corr_coeff(1)];
    corr_R = [corr_R,corr_coeff(2)];
%sum(rev4ng(:,1,1).*rev4ng(:,1,ng),1)
end
clf(figure(2));
figure(2);
hold on;
plot(corr_R);
plot(corr_L);
title('Correlation with next groove')
xlabel('Groove number')
ylabel('Correlation Coeff')
