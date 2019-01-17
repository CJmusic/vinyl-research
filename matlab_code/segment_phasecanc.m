% clr;


% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1129-18_KingGizzard/A33.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone45rpm/45-5.2.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1015_18_LiteToneTest/5.3.wav';
filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/1101_18_LiteTone78rpm/78-5.1.wav';

% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/misc/longleadout-trimmed.wav';
% filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/misc/longleadout-beginning.wav';

[data, fs] = audioread(filename);


fs
time = linspace(0,(length(data)-1)/fs,length(data));
rotation_speed = 78;%33.33333;%45;
T = 60/rotation_speed; %this is the length of one groove segment
num_segs = (floor(length(data)/fs/T))
n_sam = round(T*fs)
time_seg = time(1:n_sam);
seg_array = []; %need to 

for ng = 1:num_segs
    seg_array(:,:,ng) = data(1+(ng-1)*n_sam:ng*n_sam,:);
end

clf(figure(1))
clf(figure(2))
for ns = 2:num_segs;
    if rms(seg_array(:,1,ns),1) < 0.05
        figure(1)
        groove = seg_array(:,1,ns) - seg_array(:,2,ns);
        plot(time_seg, groove, 'DisplayName',['segment',num2str(ns)]);
        grid on; hold on; legend;

        figure(2)
        groove = seg_array(:,1,ns) + seg_array(:,2,ns);
        plot(time_seg,groove, 'DisplayName',['segment',num2str(ns)]);
        grid on; hold on; legend;

        % freq = fs*(0:(n_sam/2))/n_sam;
        % fft_seg = fft(seg_array(:,1,ns))/n_sam; 
        % fft_seg = fft_seg(1:size(fft_seg)/2+1);
        % plot(freq,20*log10(real(fft_seg)), 'DisplayName',['segment',num2str(ns)]);
        grid on; hold on; legend;
        % set(gca, 'XScale', 'log');

    end
end