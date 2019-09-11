
disp('-----------RecordProcess.m------------')
addpath('audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/');
addpath('/Volumes/AUDIOBANK/audio_files/')
addpath('/Volumes/AUDIOBANK/audio_files/pressings/')

reference = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r28a.wav');
record = audio_recordclass('/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_bin/A0000B0000/03141_A0000B0000r30a.wav');

record.process_tracks();
reference.process_tracks();

% rec_cohere = rec_cohere;
% rec_cohereL = rec_cohere(:,1);
% rec_cohereR = rec_cohere(:,1);
rec_cohere = record.tracks('leadout');
ref_cohere = reference.tracks('leadout'); 


coh_start = 1;
coh_end   = 5; 
rec_cohere = rec_cohere(coh_start*record.fs:coh_end*record.fs,:);
ref_cohere = ref_cohere(coh_start*reference.fs:coh_end*reference.fs,:); 

clicks = audio_clickdetect(rec_cohere, record.fs);
clicks_ref = audio_clickdetect(ref_cohere, reference.fs);
[click_matrix, lagdiff] = audio_clickmatrix(clicks, clicks_ref);
% record.lagdiff = lagdiff;
record.lagdiff = audio_corrlineup(rec_cohere, ref_cohere, record.fs);

record.lagcorrect();
% record.process_tracks();

% rec_cohere = rec_cohere;
rec_cohereL = rec_cohere(:,1);
rec_cohereR = rec_cohere(:,2);

ref_cohereL = ref_cohere(:,1);
ref_cohereR = ref_cohere(:,2);

figure(20); hold on;
title('post lineup L')
plot(ref_cohereL); 
plot(rec_cohereL);
grid on;
%saveas(figure(20),strcat(record.directory, '/plots/',record.filename,'postlineupL.png'))


n_sam = length(rec_cohere);
freq_fftR = record.fs*(0:(n_sam/2-1))/n_sam;
rec_fft = fft(rec_cohere)/n_sam;
rec_fft = rec_fft(1:n_sam/2);

figure(1); grid on;
plot(freq_fftR, 20.0*log10(rec_fft))  
grid on;
set(gca,'YGrid','on')
set(gca, 'XScale', 'log');
set(gca,'XGrid','on')
title(strcat(record.filename, "'s Spectrum R"))
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
%saveas(figure(1),strcat(record.directory, '/plots/',record.filename,'spectrumR.png'))

figure(2); grid on; 
n_sam = length(ref_cohere);
freq_fftL = record.fs*(0:(n_sam/2))/n_sam;
ref_fftL = fft(ref_cohereL)/n_sam;
ref_fftL = ref_fftL(1:n_sam/2+1);

grid on;
plot(freq_fftL, ref_fftL)  
set(gca, 'XScale', 'log');
set(gca,'YGrid','on')
set(gca,'XGrid','on')
title(strcat(reference.filename, "'s Spectrum L"))
xlabel('Frequency (Hz)')
ylabel('Level (dB)') 
%saveas(figure(2),strcat(record.directory, '/plots/',record.filename,'spectrumR.png'))

[ amp_coh, freq_coh ] = audio_mscohere(rec_cohere, ref_cohere, reference.fs);

figure(3); grid on; 
plot(freq_coh, amp_coh(:,1))
set(gca, 'XScale', 'log');
set(gca,'YGrid','on')
set(gca,'XGrid','on')
xlabel('Frequency (Hz)')
% title('Coherences to Reference, Left Channel')
title(strcat('Coherences, Left Channel ', record.filename,' to ', reference.filename ))
%saveas(figure(3),strcat(record.directory, '/plots/',record.filename,'cohL.png'))


figure(4); grid on;
plot(freq_coh,amp_coh(:,2))
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
title(strcat('Coherences, Right Channel '))
%saveas(figure(4),strcat(record.directory, '/plots/',record.filename,'cohR.png'))

figure(30); grid on; hold on;
plot(freq_coh,amp_coh(:,1))
set(gca, 'XScale', 'log');
set(gca,'YGrid','on')
set(gca,'XGrid','on')
xlabel('Frequency (Hz)')
% title('Coherences to Reference, Left Channel')
title(strcat('Coherences, Left Channel '))


figure(40); grid on; hold on;
plot(freq_coh,amp_coh(:,2))
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
title(strcat('Coherences, Right Channel ', record.filename,' to ', reference.filename ))
