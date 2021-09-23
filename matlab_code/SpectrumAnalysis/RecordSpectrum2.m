clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

if ismac() == true
    addpath('/Users/cz/Code/vinyl-research/matlab_code')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')    
    addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
    addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/')
    addpath('/Users/cz/Code/vinyl-research/matlab_code/Common/')

    % record1 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/300b1600.957.wav');
    record1a = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/040318_A0000B0000r032a1557.912.wav');
    % record2a = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/240a1559.867.wav');
    record2a = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/151a1559.342.wav');
    




    record1b = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0000B0000/040318_A0000B0000r032b1555.663.wav')
    record2b = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/004b1600.166.wav')
end 

csig1a = record1a('quiet');
csig2a = record2a('quiet');

csig1b = record1b('quiet');
csig2b = record2b('quiet');
% csig = csig(0.1*96000:30.5*96000,:);
% figure(2)
% plot(csig);
[freq1a, data_fft_L1a, data_fft_R1a, data_ffta_L1a, data_ffta_R1a, data_fftccir_L1a, data_fftccir_R1a] = weighted_spectrum(csig1a);

[freq, data_fft_L1b, data_fft_R1b, data_ffta_L1b, data_ffta_R1b, data_fftccir_L1b, data_fftccir_R1b] = weighted_spectrum(csig1b);

[freq2a, data_fft_L2a, data_fft_R2a, data_ffta_L2a, data_ffta_R2a, data_fftccir_L2a, data_fftccir_R2a] = weighted_spectrum(csig2a);

[freq, data_fft_L2b, data_fft_R2b, data_ffta_L2b, data_ffta_R2b, data_fftccir_L2b, data_fftccir_R2b] = weighted_spectrum(csig2b);



figure(1)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_L2a))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_R2a))
legend('CAF', 'Reinee')
title('Right Channel ')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])

figure(2)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_L2b))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_R2b))
legend('CAF', 'Reinee')
title('Right Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])


figure(3)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_L2a))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_R2a))
legend('CAF', 'Reinee')
title('Right Channel ')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])

figure(4)
subplot(2,1,1)
plot(freq, 20.0*log10(data_ffta_L1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_L2b))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_ffta_R1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_R2b))
legend('CAF', 'Reinee')
title('Right Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])



figure(5)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_L2a))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1a))
hold on; grid on; 
plot(freq, 20.0*log10(data_fft_R2a))
legend('CAF', 'Reinee')
title('Right Channel ')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])

figure(6)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fftccir_L1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_fftccir_L2b))
legend('CAF', 'Reinee')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fftccir_R1b))
hold on; grid on; 
plot(freq, 20.0*log10(data_fftccir_R2b))
legend('CAF', 'Reinee')
title('Right Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
% ylim([-96,0])
xlim([1, 20000])


figure(7)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1a),'Color',[0/255 64/255 255/255])
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_L1a),'Color',[0/255 128/255 255/255])
plot(freq, 20.0*log10(data_fftccir_L1a),'Color',[0/255 191/255 255/255])
plot(freq, 20.0*log10(data_fft_L2a),'Color',[255/255 64/255 0/255])
plot(freq, 20.0*log10(data_ffta_L2a),'Color',[255/255 128/255 0/255])
plot(freq, 20.0*log10(data_fftccir_L2a),'Color',[255/255 191/255 0/255])
legend('CAF', 'CAF A','CAF CCIR','Reinee', 'Reinee A', 'Reinee CCIR')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
ylim([-120,-30])
xlim([1, 2000000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1a),'Color',[0/255 64/255 255/255])
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_R1a),'Color',[0/255 128/255 255/255])
plot(freq, 20.0*log10(data_fftccir_R1a),'Color',[0/255 191/255 255/255])
plot(freq, 20.0*log10(data_fft_R2a),'Color',[255/255 64/255 0/255])
plot(freq, 20.0*log10(data_ffta_R2a),'Color',[255/255 128/255 0/255])
plot(freq, 20.0*log10(data_fftccir_R2a),'Color',[255/255 191/255 0/255])
legend('CAF', 'CAF A','CAF CCIR','Reinee', 'Reinee A', 'Reinee CCIR')
title('Right Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
ylim([-120,-30])
xlim([1, 2000000])


figure(8)
subplot(2,1,1)
plot(freq, 20.0*log10(data_fft_L1b),'Color',[0/255 64/255 255/255])
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_L1b),'Color',[0/255 128/255 255/255])
plot(freq, 20.0*log10(data_fftccir_L1b),'Color',[0/255 191/255 255/255])
plot(freq, 20.0*log10(data_fft_L2b),'Color',[255/255 64/255 0/255])
plot(freq, 20.0*log10(data_ffta_L2b),'Color',[255/255 128/255 0/255])
plot(freq, 20.0*log10(data_fftccir_L2b),'Color',[255/255 191/255 0/255])
legend('CAF', 'CAF A','CAF CCIR','Reinee', 'Reinee A', 'Reinee CCIR')
title('Left Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
ylim([-120,-30])
xlim([1, 2000000])
subplot(2,1,2)
plot(freq, 20.0*log10(data_fft_R1b),'Color',[0/255 64/255 255/255])
hold on; grid on; 
plot(freq, 20.0*log10(data_ffta_R1b),'Color',[0/255 128/255 255/255])
plot(freq, 20.0*log10(data_fftccir_R1b),'Color',[0/255 191/255 255/255])
plot(freq, 20.0*log10(data_fft_R2b),'Color',[255/255 64/255 0/255])
plot(freq, 20.0*log10(data_ffta_R2b),'Color',[255/255 128/255 0/255])
plot(freq, 20.0*log10(data_fftccir_R2b),'Color',[255/255 191/255 0/255])
legend('CAF', 'CAF A','CAF CCIR','Reinee', 'Reinee A', 'Reinee CCIR')
title('Right Channel')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  
ylim([-120,-30])
xlim([1, 2000000])



saveas(figure(1), 'unweighteda.png')
saveas(figure(2), 'aweighteda.png')
saveas(figure(3), 'ccirweighteda.png')
saveas(figure(4), 'unweightedb.png')
saveas(figure(5), 'aweightedb.png')
saveas(figure(6), 'ccirweightedb.png')
saveas(figure(7), 'bothpressingsspectra.png')


% size(freq)
% size(data_fft)

% figure(1)
% plot(freq, 20.0*log10(data_fft_L),'k')
% hold on; grid on;
% % plot(freq, 20.0*log10(data_ffta))  
% % plot(freq, 20.0*log10(data_fftccir)) 
% % legend('non-weighted','A-weighted','CCIR-weighted')
% title('Record Frequency Response')
% set(gca, 'XScale', 'log');
% xlabel('Frequency (Hz)')
% ylabel('Level (dB)')  



% figure(2)
% plot(freq, 20.0*log10(data_fft_R),'k')
% hold on; grid on;
% % plot(freq, 20.0*log10(data_ffta))  
% % plot(freq, 20.0*log10(data_fftccir)) 
% % legend('non-weighted','A-weighted','CCIR-weighted')
% title('Record Frequency Response')
% set(gca, 'XScale', 'log');
% xlabel('Frequency (Hz)')
% ylabel('Level (dB)')  




% figure(3)
% plot(freq, 20.0*log10(data_fft_L))
% hold on; grid on;
% plot(freq, 20.0*log10(data_ffta_L))  
% plot(freq, 20.0*log10(data_fftccir_L)) 
% legend('non-weighted','A-weighted','CCIR-weighted')
% title('Record Frequency Response')
% set(gca, 'XScale', 'log');
% xlabel('Frequency (Hz)')
% ylabel('Level (dB)')  

% figure(4)
% % plot(freq, 20.0*log10(data_fft_L))
% hold on; grid on;
% % plot(freq, 20.0*log10(data_ffta_L))  
% plot(freq, 20.0*log10(data_fftccir_L),'k') 
% % legend('non-weighted','A-weighted','CCIR-weighted')
% title('Record Frequency Response CCIR-Weighting')
% set(gca, 'XScale', 'log');
% xlabel('Frequency (Hz)')
% ylabel('Level (dB)')  

% figure(5)
% % plot(freq, 20.0*log10(data_fft_L))
% hold on; grid on;
% plot(freq, 20.0*log10(data_ffta_L),'k')  
% % plot(freq, 20.0*log10(data_fftccir_L)) 
% % legend('non-weighted','A-weighted','CCIR-weighted')
% title('Record Frequency Response A-Weighting')
% set(gca, 'XScale', 'log');
% xlabel('Frequency (Hz)')
% ylabel('Level (dB)')  




% A_L = 20.0*log10(rms_response(Aw(1,:)));
% A_R = 20.0*log10(rms_response(Aw(2,:)));

% CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
% CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));

function [freq, data_fft_L, data_fft_R, data_ffta_L, data_ffta_R, data_fftccir_L, data_fftccir_R] = weighted_spectrum(csig);


    Aw = audio_Aweighting(csig);
    CCIRw = audio_CCIRweighting(csig);

    % [data_fft, freq] = audio_spectrum(csig, 96000, 1, length(csig)-1);
    % [data_ffta, freq] = audio_spectrum(Aw, 96000, 1, length(csig)-1);
    % [data_fftccir, freq] = audio_spectrum(CCIRw, 96000, 1, length(csig)-1);

    [data_fft, freq] = audio_spectrum(csig, 96000, 1, 2^16);
    [data_ffta, freq] = audio_spectrum(Aw, 96000, 1, 2^16);
    [data_fftccir, freq] = audio_spectrum(CCIRw, 96000, 1, 2^16);


    data_fft_L = data_fft(:,1);
    data_fft_R = data_fft(:,2);
    data_ffta_L = data_ffta(:,1);
    data_ffta_R = data_ffta(:,2);
    data_fftccir_L = data_fftccir(:,1);
    data_fftccir_R = data_fftccir(:,2);


    data_fft_L = pwroctsmooth(data_fft_L,0.33);
    data_ffta_L = pwroctsmooth(data_ffta_L,0.33);
    data_fftccir_L = pwroctsmooth(data_fftccir_L,0.33);
    data_fft_R = pwroctsmooth(data_fft_R,0.33);
    data_ffta_R = pwroctsmooth(data_ffta_R,0.33);
    data_fftccir_R = pwroctsmooth(data_fftccir_R,0.33);


end