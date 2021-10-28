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
    record1 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/A0137B0137/240a1559.867.wav')
end 

csig = record1('quiet2');
% csig = csig(0.1*96000:30.5*96000,:);
% figure(2)
% plot(csig);

Aw = audio_Aweighting(csig);
CCIRw = audio_CCIRweighting(csig);


[data_fft, freq] = audio_spectrum(csig, 96000, 1, length(csig)-1);
[data_ffta, freq] = audio_spectrum(Aw, 96000, 1, length(csig)-1);
[data_fftccir, freq] = audio_spectrum(CCIRw, 96000, 1, length(csig)-1);
data_fft = pwroctsmooth_singlesided(data_fft,0.33);
data_ffta = pwroctsmooth_singlesided(data_ffta,0.33);
data_fftccir = pwroctsmooth_singlesided(data_fftccir,0.33);

size(freq)
size(data_fft)

data_fft_L = data_fft(:,1);
data_fft_R = data_fft(:,2);

data_ffta_L = data_ffta(:,1);
data_ffta_R = data_ffta(:,2);
data_fftccir_L = data_fftccir(:,1);
data_fftccir_R = data_fftccir(:,2);


figure(1)
plot(freq, 20.0*log10(data_fft_L),'k')
hold on; grid on;
% plot(freq, 20.0*log10(data_ffta))  
% plot(freq, 20.0*log10(data_fftccir)) 
% legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  



figure(2)
plot(freq, 20.0*log10(data_fft_R),'k')
hold on; grid on;
% plot(freq, 20.0*log10(data_ffta))  
% plot(freq, 20.0*log10(data_fftccir)) 
% legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  




figure(3)
plot(freq, 20.0*log10(data_fft_L))
hold on; grid on;
plot(freq, 20.0*log10(data_ffta_L))  
plot(freq, 20.0*log10(data_fftccir_L)) 
legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

figure(4)
% plot(freq, 20.0*log10(data_fft_L))
hold on; grid on;
% plot(freq, 20.0*log10(data_ffta_L))  
plot(freq, 20.0*log10(data_fftccir_L),'k') 
% legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response CCIR-Weighting')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

figure(5)
% plot(freq, 20.0*log10(data_fft_L))
hold on; grid on;
plot(freq, 20.0*log10(data_ffta_L),'k')  
% plot(freq, 20.0*log10(data_fftccir_L)) 
% legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response A-Weighting')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  




% A_L = 20.0*log10(rms_response(Aw(1,:)));
% A_R = 20.0*log10(rms_response(Aw(2,:)));

% CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
% CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));