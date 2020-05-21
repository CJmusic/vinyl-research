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

    record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');
    % record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/040318_A0000B0000r001a.wav');
    record1 = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/styluswear/042820_A0000B0000r1a.wav');

end 
if ispc() == true
    addpath('D:\Code\vinyl-research\matlab_code\A0137B0137')
    addpath('D:\Code\vinyl-research\matlab_code\')
    addpath('D:\Code\vinyl-research\matlab_code\audio_functions\')
    record1 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxsteam1a.wav');
    record2 = SeperateTracks('d:/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/testing/maxbarrelzones3a.wav');    
end



csig = record1('sweep');
% csig = csig(0.1*96000:30.5*96000,:);
% figure(2)
% plot(csig);

Aw = audio_Aweighting(csig(:,1));
CCIRw = audio_CCIRweighting(csig(:,1));


[data_fft, freq] = audio_spectrum(csig, 96000, 1, length(csig)-1);
% [data_ffta, freq] = audio_spectrum(Aw, 96000, 1, length(csig)-1);
% [data_fftccir, freq] = audio_spectrum(CCIRw, 96000, 1, length(csig)-1);
data_fft = pwroctsmooth(data_fft,0.33);
% data_ffta = pwroctsmooth(data_ffta,0.33);
% data_fftccir = pwroctsmooth(data_fftccir,0.33);

size(freq)
size(data_fft)

figure(1)
plot(freq, 20.0*log10(data_fft),'k')
hold on; grid on;
% plot(freq, 20.0*log10(data_ffta))  
% plot(freq, 20.0*log10(data_fftccir)) 
% legend('non-weighted','A-weighted','CCIR-weighted')
title('Record Frequency Response')
set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  


% A_L = 20.0*log10(rms_response(Aw(1,:)));
% A_R = 20.0*log10(rms_response(Aw(2,:)));

% CCIR_L = 20.0*log10(avg_response(CCIRw(1,:)));
% CCIR_R = 20.0*log10(avg_response(CCIRw(2,:)));