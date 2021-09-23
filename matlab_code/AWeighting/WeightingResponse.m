
fs = 96000;
Ns = 2^16;


%-----------------------calculate logsweep-----------------------------%
f1=20;% end of turnon half-hann
f2=16000;%20000;

sweep=zeros(Ns,1);% initialize

time = [0:Ns,1]/fs;
time = time.';
time(end) = Ns/fs;
dt = time(2)-time(1);

FREQ = [];
sweep = zeros(length(time),1);
for i = (1:length(sweep))
    freq = f1 + (f2-f1)/log10(Ns)*log10(i);
    FREQ = [FREQ; freq];
    sweep(i,:) = sin(2*pi*freq*time(i));
end
% dt = 1/fs;

% sweep = sin(2*pi.*FREQ.*time);
    % time = [1:Ns];
% sweep(1:Ns)=time*sin(2*pi*f_start*Ls*(exp(Ts/Ls)-1));% logsweep
% sweep(1,:)=1;
%-----------------------------------------------------------------------%

[sweep, fs] = audioread('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/reference_files/digital_signals/sweep16kHz.wav');
sweep = sweep(1*fs:35*fs,:)/max(sweep);

data_A = audio_Aweighting(sweep);
data_CCIR = audio_CCIRweighting(sweep);

figure(12)
plot(time,FREQ)
figure(11)
plot(time,time)
figure(10)
% plot(time,sweep)
plot(sweep)

figure(1)
[data_fft, freq] = audio_spectrum(data_A, fs, 1, length(data_A));
plot(freq, 20.0*log10(abs(data_fft)),'k') 
grid on 
set(gca, 'XScale', 'log');
xlim([20,16000])
title('A Weighting Response')
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  


figure(2)
[data_fftCCIR, freq] = audio_spectrum(data_CCIR, fs, 1, length(data_CCIR));
plot(freq, 20.0*log10(abs(data_fftCCIR)),'k') 
grid on 
set(gca, 'XScale', 'log');
title('CCIR Weighting Response')
xlim([20,16000])
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  

figure(1)
[data_fftA, freq] = audio_spectrum(data_A, fs, 1, length(data_A));
plot(freq, 20.0*log10(abs(data_fftA)),'k') 
grid on 
set(gca, 'XScale', 'log');
xlim([20,16000])
title('A Weighting Response')
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  


figure(3)
plot(freq, 20.0*log10(abs(data_fftA))) 
hold on;
plot(freq, 20.0*log10(abs(data_fftCCIR))) 
grid on 
set(gca, 'XScale', 'log');
% title('CCIR Weighting Response')
title('Noise Weighting Responses')
legend('A-Weighting', 'CCIR-Weighting')
xlim([10,20000])
ylim([-96,-40])
xlabel('Frequency (Hz)')
ylabel('Level (dB)')  