filename1 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00_calibration_signals/Focusrite0.6Vch2mingain.wav';
filename2 = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00_calibration_signals/dualpre0.6Vch1mingain.wav';



[data1, fs] = audioread(filename1);
[data2, fs] = audioread(filename2);


time1 = linspace(0,(length(data1)-1)/fs,length(data1));
time2 = linspace(0,(length(data2)-1)/fs,length(data2));

clf(figure(1));
figure(1);hold on;
%clf(figure(1));
plot(time1(1:1000), data1(1:1000),'blue');
%figure(2)
%clf(figure(2));
plot(time2(1:1000), data2(1:1000),'red');
legend('focusrite','dualpre');
title('0.6V with miniumum gain');


n_sam = 2^16;


freq1 = fs*(0:(n_sam/2))/n_sam;
freq2 = fs*(0:(n_sam/2))/n_sam;

fft_1 = fft(data1(1:n_sam))/n_sam;
fft_1 = fft_1(1:size(fft_1)/2+1);

fft_2 = fft(data2(1:n_sam))/n_sam;
fft_2 = fft_2(1:size(fft_2)/2+1);


clf(figure(2));
figure(2);hold on;

plot(freq1,20*log10(real(fft_1)),'blue');
plot(freq2,20*log10(real(fft_2)),'red');
grid on; hold on; legend;
set(gca, 'XScale', 'log');

legend('focusrite','dualpre');
title('0.6V with miniumum gain');


