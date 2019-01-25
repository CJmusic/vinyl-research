%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%This file takes a WAV recording of a log sweep and calculates the frequeuncy response through that. 
%Christopher Zaworski Jan 25, 2019
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%



filename = '/Users/cz/OneDrive - University of Waterloo/Vinyl_Project/audio_files/00_digital_files/sweep16kHz.wav'

[data, fs] = audioread(filename); 

n_sam = length(data) %number of samples 

freq = fs*(0:n_sam/2)/n_sam; 
fft_data = fft(data)/n_sam;
fft_data = fft_data(1:size(fft_data)/2+1);

size(freq)
size(fft_data)

plot(freq,20*log10(real(fft_data))) 
grid on; 
