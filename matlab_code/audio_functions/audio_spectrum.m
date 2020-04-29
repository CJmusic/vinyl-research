
function [data_fft, freq_fft] = audio_spectrum(data, fs, start_sam, n_sam);
    disp('inside audio_spectrum')
    freq_fft = fs*(0:(n_sam/2-2))/n_sam;
    data_fft = fft(data(start_sam:start_sam+n_sam, :))/n_sam;
    data_fft = data_fft(1:size(data_fft)/2-1);
    disp('finished audio_spectrum')
end %audio_spectrum

