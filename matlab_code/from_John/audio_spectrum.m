
function [data_fft, freq_fft] = audio_spectrum(data, fs, start_sam, n_sam)
    % disp('inside audio_spectrum')
    data_fft = fft(data(start_sam:start_sam+n_sam-1,:))/n_sam;% double sided
    data_fft = 2*data_fft(start_sam:floor(start_sam+n_sam/2));% single sided
    data_fft(1) = 0.5*data_fft(1);% fix DC
    if floor(n_sam/2)==n_sam/2
    data_fft(floor(n_sam/2+1)) = 0.5*data_fft(floor(n_sam/2+1));% fix Nyquist
    end
    freq_fft = fs*(0:floor(n_sam/2))/n_sam;
    % disp('finished audio_spectrum')
end %audio_spectrum

