function [data_fft, freq_fft] = audio_spectrum(data, fs, start_sam, n_sam);
    % data = hanning(length(data)).*data;
    disp('inside audio_spectrum')
    % start_sam
    % size(data);
    data_fft = (fft(data(start_sam:start_sam+n_sam-1, :))/n_sam);
    
    data_fft = data_fft(1:floor(n_sam/2)+1,:);
    % data_fft(2:end-1,:) = 2*data_fft(2:end-1,:); %this factor of 2 is messing up A and CCIR weighting filter responses
    data_fft(2:end-1,:) = 2*data_fft(2:end-1,:); %this factor of 2 is messing up A and CCIR weighting filter responses
    
    %fix DC
    % data_fft(1,:) = 0.5*data_fft(1,:);
    % % fix Nyquist
    % if floor(n_sam/2)==n_sam/2
    %     data_fft(floor(n_sam/2+1)) = 0.5*data_fft(floor(n_sam/2+1));% fix Nyquist
    % end

    freq_fft = fs*(0:(n_sam/2))/n_sam;
    % disp('finished audio_spectrum')
end %audio_spectrum

