%This function cleans up mscohere's inputs and outputs to produce coherence per frequency 
%arrays.  
%

function [amp_coh, freq_coh] = audio_mscohere(data1, data2, fs)
    nfft=2^14;
    %freq_coh=([0:nfft/2])*fs/nfft;
    [amp_coh, freq_coh] = mscohere(data1,data2,hanning(nfft),[],nfft/2,fs);
end

