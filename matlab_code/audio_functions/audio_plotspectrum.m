
function audio_plotspectrum(freq, data_fft, title_string) 
    plot(freq, 20.0*log10(data_fft), 'k') 
    grid on 
    set(gca, 'XScale', 'log');
    title(title_string)
    xlabel('Frequency (Hz)')
    ylabel('Level (dB)')  
end
