

function stereo_bleed = StereoBleed(sig, chan)
    %%% fft based 
    L = 2^16;
    fs = 96000;

    seg = sig(floor(length(sig)/2) - L/2:floor(length(sig)/2) + L/2 - 1,:);

    [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
    seg(:,1) = filter(b,a,seg(:,1));
    seg(:,2) = filter(b,a,seg(:,2));

    win = flattopwin(L);
    seg = seg.*win;

    fftsigL = fft(seg(:,1))/L;
    fftsigL = fftsigL(1:L/2+1);

    fftsigR = fft(seg(:,2))/L;
    fftsigR = fftsigR(1:L/2+1);

    fftfreq = fs*(0:(L/2))/L;


    % figure(1+t); hold on; grid on;
    % plot(fftfreq, 20*log10(fftsigL))
    % plot(fftfreq, 20*log10(fftsigR))
    % set(gca, 'XScale', 'log')

    peakL = max(real(fftsigL));
    peakR = max(real(fftsigR));

    
    if chan == 1;
        stereo_bleed = 20*log10(peakR/peakL);
    elseif chan == 2;
        stereo_bleed = 20*log10(peakL/peakR);
    end
    % title(stereo_bleed)

    % disp('divide')
    % 20*log10(peakL/peakR)
    % disp('subtract')
    % 20*log10(peakL) - 20*log10(peakR)
end