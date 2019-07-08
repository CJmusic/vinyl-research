
        record.tracks('transition') = record.tracks('transition');
        % record.tracks('transition')L = record.tracks('transition')(:,1);
        % record.tracks('transition')R = record.tracks('transition')(:,1);

        clicks = audio_clickdetect(record.tracks('transition'), record.fs);
        [click_matrix, lagdiff] = audio_clickmatrix(clicks, clicks_ref);
        % record.lagdiff = lagdiff;
        record.lagdiff = audio_corrlineup(record.tracks('transition'), record.tracks('transition')_ref, record.fs);

        record.lagcorrect();
        % record.process_tracks();

        % record.tracks('transition') = record.tracks('transition');
        % record.tracks('transition')L = record.tracks('transition')(:,1);
        % record.tracks('transition')R = record.tracks('transition')(:,1);

        % figure(20);  grid on;
        % title('post lineup L')
        % plot(reference.track_times('transition'), record.tracks('transition')_refL); 
        % plot(record.track_times('transition'), record.tracks('transition')L);
        % saveas(figure(20),strcat(record.directory, '/plots/',record.filename,'postlineupL.png'))

        coh_start = 20;
        coh_end   = 80; 
        rec_cohere = record.tracks('transition')(coh_start*record.fs:coh_end*record.fs,:);
        ref_cohere = record.tracks('transition')_ref(coh_start*reference.fs:coh_end*reference.fs,:); 

        figure(1); grid on;
        n_sam = length(record.tracks('transition')_ref);
        freq_fftR = record.fs*(0:(n_sam/2-1))/n_sam;
        record.tracks('transition')_fftR = fft(record.tracks('transition')R)/n_sam;
        record.tracks('transition')_fftR = record.tracks('transition')_fftR(1:n_sam/2);

        plot(freq_fftR, 20.0*log10(record.tracks('transition')_fftR))  
        grid on;
        set(gca,'YGrid','on')
        set(gca, 'XScale', 'log');
        set(gca,'XGrid','on')
        title(strcat(record.filename, "'s Spectrum R"))
        xlabel('Frequency (Hz)')
        ylabel('Level (dB)')  
        saveas(figure(1),strcat(record.directory, '/plots/',record.filename,'spectrumR.png'))

        figure(2); grid on; 
        n_sam = length(record.tracks('transition')L);
        freq_fftL = record.fs*(0:(n_sam/2))/n_sam;
        record.tracks('transition')_fftL = fft(record.tracks('transition')L)/n_sam;
        record.tracks('transition')_fftL = record.tracks('transition')_fftL(1:n_sam/2+1);

        grid on;
        plot(freq_fftL, 20.0*log10(record.tracks('transition')_fftL))  
        set(gca, 'XScale', 'log');
        set(gca,'YGrid','on')
        set(gca,'XGrid','on')
        title(strcat(record.filename, "'s Spectrum L"))
        xlabel('Frequency (Hz)')
        ylabel('Level (dB)') 
        saveas(figure(2),strcat(record.directory, '/plots/',record.filename,'spectrumR.png'))

        [ amp_coh, freq_coh ] = audio_mscohere(record.tracks('transition'), record.tracks('transition')_ref, reference.fs);

        figure(3); grid on; 
        plot(freq_coh,amp_coh(:,1))
        set(gca, 'XScale', 'log');
        set(gca,'YGrid','on')
        set(gca,'XGrid','on')
        xlabel('Frequency (Hz)')
        % title('Coherences to Reference, Left Channel')
        title(strcat('Coherences, Left Channel ', record.filename,' to ', reference.filename ))
        saveas(figure(3),strcat(record.directory, '/plots/',record.filename,'cohL.png'))


        figure(4); grid on;
        plot(freq_coh,amp_coh(:,2))
        set(gca,'YGrid','on')
        set(gca,'XGrid','on')
        set(gca, 'XScale', 'log');
        xlabel('Frequency (Hz)')
        title(strcat('Coherences, Right Channel '))
        saveas(figure(4),strcat(record.directory, '/plots/',record.filename,'cohR.png'))

        figure(30); grid on; hold on;
        plot(freq_coh,amp_coh(:,1))
        set(gca, 'XScale', 'log');
        set(gca,'YGrid','on')
        set(gca,'XGrid','on')
        xlabel('Frequency (Hz)')
        % title('Coherences to Reference, Left Channel')
        title(strcat('Coherences, Left Channel '))


        figure(40); grid on; hold on;
        plot(freq_coh,amp_coh(:,2))
        set(gca,'YGrid','on')
        set(gca,'XGrid','on')
        set(gca, 'XScale', 'log');
        xlabel('Frequency (Hz)')
        title(strcat('Coherences, Right Channel ', record.filename,' to ', reference.filename ))
        