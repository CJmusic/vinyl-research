% Record Wow Measurement
% using 3150 track
% John Vanderkooy
% April 2019
%  


function wow = WowFlutter(rev4_1);
    % disp('--------------WOW AND FLUTTER--------------')
    filename = 'wow';
    fs = 96000;
    ts = length(rev4_1)*(1/4)/fs;
    tf = length(rev4_1)*(3/4)/fs;
    % [rev4_1,fs]=audioread(filename);
    %lr=1; %1=left, 2=right
    % disp(['lr: ' num2str(lr)])
    Nt=length(rev4_1);
    % disp(['fs: ' num2str(fs) '  N_total: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    t=linspace(0,(Nt-1)/fs,Nt)';%column vector
    %-----------------prefiltering-----------
    [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
    %rev4_1=filter(b,a,rev4_1);
    %-------------------plot all data------------
    % figure(20)
    % plot(t,rev4_1(1:Nt),'b') 
    % grid on;
    % % axis([0 Nt/fs ylim])  
    % xlabel('time[s]')
    % legend('left or right')
    % title('untrimmed file data')
    %----------select useful portion------------
    ns=round(ts*fs)+1;nf=round(tf*fs);
    rev4_1=rev4_1(ns:nf);% after this the lr dimension is gone
    Nt=length(rev4_1);
    % disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    %-----------------plots------------
    f=[0:Nt/2]*fs/Nt;
    %Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
    Rev4=abs(fft(rev4_1));
    % figure(30)
    % plot(f,20*log10(Rev4(1:Nt/2+1)),'b');
    % grid on;
    % xlabel('freq[Hz]')
    % ylabel('Power Spectrum [dB]')
    % legend('N_t data','Location','Best');
    % title('PSD')
    %-------------------get rough estimate of test freq--------------------
    [M,I]=max(Rev4(1:Nt/2+1));
    freq_ref=(I-1)*fs/Nt;
    %------------------section spectra, get weighted line freq-------------
    nfft=2^12 ;%disp(['nfft: ' num2str(nfft)]) %not critical
    fractional_bin=1+freq_ref*nfft/fs;
    nref=1+round((freq_ref/fs)*nfft);%freq bin nearest reference
    nseg=floor(2*Nt/nfft-1);
    w=window(@blackmanharris,nfft,'periodic');
    n_sum=7;% a blackmanharris window allows smaller range
    for k=1:nseg
        rev=w.*rev4_1((k-1)*nfft/2+1:(k+1)*nfft/2);
        Prev=abs(fft(rev)).^2;% power in each bin
        P(k)=0;Pw(k)=0;%initialize power and weighted power sums
        for p=-n_sum:n_sum
            % k
            % nref
            % p
            P(k)=P(k)+Prev(nref+p);
            Pw(k)=Pw(k)+Prev(nref+p)*(nref-1+p)*fs/nfft;
        end
        freq(k)=Pw(k)/P(k);%power weighted frequency average
    end
    freq(1)=freq(3);freq(2)=freq(3);%%%%%% 2i2 seems to need this %%%%%%%%%%
    tseg=[0:nseg-1]*(nfft/2)/fs;
    % figure(40)
    % plot(tseg,freq)
    % grid on;
    % xlabel('Time[sec]')
    % ylabel('Freq[Hz]')
    % axis([xlim 3149.9932 3149.9938])
    % axis([xlim ylim])
    % title([filename ' ref:' num2str(freq_ref) ' nsum:' num2str(n_sum)])
    % ---------- ---closeup plot--------
    % figure(50)
    % plot(tseg,freq)
    % grid on;
    % axis([0 5 ylim])
    % xlabel('Time[sec]')
    % ylabel('Freq[Hz]')
    % title([filename ' ref:' num2str(freq_ref) ' nsum:' num2str(n_sum)])


    wow = max(freq) - min(freq); %%peak to peak measurement

    % disp('-------------------finished--------------------') 
end