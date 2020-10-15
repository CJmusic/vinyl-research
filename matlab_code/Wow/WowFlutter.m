% Record Wow Measurement
% using 3150 track
% John Vanderkooy
% April 2019
%  
filename = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/wow/maxbarrelzones1a.wav'; ts = 200;tf = 230;
[rev4_1,fs]=audioread(filename);

[wow, cho] = WowFlutterTest(rev4_1)
wow 
cho

% function [wow, centreholeoffset] = WowFlutter(rev4_1 );
function [wow, centreholeoffset] = WowFlutterTest(rev4_1);
    disp('--------------WOW AND FLUTTER--------------')
    filename = 'wow';
    fs = 96000;
    ts = length(rev4_1)*(1/4)/fs;
    tf = length(rev4_1)*(3/4)/fs;
    Nt=length(rev4_1);
    t=linspace(0,(Nt-1)/fs,Nt)';%column vector
    %-----------------prefiltering-----------
    disp('prefiltering')
    [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
    %rev4_1=filter(b,a,rev4_1);

    %----------select useful portion------------
    ns=round(ts*fs)+1;
    nf=round(tf*fs);
    rev4_1=rev4_1(ns:nf);% after this the lr dimension is gone
    Nt=length(rev4_1);

    f=[0:Nt/2]*fs/Nt;
    %Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
    disp('taking fft')

    Rev4=abs(fft(rev4_1));

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
    nseg
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
    disp('finished for loop')
    freq(1)=freq(3);freq(2)=freq(3);%%%%%% 2i2 seems to need this %%%%%%%%%%
    tseg=[0:nseg-1]*(nfft/2)/fs;

    disp('calculating wow')


    wow = max(freq) - min(freq); %%peak to peak measurement

    disp('-------------------finished--------------------') 

    %% NOW CALCULATE THE CENTRE HOLE OFFSET 
    NUM = X(5)*[(2*pi)^3*X(2)*X(3)*X(4) 0 0 0];% s^3 character
    DEN = conv(conv(conv([1 2*pi*X(2)],[1 2*pi*X(3)]) ,[1 2*pi*X(4)]), [1 2*pi*X(1)]); 
    % Bilinear transformation of analog design to get the digital filter.
    fsn=fs/(nfft/2);%rate of freq_deviation calculations
    disp(['WF sampling freq[Hz]: ' num2str(fsn)]);disp(' ')
    [b,a] = bilinear(NUM,DEN,fsn);
    %-----------------------------------------
    figure(55)
    W=logspace(-1,0,100);
    [H,F]=freqz(b,a,W);
    semilogx(W,(abs(H)))
    grid on
    %-----------------------------------
    freq=freq-sum(freq)/nseg;% remove DC from WF sampled deviation array

    figure(60)
    plot(tseg,freq)
    grid on;
    axis([xlim ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('freq(no DC)')

    figure(80)
    plot(tseg,freq)
    grid on;
    axis([1 6 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom freq(no DC)')
    %-----------------get fundamental wow amplitude---------
    % 'freq': demodulated frequencies, sampled at fsn, no DC, nseg points
    % freq=sin(2*pi*1.0331*(1:nseg)/fsn);%test window amplitude
    FREQ=(1/0.2155)*(2/nseg)*fft(freq'.*flattopwin(nseg));%extra factor for flattop
    wamplitude=max(abs(FREQ))
    R=100;%wow track radius in mm$$$$$$$$$$$$$$$$$$$$$$
    centreholeoffset=R*wamplitude/test_freq %offset in mm
    fw=(0:floor(nseg/2))*fsn/nseg;%wow frequencies
    figure(85)
    plot(fw,abs(FREQ(1:floor(nseg/2+1))))
    grid on
    axis([0 3 ylim])
    xlabel('wow freq [Hz]')
    ylabel('wow spectrum')
    %--------------apply W&F weighting---------------
    % WFfreq=filter(b,a,freq);

    figure(90)
    plot(tseg,WFfreq)
    grid on;
    axis([1 6 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom weighted WFfreq')
    %----------------characterize W&F result-----------------
    freqrms=rms_response(freq);
    disp(['rms unweighted freq variation [Hz]: ' num2str(freqrms)])
    freq_unweighted=freqrms/test_freq;
    disp(['rms unweighted W&F: ' num2str(freq_unweighted)])

    WFrms=rms_response(WFfreq);
    disp(['rms weighted freq variation [Hz]: ' num2str(WFrms)])
    WF_weighted=WFrms/test_freq;
    disp(['rms weighted W&F: ' num2str(WF_weighted)])
    disp('-------------------finished--------------------') 

end