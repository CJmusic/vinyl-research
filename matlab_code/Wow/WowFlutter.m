% Record Wow Measurement
% using 3150 track
% John Vanderkooy
% April 2019
%  

clear all; clc;close all;
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')

tracks = SeperateTracks('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/031419_A0000B0000r029a.wav');
fs = 96000;
rev4_1 = tracks('3150Hz');

[wow, cho] = WowFlutterTest(rev4_1)
wow 
cho

% function [wow, centreholeoffset] = WowFlutter(rev4_1 );
function [wow, centreholeoffset] = WowFlutterTest(rev4_1);
    L = length(rev4_1);
    rev4_1 = rev4_1(L/3:2*L/3,1);
    fs = 96000;
    % [rev4_1,fs]=audioread(filename);
    lr=1; %1=left, 2=right
    disp(['left_right: ' num2str(lr)])
    Nt=length(rev4_1);
    disp(['fs: ' num2str(fs) '  N_total: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    t=linspace(0,(Nt-1)/fs,Nt)';%column vector
    %-----------------prefiltering-----------
    [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
    %rev4_1=filter(b,a,rev4_1);
    %-------------------plot all data------------
    figure(20)
    plot(t,rev4_1(1:Nt,lr),'b') 
    grid on;
    % axis([0 Nt/fs ylim])
    xlabel('time[s]')
    legend('left or right')
    title('untrimmed file data')
    %----------select useful portion------------
    % ns=round(ts*fs)+1;nf=round(tf*fs);
    ns = length(rev4_1)
    % rev4_1=rev4_1(ns:nf,lr);% after this the lr dimension is gone
    rev4_1 = rev4_1(:,1);
    Nt=length(rev4_1);
    disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    figure(21)
    plot(rev4_1(1:Nt,lr),'b') 
    grid on;
    % axis([0 Nt/fs ylim])
    xlabel('time[s]')
    legend('left or right')
    title('trimmed file data')

    %-----------------plots------------
    f=[0:Nt/2]*fs/Nt;
    %Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
    Rev4=abs(fft(rev4_1));

    figure(30)
    plot(f,20*log10(Rev4(1:floor(Nt/2+1))),'b');
    grid on;
    xlabel('freq[Hz]')
    ylabel('Power Spectrum [dB]')
    %axis([3150-50 3150+50 80 120])
    legend('N_t data','Location','Best');
    title('PSD')
    %-------------------get rough estimate of test freq--------------------
    [M,I]=max(Rev4(1:floor(Nt/2+1)));
    test_freq=(I-1)*fs/Nt;
    disp(['rough test freq [Hz]: ' num2str(test_freq)])
    % -----------it might be wise to calculate test freq better
    %test_freq=3150 %this doesn't change things at all
    %------------------section spectra, get weighted line freq-------------
    nfft=2^12;% not critical
    disp(['nfft: ' num2str(nfft)]);disp(' ')
    %fractional_bin=1+test_freq*nfft/fs;
    nref=1+round(test_freq*nfft/fs);%freq bin nearest reference
    nseg=floor(2*Nt/nfft-1);% prepare for 50% overlap
    w=window(@hann,nfft,'periodic');
    % w=window(@blackmanharris,nfft,'periodic');
    % w=window(@flattopwin,nfft,'periodic');
    % w=window(@rectwin,nfft);
    % w=window(@triang,nfft);
    n_sum=10;% a better window allows smaller range
    for k=1:nseg
        rev=w.*rev4_1((k-1)*nfft/2+1:(k+1)*nfft/2);% 50% overlap
        Prev=abs(fft(rev)).^2;% power in each bin
        P(k)=0;Pw(k)=0;%initialize power and weighted power sums
        for p=-n_sum:n_sum
            k
            p
            nref
            P(k)=P(k)+Prev(nref+p);
            Pw(k)=Pw(k)+Prev(nref+p)*(nref-1+p)*fs/nfft;
        end
        freq(k)=Pw(k)/P(k);%power weighted frequency average
    end
    freq(1)=freq(2);%freq(2)=freq(3);%%%%%% 2i2 seems to need this %%%%%%%%%%
    tseg=[0:nseg-1]*(nfft/2)/fs;%time variable for ts to tf points

    figure(40)
    plot(tseg,freq)
    grid on;
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    %axis([xlim 3149.9932 3149.9938])
    axis([xlim ylim])
    title('freq(t)')
    %-------------closeup plot--------
    figure(50)
    plot(tseg,freq)
    grid on;
    axis([0 5 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom freq(t)')
    %% -------------------------WF-wtg table-----------------------------------
    % fr=[0.1 0.19 0.43 0.77 1.0 2.0 5.0 10.0 20.0 50.0 165 1000];
    % dBWFtable=[-57 -40 -20 -10 -7.25 -1.52 0 -1 -4 -10 -20 -36];
    %---------------------------------------
    f1 = 15.0;%HF rolloff
    f2 = 0.65;%LF rollup
    f3 = 0.9;%LF rollup
    f4 = 1.;%LF rollup
    WF4 = 0.71;%sets dB gain
    X=[f1 f2 f3 f4 WF4];
    %---------Analog W&F-weighting filter from filter convolution---------
    NUM = X(5)*[(2*pi)^3*X(2)*X(3)*X(4) 0 0 0];% s^3 character
    DEN = conv(conv(conv([1 2*pi*X(2)],[1 2*pi*X(3)]) ,[1 2*pi*X(4)]), [1 2*pi*X(1)]); 
    % Bilinear transformation of analog design to get the digital filter.
    fsn=fs/(nfft/2);%rate of freq_deviation calculations
    disp(['WF sampling freq[Hz]: ' num2str(fsn)]);disp(' ')
    [b,a] = bilinear(NUM,DEN,fsn);
    %-----------------------------------------
    % figure(55)
    % W=logspace(-1,0,100);
    % [H,F]=freqz(b,a,W);
    % semilogx(W,(abs(H)))
    % grid on
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
    WFfreq=filter(b,a,freq);

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

    wow = max(freq) - min(freq); %%peak to peak measurement

end