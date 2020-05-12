%  
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
addpath('/Users/cz/Code/vinyl-research/matlab_code')
addpath('/Users/cz/Code/vinyl-research/matlab_code/audio_functions')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/from_John/Sept 26')
addpath('/Users/cz/Code/vinyl-research/matlab_code/Lacquer')
addpath('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lacquer_recordings/')
addpath('/Users/cz/Code/vinyl-research/matlab_code/A0137B0137')
% file = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/003141_A0000B0000r30a.wav'

% folder = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_bin/A0000B0000/';
% addpath('/Volumes/AUDIOBANK/audio_files/A0137B0137/')
folder = ('/Volumes/AUDIOBANK/audio_files/')
% addpath('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
% folder = ('/Volumes/AUDIOBANK/audio_files/duplicaterecordingtest/')
% addpath('/Volumes/AUDIOBANK/audio_files/duplicatefiletest/')
% folder = ('/Volumes/AUDIOBANK/audio_files/multiplerecordstest/')
% folder = ('/Volumes/AUDIOBANK/audio_files/multiplerecordstest/')
% folder = ('/Volumes/AUDIOBANK/audio_files/A01B01test/')
% folder = ('/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/A0000B0000_misc/wow/')
folder = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lacquer_recordings/lacquerwow/';



% ~~~~ MAC ENDS ~~~~ %

pressingID = 'A0137B0137'


disp(['loading folder...:', folder])
% files = dir(strcat(folder,'*.wav'))
files = dir(fullfile(folder,'*.wav'))


for i = (1:length(files)) %%loop through records
    % filename = files(i).name
    % file= strcat(files(i).folder,'/',filename)
    % tracks=SeperateTracks(file);
    % rev4_1 = tracks('3150Hz');


    filename = files(i).name
    file= strcat(files(i).folder,'/',filename)
    [rev4_1, fs] = audioread(file);
    rev4_1 = Aweighting()
    rev4_1 = rev4_1(30*fs:45*fs,:);

    % if contains(filename,'lacquerpartone-offset126.wav');
    %     lacquer = LacquerProcess(filename, 12.6);
    %     rev4_1 = lacquer{'3150Hz'};
    % end

    % if contains(filename, 'lacquerpartthree-offset000.wav');
    %     lacquer = LacquerProcess(filename, 0);
    %     rev4_1 = lacquer{'3150Hz'}  ;  
    % end


    fs = 96000;
    ts = length(rev4_1)*(1/4)/fs;
    tf = length(rev4_1)*(3/4)/fs;

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
    grid on; hold on;
    % axis([0 Nt/fs ylim])
    xlabel('time[s]')
    legend('left or right')
    title('untrimmed file data')
    %----------select useful portion------------
    ns=round(ts*fs)+1;nf=round(tf*fs);
    rev4_1=rev4_1(ns:nf,lr);% after this the lr dimension is gone ~~~~~~~~~ THIS THE TRIMMING LINE ~~~~~~~~~~~~
    Nt=length(rev4_1);
    disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    %-----------------plots------------
    f=[0:Nt/2]*fs/Nt;
    %Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
    Rev4=abs(fft(rev4_1));
    
    figure(30)
    plot(f,20*log10(Rev4(1:floor(Nt/2+1))),'b');
    grid on; hold on;
    xlabel('freq[Hz]')
    ylabel('Power Spectrum [dB]')
    legend('N_t data','Location','Best');
    title('PSD')
    %-------------------get rough estimate of test freq--------------------
    [M,I]=max(Rev4(1:floor(Nt/2+1)));
    test_freq=(I-1)*fs/Nt;
    disp(['test freq [Hz]: ' num2str(test_freq)])
    %------------------section spectra, get weighted line freq-------------
    nfft=2^12;% not critical
    disp(['nfft: ' num2str(nfft)]);disp(' ')
    fractional_bin=1+test_freq*nfft/fs;
    nref=1+round((test_freq/fs)*nfft);%freq bin nearest reference
    nseg=floor(2*Nt/nfft-1);% prepare for 50% overlap
    w=window(@blackmanharris,nfft,'periodic');
    n_sum=7;% a blackmanharris window allows smaller range
    for k=1:nseg
        rev=w.*rev4_1((k-1)*nfft/2+1:(k+1)*nfft/2);% 50% overlap
        Prev=abs(fft(rev)).^2;% power in each bin
        P(k)=0;Pw(k)=0;%initialize power and weighted power sums
        for p=-n_sum:n_sum
            P(k)=P(k)+Prev(nref+p);
            Pw(k)=Pw(k)+Prev(nref+p)*(nref-1+p)*fs/nfft;
        end
        freq(k)=Pw(k)/P(k);%power weighted frequency average
    end
    freq(1)=freq(2);%freq(2)=freq(3);%%%%%% 2i2 seems to need this %%%%%%%%%%
    tseg=[0:nseg-1]*(nfft/2)/fs;
    
    figure(40)
    plot(tseg,freq)
    grid on; hold on;
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    %axis([xlim 3149.9932 3149.9938])
    axis([xlim ylim])
    title('freq(t)')
    %-------------closeup plot--------
    figure(50)
    plot(tseg,freq)
    grid on; hold on;
    axis([0 5 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom freq(t)')

    figure(100)
    subplot(2,2,i)
    plot(tseg,freq)
    grid on; hold on;
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
    fsn=fs/(nfft/2);
    disp(['WF sampling freq[Hz]: ' num2str(fsn)]);disp(' ')
    [b,a] = bilinear(NUM,DEN,fsn);
    %-----------------------------------------
    freq=freq-sum(freq)/nseg;% remove most of DC
    
    figure(60)
    plot(tseg,freq)
    grid on; hold on;
    axis([xlim ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('freq(no DC)')
    
    figure(80)
    plot(tseg,freq)
    grid on; hold on;
    axis([0 5 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom freq(no DC)')
    


    %--------------apply W&F weighting---------------
    WFfreq=filter(b,a,freq);
    % WFfreq=freq;% no WF filter
    
    figure(90)
    plot(tseg,WFfreq)
    grid on; hold on;
    axis([0 5 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    title('zoom weighted WFfreq')

    figure(90 + i)
    plot(tseg,WFfreq,'k')
    grid on; hold on;
    axis([0 5 ylim])
    xlabel('Time[sec]')
    ylabel('Freq[Hz]')
    % title('zoom weighted WFfreq')

    %----------------characterize W&F result-----------------
    freqrms=rms_response(freq);
    disp(['rms unweighted freq variation: ' num2str(freqrms)])
    freq_unweighted=freqrms/test_freq;
    disp(['rms unweighted W&F: ' num2str(freq_unweighted)])
    
    WFrms=rms_response(WFfreq);
    disp(['rms weighted freq variation: ' num2str(WFrms)])
    WF_weighted=WFrms/test_freq;
    disp(['rms weighted W&F: ' num2str(WF_weighted)])
    disp('-------------------finished--------------------') 
    


    %OLD CODE BELOW
% disp('--------------WOW AND FLUTTER--------------')
    % filename = 'wow';
    % fs = 96000;
    % tracks=SeperateTracks(file);
    % rev4_1 = tracks('3150Hz');
    % ts = length(rev4_1)*(1/4)/fs;
    % tf = length(rev4_1)*(3/4)/fs;
    % %lr=1; %1=left, 2=right
    % % disp(['lr: ' num2str(lr)])
    % Nt=length(rev4_1);
    % % disp(['fs: ' num2str(fs) '  N_total: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    % t=linspace(0,(Nt-1)/fs,Nt)';%column vector
    % %-----------------prefiltering-----------
    % [b,a]=butter(2,2*100/fs,'high');% not really necessary with fft filter
    % rev4_1=filter(b,a,rev4_1);
    % %-------------------plot all data------------
    % figure(20)
    % plot(t,rev4_1(1:Nt),'b') 
    % hold on;
    % grid on;
    % % axis([0 Nt/fs ylim])  
    % xlabel('time[s]')
    % legend('left or right')
    % title('untrimmed file data')
    % %----------select useful portion------------
    % ns=round(ts*fs)+1;nf=round(tf*fs);
    % rev4_1=rev4_1(ns:nf);% after this the lr dimension is gone
    % Nt=length(rev4_1);
    % % disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
    % %-----------------plots------------
    % f=[0:Nt/2]*fs/Nt;
    % %Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
    % Rev4=abs(fft(rev4_1));
    % figure(30)
    % plot(f,20*log10(Rev4(1:Nt/2+1)),'b');
    % hold on;
    % grid on;
    % xlabel('freq[Hz]')
    % ylabel('Power Spectrum [dB]')
    % legend('N_t data','Location','Best');
    % title('PSD')
    % %-------------------get rough estimate of test freq--------------------
    % [M,I]=max(Rev4(1:Nt/2+1));
    % freq_ref=(I-1)*fs/Nt;
    % %------------------section spectra, get weighted line freq-------------
    % nfft=2^12 ;%disp(['nfft: ' num2str(nfft)]) %not critical
    % fractional_bin=1+freq_ref*nfft/fs;
    % nref=1+round((freq_ref/fs)*nfft);%freq bin nearest reference
    % nseg=floor(2*Nt/nfft-1);
    % w=window(@blackmanharris,nfft,'periodic');
    % n_sum=7;% a blackmanharris window allows smaller range
    % for k=1:nseg
    %     rev=w.*rev4_1((k-1)*nfft/2+1:(k+1)*nfft/2);
    %     Prev=abs(fft(rev)).^2;% power in each bin
    %     P(k)=0;Pw(k)=0;%initialize power and weighted power sums
    %     for p=-n_sum:n_sum
    %         % k
    %         % nref
    %         % p
    %         P(k)=P(k)+Prev(nref+p);
    %         Pw(k)=Pw(k)+Prev(nref+p)*(nref-1+p)*fs/nfft;
    %     end
    %     freq(k)=Pw(k)/P(k);%power weighted frequency average
    % end
    % freq(1)=freq(3);freq(2)=freq(3);%%%%%% 2i2 seems to need this %%%%%%%%%%
    % tseg=[0:nseg-1]*(nfft/2)/fs;
    % figure(40)
    % plot(tseg,freq)
    % grid on; hold on;
    % xlabel('Time[sec]')
    % ylabel('Freq[Hz]')
    % axis([xlim 3149.9932 3149.9938])
    % axis([xlim ylim])
    % title([filename ' ref:' num2str(freq_ref) ' nsum:' num2str(n_sum)])
    % % ---------- ---closeup plot--------
    % figure(50)
    % plot(tseg,freq)
    % grid on; hold on;
    % axis([0 5 ylim])
    % xlabel('Time[sec]')
    % ylabel('Freq[Hz]')
    % title([filename ' ref:' num2str(freq_ref) ' nsum:' num2str(n_sum)])


    % wow = max(freq) - min(freq); %%peak to peak measurement

    % disp('-------------------finished--------------------') 
end