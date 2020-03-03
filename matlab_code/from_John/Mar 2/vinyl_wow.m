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
%
try
    pkg load signal %for Octave
catch
end

% filename='JV3150-96-AG33250A.wav';ts=0;tf=60;
% filename='Mann3150-96-AG33250A.wav';ts=0;tf=60;
% filename='JV3150-96-DS345.wav';ts=0;tf=60;
% filename='Mann3150-96-DS345.wav';ts=0;tf=60;
% 
% filename='wow3150oscillator.wav';ts=10;tf=60;
% filename='wow3150funcgen.wav';ts=10;tf=60;
% filename='1kHz_DualPre.wav';ts=0;tf=50;
% 
% filename='lacquerA_sil_1kHz_beginning.wav';ts=17;tf=24.7;
%filename='lacquerA_silence_1kHz.wav';ts=50;tf=90;
% 
% filename='052319_pioneer3150.wav';ts=30;tf=70;% tts=30;ttf=70;
% filename='wow3150innerJ15aDual.wav';ts=8;tf=58;

%filename='r27wow0deg.wav';ts=10.0;tf=60.0;
%filename='r27wow90deg.wav';ts=25.0;tf=75.0;
filename='r27wow180deg.wav';ts=10.0;tf=60.0;
% filename='r27wow270deg.wav';ts=20.0;tf=70.0;
% filename='r27wowinner270deg.wav';ts=10.0;tf=60.0;

% filename='analogmagik_silent_1kHz.wav';ts=8.61;tf=15.7;%44.1 kHz

[rev4_1,fs]=audioread(filename);
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
ns=round(ts*fs)+1;nf=round(tf*fs);
rev4_1=rev4_1(ns:nf,lr);% after this the lr dimension is gone
Nt=length(rev4_1);
disp(['N_analyze: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
%-----------------plots------------
f=[0:Nt/2]*fs/Nt;
%Rev4=abs(fft(window(@blackmanharris,Nt,'periodic').*rev4_1));
Rev4=abs(fft(rev4_1));

figure(30)
plot(f,20*log10(Rev4(1:floor(Nt/2+1))),'b');
grid on;
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
fsn=fs/(nfft/2);
disp(['WF sampling freq[Hz]: ' num2str(fsn)]);disp(' ')
[b,a] = bilinear(NUM,DEN,fsn);
%-----------------------------------------
freq=freq-sum(freq)/nseg;% remove most of DC

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
axis([0 5 ylim])
xlabel('Time[sec]')
ylabel('Freq[Hz]')
title('zoom freq(no DC)')

%--------------apply W&F weighting---------------
WFfreq=filter(b,a,freq);
% WFfreq=freq;% no WF filter

figure(90)
plot(tseg,WFfreq)
grid on;
axis([0 5 ylim])
xlabel('Time[sec]')
ylabel('Freq[Hz]')
title('zoom weighted WFfreq')
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
