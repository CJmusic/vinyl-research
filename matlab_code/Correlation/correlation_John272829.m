% Record Surface Noise analysis
% using correlation & coherence between different recordings same groove
% John Vanderkooy
% Feb. 2019
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

addpath('/Users/cz/Code/vinyl-research/matlab_code/Common')


trackname = 'sweep2'
ts = 2; tf = 22;
ts2 = ts;
ts3 = ts;

ts = 1;
ts2 =  1;
ts3 = 1;
ts4 = 1;
ts5 = 1;


tracks1 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1550.844.wav');
tracks2 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1553.621.wav');
tracks3 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1558.697.wav');
tracks4 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1559.267.wav');
tracks5 = SeperateTracks('/Volumes/AUDIOBANK/audio_files/SameRecordTest2/r300a1559.287.wav');

rev4_1=tracks1(trackname);
rev4_2=tracks2(trackname);
rev4_3=tracks3(trackname);
rev4_4=tracks4(trackname);
rev4_5=tracks5(trackname);
fs = 96000;
lr=1; %1=left, 2=right
disp(['lr: ' num2str(lr)])
Nt=length(rev4_1);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration_1: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-------------------plot all data------------
figure(20)
plot(t,rev4_1(1:Nt,lr),'b') 
grid on;
% plot(t,rev4_2(:,lr),'r')
% plot(t,rev4_3(:,lr),'g')
% axis([0 Nt/fs ylim])
% xlabel('time[s]')
% legend('left','right')
title('untrimmed file data')
%----------------choose equal length files----------------
% ns and nt should be defined with filename above
ns=round(ts*fs);nf=round(tf*fs);
rev4_1=rev4_1(ns:nf,:);
Nt=length(rev4_1);
ns2=round(ts2*fs);
rev4_2=rev4_2(ns2:ns2+Nt-1,:);
ns3=round(ts3*fs);
rev4_3=rev4_3(ns3:ns3+Nt-1,:);
ns4=round(ts4*fs);
rev4_4=rev4_4(ns4:ns4+Nt-1,:);
ns5=round(ts5*fs);
rev4_5=rev4_5(ns5:ns5+Nt-1,:);

t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------xcorr plots------------
[b,a]=butter(2,2*200/fs,'high');% remove LF arm resonance
rev4_1f=filter(b,a,rev4_1);
rev4_2f=filter(b,a,rev4_2);
rev4_3f=filter(b,a,rev4_3);
rev4_4f=filter(b,a,rev4_4);
rev4_5f=filter(b,a,rev4_5);
maxlag=2^16;
C0=xcorr(rev4_1f(:,1),rev4_1f(:,2),maxlag);
C1=xcorr(rev4_1f(:,lr),rev4_2f(:,lr),maxlag);
C2=xcorr(rev4_1f(:,lr),rev4_3f(:,lr),maxlag);
C3=xcorr(rev4_1f(:,lr),rev4_4f(:,lr),maxlag);
C4=xcorr(rev4_1f(:,lr),rev4_5f(:,lr),maxlag);
figure(21)
% plot(C0,'k')
hold on;grid on;
plot(C1,'b')
plot(C2,'r')
plot(C3,'g')
plot(C4,'o')
% legend()
%axis([ fs/100 fs/2 0 1.05])
title('cross correlations')
%-----------------------------
figure(23)
plot(t,rev4_1(:,lr),'b') 
grid on;hold on;
plot(t,rev4_2(:,lr),'r')
plot(t,rev4_3(:,lr),'g')
axis([0 Nt/fs ylim])
xlabel('time[s]')
legend('1','2','3')
title('all recorded data')
%------------fft of a click--------------
% nft=2^10;
% f=([1:nft/2+1]-1)*fs/nft;
% d=rev4_1(1:nft);
% g=fft(d);
% figure(24)
% plot(f(3:nft/2+1),abs(g(3:nft/2+1)))
% xlabel('Freq[Hz]')
% title('click spectrum')
% grid on;
%---------------histogram----------------
figure(25)
hist(rev4_1,2048)
grid on;
title('histogram')
% %-----------calculate # of clicks-----------
% [b,a]=butter(2,2*1000/fs,'high');%1kHz is good
% rev4_hf=filter(b,a,rev4_1);
% figure(30)
% plot(rev4_hf)
% title('HP filtered silent grooves')
% grid on
% n_stft=2^9;
% n_seg=floor(Nt/n_stft-1)
% rms_tot=sqrt(sum(rev4_hf(:,lr).^2)/Nt)
% factor=zeros(50,1);
% num=zeros(50,1);
% for m=1:50
%   factor(m)=2+0.2*m;%about 3.4 is optimal
%   n=0;
%    for k=1:n_seg
%      rev4_seg(:,k)=rev4_hf((k-1)*n_stft+1:k*n_stft,lr);
%      rms_seg(k)=sqrt(sum(rev4_seg(:,k).^2)/n_seg);
%        if rms_seg(k)>factor*rms_tot
%            n=n+1;
%            disp(['factor: ' num2str(factor(m)) '  click at: ' num2str(k)])
%            num(m)=n;
%        end
%     end
% end
% figure(40)
% plot(factor,num)
% xlabel('factor')
% ylabel('number')
% grid on
%--------------psd of whole segment----------------
nfft=2^14;
window=hanning(nfft,'periodic');
[Pxx,f]=pwelch(rev4_1(:,lr),window,nfft/2,nfft,fs,'psd');
[Pxx2,f]=pwelch(rev4_2(:,lr),window,nfft/2,nfft,fs,'psd');
[Pxx3,f]=pwelch(rev4_3(:,lr),window,nfft/2,nfft,fs,'psd');
psdrms=sqrt(sum(Pxx)*fs/nfft)
%-----------------plot psd for each test-----------------------------
figure(50)
semilogx(f,10*log10(Pxx(1:nfft/2+1)),'b');
grid on;hold on;
semilogx(f,10*log10(Pxx2),'r');
semilogx(f,10*log10(Pxx3),'g');
axis([xlim -140 -30])
xlabel('freq[Hz]')
ylabel('Power Spectrum')
legend('1','2','3','Location','Best');
title('PSD')
%----------------------coherence-----------------
[Cc0,fc]=mscohere(rev4_1(:,1),rev4_1(:,2),window,nfft/2,nfft,fs);
[Cc1(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_2(:,lr),window,nfft/2,nfft,fs);
[Cc2(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_3(:,lr),window,nfft/2,nfft,fs);
[Cc3(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_3(:,lr),window,nfft/2,nfft,fs);
[Cc4(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_4(:,lr),window,nfft/2,nfft,fs);
[Cc5(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_5(:,lr),window,nfft/2,nfft,fs);
%---------------plot coherences---------
figure(70)
% semilogx(fc,Cc0,'k')
grid on;hold on;
semilogx(fc,Cc1(:,lr),'r')
semilogx(fc,Cc2(:,lr),'g')
semilogx(fc,Cc3(:,lr),'c')
semilogx(fc,Cc4(:,lr),'p')
semilogx(fc,Cc5(:,lr),'b')
axis([0,fs/2,0,1])
% legend('27L&R','27&28','27&29','28&29','Location','Best')
xlabel('frequency [Hz]')
title('Run Coherences')
disp('-------------------finished--------------------') 
