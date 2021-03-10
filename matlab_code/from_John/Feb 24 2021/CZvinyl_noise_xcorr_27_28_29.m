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
filename='27a-cutsilence.wav';ts=2.000;tf=22.000;
[rev4_1,fs]=audioread(filename);
filename='28a-cutsilence.wav';ts2=2.000+(0+34282)/fs;
[rev4_2,fs]=audioread(filename);
filename='29a-cutsilence.wav';ts3=2.400-(0+23770)/fs;
[rev4_3,fs]=audioread(filename);
filename='lacquer-cutsilence.wav';ts4=1.000;tf4 = 11;
[rev4_4,fs]=audioread(filename);
filename='lacquer-cutsilence2.wav';ts5=1.000;tf5 = 3;
[rev4_5,fs]=audioread(filename);
filename='lacquer-cutsilence3.wav';ts6=1.000;tf6 = 4;
[rev4_6,fs]=audioread(filename);
filename='lacquer-cutsilence4.wav';ts7=1.000;tf7 = 9;
[rev4_7,fs]=audioread(filename);


lr=1; %1=left, 2=right
disp(['lr: ' num2str(lr)])
Nt=length(rev4_1);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration_1: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-------------------plot all data------------
plotnum = 0;
plotnum = plotnum + 1;
figure(plotnum)
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
ns4=round(ts4*fs); Nt4 = 10*fs;
ns5=round(ts5*fs); Nt5 = 2*fs;
ns6=round(ts6*fs); Nt6 = 3*fs;
ns7=round(ts7*fs); Nt7 = 8*fs;
rev4_1=rev4_1(ns:nf,:);
Nt=length(rev4_1);
ns2=round(ts2*fs);
rev4_2=rev4_2(ns2:ns2+Nt-1,:);
ns3=round(ts3*fs);
rev4_3=rev4_3(ns3:ns3+Nt-1,:);
ns4=round(ts4*fs);
rev4_4=rev4_4(ns4:ns4+Nt4-1,:);
ns5=round(ts5*fs);
rev4_5=rev4_5(ns5:ns5+Nt5-1,:);
ns6=round(ts6*fs);
rev4_6=rev4_6(ns6:ns6+Nt6-1,:);
ns7=round(ts7*fs);
rev4_7=rev4_7(ns7:ns7+Nt6-1,:);





t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------xcorr plots------------
[b,a]=butter(2,2*200/fs,'high');% remove LF arm resonance
rev4_1f=filter(b,a,rev4_1);
rev4_2f=filter(b,a,rev4_2);
rev4_3f=filter(b,a,rev4_3);
rev4_4f=filter(b,a,rev4_4);
rev4_5f=filter(b,a,rev4_5);
rev4_6f=filter(b,a,rev4_6);
rev4_7f=filter(b,a,rev4_7);
maxlag=2^16;
C0=xcorr(rev4_1f(:,1),rev4_1f(:,2),maxlag);
C1=xcorr(rev4_1f(:,lr),rev4_2f(:,lr),maxlag);
C2=xcorr(rev4_1f(:,lr),rev4_3f(:,lr),maxlag);

plotnum = plotnum + 1;
figure(plotnum)
plot(C0,'k')
hold on;grid on;
plot(C1,'b')
plot(C2,'r')
legend('27 L&R','27-28','27-29')
%axis([ fs/100 fs/2 0 1.05])
title('cross correlations')
%-----------------------------
plotnum = plotnum + 1;
figure(plotnum)
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
plotnum = plotnum + 1;
figure(plotnum)
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
plotnum = plotnum + 1;
figure(plotnum)
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

%~~~~~~~~~~~~coherences between records~~~~~~~~~~%
[Cc1(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_2(:,lr),window,nfft/2,nfft,fs);
[Cc2(:,lr),fc]=mscohere(rev4_1(:,lr),rev4_3(:,lr),window,nfft/2,nfft,fs);
[Cc3(:,lr),fc]=mscohere(rev4_2(:,lr),rev4_3(:,lr),window,nfft/2,nfft,fs);

[Cc4(:,lr),fc]=mscohere(rev4_4(:,1),rev4_4(:,2),window,nfft/2,nfft,fs);
[Cc5(:,lr),fc]=mscohere(rev4_5(:,1),rev4_5(:,2),window,nfft/2,nfft,fs);
[Cc6(:,lr),fc]=mscohere(rev4_6(:,1),rev4_6(:,2),window,nfft/2,nfft,fs);
[Cc7(:,lr),fc]=mscohere(rev4_7(:,1),rev4_7(:,2),window,nfft/2,nfft,fs);

%---------------plot coherences---------
plotnum = plotnum + 1;
figure(plotnum)
semilogx(fc,Cc1(:,lr),'r')
grid on;hold on;
semilogx(fc,Cc2(:,lr),'g')
semilogx(fc,Cc3(:,lr),'c')
axis([0,fs/2,0,1])
legend('27&28','27&29','28&29', 'Location','Best')
xlabel('frequency [Hz]')
title('Left channel coherence between records 27-29')
saveas(figure(plotnum),'recordCoherence.png')

plotnum = plotnum + 1;
figure(plotnum)
semilogx(fc,Cc4(:,lr),'r')
grid on; hold on
semilogx(fc,Cc5(:,lr),'g')
semilogx(fc,Cc6(:,lr),'c')
axis([0,fs/2,0,1])
legend('27&28 L&R','27&29 L&R','28&29 L&R', 'Location','Best')
xlabel('frequency [Hz]')
title('Left and right channel coherence between records 27-29')
saveas(figure(plotnum),'LRrecordCoherence.png')




%~~~~~~~~~~~~coherences between channels ~~~~~~~~~~%
[Cc1(:,lr),fc]=mscohere(rev4_1(:,1),rev4_1(:,2),window,nfft/2,nfft,fs);
[Cc2(:,lr),fc]=mscohere(rev4_2(:,1),rev4_2(:,2),window,nfft/2,nfft,fs);
[Cc3(:,lr),fc]=mscohere(rev4_3(:,1),rev4_3(:,2),window,nfft/2,nfft,fs);
[Cc4(:,lr),fc]=mscohere(rev4_4(:,1),rev4_4(:,2),window,nfft/2,nfft,fs);
[Cc5(:,lr),fc]=mscohere(rev4_5(:,1),rev4_5(:,2),window,nfft/2,nfft,fs);
[Cc6(:,lr),fc]=mscohere(rev4_6(:,1),rev4_6(:,2),window,nfft/2,nfft,fs);
[Cc7(:,lr),fc]=mscohere(rev4_7(:,1),rev4_7(:,2),window,nfft/2,nfft,fs);
%---------------plot coherences---------
plotnum = plotnum + 1;
figure(plotnum)
semilogx(fc,Cc1(:,lr),'r')
grid on;hold on;
semilogx(fc,Cc2(:,lr),'g')
semilogx(fc,Cc3(:,lr),'c')
axis([0,fs/2,0,1])
legend('27 L&R','28 L&R','29 L&R', 'Location','Best')
xlabel('frequency [Hz]')
title('Coherence between left and right channels')
saveas(figure(plotnum),'LRCoherence.png')


plotnum = 0;
plotnum = plotnum + 1;
figure(plotnum)
semilogx(fc,Cc4(:,lr),'k')
% grid on;hold on;
% semilogx(fc,Cc5(:,lr),'g')
% semilogx(fc,Cc6(:,lr),'c')
% semilogx(fc,Cc7(:,lr),'c')
grid on;
axis([0,fs/2,0,1])
legend('lacquer L&R', 'Location','northeast')
xlabel('frequency [Hz]')
title('Lacquer coherence between left and right channels')
saveas(figure(plotnum),'lacqeurCoherence.png')



disp('-------------------finished--------------------') 
