% Record Surface Noise analysis
% using xcorrelation and coherence between groove signals
% July. 2019% John Vanderkooy

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
offset=-20;% default %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename='lacquerA_silence_edited_onto_1kHz_declicked.wav';
% ts=0;tf=43.9;

% filename='043019_lacquerA_declicked_silence.wav';
% ts=10;tf=53;

filename='/Volumes/AUDIOBANK/audio_files/A0000B0000/031418_A0000B0000r027a1553.770/transition.wav';
ts=0.0;tf=84.6;% silent closely packed 30 grooves
t2s=71.0;t2f=101.6;% 1mm spaced 17 grooves

% filename='27b_declick_full_transition.wav';
% seems to have recording gap after ~20 grooves!!!
% ts=12.0;tf=64;% silent closely packed 30 grooves
%ts=68;tf=100;% silent closely packed 30 grooves
% t2s=71.0;t2f=101.6;% 1mm spaced 17 grooves

% filename='28b_declick_full_transition.wav';
% % seems to have recording gap after ~8 grooves!!!
% ts=80.0;tf=106.6;%13.4 70 silent closely packed 30 grooves
% t2s=71.0;t2f=101.6;% 1mm spaced 17 grooves

[sig,fs]=audioread(filename);
Nt=length(sig);
%-------------------plot all data------------
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
figure(10)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([xlim ylim])
xlabel('time[s]')
legend('left','right')
title('all raw data')
%------------------- filter--------------
[b,a]=butter(2,2*30/fs,'high');%remove most warp
%sig=filtfilt(b,a,sig);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------select closely spaced data----------------
ns=1+floor(ts*fs);nf=floor(tf*fs);
sig=sig(ns:nf,:);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lr=2 %1=lateral, 2=vertical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=length(sig);
disp([filename ' fs: ' num2str(fs) '  Nt:' num2str(Nt) '  close grooves duration:' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------------data plots------------
figure(20)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs -.1 .1])
xlabel('time[s]')
legend('left','right')
title('All closely packed groove data')

%------------------generate separate groove results-----------
% we'll make n_rev files
n_per_rev=round(fs*60/33.33333333)+offset %best match %%%%%%%%%%%%%%%
n_rev=floor(Nt/n_per_rev) % number of grooves

% select the groove responses
for ng=1:n_rev;% groove number
    sig_gr(:,:,ng)=sig(1+(ng-1)*n_per_rev:ng*n_per_rev,:);
%     sig_gr(:,1,ng)=sig(1+(ng-1)*n_per_rev:ng*n_per_rev,1)+sig(1+(ng-1)*n_per_rev:ng*n_per_rev,2);% lateral
%     sig_gr(:,2,ng)=sig(1+(ng-1)*n_per_rev:ng*n_per_rev,1)-sig(1+(ng-1)*n_per_rev:ng*n_per_rev,2);% vertical
end
%---------------edit groove responses to remove problems------------
% for k=round(1.3836*fs):round(1.3845*fs)
% sig_gr(k,:,:)=sig_gr(round(1.3836*fs),:,:);
% end
% sig_gr(1:0.41*fs,lr,:)=0;
% sig_gr(0.76*fs:n_per_rev,lr,:)=0;
%-----------check revolution accuracy------------
nlags=10000;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(25)
for k=2:n_rev
[C,lags]=xcorr(sig_gr(:,lr,1),sig_gr(:,lr,k),nlags);
plot(lags,C)
grid on;hold on;
end
title('xcorr wrt groove 1')
%-------------plot first bad groove----------
test_gr=22
figure(26)
plot(sig_gr(:,lr,test_gr),'b')
grid on;hold on;
plot(sig_gr(:,lr,test_gr+1),'r')
title('comparing grooves')
legend(['groove ' num2str(test_gr)],'next groove')
%---------------correlation between grooves--------------
for ng=1:n_rev %whole grooves
    cf(ng)=sum(sig_gr(:,lr,1).*sig_gr(:,lr,ng),1);% corr over samples
end
figure(30)
plot(cf/cf(1),'b')
grid on;
title('correlation between ref and each groove')
axis([0 28 -0.2 1.2])
xlabel('groove number')
ylabel('normalized corr coeff')
title('Correlation decay')
%---------time variable in plots------------
tg=t(1:n_per_rev);
minplot=-0.2;
maxplot=+0.2;
%-----------------plot signal for each groove-----------------------------
figure(40)
grid on;hold on;
for np=1:n_rev
plot(tg,sig_gr(:,lr,np),'b');% variables are sample#, left/right, groove#
end
axis([xlim ylim])
xlabel('time[s]')
title('groove signals')
% %-----------------plot L/R signal for a groove-----------------------------
gr_num=1
figure(50)
plot(tg,sig_gr(:,1,gr_num),'b') 
grid on;hold on;
plot(tg,sig_gr(:,2,gr_num),'r') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('lateral','vertical')
title(['stereo groove ' num2str(gr_num)])
%----------------------coherence-----------------
nfft=2^11;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nfft=n_per_rev/2;
smoothwidth=10;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=(0:nfft/2)*fs/nfft;

for k=1:n_rev-1;
C=mscohere(sig_gr(:,lr,k),sig_gr(:,lr,k+1),hanning(nfft,'periodic'),[],nfft,fs);
%C=boxsmooth(C,smoothwidth);
figure(70)
semilogx(fc,C,'b')
axis([xlim,0,1.1])
grid on;hold on;
xlabel('frequency [Hz]')
title('Coherence lr between adjacent grooves')
end

ref_gr=15 % reference groove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(75)
for k=ref_gr:n_rev
C=mscohere(sig_gr(:,lr,ref_gr),sig_gr(:,lr,k),hanning(nfft),[],nfft,fs);
%C=boxsmooth(C,smoothwidth);
semilogx(fc,C,'b')
hold on;grid on;
end
axis([0,fs/2,0,1.1])
xlabel('frequency [Hz]')
title('Coherence lat or vert all ff grooves wrt reference groove')

figure(90)
for k=1:n_rev
C=mscohere(sig_gr(:,1,k),sig_gr(:,2,k),hanning(nfft),[],nfft,fs);
% C=boxsmooth(C,smoothwidth);
plot(fc,C)
axis([0,fs/2,0,1.1])
grid on;hold on;
xlabel('frequency [Hz]')
title('Coherence between channels all grooves')
end

%-----------------spectra------------------
SIG=fft(sig(:,lr));
f=[0:Nt/2]'*fs/Nt;
figure(110)
% SIG=boxsmooth(abs(SIG(1:floor(Nt/2+1))),35);
semilogx(f,20*log10(SIG),'b');
grid on;
axis([1 fs/2 ylim])
xlabel('freq[Hz]');
ylabel('[dB]');
title('Spectrum all data filtered')

figure(120)
semilogx(f,20*log10(SIG.*f(1:floor(Nt/2+1)))-60,'b');
grid on;
axis([1 fs/2 ylim])
xlabel('freq[Hz]');
ylabel('[dB]');
title('Stylus amplitude spectrum filtered')

% disp('-------------------finished--------------------') 
