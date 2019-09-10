% Record Surface Noise analysis
% using xcorrelation between groove signals
% John Vanderkooy
% Jan. 2019
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
%filename='ChirpSilentLeadin.wav';
%filename='leadin35sec.wav';
filename='leadin35sec_declick.wav';
[sig,fs]=audioread(filename);
ns=round(35.0*fs);%this is where we set the analysis portion
nf=round(72.0*fs);
% ns=round(10.0*fs);%this is where we set the analysis portion
% nf=round(30.0*fs);
sig=sig(ns:nf,:);
lr=1;%set to 1 for left channel, 2 for right channel
Nt=length(sig);
disp(['fs: ' num2str(fs) '  Nt: ' num2str(Nt) '  duration: ' num2str(Nt/fs)])
t=linspace(0,(Nt-1)/fs,Nt)';%column vector
%-----------------------data plots------------
figure(10)
plot(t,sig(:,1),'b') 
grid on;hold on;
plot(t,sig(:,2),'r')
axis([0 Nt/fs -.1 .1])
xlabel('time[s]')
legend('left','right')
title('All data')
%------------------generate separate groove results-----------
% we'll make n_rev files
n_per_rev=round(fs*60/33.3333)
n_rev=floor(Nt/n_per_rev) % number of grooves

% choose the data to plot
for ng=1:n_rev;% groove number
    rev4ng(:,:,ng)=sig(1+(ng-1)*n_per_rev:ng*n_per_rev,:);
end
%---------------correlation between grooves--------------
for ng=1:n_rev %whole grooves
    cf(ng)=sum(rev4ng(:,lr,10).*rev4ng(:,lr,ng),1);% correlation over samples only
end
figure(20)
plot(cf/cf(1),'b')
grid on;hold on
% analyze just a portion of each groove
ns=round(1.0*fs);%this is where we set the local analysis portion
nf=round(1.2*fs);
% ns=round(1.2*fs);
% nf=round(1.8*fs);
% ns=round(0.4*fs);
% nf=round(0.9*fs);
% ns=round(0.8*fs);
% nf=round(1.2*fs);
for ng=1:n_rev
    cp(ng)=sum(rev4ng(ns:nf,lr,10).*rev4ng(ns:nf,lr,ng),1);
end
plot(cp/cp(1),'r')
xlabel('groove number')
ylabel('normalized correlation coefficient')
legend('whole revolution','selected portion')
title('Correlation decay')
grid on
%---------time variable in plots------------
tg=t(1:n_per_rev);
minplot=-0.1;
maxplot=+0.1;
%-----------------plot signal for each groove-----------------------------
figure(30)
grid on;hold on;
for np=1:n_rev
plot(tg,rev4ng(:,lr,np),'b');% variables are sample#, left/right, groove#
end
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
title('all grooves')

%-----------------plot signal for each groove-----------------------------
figure(35)
plot(tg,rev4ng(:,lr,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,lr,2),'r'); % variables are sample#, left/right, groove#
plot(tg,rev4ng(:,lr,3),'g') 
plot(tg,rev4ng(:,lr,4),'k') 
plot(tg,rev4ng(:,lr,5),'c') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('groove1','groove2','groove3','groove4','groove5')
title('first grooves')
%-----------------plot signal for each groove-----------------------------
figure(40)
plot(tg,rev4ng(:,1,1),'b') 
grid on;hold on;
plot(tg,rev4ng(:,2,1),'r') 
axis([0 n_per_rev/fs minplot maxplot])
xlabel('time[s]')
legend('groove1 left','groove1 right')
title('stereo groove')
% 
% %---------------cross correlation properties-----------
% ex1=rev4ng(:,1,1); % groove1 left
% ex2=rev4ng(:,1,2); % groove2 left
% wy1=rev4ng(:,2,1); % groove1 right
% wy2=rev4ng(:,2,2); % groove2 right
% 
% Cxy=xcorr(ex1,wy1);
% Cxx=xcorr(ex1);
% Cyy=xcorr(wy1);
% RHOxy=Cxy./(sqrt(Cxx(n_per_rev+1)).*Cyy(n_per_rev+1));
% 
% tau=t(1:2*n_per_rev-1);
% figure(50)
% plot(tau-tau(n_per_rev+1),RHOxy,'b') 
% grid on
% %axis([0 n_per_rev/fs minplot maxplot])
% xlabel('time[s]')
% legend('correlation coefficient')
% title('Declicked stereo groove')
% %----------------------coherence-----------------
% nfft=1024;
% fc=([0:nfft/2])*fs/nfft;
% 
% % rev4ng(:,1,1)=rand(n_per_rev,1);%testing coherence
% % rev4ng(:,1,4)=rand(n_per_rev,1);
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,1,4),hanning(nfft),[],nfft,fs);
% %C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),[],[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(60)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 L groove4')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,1,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(70)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence hann L groove1 & L groove2')
% 
% C=mscohere(rev4ng(:,1,1),rev4ng(:,2,1),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(80)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence L groove1 & R groove1')
% 
% C=mscohere(rev4ng(:,2,1),rev4ng(:,2,2),hanning(nfft),[],nfft,fs);
% Cs=boxsmooth(C,2);
% figure(90)
% plot(fc,C,'b')
% axis([0,fs/2,0,1.1])
% grid on;
% xlabel('frequency [Hz]')
% title('Coherence R groove1 & R groove2')

% disp('-------------------finished--------------------') 
