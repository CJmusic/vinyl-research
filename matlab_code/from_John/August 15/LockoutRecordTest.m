% Test for USB recording integrity using repeated lockout groove clicks
% then inspection of click separations will indicate sample integrity
% Aug. 2019 John Vanderkooy
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

%large files take time to load, especially first time!!!!

%----------------prepare to sync sig into revolutions-----------
offset=-18;%close estimate for Technics turntable
lr=1;%1=left 2=right channel
n_start=100;% start sample of revolution decomposition

%filename='r27a-ableton64.wav';
%filename='r27a-ableton1024.wav';
%filename='r27a-ableton1024-192kHz.wav';
%filename='r27a-ableton2048.wav';
%filename='r28a-ableton1024.wav';

%filename='shure-ttr103-dual.wav';offset=50;

%filename='r28a-aceraudacity.wav';%has missing sample glitches after reboot!
%filename='r28a-macaudacity.wav';

%filename='lockout-acer2.wav'
%filename='lockout-acer.wav'

filename='lockout-groove10b.wav';

[sig,fs]=audioread(filename);
Nt=length(sig);
disp([filename ' Nt: ' num2str(Nt) ' fs; ' num2str(fs)])
n_per_rev=1.8*fs+offset
n_rev=floor((Nt-n_start)/n_per_rev)
%n_rev=10;

%----------decompose sig into segments----------
seg=zeros(n_per_rev,n_rev);
delta_n=zeros(n_rev,1);
n=zeros(n_rev,1);
for k=1:n_rev
   seg(:,k)=sig(1+(k-1)*n_per_rev+n_start:k*n_per_rev+n_start,lr)';
end
[x,n(1)]=max(seg(:,1));% first lockout pulse
% n(k) are the indices of the clicks
% delta_n(k) are the click separations
for k=2:n_rev
[x,n(k)]=max(seg(:,k));
delta_n(k)=n(k)-n(k-1)+n_per_rev;
end
delta_n(1)=n_per_rev;% let's start with a correct spacing
Nr=n_per_rev;
t=linspace(0,(Nr-1)/fs,Nr)';%column vector

%-------------------plot all data------------
figure(10)
plot(delta_n,'b') 
grid on;
axis([xlim 0 1.2*max(delta_n)])
xlabel('rev number')
ylabel('delta-n')
legend('ref file1')
title('segment data')

figure(20)
for k=1:n_rev
plot(seg(:,k)) 
grid on;hold on;
end
%axis([1.33e5 1.36e5 ylim])
axis([xlim ylim])
xlabel('sample number')
ylabel('signal')
legend('ref file1')
title('segment data')

% disp('-------------------finished--------------------') 
