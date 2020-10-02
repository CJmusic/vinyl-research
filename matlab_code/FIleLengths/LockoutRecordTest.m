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
offset=-9;%typical value
n_start=0;% start sample of revolution decomposition

%filename='r27a-ableton1024-192kHz.wav';offset=-23;%takes long time...

%filename='lockout-groove10b.wav';offset=-19;

%filename='D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\lockoutgroove\lockout2.wav';offset=-9;
filename='D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\lockoutgroove\lockout1.wav';offset=-9;
filename = '/Users/cz/OneDrive - University of Waterloo/School/Vinyl_Project/audio_files/lockoutgroove/lockout1.wav'; offset = -9;
[sig,fs]=audioread(filename);

lr=1;%---left or right channel
sig=sig(:,lr);%reduce memory space
Nt=length(sig);
%sig=sig(1:round(Nt/2));%for really long files
%sig=sig(round(Nt/2+1):Nt);%for really long files
Nt=length(sig)
disp([filename ' Nt: ' num2str(Nt) ' fs; ' num2str(fs)])
n_per_rev=round(1.8*fs+offset)
n_rev=floor((Nt-n_start)/n_per_rev)

%------------------plot raw file-----------------------
figure(1)
plot(sig)
grid on;
%----------decompose sig into segments----------
seg=zeros(n_per_rev,n_rev);
for k=1:n_rev
   seg(:,k)=sig(1+(k-1)*n_per_rev+n_start:k*n_per_rev+n_start)';
end
delta_n=zeros(n_rev,1);
n=zeros(n_rev,1);
[x,n(1)]=max(seg(:,1));% first lockout pulse
% n(k) are the indices of the clicks
% delta_n(k) are the click separations
for k=2:n_rev
[x,n(k)]=max(seg(:,k));
delta_n(k)=n(k)-n(k-1)+n_per_rev;
end
delta_n(1)=n_per_rev;% let's start with a correct spacing

%t=(0:Nr-1)/fs';%column vector
%-------------------plot all data------------
figure(10)
plot(delta_n-1.8*fs,'b') 
grid on;
axis([xlim ylim])
xlabel('rev number')
ylabel('delta-n offset')
legend('ref file1')
title('segment data')

figure(20)
for k=1:n_rev
plot(seg(:,k)) 
grid on;hold on;
end
axis([xlim ylim])
xlabel('sample number')
ylabel('signal')
legend('ref file1')
title('segment data')
disp('-------------------finished--------------------') 
