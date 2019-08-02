% RIAA digital filter
% John Vanderkooy Feb 2019
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)

Po=2e-5;% SPL ref
try
    pkg load signal; %for Octave
catch
end
%% -----------------------------------------------------------------
fs=44100; 
N=2^15;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(1)=1;%unit impulse

figure(10)
plot(time,signal,'b');
grid on;
title('input signal')
xlabel('Time [s]')

% RIAA filter according to specs
% the bilinear transformation leaves the response poor above approx 3/4 fs
f1 = 50; 
f2 = 500;
f3 = 2122;
NUM = [1 2*pi*f2];%reversed coefficients for inverted fraction 
DEN = conv([1 2*pi*f1],[1 2*pi*f3]); 
%Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);

%-------------------------------------------
G=2*pi*f1*f3/f2;
G2=10^(19.9183/20);
output=G*G2*filter(b,a,signal);

figure(20)
plot(time,output,'b');
grid on;
title('RIAA impulse response')
xlabel('Time [s]')

figure(30)
plot(time,output,'b');
grid on;
axis([0 .01 ylim]);
title('RIAA impulse signal')
xlabel('Time [s]')

%---------------frequency domain----------------------
%fftfactor=sqrt(2/(N*fs)); %Pa/root_Hz single sided
f=([1:N/2+1]'-1)*fs/N;
%OUTPUT=fftfactor*fft(output);
OUTPUT=fft(output);

figure(60);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;
pk=10*ceil(max(20*log10(abs(OUTPUT)))/10);
axis([fs/N,fs/2,pk-60,pk])
legend('RIAA response','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('magnitude response');

% figure(70);
% semilogx(f,(180/pi)*angle(OUTPUT(1:floor(N/2+1))),'b');
% grid on;
% axis([fs/N,fs/2,-180,180])
% xlabel('Frequency [Hz]')
% ylabel('phase')
% title('net phase response');

disp('-----------------------------finished------------------------')


