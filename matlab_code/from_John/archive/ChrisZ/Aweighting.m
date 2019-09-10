%A-weighting digital filter
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
fs=96000; 
N=2^14;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(1)=1;%unit impulse

figure(10)
plot(time,signal,'b');
grid on;
title('input signal')
xlabel('Time [s]')

%Analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);

%Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);

%-------------------------------------------
output=filter(b,a,signal);

figure(20)
plot(time,output,'b');
grid on;
title('A-wtd impulse signal')
xlabel('Time [s]')

figure(30)
plot(time,output,'b');
grid on;
axis([0 .01 -1 1]);
title('A-wtd impulse signal')
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
legend('A-wtg response','Location','Best')
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
title('magnitude response');

figure(70);
semilogx(f,(180/pi)*angle(OUTPUT(1:floor(N/2+1))),'b');
grid on;
axis([fs/N,fs/2,-180,180])
xlabel('Frequency [Hz]')
ylabel('phase')
title('net phase response');

disp('-----------------------------finished------------------------')


