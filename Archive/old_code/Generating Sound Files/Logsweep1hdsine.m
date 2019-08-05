% Logsweep1hd.m *** Logsweep to get TF, h(t) & harmonic spectra
% the left channel (ref) can be looped to determine the record latency,
% while the right channel (dut) measures the system.
% for most sound cards, the player is delayed wrt the recorder.
% the reference, if not detected, will be the sweep file (sweep).
%
% sweep is padded with zeros, giving 1/4 length clearance at sweep end. 
% the excitation sweep is windowed at both LF and HF (windosweep).
% TF is obtained from DUT/SWEEP, using unwindowed sweep.
% this allows sensible behaviour at LF and HF.
% TF2 is calculated from DUT/REF, nicely normalized but noisy near Nyquist.
% John Vanderkooy, Sept 2016, Audio Research Group, University of Waterloo
%% 0-----------------preliminaries and settings---------------------------%
clear all; clc; close all;
try
  pkg load signal;% for Octave signal processing toolbox
catch
end
disp('----------------start of program--------------------')
set(0,'DefaultLineLinewidth',1.5)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
Po=2e-5;% SPL reference pressure
c=343;% speed of sound
dBclearance=10;% white space above TF plots
dBspan=80;% total scale of TF plots
%----------------enter important parameters-----------------
sig_frac=0.5;% fraction of full scale
fs=44100;% must ensure that Windows settings are the same!
N=2^21;%19;% make this larger if there is insufficient time clearance
% S_dac=-1.22;% UCA202
% S_adc=-1.40;% UCA202 (gain=3 Windows 7). This avoids ADC overload
% S_dac=0.98;% UA-1ex
% S_adc=0.5;% UA-1ex 
S_dac=-1.59;% Focusrite 2i2 Yes, it inverts its monitor output!
S_adc=+1.16;% Focusrite 2i2 line input gain @ 12:00 o'clock
% S_dac=+1.148;% USB Dual Pre peak volts out for digital Full Scale 
% S_adc=+1.49;% USB Dual Pre JV pot minimum (gain=3, Windows 7)
%---------------system calibration factor CF for dut--------------------
system_type=0;% 1:acoustic, 0:electronic, -1:level-dependent, -2:custom
power_amp_gain=1;% V/V acoustic
mic_cal=0.012;% V/Pa acoustic
mic_preamp_gain=10.0;% V/V acoustic
Vw=1;%2.83;% V/watt8ohms. acoustic Leave this=1 if not wanted
electronic_gain=1;% total gain in series with electronic system
if system_type==1 % acoustic
  CF=Vw/(power_amp_gain*mic_cal*mic_preamp_gain*Po);% SPL calibration f===actor
elseif system_type==0 % electronic system
  CF=1/electronic_gain;% electronic
elseif system_type==-1
  CF=1/(sig_frac*mic_cal*mic_preamp_gain*Po);% SPL level-dependent
else %custom
  CF=1;  
end 
%------------- system parameters ----------------------------
Npad=N/4;% this is the total zeropad, added to end of play file
Ns=N-Npad;% most of array is used for sweep
disp(['fs: ' num2str(fs) '  N: ' num2str(N) '  duration: ' num2str(N/fs)])
disp(['signal fraction: ' num2str(sig_frac) '  time clearance: ' num2str(N/(4*fs))])
t=linspace(0,(N-1)/fs,N)';% column vector
ts=linspace(0,(Ns-1)/fs,Ns)';% to calculate sweep

%% 1-----------------------calculate logsweep-----------------------------%
f_start=10;% beginning of turnon half-Hann
f1=20;% end of turnon half-hann
%f2=0.91*fs/2;% beginning of turnoff half-Hann
f2=20000;
f_stop=fs/2;% end of turnoff half-Hann
Ts=Ns/fs;% sweep duration.  This is N-Npad samples in length.
Ls=Ts/log(f_stop/f_start);% time for frequency to increase by factor e
sweep=zeros(N,1);% initialize

%sweep(1:Ns)=sin(2*pi*f_start*Ls*(exp(ts/Ls)-1));% logsweep
sweep(1:Ns)=sin(2*pi*3150*ts);

%------------------tapered sweep window------------------------------
indexf1=round(fs*Ls*log(f1/f_start))+1;% end of starting taper
indexf2=round(fs*Ls*log(f2/f_start))+1;% beginning of ending taper
windo=ones(N,1);
windo(1:indexf1)=0.5*(1-cos(pi*(1:indexf1)/indexf1));% pre-taper
windo(indexf2+1:Ns)=0.5*(1+cos(pi*(1:Ns-indexf2)/(Ns-indexf2)))';% post-taper
windo(Ns+1:N)=0; % zeropad end of sweep
windosweep=windo.*sweep;% tapered at each end for output to DAC

figure(10)
plot(t,sweep,'b')
grid on;hold on
plot(t,windosweep,'r')
xlabel('time[s]')
axis([0 N/fs -1.5 1.5])
legend('raw sweep','windowed sweep')
title('Sweep')
disp(['f_start*Ls: ',num2str(f_start*Ls) '  Ls: ' num2str(Ls)])
disp('finished sweep generation...')
%% 2------------data gathering: send out sweep, record system output------%
% y=sig_frac*[windosweep -windosweep];% antiphase avoids codec midtap modulation
y = sig_frac*[sweep sweep]; %%can't have that 
% ysine315 = sig_frac*[sine315]
% wavwrite(y, fs, 16, "sweep21s.wav")

%wavwrite(y, fs, 16, "sine315.wav");
wavwrite(y, fs, 32, "sine315.wav");
% wavwrite(y, fs, 16, "sine315.wav");




rec = audiorecorder(fs,16,2);
% ply = audioplayer(y,fs,16);
% play(ply); %this should play and allow immediate further execution
recordblocking(rec,N/fs);% program waits until all samples recorded
z = getaudiodata(rec,'double');% this retrieves the recorded data
disp('finished recording...');

disp('Some sound cards act strangely. Check carefully!')
if max(abs(z(:,1)))>0.999
   disp('Check for ch0 REF record saturation!');
end
if max(abs(z(:,2)))>0.999
   disp('Check for ch1 DUT record saturation!');
end
%----------------------- end of data gathering ---------------------
%load('Logsweep1data.mat');
ref=z(1:N,1);%if we use left channel as reference
dut=z(1:N,2);
N=length(dut);
%save('Logsweep1data.mat','sweep','z','N','fs','sig_frac','CF');
clear y z;

disp(['refrmsLSBs: ' num2str(sqrt(2^30*sum(ref.^2)/N))])
disp(['dutrmsLSBs: ' num2str(sqrt(2^30*sum(dut.^2)/N))])

figure(20);
plot(t,ref,'b');
grid on;hold on
plot(t,dut,'r');
axis([0 N/fs -1.5 1.5])
legend('raw ref','raw dut');
xlabel('Time [s]')
title('Recorded responses');

%% 3----------determine record/play delay using crosscorrelation----------%
lags=N/2;% large enough to catch most delays
if max(ref) < 0.1*max(dut) % automatic reference selection
    X=xcorr(sweep,dut,lags);% in case reference is low, use data itself
else
    X=xcorr(sweep,ref,lags);% this uses recorded reference
end
[~,nmax]=max(abs(X));

figure(30)
tl=t-lags/fs;% centred plot
plot(tl(1:N),X(1:N))
grid on
legend('crosscorrelation')
xlabel('time [s]')
title('recorder leads player-----|------recorder lags player')

offset=nmax-lags-1;
disp(['record offset: ' num2str(offset) '  time: ' num2str(offset/fs)])
if abs(offset) > Npad 
    disp('******INSUFFICIENT TIME CLEARANCE!******')
    disp('******INSUFFICIENT TIME CLEARANCE!******')
    disp('******INSUFFICIENT TIME CLEARANCE!******')
end
disp('negative offset means player lags recorder!')

lwindo=ones(N,1);
lwindo(1:indexf1)=0.5*(1-cos(pi*(1:indexf1)/indexf1));% LF pre-taper
lwindosweep=lwindo.*sweep;
% remove play-record delay by shifting computer sweep array
%sweep=circshift(sweep,-offset);
lwindosweep=circshift(lwindosweep,-offset);
%% 4-----------calculate TFs using Frequency Domain ratios----------------%
% all frequency variables are meant to be voltage spectra
%SWEEP=sig_frac*S_dac*fft(sweep);
LWINDOSWEEP=sig_frac*S_dac*fft(lwindosweep);
REF=S_adc*fft(ref);
DUT=CF*S_adc*fft(dut);

TF=DUT./LWINDOSWEEP;% this has good Nyquist behaviour
%TF=DUT./SWEEP;
TF2=DUT./REF;% this has good delay and normalization
ft=((1:N/2+1)'-1)*fs/N;
ft(1)=NaN;% this avoids non-positive data in semilogx warning from Octave

top=20*log10(max(abs(TF(2:N/2+1))));% highest dB in plot
figure(100);
semilogx(ft,20*log10(abs(TF(1:N/2+1))),'b');
grid on;
axis([fs/N fs/2 top+dBclearance-dBspan top+dBclearance])
legend('DUT/SWP')
xlabel('frequency [Hz]');ylabel('dB')
title('Transfer Function');

top=20*log10(max(abs(TF2(2:N/2+1))));% highest dB in plot
figure(110);
semilogx(ft,20*log10(abs(TF2(1:N/2+1))),'b');
grid on;
axis([fs/N fs/2 top+dBclearance-dBspan top+dBclearance])
legend('DUT/REF')
xlabel('frequency [Hz]');ylabel('dB')
title('Transfer Function 2');
%% 5---------------load previous file and compare TF ---------------------%
try
    load('old_TF.mat')% loads old_TF
    % this figure will appear when file exists
    top=20*log10(max(abs(TF(2:N/2+1))));% highest dB in plot
    figure(120);
    semilogx(ft,20*log10(abs(TF(1:N/2+1))),'b');
    grid on;hold on
    if length(old_TF) == length(TF)
     semilogx(ft,20*log10(abs(old_TF(1:N/2+1))),'r');
    else
     disp('Different file lengths: do another run')   
    end %if length
    axis([1 fs/2 top+dBclearance-dBspan top+dBclearance])
    legend('new DUT/SWP','old DUT/SWP')
    xlabel('frequency [Hz]');ylabel('dB')
    title('old & new Transfer Functions');
catch
    disp('continuing...')
end %try
old_TF=TF;
save('old_TF.mat','old_TF');% saves current TF file
% all the programs are the same to this point--------------------
%% 6----------------obtain impulse response ------------------------------%
h=ifft(TF);% this should be real
impulse=h;
figure(300);% has spike at t near zero, so acausal bits go to negative time
plot(t,h,'b')
grid on;
axis([-.05*N/fs 1.05*N/fs min(h)-0.05*max(h) 1.2*max(h)])
xlabel('time [s]')
title('h(t) from DUT/SWP');

hpre=circshift(h,floor(0.01*fs));% 10 ms precursor view
figure(310);% zoomed, shifted to show precursor
plot(t,hpre,'b')
grid on;
axis([0 .05 min(h)-0.05*max(h) 1.2*max(h)]) % show 50ms for editing
xlabel('time [s]')
title('zoomed h(t) with precursor');

%% 7-------select the various h(t)s for the fundamental and harmonics-----%
Tpre=abs(Ls/50);Npre=round(Tpre*fs);%wraparound tolerance
T2=abs(Ls)*log(2);T2pre=round(T2*fs);% harmonic positions in bins
T3=abs(Ls)*log(3);T3pre=round(T3*fs);
h1=zeros(N,1);h2=zeros(N,1);h3=zeros(N,1);%initialize

% h1 is already at origin
h1(1:N/2+1)=impulse(1:N/2+1);
h1(N-Npre:N)=impulse(N-Npre:N);%wraparound

h=circshift(impulse,round(T2pre));% sets h2 to origin
h2(1:T2pre-Npre)=h(1:T2pre-Npre);
h2(N-Npre:N)=h(N-Npre:N);%wraparound

h=circshift(impulse,round(T3pre));% sets h3 to origin
h3(1:T3pre-T2pre-Npre)=h(1:T3pre-T2pre-Npre);
h3(N-Npre:N)=h(N-Npre:N);%wraparound

h1=circshift(h1,Npre);% use Npre to display acausal segment
h2=circshift(h2,Npre);
h3=circshift(h3,Npre);

% figure(320);
% plot(t,h1(1:N),'b')
% grid on;hold on
% plot(t,h2(1:N),'r')
% plot(t,h3(1:N),'g')
% legend('h1','h2','h3','Location','Best')
% xlabel('time [s]')
% title('DUT/SWP Lin&dist IRs t-corr, t-shift');
% 
% figure(325);
% plot(t-Tpre,h2(1:N),'b')
% grid on;
% legend('h2','Location','Best')
% axis([-.002 0.006 -4e-4 4e-4])
% axis([-.001*N/fs 0.01*N/fs -1.2*max(h2) 1.2*max(h2)])
% xlabel('time [s]')
% title('2nd harmonic IR time corrected');

figure(330);
plot(t,h1(1:N),'b')
grid on;hold on
plot(t,h2(1:N),'r')
plot(t,h3(1:N),'g')
legend('h1','h2','h3','Location','Best')
xlabel('time [s]')
title('DUT/SWP Lin & dist IRs time corrected');

%% 8-------------------plot TF magnitude and phases-----------------------%
h1=circshift(h1,-Npre);% unshift h's
h2=circshift(h2,-Npre);
h3=circshift(h3,-Npre);

TFh1=fft(h1);
TFh2=fft(h2);
TFh3=fft(h3);

oct_width=0.1;% do some smoothing.  This takes time in Octave
TFh1=pwroctsmooth(TFh1,oct_width);
TFh2=pwroctsmooth(TFh2,oct_width);
TFh3=pwroctsmooth(TFh3,oct_width);

figure(400);
top=20*log10(max(abs(TF(2:N/2+1))));% highest dB in plot
semilogx(ft,20*log10(abs(TFh1(1:N/2+1))),'b');
grid on;hold on
axis([10*fs/N fs/2 top+dBclearance-dBspan-40 top+dBclearance])
semilogx(ft/2,20*log10(abs(TFh2(1:N/2+1))),'r');% note scale ft/2
semilogx(ft/3,20*log10(abs(TFh3(1:N/2+1))),'g');% note scale ft/3
xlabel('frequency [Hz]')
ylabel('SPL [dB]')
legend('fundamental','2nd harm','3rd harm','Location','Best')
title('DUT/SWP Transfer Functions');
disp('-------------------finished--------------------') 
