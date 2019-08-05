dBfs = -2.0e32%-14.0
time_len = 60.0 
tone_f = 1000.0

fs=44100;% must ensure that Windows settings are the same!
N = (time_len*fs)/0.9;%/0.75;
sig_frac = 10^(dBfs/20.0);

% N=2^21;%19;% make this larger if there is insufficient time clearance

S_dac=-1.59;% Focusrite 2i2 Yes, it inverts its monitor output!
S_adc=+1.16;% Focusrite 2i2 line input gain @ 12:00 o'clock

Npad=N/10;% this is the total zeropad, added to end of play file
Ns=N-Npad;% most of array is used for sweep

t=linspace(0,(N-1)/fs,N)';% column vector
ts=linspace(0,(Ns-1)/fs,Ns)';% to calculate sweep

%-----------------------calculate logsweep-----------------------------%
f_start=10;% beginning of turnon half-Hann
f1=20;% end of turnon half-hann
%f2=0.91*fs/2;% beginning of turnoff half-Hann
f2=20000;
f_stop=fs/2;% end of turnoff half-Hann
Ts=Ns/fs;% sweep duration.  This is N-Npad samples in length.
Ls=Ts/log(f_stop/f_start);% time for frequency to increase by factor e
sweep=zeros(N,1);% initialize

%sweep(1:Ns)=sin(2*pi*f_start*Ls*(exp(ts/Ls)-1));% logsweep
sweep(1:Ns)=sin(2*pi*tone_f*ts);
%-----------------------------------------------------------------------%

%------------------tapered sweep window------------------------------
indexf1=round(fs*Ls*log(f1/f_start))+1;% end of starting taper
indexf2=round(fs*Ls*log(f2/f_start))+1;% beginning of ending taper
windo=ones(N,1);
windo(1:indexf1)=0.5*(1-cos(pi*(1:indexf1)/indexf1));% pre-taper
windo(indexf2+1:Ns)=0.5*(1+cos(pi*(1:Ns-indexf2)/(Ns-indexf2)))';% post-taper
windo(Ns+1:N)=0; % zeropad end of sweep
windosweep=windo.*sweep;% tapered at each end for output to DAC
%-----------------------------------------------------------------------%

y=sig_frac*[windosweep windosweep];% in phase for vinyl pressing

fname = sprintf ( 'toneMIN%i.wav', tone_f );
wavwrite(y, fs, 32, fname)
