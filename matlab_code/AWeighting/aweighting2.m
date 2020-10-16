% A-weighting digital filter
clear all;close all;clc
set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','bold')
set(0,'DefaultAxesLineWidth',1.5)
% dBfs = 0;%-34.0%-6.0


% sig_frac=0.5;% fraction of full scale
% sig_frac = 10^(dBfs/20.0)
sig_frac = 1;


fs=96000;% must ensure that Windows settings are the same!
N=2^21;%19;% make this larger if there is insufficient time clearance

% S_dac=-1.59;% Focusrite 2i2 Yes, it inverts its monitor output!
% S_adc=+1.16;% Focusrite 2i2 line input gain @ 12:00 o'clock

Npad=N/4;% this is the total zeropad, added to end of play file
Ns=N-Npad;% most of array is used for sweep

t=linspace(0,(N-1)/fs,N)';% column vector
ts=linspace(0,(Ns-1)/fs,Ns)';% to calculate sweep

%-----------------------calculate logsweep-----------------------------%
f_start=10;% beginning of turnon half-Hann
f1=20;% end of turnon half-hann
%f2=0.91*fs/2;% beginning of turnoff half-Hann
f2=16000%20000;
f_stop=fs/2;% end of turnoff half-Hann
Ts=Ns/fs;% sweep duration.  This is N-Npad samples in length.
Ls=Ts/log(f_stop/f_start);% time for frequency to increase by factor e
sweep=zeros(N,1);% initialize
sweep(1:Ns)=sin(2*pi*f_start*Ls*(exp(ts/Ls)-1));% logsweep
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

data = y;


fs = 96000;
figure(10);
plot(data)
grid on;


[data_A]=Aweighting_filterTest(data,96000);
% data_A = audio_AweightingTest(data);


% function[weighted_array]=Aweighting_filter(array,fs)
function[weighted_array]=Aweighting_filterTest(array,fs)
    % function[weighted_array]=Aweighting_filter(array,fs)
    %----------------------------------------------------
    % A weighting filter for 44.1 or 96 kHz data only, otherwise bypass.
    % The response above 16 kHz is not specified so a smooth rolloff is chosen.
    % ---------------------------------------------------
    if fs==44100
    %------------------44.1 jule walker implementation-----------------
    fr=[0 6.3 10 16 25 40 63 100 160 250 400 630 1000 1600 2500 4000 6300 10000 16000 20000 22050];
    dBAtable=[-inf -85.4 -70.4 -56.7 -44.7 -34.6 -26.2 -19.1 -13.4 -8.6 -4.8 -1.9 0 1.0 1.3 1.0 -0.1 -2.5 -6.6 -9.3 -10];
    Wn=2*fr/fs;
    Amag=10.^(dBAtable/20);
    [b,a]=yulewalk(10,Wn,Amag);
    [d,c]=butter(1,2*240/fs,'high');% this corrects DC-LF with highpass
    [k,h]=butter(1,2*150/fs,'high');% this corrects DC-LF with highpass
    [m,l]=butter(1,2*40/fs,'high');% this corrects DC-LF with highpass
    %------------------compactify into single filter----------------
    bk=conv(b,d);ah=conv(a,c);% YW+HP1
    bm=conv(bk,k);al=conv(ah,h);% YW1+HP2
    bo=conv(bm,m);an=conv(al,l);% YW2+HP3 final coefficients
    weighted_array=filter(bo,an,array);
    %--------------------------------------------------------------------
    elseif fs==96000
    %------------------96kHz jule walker implementation-----------------
    fr=[0 6.3 10 16 25 40 63 100 160 250 400 630 1000 1600 2500 4000 6300 10000 16000 20000 25000 35000 48000];
    dBAtable=[-inf -85.4 -70.4 -56.7 -44.7 -34.6 -26.2 -19.1 -13.4 -8.6 -4.8 -1.9 0 1.0 1.3 1.0 -0.1 -2.5 -6.6 -9.3 -15 -25 -40];
    Wn=2*fr/fs;
    Amag=10.^(dBAtable/20);
    [b,a]=yulewalk(12,Wn,Amag);
    [d,c]=butter(1,2*280/fs,'high');% this corrects DC-LF with highpass
    [k,h]=butter(1,2*200/fs,'high');% this corrects DC-LF with highpass
    [m,l]=butter(1,2*40/fs,'high');% this corrects DC-LF with highpass
    %------------------compactify into single filter----------------
    bk=conv(b,d);ah=conv(a,c);% YW+HP1
    bm=conv(bk,k);al=conv(ah,h);% YW1+HP2
    bo=conv(bm,m);an=conv(al,l);% YW2+HP3 final coefficients
    weighted_array=filter(bo,an,array);
    
    else
        disp('%%%%%%%%%%%%%%%% Wrong fs, filter bypassed! %%%%%%%%%%%%%%%%')
        weighted_array=array;
    end
    
    
end
    