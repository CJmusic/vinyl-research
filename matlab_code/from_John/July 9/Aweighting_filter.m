function[weighted_array]=Aweighting_filter(array,fs)
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



