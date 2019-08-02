function[filtered_array]=CCIRarm_filter(array,fs)
% [filtered_array]=CCIRarm_filter(array,fs)
% ----------------------------------------------
% This function applies a CCIRarm filter with unity gain at 2kHz.
% If the sampling frequency is not 44.1 or 96 kHz, it will fail.
% A jule-walker design is used, with an extra highpass for LF accuracy

if fs==44100
%------------------yule-walker 44.1kHz filter design----------------
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000  fs/2];
CCIR=[-100 -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -32];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(12,Wn,CCIRmag);
[d,c]=butter(1,2*370/fs,'high');% this corrects DC-LF with highpass
%------------------compactify into single filter----------------
fb=conv(b,d);ea=conv(a,c);
filtered_array=filter(fb,ea,array);

elseif fs==96000
%------------------yule-walker 96kHz filter design----------------
frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
Wn=2*frdc/fs;
CCIRmag=10.^(CCIR/20);
[b,a]=yulewalk(12,Wn,CCIRmag);
[d,c]=butter(1,2*750/fs,'high');% this corrects DC-LF with highpass
%------------------compactify into single filter----------------
fb=conv(b,d);ea=conv(a,c);
filtered_array=filter(fb,ea,array);

else
    disp('Wrong sampling frequency!  Use 44.1 or 96 kHz.')
end