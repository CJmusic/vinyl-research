

function data_CCIR = audio_CCIRweighting(data)
   
   fs=96000;
   frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
   CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
   Wn=2*frdc/fs;
   CCIRmag=10.^(CCIR/20);
   [b,a]=yulewalk(12,Wn,CCIRmag);%%%%%%%%%%%%%%%%%%%%%
   [d,c]=butter(1,2*750/fs,'high');% this corrects DC-LF with highpass
   fb=conv(b,d);ea=conv(a,c);
   data_CCIR=filter(fb,ea,data);
   %----------------------plot yulewalk freq------------
   % f=[0:N/2]*fs/N;
   % OUTPUT=fft(signal);
   %------------------------------------------------------------------
   
end