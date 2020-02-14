

function data_CCIR = audio_CCIRweighting(data)
   
   % Po=2e-5;% SPL ref
   % -----------------------------------------------------------------
   N=2^14;
   %----------------------CCIR/ARM dB table-------------------------------
   fr=[31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000];
   CCIR=[-35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8];
   %------------------this 44.1kHz fit needs optimization---------------------
   fs=96000;
   time=([1:N]'-1)/fs; % make time positive column vector starting at zero
   signal=zeros(N,1);
   impulse=zeros(N,1);
   impulse(1)=1;%unit impulse
   %--------
   Wn=2*(8400)/fs;
   % fair design Wn=8700, 2nd LP, 1st HP factor 1.05*4.7544
   % good 44.1 peak design Wn=8400, 2 2nd LP W=1.1Wn, 2 (1st LP W=1.2Wn ?), 1st HP W=0.8Wn, factor=0.77*4.7544
   %-------3 2nd-order LP + one 1st-order LP-------------
   [b,a]=butter(2,1.3*Wn);
   signal=impulse;
   for k=1:1
      signal=filter(b,a,signal);
   end
   [b,a]=butter(1,1.4*Wn);
   for k=1:2
      signal=filter(b,a,signal);
   end
   %---------1st order highpass---------
   [b,a]=butter(1,0.8*Wn,'high');
   signal=filter(b,a,signal);
   multfactor=0.77*4.7544;% to bring gain=1 at 2kHz
   signal=multfactor*signal;
   data_CCIR = signal;
end