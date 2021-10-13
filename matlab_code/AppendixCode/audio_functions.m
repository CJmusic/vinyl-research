
function data_A = audio_Aweighting(data)
        fs = 96000;
        f1 = 20.598997; 
        f2 = 107.65265;
        f3 = 737.86223;
        f4 = 12194.217;
        f4 = 14100;
        A1000 = 1.9997;
        NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
        DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
        DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
        % % Bilinear transformation of analog design to get the digital filter. 
        [b,a] = bilinear(NUM,DEN,fs);
        data_A = filter(b, a, data);
    end


function data_CCIR = audio_CCIRweighting(data)
       fs=96000;
       frdc=[0 31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10000 12500 14000 16000 20000 25000 30000 fs/2];
       CCIR=[-inf -35.5 -29.5 -25.4 -19.4 -13.4 -7.5 -5.6 0.0 3.4 4.9 6.1 6.6 6.4 5.8 4.5 2.5 -5.6 -10.9 -17.3 -27.8 -35 -50 -inf];
       Wn=2*frdc/fs;
       CCIRmag=10.^(CCIR/20);
       [b,a]=yulewalk(12,Wn,CCIRmag);
       [d,c]=butter(1,2*750/fs,'high');
       fb=conv(b,d);ea=conv(a,c);
       data_CCIR=filter(fb,ea,data);
    end


function [amp_coh, freq_coh] = audio_mscohere(data1, data2, fs)
    nfft=2^14;
    [amp_coh, freq_coh] = mscohere(data1,data2,hanning(nfft),nfft/2,nfft,fs);
end


function [data_fft, freq_fft] = audio_spectrum(data, fs, start_sam, n_sam);
    data_fft = (fft(data(start_sam:start_sam+n_sam-1, :))/n_sam);
    data_fft = data_fft(1:floor(n_sam/2)+1,:);
    data_fft(2:end-1,:) = 2*data_fft(2:end-1,:);
    freq_fft = fs*(0:(n_sam/2))/n_sam;
end 




function [smoothed_tf]=pwroctsmooth(freq_tf,octave_width)
    % function[smoothed_tf]=pwroctsmooth(freq_tf,octave_width)
    % fractional-octave pwr preserving smoothing. Max width<2.0_
    % freq_tf is complex double sided spectrum, length N, conjugate even
    % smoothed_tf (amplitude) is returned, real only
    
    if octave_width>=2.0 
       octave_width=0.33;
       disp(['oct_width set to ' num2str(octave_width)])
    end
    an=2^(octave_width/2);
    
    N=length(freq_tf); % should work for even or odd
    smoothed_tf=zeros(size(freq_tf));% making this complex causes long execution!
    old_lo=0;old_hi=0;
    for I=1:floor(N/2)+1 % fast smoothing algorithm
       lo_bin=round((I-1)/an)+1;hi_bin=round((I-1)*an)+1;
       if hi_bin==lo_bin 
          pwr_sum=abs(freq_tf(I))^2;% no change
       else   
          for J=old_lo:(lo_bin-1)
            pwr_sum = pwr_sum - abs(freq_tf(J))^2;
          end   
          for K=(old_hi+1):hi_bin
            pwr_sum = pwr_sum + abs(freq_tf(K))^2;
          end
       end
        smoothed_tf(I) = sqrt(pwr_sum/(hi_bin-lo_bin+1));
        old_lo=lo_bin;old_hi=hi_bin;
    end
    smoothed_tf(floor(N/2)+2:N)=conj(smoothed_tf(ceil(N/2):-1:2));% ensure conjugate even
    if N/2==floor(N/2) % even N
      smoothed_tf(N/2+1)=abs(smoothed_tf(N/2+1));% Nyquist should be real
    end
    
    % smoothed_tf=complex(smoothed_tf);
    