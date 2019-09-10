function[smoothed_power]=poweroctsmooth(freq_nyq,octave_width)
% fractional-octave power smoothing Max width=2.0_
% freq_nyquist is complex single sided spectrum, length N/2+1
% smoothed_power is returned
if octave_width>2.0 octave_width=0.33;end
an=2^(octave_width/2);

N=2*length(freq_nyq)-2;
smoothed_power=zeros(N/2+1,1);
old_lo=0;old_hi=0;
for I=1:N/2+1 % fast smoothing algorithm
   lo_bin=round((I-1)/an)+1;hi_bin=round((I-1)*an)+1;
   if hi_bin==lo_bin 
      power_sum=(abs(freq_nyq(I)))^2; % no change
   else   
    for J=old_lo:(lo_bin-1)
      power_sum = power_sum - (abs(freq_nyq(J)))^2;
    end   
    for K=(old_hi+1):hi_bin
        if K > (N/2+1)
            J=N+2-K;
        else
            J=K;
        end    
      power_sum = power_sum + (abs(freq_nyq(J)))^2;
    end
   end
   smoothed_power(I)=power_sum/(hi_bin-lo_bin+1);
   old_lo=lo_bin;old_hi=hi_bin;
end 
