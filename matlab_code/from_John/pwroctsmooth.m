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
