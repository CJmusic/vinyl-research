function[smoothed_tf]=pwroctsmooth_singlesided(freq_tf,octave_width)
% function[smoothed_tf]=pwroctsmooth(freq_tf,octave_width)
%
% fractional-octave pwr preserving smoothing. Max width<2.0 octaves
% freq_tf is complex single- sided spectrum, length N, end @ Nyquist.
% smoothed_tf is returned, real single-sided, size N.
% We assume that the double-sided array is even, so Nyquist is N.
%
% John Vanderkooy Oct, 2021
if octave_width>=2.0 
   octave_width=0.33;
   disp(['oct_width set to ' num2str(octave_width)])
end
an=2^(octave_width/2);
N=length(freq_tf); % should work for even or odd N
% Nyquist will be at Matlab index N, DC is index 1.
smoothed_tf=zeros(size(freq_tf));% if complex causes long execution?
old_lo=0;old_hi=0;
for I=1:N % fast smoothing algorithm
   lo_bin=round((I-1)/an)+1;
   hi_bin=round((I-1)*an)+1;
   if hi_bin==lo_bin 
      pwr_sum=abs(freq_tf(I))^2;% no change
   else   
      for J=old_lo:(lo_bin-1)
        pwr_sum = pwr_sum - abs(freq_tf(J))^2;%forget these
      end   
      for K=(old_hi+1):hi_bin
        index=K+(K>N)*2*(N-K);%makes it 2N-K if K>N
        pwr_sum = pwr_sum + abs(freq_tf(index)^2);%add these
      end
   end
    smoothed_tf(I) = sqrt(pwr_sum/(hi_bin-lo_bin+1));
    old_lo=lo_bin;old_hi=hi_bin;
end

