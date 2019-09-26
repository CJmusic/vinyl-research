function[smoothed_array]=boxsmooth(array,box_halfwidth)
% function[smoothed_array]=boxsmooth(array,box_halfwidth)
%
% box smoothing for arrays, real or complex
% assumes mirror each end boundary condition
% box_halfwidth is half of total sample spread
% JV January 2016
N=length(array); % should work even or odd, complex or real
box_halfwidth=ceil(abs(box_halfwidth));% force integer
smoothed_array=zeros(size(array));% keeps dimensions
sum=array(1);% start calculation of smoothed_array(1)
for J=2:box_halfwidth+1 % set up smoothed_array(1)
    %sum = sum + array(J)+array(1);% flat end condition
    sum = sum + 2*array(J);% mirror end condition
end
  smoothed_array(1)=sum/(2*box_halfwidth+1);

  for K=2:N; % fast smoothing, subtract array(I-box-1) add array(I+box)
    if K > box_halfwidth+1 % normal index
        index=K-box_halfwidth-1;% in midst of array, min=1
    else
        %index=1;% flat end condition
         index=3-K+box_halfwidth;% 3 is correct!
    end
        s=array(index);% subtract from sum
    if K <= N-box_halfwidth % normal index
        index=K+box_halfwidth;% in midst of array, max=N
    else
        %index=N;% flat end condition
         index=2*N-K-box_halfwidth;% looks OK
    end
        a=array(index);% add to sum
    sum = sum - s + a;
  smoothed_array(K) = sum/(2*box_halfwidth+1);
end
