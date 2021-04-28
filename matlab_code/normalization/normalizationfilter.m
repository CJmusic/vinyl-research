

fs=9600;%decimated down from 96000
N=2^20;disp(['Total duration [s] ' num2str(N/fs)])
time=([1:N]'-1)/fs; % make time positive column vector starting at zero
signal=zeros(N,1);
signal(10)=1;%unit impulse



w1 = 2*707/fs; w2 = 2*1404/fs;
[b,a] = butter(4, [w1 w2]);
% segfilt = filter(b,a,seg);
output=filter(b,a,signal);
OUTPUT=fft(output);
f=[0:N/2]'*fs/N;

figure(400);
semilogx(f,20*log10(abs(OUTPUT(1:floor(N/2+1)))),'b');
grid on;hold on;
