%A-weighting digital filter

% [data, fs] = audioread('D:\OneDrive - University of Waterloo\School\Vinyl_Project\audio_files\A0000B0000\03141_A0000B0000r030b.wav');
% data = data(1:10*fs);
% data_A = audio_Aweighting(data);

% N = length(data);
% freq=([1:N/2+1]'-1)*fs/N;

% data_A_fft = fft(data_A);
% data_A_fft = abs(data_A_fft(1:floor(N/2+1)));

% data_fft = fft(data);
% data_fft = abs(data_fft(1:floor(N/2+1)));
% figure(1); grid on;
% semilogx(freq,20*log10(data_A_fft),'b');
% grid on;
% xlabel('Frequency [Hz]')
% ylabel('SPL [dB]')
% title('A weighted noise');


% data_fft = fft(data);
% data_fft = abs(data_fft(1:floor(N/2+1)));

% figure(2); 
% semilogx(freq,20*log10(data_fft),'b');
% grid on;
% xlabel('Frequency [Hz]')
% ylabel('SPL [dB]')
% title('Unweighted noise');

function data_A = audio_Aweighting(data)
    fs = 96000;
    f1 = 20.598997; 
    f2 = 107.65265;
    f3 = 737.86223;
    f4 = 12194.217;%for infinite fs
    f4 = 14100;%for finite fs, makes magnitude -9.5dB at 20kHz
    A1000 = 1.9997;
    NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
    DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
    DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);
    % % Bilinear transformation of analog design to get the digital filter. 
    [b,a] = bilinear(NUM,DEN,fs);
    time = (0:length(data)-1)*fs;
    data_A = filter(b, a, data);
end