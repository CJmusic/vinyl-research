% addpath('audio_functions')

% function RIAA_eq();

T1 = 318e-6
T2 = 75e-6
T3 = 3180e-6

b = [1+T1,1+T2];
a = [1+T3];
fs = 96000
data = rand(1*fs,1);
data = real(data)/max(data);

RIAA = filter(b, a, data);
RIAA = real(RIAA);
RIAA_f, freq = audio_spectrum(RIAA, fs, 1, length(RIAA)-1);
size(RIAA_f)
size(freq)
plot(freq,RIAA_f) 