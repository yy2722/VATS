function [fps,freq] = makeFPS(stimulus,fs)

% Do the FFT
nFFT = 2^nextpow2(length(stimulus));
fft1 = fft(stimulus,nFFT)/length(stimulus);
fft2 = 2*abs(fft1(1:nFFT/2+1));
fft3 = resample(fft2,1,round(50*(nFFT/2^17))); % smooth the fft
fps  = fft3/max(fft3);

% Get freq values
fq   = fs/2*linspace(0,1,nFFT/2+1);
freq = resample(fq,1,round(50*(nFFT/2^17)));
