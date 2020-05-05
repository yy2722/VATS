function [outSpectrum,outFreqs] = makeSpectrogram(stimulus,fs,dBNoise,dBMax,ampSampRate,fWidth,minFreq,maxFreq)
% Computes the spectrogram of an input wave file.
% Stimulus must be resampled to 48 kHz in the calling program.

%% These parameters should be defined within the calling script
% dBNoise     = 50;           % noise floor, low = more blue
% dBMax       = 35;           % low = more red
% ampSampRate = 1000;         % must be a factor of 48000 so 'increment' is an integer (originally 1000)
% fWidth      = 125;          % must be a # so binSize is even (80,100,125,145; originally 125)
% minFreq     = 0;            % min freq
% maxFreq     = 8000;         % max freq

%% Define more parameters
increment = floor(fs/ampSampRate);                      %# samples to shift for each new bin start point
binSize = floor(1/(2*pi*fWidth)*6*fs);                  %size of time bin over which fft is computed (in samples)
nFTfreqs = binSize;                                     %# of frequencies at which fft is computed

%% Zero-pad the stimulus
stimulus = stimulus(:)';
padSize = binSize/2;
stimulus = [zeros(1,padSize), stimulus, zeros(1,padSize)];
stimSize = length(stimulus);
frameCount = floor((stimSize-binSize)/increment)+1;     %# of bins

%% Compute the gaussian filter
wx2 = ((1:binSize)-binSize/2).^2;
wvar = (binSize/6)^2;
ws = exp(-0.5*(wx2./wvar)); 

%% Compute spectrogram of entire stimulus (using a sliding Fourier transform)
s = zeros(binSize/2+1, frameCount);
for bin = 1:frameCount
    frst = (bin-1)*increment + 1;
    last = frst + binSize - 1;
    f = zeros(binSize, 1);
    f(1:binSize) = ws.*stimulus(frst:last);
    binspec = fft(f);
    s(:,bin) = binspec(1:(nFTfreqs/2+1));
end

%% Translate to dB, rectify
tmp = 20*log(abs(s)/max(abs(s(:))));
tmp = max(0, 20*log(abs(s)/max(abs(s(:)))) + dBNoise*range(tmp(isfinite(tmp(:,:)))));
tmp = min(dBMax*max(tmp(:)),tmp);

% % original, dBNoise & dBMax in dB
% tmp = max(0, 20*log(abs(s)./max(abs(s(:))))+dBNoise);
% tmp = min(dBMax,tmp);

%% Edit out the appropriate range of frequencies
select = 1:nFTfreqs/2+1;
fo = (select-1)'*fs/nFTfreqs;
freq_range = fo>=minFreq & fo<=maxFreq;
outSpectrum = tmp(freq_range,:);
outFreqs = fo(freq_range);

%% Align spectrogram & stimulus waveform
stimulus = stimulus(1,padSize+1:end-padSize);           %remove zero-pads
stimDur = length(stimulus)/fs;                          %duration (sec) of stimulus
specDur = size(outSpectrum,2)/ampSampRate;              %duration (sec) of spectrogram
diffStimSpec = stimDur - specDur;
diffSpecSamps = round(diffStimSpec*ampSampRate);
if diffSpecSamps > 0
    outSpectrum = [outSpectrum, zeros(size(outSpectrum,1),diffSpecSamps)];
elseif diffSpecSamps < 0
    outSpectrum = outSpectrum(:,1:(end+diffSpecSamps));
end
