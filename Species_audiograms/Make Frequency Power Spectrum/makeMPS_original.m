function [mps,mps2,mpsCV,tmf,smf] = makeMPS(stimulus,fs)
%% Define some parameters

sampRate    = 48000;
fBand       = 125;                          % freq band (Hz)
nStd        = 6;
winTime     = (1000*nStd)/(fBand*2*pi);     % window size (msec), nStd Gaussian window
winLength   = fix(winTime*sampRate/1000);   % window size (nSamps); must be even
segSize     = sampRate;                     % max size (nSamps) for estimating mps
increment   = fix(0.001*segSize);           % sampRate of spectrogram
fLow        = 40;                           % Lower fq bound to get avg amplitude in spectrogram
fHigh       = 10000;                        % Upper fq bound to get avg amplitude in spectrogram

if mod(winLength, 2)
    winLength = winLength + 1;
end

frameCount = floor((segSize-winLength)/increment)+1;
nTMF = floor((segSize-winLength)/increment)+1;
nSMF = winLength/2+1;

%%
% filter stimulus segments to avoid edge artifacts
ham1 = hamming(segSize)';
ham2 = ham1(1:(segSize*0.0005):segSize/2);
filter = [ham2,ones(1,segSize-length(ham2)*2),fliplr(ham2)];

% check stimulus properties
if fs ~= sampRate
    stimulus = resample(stimulus,sampRate,fs);
    fs = sampRate;
end
if iscolumn(stimulus)
    stimulus = stimulus';
end
stimSize = length(stimulus);
nSegs = ceil(stimSize/segSize);
if nSegs == 1
    segIdx = [1, stimSize];
else
    segIdx = [(1:segSize:stimSize)', [segSize:segSize:stimSize, stimSize]'];
end
mps = zeros(nSMF,nTMF,nSegs);
mps2 = zeros(nSMF,nTMF,nSegs);

for seg = 1:nSegs
    currSeg = zeros(1,segSize);
    if seg < nSegs
        currSeg(1:diff(segIdx(seg,:))+1) = stimulus(segIdx(seg,1):segIdx(seg,2));
    else
        buffer = fix((segSize - diff(segIdx(seg,:)))/2);
        currSeg(buffer+(1:diff(segIdx(seg,:))+1)) = stimulus(segIdx(seg,1):segIdx(seg,2));
    end
    currSeg = currSeg .* filter;
    
    %----- GaussianSpectrum ---------------------------------------
    fftLength = winLength;
    
    % create Gaussian window
    wx2  = ((1:winLength)-((winLength+1)/2)).^2;
    wvar = (winLength/nStd)^2;
    ws   = exp(-0.5*(wx2./wvar));
    s = zeros(fftLength/2+1, frameCount);
    
    pg = zeros(1,frameCount);
    for win = 1:frameCount
        currIdx = (win-1)*increment + (1:winLength);
        f = ws.*currSeg(currIdx);
        
        specSlice = fft(f);
        s(:,win) = specSlice(1:(fftLength/2+1));
        pg(1,win) = std(f);
    end
    
    % Assign frequency & time labels
    freqLbl = (1:nSMF-1)' * (sampRate/fftLength);
    %timeLbl = (1:size(mps,2))';
    %--------------------------------------------------------------
    
    s_abs = log(abs(s)+1);
    
    % calculate the 2D fft
    f_abs = fft2(s_abs);
    f_pwr = real(f_abs.*conj(f_abs));
    f_pwr = fftshift(f_pwr);
    
    mps(:,:,seg)  = f_pwr;
    mps2(:,:,seg) = f_pwr.^2;
    
    % find axis labels
    fStep = freqLbl(2) - freqLbl(1);
    tmf = zeros(1,nTMF);
    smf = zeros(1,nSMF);
    
    if mod(nTMF,2)
        for tmp = 1:nTMF
            tmf(tmp) = (tmp-(nTMF+1)/2)*(1000/nTMF);
        end
    else
        for tmp = 1:nTMF
            tmf(tmp) = (tmp-nTMF/2-1)*(1000/nTMF);
        end
    end
    
    if mod(nSMF,2)
        for spc = 1:nSMF
            smf(spc) = (spc-(nSMF+1)/2)*(1/(fStep*nSMF))*1000;
        end
    else
        for spc = 1:nSMF
            smf(spc) = (spc-nSMF/2-1)*(1/(fStep*nSMF))*1000;
        end
    end
end

mps = sum(mps,3);
mps = mps/nSegs;

mps2 = sum(mps2,3);
mps2 = (mps2/(nSegs-1)) - (nSegs/(nSegs-1))*mps.^2;
mps2 = sqrt(mps2);

mpsCV = mps2 ./ mps;
