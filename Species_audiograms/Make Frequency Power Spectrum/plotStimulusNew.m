%% Tim - modified from Moises'script 02/13/20
% The script plots FPS of a sound file

% fq range changed to [0:1000:10000]
% y axis FPS changed to relative dB 

function plotStimulusNew
% displays the waveform, spectrogram, frequency power spectrum, & modulation power spectrum of user-selected sounds.
% ifdv2.m not functional

handles = struct;
if ispc
    handles.stimDir = 'C:\Users\Woolley Lab-MB\Desktop\TEST\';
elseif ismac
    handles.stimDir = '/Users/timyeh/Documents/Experimental Data/';
end
handles.figDir = handles.stimDir;

%% Define some parameters
handles.spec = 'fft';
handles.matchDur = false;       % 'fft','sparse'; sparse still has bugs
handles.prepost = false;        % makes all stimuli the same duration by adding silence to the end

switch handles.spec
    case 'fft'
        handles.tFs  = 500;                 % spec time sampling rate (bins/sec)
        handles.fFs  = 10;                  % spec freq sampling rate (bins/kHz)
        %handles.fLim = [500,16000];            % (MR: Originally [0,8000])frequency range (Hz)
        handles.fLim = [0, 10000];
        handles.dbRg = [0.25,0.75];         % sound floor/ceiling for computing spectrogram
        
    case 'sparse'
        handles.nFilters = 2000;
        handles.fLim     = [500,16000]; %(MR: Originally [1,8000])
end

%fBin = (1/handles.fFs)*1000;
%handles.FRQ = (handles.fLim(1)+fBin/2:fBin:handles.fLim(2));
handles.FRQ = (0:1000:10000);

handles.TMF = (-50:50); %(MR: Originally (-50:5:50))
handles.SMF = (0.1:0.05:2.3); %(MR: Originally (0:0.1:2.5))
handles.CON = 0.9;

stimInfo = struct;

%% Load stimuli
if ispc
    [fileNames, pathName] = uigetfile([handles.stimDir,'\*.wav'],'Choose stimulus files', 'MultiSelect','on');
elseif ismac
    [fileNames, pathName] = uigetfile([handles.stimDir,'/*.wav'],'Choose stimulus files', 'MultiSelect','on');
end

if ~iscell(fileNames)
    fileNames = {fileNames};
end
nFiles = length(fileNames);

for s = 1:nFiles
    filename = fileNames{s};
    if ispc
        [wave,fs] = audioread([pathName,'\',filename]);
    elseif ismac
        [wave,fs] = audioread([pathName,'/',filename]);
    end
    wave = wave(:,1);
    
    if fs ~= 48000
        wave = resample(wave,48000,fs);
        fs = 48000;
    end
    wave = wave(:)';
    
    stimInfo(s).name = filename(1:end-4);
    stimInfo(s).wave = wave;
    stimInfo(s).fs   = fs;
    stimInfo(s).size = length(wave);
end

% add prepost for plotting with spike data
if handles.prepost
    for s = 1:nFiles
        zPad = stimInfo(s).fs * 0.2;
        stimInfo(s).wave = [zeros(1,zPad), stimInfo(s).wave, zeros(1,zPad)];
        stimInfo(s).size = length(stimInfo(s).wave);
    end
end

% match stimulus durations
if handles.matchDur
    maxSize = max([stimInfo.size]);
    for s = 1:nFiles
        addSamps = maxSize - stimInfo(s).size;
        stimInfo(s).wave = [stimInfo(s).wave,zeros(1,addSamps)];
    end
end
stimInfo = rmfield(stimInfo,'size');

%% Process the sound
stimInfo = getEnvelope(handles,stimInfo);
stimInfo = getSpectrogram(handles,stimInfo);
stimInfo = getFPS(handles,stimInfo);
%stimInfo = getMPS(handles,stimInfo);

%% Save the stimulus data
if ispc
    save([handles.figDir,'\stimInfo.mat'],'stimInfo','-v7.3');
    
    % songInfo = stimInfo;
    % toneInfo = stimInfo;
    % save([figDir,'\toneInfo.mat'],'toneInfo','-v7.3');
elseif ismac
    
end

%% Print a figure
for s = 1:length(stimInfo)
    fig = drawFig(handles,stimInfo(s));
    print(fig,'-depsc','-loose','-painters', ...
        [handles.figDir,'\',stimInfo(s).name,'.eps']);
    %close(fig);
end


function [stimInfo] = getEnvelope(handles,stimInfo)
%% Compute the amplitude envelope
for s = 1:length(stimInfo)
    fs = (handles.fLim(2)*2);
    amp = resample(stimInfo(s).wave,fs,stimInfo(s).fs);
    time = (1:length(amp))/fs;
    
    % symmetric envelope --------------------------------------------------
    % env = abs(hilbert(stimInfo(s).wave));
    
    % asymmetric envelope -------------------------------------------------
    [posPks,posIdx] = findpeaks(amp);
    [negPks,negIdx] = findpeaks(-amp);
    negPks = -negPks;
    
    thr = max(abs(amp)) * 0.01;
    posPks(posPks<thr) = 0;
    negPks(negPks>-thr) = 0;
    
    posX = [1,posIdx(1),posIdx,posIdx(end),length(time)];
    posY = [0,0,posPks,0,0];
    negX = [1,negIdx(1),negIdx,negIdx(end),length(time)];
    negY = [0,0,negPks,0,0];
    
    stimInfo(s).time = time;
    stimInfo(s).amp  = amp;
    stimInfo(s).env_pos = [time(posX); posY];
    stimInfo(s).env_neg = [time(negX); negY];
end

function [stimInfo] = getSpectrogram(handles,stimInfo)
%% Compute the spectrogram
for s = 1:length(stimInfo)
    switch handles.spec
        case 'fft'
            Spec = makeSpectrogram(stimInfo(s).wave,stimInfo(s).fs, ...
                handles);
            
        case 'sparse'
%             [ifdgram,sonogram,fo] = ifdv2(stimInfo(s).wave,stimInfo(s).fs, ...
%                 handles);
%             spec = log(ifdgram+1);
%             freq = fo;
    end
    
    % adjust sampling intervals
    dur = ceil((length(stimInfo(s).wave)/stimInfo(s).fs)*1000)/1000;
    TT = ((1/handles.tFs)/2:(1/handles.tFs):dur);
    TT = round(TT*10000)/10000;
    FF = (handles.fLim(1)+(1000/handles.fFs):(1000/handles.fFs):handles.fLim(2))' - (1000/handles.fFs)/2;
    
    [X,Y] = meshgrid(Spec.T,Spec.F);
    [Xq,Yq] = meshgrid(TT,FF);
    
    Vq1 = interp2(X,Y,Spec.S,Xq,Yq,'linear');
    finIdx = isfinite(Vq1);
    Vq1(~finIdx) = min(Vq1(finIdx));
    Vq1 = (Vq1-nanmin(Vq1(:)))/range(Vq1(:));
    
    stimInfo(s).T = TT;
    stimInfo(s).F = FF;
    stimInfo(s).S = Vq1;
end

function [stimInfo] = getFPS(handles,stimInfo)
%% Compute frequency power spectrum
for s = 1:length(stimInfo)
    Freqs = makeFPS(stimInfo(s).wave,stimInfo(s).fs);
    
    % adjust sampling intervals
    smooth = hanning(5);
    smoother = smooth'/sum(smooth);
    fpsSm = convn(Freqs.fps,smoother,'same');
    FPS = interp1(Freqs.fq,fpsSm,handles.FRQ,'pchip');
    FPS = FPS/nansum(FPS(:));
    FPS = (FPS/max(FPS))*80; % Tim - normalize and convert to dB using 80 as highest value default
    
    % get power contour
    XS = sort(FPS(:),'descend');
    cumXS = cumsum(XS)/sum(XS);
    idx = find(cumXS>handles.CON,1,'first');
    CON = (FPS>=XS(idx));
    fW = diff(handles.FRQ(1:2));
    bw90 = sum(CON)*fW; % BW of primary frequencies constituting CON 5 of song
    
    stimInfo(s).FRQ         = handles.FRQ;
    stimInfo(s).FPS         = FPS;
    stimInfo(s).FPS_cMass   = sum(handles.FRQ.*FPS)/sum(FPS);         % center of mass, ie weighted average
    stimInfo(s).FPS_bw90    = bw90;
    
%     stimInfo(s).FPS_noiseIdx  = sum(FPS/max(FPS))/length(handles.FRQ);  % normalized ((1/length(freq))=tone, 1=noise)
%     stimInfo(s).FPS_bwRg_90   = [min(handles.FRQ(FPS>0.1)), max(freq(fps>0.1))]; %bw range
%     stimInfo(s).FPS_bwDist_90 = length(freq(fps>0.1))*fqBin; %bw distribution
end

%function [stimInfo] = getDynamicBW(handles,stimInfo)
%% Compute dynamic bandwidth
% for stim = 1:nFiles
%     pwrRg = [min(stimInfo(stim).spec(:)), max(stimInfo(stim).spec(:))];
%     lvlThr = pwrRg(2)-0.8*diff(pwrRg);
%     fqInt = mean(diff(stimInfo(stim).spec_fq));
%     
%     sigFrq = sum(stimInfo(stim).spec>lvlThr,1);
%     sigIdx = sigFrq > (200/fqInt);      % sig time bins must have bw >= 200 Hz
%     
%     [barY,barX] = hist(sigFrq(sigIdx),1:length(stimInfo(stim).spec_fq));
%     bw = (barX*fqInt)/1000;
%     bwPct = barY/sum(barY);
%     
%     stimInfo(stim).dynBW_X       = bw;
%     stimInfo(stim).dynBW_Y       = bwPct; % pct song that has a particular bandwidth
%     stimInfo(stim).dynBW_cMass   = sum(bw.*bwPct);
%     stimInfo(stim).dynBW_dist90  = sum(bwPct>(max(bwPct)*0.1))*mean(diff(bw));
% end

% function [stimInfo] = getMPS(handles,stimInfo)
% %% Compute modulation power spectrum
% for s = 1:length(stimInfo)
%     Mods = makeMPS(stimInfo(s).wave,stimInfo(s).fs,handles);
%     
%     % adjust sampling intervals
%     smooth = hanning(3);
%     smoother = repmat(smooth,1,length(smooth)) + repmat(smooth',length(smooth),1);
%     smoother = smoother/sum(smoother(:));
%     mpsSm = convn(Mods.mps,smoother,'same');
%     
%     [X,Y] = meshgrid(Mods.tmf,Mods.smf);
%     [Xq,Yq] = meshgrid(handles.TMF,handles.SMF);
%     MPS = interp2(X,Y,mpsSm,Xq,Yq,'linear');
%     MPS = MPS/sum(MPS(:));
%     MPSlog = log(MPS);
%     MPSrel = (MPSlog-min(MPSlog(:)))/range(MPSlog(:));
%     
%     TMC = nansum(MPS,1);
%     SMC = nansum(MPS,2);
%     
%     % get power contour
%     XS = sort(MPS(:),'descend');
%     cumXS = cumsum(XS)/sum(XS);
%     idx = find(cumXS>handles.CON,1,'first');
%     CON = (MPS>=XS(idx));
%     
%     stimInfo(s).TMF = handles.TMF;
%     stimInfo(s).SMF = handles.SMF';
%     stimInfo(s).MPS = MPSlog;
%     stimInfo(s).CON = CON;
%     stimInfo(s).TMC = TMC;
%     stimInfo(s).SMC = SMC;
%end

function [fig] = drawFig(handles,stimInfo)
%% Plot the figure
fig = figure;
set(fig,'PaperPositionMode','auto', ...
    'Units','normalized', ...
    'Position',[0 0 1 1]);

handles.cmap = [0,0,0; 0,0,0.5; jet];

% % stimulus waveform/envelope ----------------------------------------------
% ax(1) = subplot(3,5,[1,5]);
% %hold on;
% 
% ampPlot = 'wav'; % {'env','wav'}
% switch ampPlot
%     case 'wav'
%         line(stimInfo.time',stimInfo.amp','Color','k');
%         yRg = [min(stimInfo.amp),max(stimInfo.amp)];
%     case 'env'
%         line(stimInfo.env_pos(1,:)',stimInfo.env_pos(2,:)','Color','r');
%         line(stimInfo.env_neg(1,:)',stimInfo.env_neg(2,:)','Color','r');
%         yRg = [min(stimInfo.env_neg(2,:)),max(stimInfo.env_pos(2,:))];
% end
% yLim = [-1.05,1.05] * max(abs(yRg));
% 
% dur = ceil((length(stimInfo.wave)/stimInfo.fs)*handles.tFs)/handles.tFs;
% set(ax(1),'XLim',[0,dur], ...
%     'XTick',(0:1:dur), ...
%     'XTickLabel',(0:1:dur), ...
%     'YLim',yLim, ...
%     'YTick',[yRg(1), 0, yRg(2)], ...
%     'YTickLabel',[round(yRg(1)/0.01)*0.01, 0, round(yRg(2)/0.01)*0.01], ...
%     'TickDir','Out', ...
%     'TickLength',[0.005,0.005], ...
%     'FontSize',9);
% 
% ylabel('Amplitude','FontSize',9);
% title(stimInfo.name,'Interpreter','none', ...
%     'FontSize',12, ...
%     'FontWeight','bold');

% % spectrogram ---------------------------------------------------------
% ax(2) = subplot(3,5,[6,10]);
% %hold on;
% 
% colormap(ax(2),handles.cmap);
% 
% T = stimInfo.T;
% F = stimInfo.F;
% S = stimInfo.S;
% cRg  = [min(S(:)), max(S(:))];
% cLim = [cRg(1)+range(cRg)*0.01, cRg(2)-range(cRg)*0.1];
% 
% imagesc(T,F,S,cLim);
% 
% tW = 1/handles.tFs;
% set(ax(2),'XLim',[0,T(end)+(tW/2)], ...
%     'XTick',0:1:floor(T(end)), ...
%     'XTickLabel',0:1:floor(T(end)), ...
%     'YLim',handles.fLim, ...
%     'YTick',(0:2000:handles.fLim(end)), ...
%     'YTickLabel',(0:2000:handles.fLim(end))/1000, ...
%     'TickDir','Out', ...
%     'TickLength',[0.005 0.005], ...
%     'FontSize',9);
% xlabel('Time (sec)','FontSize',9);
% ylabel('Frequency (kHz)','FontSize',9);
% 
% linkaxes([ax(1),ax(2)],'x');

% frequency power spectrum --------------------------------------------
%ax(3) = subplot(3,5,11);
%hold on;

FRQ = [stimInfo.FRQ];
FPS = [stimInfo.FPS];
%yy1 = smooth(FPS, FRQ, 0.1);
plot(FRQ,FPS,'--k');
hold on

%line(FRQ, yy1, '-')
%yRg = [0,max(FPS)];
%yLim = ceil(yRg*100)/100;
%axis square;
% set(ax(3),'XLim',handles.fLim, ...
%     'XTick',(handles.fLim(1):1000:handles.fLim(2)), ...
%     'XTickLabel',(handles.fLim(1):1000:handles.fLim(2))/1000, ...
%     'YLim',yLim);
%xlabel('Frequency (kHz)');
%ylabel('Rel. Power');
%ylabel('Relative dB');
%title('LF FPS', ...
%    'FontWeight','normal');

%text(handles.fLim(2)*0.98,yLim(2)*0.95,sprintf('%s%4.2f','CtrMass = ',stimInfo.FPS_cMass/1000), ...
%   'FontSize',9, ...
%    'HorizontalAlignment','Right', ...
%    'VerticalAlignment','Middle');
%text(handles.fLim(2)*0.98,yLim(2)*0.85,sprintf('%s%4.2f','BW = ',stimInfo.FPS_bw90/1000), ...
%    'FontSize',9, ...
%    'HorizontalAlignment','Right', ...
%    'VerticalAlignment','Middle');
%
% dynamic bandwidth ------------------------------------------
%     ax(4) = subplot(3,5,12);
%     hold on;
%
%     cmap = [0,0,0.5; jet];
%     colormap(ax(4),cmap);
%
%     frq = [0,stimInfo(stim).dynBW_X,stimInfo(stim).dynBW_X(end)];
%     bw  = [0,stimInfo(stim).dynBW_Y,0];
%     line(frq,bw, ...
%         'Color','k');
%
%     xLim = [0,frq(end)];
%     yLim = [0,max(bw)];
%     axis square;
%     set(ax(4),'XLim',xLim, ...
%         'XTick',(0:1:frq(end)), ...
%         'XTickLabel',(0:1:frq(end)));
%     xlabel('Bandwidth (kHz)');
%     ylabel('Pct Song');
%     title('Dynamic Bandwidth', ...
%         'FontWeight','normal');
%
%     text(xLim(2)*0.98,yLim(2)*0.95,sprintf('%s%4.2f','ctrMass = ',stimInfo(stim).dynBW_cMass), ...
%         'FontSize',8, ...
%         'HorizontalAlignment','Right');
%     text(xLim(2)*0.98,yLim(2)*0.85,sprintf('%s%4.2f','dist90 = ',stimInfo(stim).dynBW_dist90), ...
%         'FontSize',8, ...
%         'HorizontalAlignment','Right');

% % modulation power spectrum -------------------------------------------
% ax(5) = subplot(3,5,13);
% hold on;
% 
% cmap = [0,0,0.5; jet];
% colormap(ax(5),cmap);
% 
% TMF = stimInfo.TMF;
% SMF = stimInfo.SMF;
% MPS = stimInfo.MPS;
% CON = stimInfo.CON;
% 
% cMargin = 0.05;
% 
% % mpsLog(1,:) = min(min(mpsLog(2:end,:)));
% mpsRel = (MPS-min(MPS(:)))/range(MPS(:));
% %disp(mpsRel)
% imagesc(TMF,SMF,mpsRel,[cMargin*2,1-cMargin*4]);
% contour(TMF,SMF,CON,1,'k-', ...
%     'LineWidth',2);
% 
% tW = mean(diff(TMF));
% sW = mean(diff(SMF));
% axis square;
% set(ax(5),'XLim',[TMF(1)-(tW/2),TMF(end)+(tW/2)], ...
%     'XTick',[TMF(1),0,TMF(end)], ...
%     'XTickLabel',[TMF(1),0,TMF(end)], ...
%     'YLim',[SMF(1)-(sW/2),SMF(end)+(sW/2)], ...
%     'YTick',[0,1,2], ...
%     'YTickLabel',[0,1,2]);
% xlabel('{\omega}_{t} (Hz)');
% ylabel('{\Omega}_{s} (cycles/kHz)');
% title('MPS', ...
%     'FontWeight','normal');
% 
% % temporal modulation curve -------------------------------------------
% ax(6) = subplot(3,5,14);
% hold on;
% 
% TMF = stimInfo.TMF;
% TMC = stimInfo.TMC;
% line(TMF,TMC, ...
%     'Color','k');
% 
% tW = mean(diff(TMF));
% axis square;
% set(ax(6),'XLim',[TMF(1)-(tW/2),TMF(end)+(tW/2)], ...
%     'XTick',(-50:25:50), ...
%     'XTickLabel',(-50:25:50), ...
%     'YScale','log');
% xlabel('{\omega}_{t} (Hz)');
% ylabel('Rel. Power');
% title('TMC', ...
%     'FontWeight','normal');
% 
% % spectral modulation curve -------------------------------------------
% ax(7) = subplot(3,5,15);
% hold on;
% 
% SMF = stimInfo.SMF;
% SMC = stimInfo.SMC;
% line(SMC,SMF, ...
%     'Color','k');
% 
% sW = mean(diff(SMF));
% axis square;
% set(ax(7),'XScale','log', ...
%     'YLim',[SMF(1)-(sW/2),SMF(end)+(sW/2)], ...
%     'YTick',(0:0.5:2), ...
%     'YTickLabel',(0:0.5:2));
% xlabel('Rel. Power');
% ylabel('{\Omega}_{s} (cycles/kHz)');
% title('SMC', ...
%     'FontWeight','normal');


function [Spec] = makeSpectrogram(wave,fs,handles)
%%
Spec = struct;

dBNoise     = handles.dbRg(1);
dBMax       = handles.dbRg(2);
ampSampRate = handles.tFs;
fWidth      = (1000/handles.fFs);
minFreq     = handles.fLim(1);
maxFreq     = handles.fLim(2);

% Define more parameters
increment = floor(fs/ampSampRate);                      %# samples to shift for each new bin start point
binSize = floor(1/(2*pi*fWidth)*6*fs);                  %size of time bin over which fft is computed (in samples)
nFTfreqs = binSize;                                     %# of frequencies at which fft is computed

% Zero-pad the stimulus
padSize = binSize;
stimulus = [zeros(1,padSize), wave, zeros(1,padSize)];
stimSize = length(stimulus);
frameCount = floor((stimSize-binSize*2)/increment)+1;   %# of bins
% frameCount = ceil((length(wave)/fs)*1000);

% Compute the gaussian filter
wx2 = ((1:binSize)-binSize/2).^2;
wvar = (binSize/6)^2;
ws = exp(-0.5*(wx2./wvar)); 

% Compute spectrogram of entire stimulus (using a sliding Fourier transform)
s = zeros(binSize/2+1, frameCount);
for bin = 1:frameCount
    frst = padSize -binSize/2 + (bin-1)*increment + round(increment/2);
%     frst = (bin-1)*increment + (increment/2);
    last = frst + binSize - 1;
    f = zeros(binSize, 1);
    f(1:binSize) = ws.*stimulus(frst:last);
    binspec = fft(f);
    s(:,bin) = binspec(1:(nFTfreqs/2+1));
end

% Translate to dB, rectify
tmp = 20*log(abs(s)/max(abs(s(:))));
tmp = max(0, 20*log(abs(s)/max(abs(s(:)))) + dBNoise*range(tmp(isfinite(tmp(:,:)))));
tmp = min(dBMax*max(tmp(:)),tmp);

% % original, dBNoise & dBMax in dB
% tmp = max(0, 20*log(abs(s)./max(abs(s(:))))+dBNoise);
% tmp = min(dBMax,tmp);

% Edit out the appropriate range of frequencies
select = 1:nFTfreqs/2+1;
fo = (select-1)'*fs/nFTfreqs;
freq_range = fo>=minFreq & fo<=maxFreq;
outSpectrum = tmp(freq_range,:);
outFreqs = fo(freq_range);

Spec.T = linspace((increment/2),(increment/2)+increment*frameCount,frameCount)/fs;
Spec.F = outFreqs;
Spec.S = outSpectrum;

function [Freqs] = makeFPS(stimulus,fs)
Freqs = struct;

%% compute the frequency power spectrum
% Do the FFT
nFFT = 2^nextpow2(length(stimulus));
fft1 = fft(stimulus,nFFT)/length(stimulus);
fft2 = 2*abs(fft1(1:nFFT/2+1));
if round(50*(nFFT/2^17))>0
    fft3 = resample(fft2,1,round(50*(nFFT/2^17))); % smooth the fft
else
    fft3 = fft2;
end
fps  = fft3/max(fft3);

% Get freq values
fq = fs/2*linspace(0,1,nFFT/2+1);
if round(50*(nFFT/2^17))>0
    freq = resample(fq,1,round(50*(nFFT/2^17)));
else
    freq = fq;
end

Freqs.fq  = freq;
Freqs.fps = fps;

% function [Mods] = makeMPS(stimulus,fs,handles)
% Mods = struct;
% 
% %% Define some parameters
% sampRate    = 48000;
% fBand       = 1000/handles.fFs;             % freq band (Hz)
% nStd        = 6;
% winTime     = (1000*nStd)/(fBand*2*pi);     % window size (msec), nStd Gaussian window
% winLength   = fix(winTime*sampRate/1000);   % window size (nSamps); must be even
% segSize     = sampRate;                     % max size (nSamps) for estimating mps
% increment   = fix(0.001*segSize);           % sampRate of spectrogram
% % fLow        = 40;                           % Lower fq bound to get avg amplitude in spectrogram
% % fHigh       = 10000;                        % Upper fq bound to get avg amplitude in spectrogram
% 
% if mod(winLength, 2)
%     winLength = winLength + 1;
% end
% 
% frameCount = floor((segSize-winLength)/increment)+1;
% nTMF = floor((segSize-winLength)/increment)+1;
% nSMF = winLength/2+1;
% 
% % filter stimulus segments to avoid edge artifacts
% ham1 = hamming(segSize)';
% ham2 = ham1(1:(segSize*0.0005):segSize/2);
% filter = [ham2,ones(1,segSize-length(ham2)*2),fliplr(ham2)];
% 
% % check stimulus properties
% if fs ~= sampRate
%     stimulus = resample(stimulus,sampRate,fs);
%     fs = sampRate;
% end
% if iscolumn(stimulus)
%     stimulus = stimulus';
% end
% stimSize = length(stimulus);
% nSegs = ceil(stimSize/segSize);
% if nSegs == 1
%     segIdx = [1, stimSize];
% else
%     segIdx = [(1:segSize:stimSize)', [segSize:segSize:stimSize, stimSize]'];
% end
% mps = zeros(nSMF,nTMF,nSegs);
% mps2 = zeros(nSMF,nTMF,nSegs);
% 
% for seg = 1:nSegs
%     currSeg = zeros(1,segSize);
%     if seg < nSegs
%         currSeg(1:diff(segIdx(seg,:))+1) = stimulus(segIdx(seg,1):segIdx(seg,2));
%     else
%         buffer = fix((segSize - diff(segIdx(seg,:)))/2);
%         currSeg(buffer+(1:diff(segIdx(seg,:))+1)) = stimulus(segIdx(seg,1):segIdx(seg,2));
%     end
%     currSeg = currSeg .* filter;
%     
%     %----- GaussianSpectrum ---------------------------------------
%     fftLength = winLength;
%     
%     % create Gaussian window
%     wx2  = ((1:winLength)-((winLength+1)/2)).^2;
%     wvar = (winLength/nStd)^2;
%     ws   = exp(-0.5*(wx2./wvar));
%     s = zeros(fftLength/2+1, frameCount);
%     
%     pg = zeros(1,frameCount);
%     for win = 1:frameCount
%         currIdx = (win-1)*increment + (1:winLength);
%         f = ws.*currSeg(currIdx);
%         
%         specSlice = fft(f);
%         s(:,win) = specSlice(1:(fftLength/2+1));
%         pg(1,win) = std(f);
%     end
%     
%     % Assign frequency & time labels
%     freqLbl = (1:nSMF-1)' * (sampRate/fftLength);
%     %timeLbl = (1:size(mps,2))';
%     %--------------------------------------------------------------
%     
%     s_abs = log(abs(s)+1);
%     
%     % calculate the 2D fft
%     f_abs = fft2(s_abs);
%     f_pwr = real(f_abs.*conj(f_abs));
%     f_pwr = fftshift(f_pwr);
%     
%     mps(:,:,seg)  = f_pwr;
%     mps2(:,:,seg) = f_pwr.^2;
%     
%     % find axis labels
%     fStep = freqLbl(2) - freqLbl(1);
%     tmf = zeros(1,nTMF);
%     smf = zeros(1,nSMF);
%     
%     if mod(nTMF,2)
%         for tmp = 1:nTMF
%             tmf(tmp) = (tmp-(nTMF+1)/2)*(1000/nTMF);
%         end
%     else
%         for tmp = 1:nTMF
%             tmf(tmp) = (tmp-nTMF/2-1)*(1000/nTMF);
%         end
%     end
%     
%     if mod(nSMF,2)
%         for spc = 1:nSMF
%             smf(spc) = (spc-(nSMF+1)/2)*(1/(fStep*nSMF))*1000;
%         end
%     else
%         for spc = 1:nSMF
%             smf(spc) = (spc-nSMF/2-1)*(1/(fStep*nSMF))*1000;
%         end
%     end
% end
% 
% mps = sum(mps,3);
% mps = mps/nSegs;
% 
% mps2 = sum(mps2,3);
% mps2 = (mps2/(nSegs-1)) - (nSegs/(nSegs-1))*mps.^2;
% mps2 = sqrt(mps2);
% mpsCV = mps2 ./ mps;
% 
% Mods.tmf = tmf;
% Mods.smf = smf;
% Mods.mps = mps;
% Mods.cv  = mpsCV;
% 
