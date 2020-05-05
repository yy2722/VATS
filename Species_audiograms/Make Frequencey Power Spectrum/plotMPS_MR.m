function plotStimulusNew
% displays the waveform, spectrogram, frequency power spectrum, & modulation power spectrum of user-selected sounds.
% ifdv2.m not functional
% Edit, Moises (4/24/19): saves MPS fig and RelMPS data in figDir for
% further analysis using analyzeMPS.m

handles = struct;
if ispc
    handles.stimDir = 'Z:\undergrad\User\Pooja\MPS-FPS Song Selection\Multispecies\MPS\Sound Files\'; %where to select your wav files
elseif ismac
    handles.stimDir = 'Z:\mr3683\Song Recordings\Sample Files\';
end
handles.figDir = 'Z:\undergrad\User\Pooja\MPS-FPS Song Selection\Multispecies\MPS\figs\'; %where MPS will be saved

%% Define some parameters
handles.spec = 'fft';
handles.matchDur = false;       % 'fft','sparse'; sparse still has bugs
handles.prepost = false;        % makes all stimuli the same duration by adding silence to the end

switch handles.spec
    case 'fft'
        handles.tFs  = 500;                 % spec time sampling rate (bins/sec)
        handles.fFs  = 10;                  % spec freq sampling rate (bins/kHz)
        handles.fLim = [500,16000];            % (MR: Originally [0,8000])frequency range (Hz)
        handles.dbRg = [0.25,0.75];         % sound floor/ceiling for computing spectrogram
        
    case 'sparse'
        handles.nFilters = 2000;
        handles.fLim     = [500,16000]; %(MR: Originally [1,8000])
end

fBin = (1/handles.fFs)*1000;
handles.FRQ = (handles.fLim(1)+fBin/2:fBin:handles.fLim(2));

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

handles.input = fileNames;

stimInfo = rmfield(stimInfo,'size');

%% Process the sound
stimInfo = getMPS(handles,stimInfo);

%% Save the stimulus data
if ispc
    % save([figDir,'\stimInfo.mat'],'stimInfo','-v7.3');
    
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
    close(fig);
end

function [stimInfo] = getMPS(handles,stimInfo)
%% Compute modulation power spectrum
for s = 1:length(stimInfo)
    Mods = makeMPS(stimInfo(s).wave,stimInfo(s).fs,handles);
    
    % adjust sampling intervals
    smooth = hanning(3);
    smoother = repmat(smooth,1,length(smooth)) + repmat(smooth',length(smooth),1);
    smoother = smoother/sum(smoother(:));
    mpsSm = convn(Mods.mps,smoother,'same');
    
    [X,Y] = meshgrid(Mods.tmf,Mods.smf);
    [Xq,Yq] = meshgrid(handles.TMF,handles.SMF);
    MPS = interp2(X,Y,mpsSm,Xq,Yq,'linear');
    MPS = MPS/sum(MPS(:));
    MPSlog = log(MPS);
    MPSrel = (MPSlog-min(MPSlog(:)))/range(MPSlog(:));
    
    TMC = nansum(MPS,1);
    SMC = nansum(MPS,2);
    
    % get power contour
    XS = sort(MPS(:),'descend');
    cumXS = cumsum(XS)/sum(XS);
    idx = find(cumXS>handles.CON,1,'first');
    CON = (MPS>=XS(idx));
    
    stimInfo(s).TMF = handles.TMF;
    stimInfo(s).SMF = handles.SMF';
    stimInfo(s).MPS = MPSlog;
    stimInfo(s).CON = CON;
    stimInfo(s).TMC = TMC;
    stimInfo(s).SMC = SMC;
end

function [fig] = drawFig(handles,stimInfo)
%% Plot the figure
fig = figure;
set(fig,'PaperPositionMode','auto', ...
    'Units','normalized', ...
    'Position',[0 0 1 1]);

handles.cmap = [0,0,0; 0,0,0.5; jet];

% modulation power spectrum -------------------------------------------
fileChar = char(handles.input);

ax(5) = subplot(1,1,1);
hold on;

cmap = [0,0,0.5; jet];
colormap(ax(5),cmap);

TMF = stimInfo.TMF;
SMF = stimInfo.SMF;
MPS = stimInfo.MPS;
CON = stimInfo.CON;

cMargin = 0.05;

% mpsLog(1,:) = min(min(mpsLog(2:end,:)));
mpsRel = (MPS-min(MPS(:)))/range(MPS(:));
imagesc(TMF,SMF,mpsRel,[cMargin*2,1-cMargin*4]);
contour(TMF,SMF,CON,1,'k-', ...
    'LineWidth',2);

tW = mean(diff(TMF));
sW = mean(diff(SMF));
axis square;
set(ax(5),'XLim',[TMF(1)-(tW/2),TMF(end)+(tW/2)], ...
    'XTick',[TMF(1),0,TMF(end)], ...
    'XTickLabel',[TMF(1),0,TMF(end)], ...
    'YLim',[SMF(1)-(sW/2),SMF(end)+(sW/2)], ...
    'YTick',[0,1,2], ...
    'YTickLabel',[0,1,2]);
xlabel('{\omega}_{t} (Hz)');
ylabel('{\Omega}_{s} (cycles/kHz)');
title(['MPS_' fileChar],'Interpreter', 'none', ...
    'FontWeight','normal');
%assignin('base','RelMPS',mpsRel);
save([handles.figDir [fileChar '_RelMPS.mat']]);


function [Mods] = makeMPS(stimulus,fs,handles)
Mods = struct;

%% Define some parameters
sampRate    = 48000;
fBand       = 1000/handles.fFs;             % freq band (Hz)
nStd        = 6;
winTime     = (1000*nStd)/(fBand*2*pi);     % window size (msec), nStd Gaussian window
winLength   = fix(winTime*sampRate/1000);   % window size (nSamps); must be even
segSize     = sampRate;                     % max size (nSamps) for estimating mps
increment   = fix(0.001*segSize);           % sampRate of spectrogram
% fLow        = 40;                           % Lower fq bound to get avg amplitude in spectrogram
% fHigh       = 10000;                        % Upper fq bound to get avg amplitude in spectrogram

if mod(winLength, 2)
    winLength = winLength + 1;
end

frameCount = floor((segSize-winLength)/increment)+1;
nTMF = floor((segSize-winLength)/increment)+1;
nSMF = winLength/2+1;

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

Mods.tmf = tmf;
Mods.smf = smf;
Mods.mps = mps;
Mods.cv  = mpsCV;