function plotStimulus
% displays the waveform, spectrogram, frequency power spectrum, & modulation power spectrum of user-selected sounds.
% calls makeSpectrogram.m, ifdv2.m, makeFPS.m, makeMPS.m
% Modified by Moises 3/12/19 to pull from up to 16kHz

specMeth = 'fft';           % 'fft','sparse'; sparse still has bugs
matchDur = true;            % adds silence to make all stimuli same duration before plotting
cmap     = 'jet';
pause on;

%% Set up directories
% stimDir = 'C:\Users\Woolley Lab\Documents\MATLAB\jm';
% figDir  = 'C:\Users\Woolley Lab\Documents\MATLAB\jm';
% stimDir = 'Z:\User\Jordan\Projects\estrildids\songs\filtered';
% figDir  = 'Z:\User\Jordan\Projects\estrildids\songs\specs_kjet';
stimDir = 'C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\MPS\Sound Files';
figDir  = 'C:\Users\Woolley Lab-MB\Desktop\Multispecies Vocalizations\MPS';


% load([figDir,'\stimInfo.mat']);
% nFiles = length(stimInfo);

%% Define some parameters
switch specMeth
    case 'fft'
        switch cmap
            case 'bone'
                dBNoise     = 0.20;          % noise floor {tone,ripp,song} (, 0.80, 0.20)
                dBMax       = 0.50;          % max plateau {tone,ripp,song} (, 0.91, 0.50)
            case 'copper'
                dBNoise    = 0.25;          % noise floor {tone,ripp,song} (, 0.80, 0.20)
                dBMax      = 0.70;          % max plateau {tone,ripp,song} (, 0.91, 0.50)
            case 'hot'
                dBNoise    = 0.20;          % noise floor {tone,ripp,song} (, 0.80, 0.20)
                dBMax      = 0.60;          % max plateau {tone,ripp,song} (, 0.91, 0.50)
            case 'jet'
                dBNoise     = 0.2;          % noise floor {tone,ripp,song} (0.9, 0.80, 0.25)
                dBMax       = 0.8;          % max plateau {tone,ripp,song} (0.92, 0.91, 0.75)
        end
        
        ampSampRate = 1000;         % must be a factor of 48000 so 'increment' is an integer (was originally 1000)
        fWidth      = 100;          % must be a # so binSize is even (80,100,125,145; was originally 125)
        minFreq     = 500;            % minimum frequency (MR: Originally 0) 
        maxFreq     = 16000;         % maximum frequency %(MR: originally 8000)
        
    case 'sparse'
        nFilters    = 2000;
        minFreq     = 500;            % min freq (MR: originally 1)
        maxFreq     = 16000;         % max freq %(MR: originally 8000)
        
end

stimInfo = struct;

%% Load stimuli
[fileNames, pathName] = uigetfile([stimDir,'\*.wav'],'Choose stimulus files', 'MultiSelect','on');
if ~iscell(fileNames)
    fileNames = {fileNames};
end
nFiles = length(fileNames);

for stim = 1:nFiles
    filename = fileNames{stim};
    [stimulus, fs] = audioread([pathName,'\',filename]);
    
    if fs ~= 48000
        stimulus = resample(stimulus,48000,fs);
        fs = 48000;
    end
    stimulus = stimulus(:)';
    
    stimInfo(stim).name = filename(1:end-4);
    stimInfo(stim).wave = stimulus;
    stimInfo(stim).fs   = fs;
    stimInfo(stim).size = length(stimulus);
end

% match stimulus durations
if matchDur
    maxSize = max([stimInfo.size]);
    for stim = 1:nFiles
        addSamps = maxSize - stimInfo(stim).size;
        stimInfo(stim).wave = [stimInfo(stim).wave,zeros(1,addSamps)];
    end
end

%% Compute the amplitude envelope
for stim = 1:nFiles
    % symmetric envelope --------------------------------------------------
    %envelope = abs(hilbert(stimInfo(stim).wave));
    
    % asymmetric envelope -------------------------------------------------
    %smoother = hanning(10);
    smoother = ones(1,10);
    smthStim = convn(stimInfo(stim).wave,smoother,'same');
    [posPks,posIdx] = findpeaks(smthStim);
    [negPks,negIdx] = findpeaks(-smthStim);
    negPks = -negPks;
    
    % remove spurious peaks
    posIdx(posPks<=0) = [];
    posPks(posPks<=0) = [];
    negIdx(negPks>=0) = [];
    negPks(negPks>=0) = [];
    
    thresh = 0.1;
    rmIdx = false(1,length(posIdx));
    for smp = 3:length(posIdx)-3
        maxPk = max(posPks(smp-2:smp+2));
        if posPks(smp) < maxPk*thresh
            rmIdx(smp) = 1;
        end
    end
    posIdx(rmIdx) = [];
    posPks(rmIdx) = [];
    
    rmIdx = false(1,length(negIdx));
    for smp = 3:length(negIdx)-3
        minPk = max(negPks(smp-2:smp+2));
        if negPks(smp) > minPk*thresh
            rmIdx(smp) = 1;
        end
    end
    negIdx(rmIdx) = [];
    negPks(rmIdx) = [];
    
    stimInfo(stim).env_posX = posIdx;
    stimInfo(stim).env_posY = posPks;
    stimInfo(stim).env_negX = negIdx;
    stimInfo(stim).env_negY = negPks;
end

%% Compute the spectrograms
for stim = 1:nFiles
    switch specMeth
        case 'fft'
            [outSpectrum,outFreqs] = makeSpectrogram(stimInfo(stim).wave,stimInfo(stim).fs, ...
                dBNoise,dBMax,ampSampRate,fWidth,minFreq,maxFreq);
            stimInfo(stim).spec     = outSpectrum;
            stimInfo(stim).spec_fq  = outFreqs;
            stimInfo(stim).spec_fs  = ampSampRate;
            
        case 'sparse'
            [ifdgram,sonogram,fo] = ifdv2(stimInfo(stim).wave,stimInfo(stim).fs, ...
                nFilters,minFreq,maxFreq);
            stimInfo(stim).spec     = log(ifdgram+1);
            stimInfo(stim).spec_fq  = fo;
            stimInfo(stim).spec_fs  = [];
            
    end
end

%% Compute frequency power spectrum
for stim = 1:nFiles
    [fps,frq] = makeFPS(stimInfo(stim).wave,stimInfo(stim).fs);
    
    smooth = hanning(5);
    smoother = smooth'/sum(smooth);
    fpsSm = convn(fps,smoother,'same');
    FRQ = (0.025:0.025:(maxFreq/1000));
    FPS = interp1((frq/1000),fpsSm,FRQ,'pchip');
    FPS = FPS/sum(FPS(:));
    
    stimInfo(stim).frq           = FRQ;
    stimInfo(stim).fps           = FPS;
    stimInfo(stim).fps_cMass     = sum(FRQ.*FPS)/sum(FPS); % center of mass, ie weighted average
    stimInfo(stim).fps_noiseIdx  = sum(FPS)/length(FRQ); %normalized ((1/length(freq))=tone, 1=noise)
    
%     stimInfo(stim).frq           = FRQ;
%     stimInfo(stim).fps           = FPS;
%     stimInfo(stim).fps_pk        = freq(fps==max(fps));
%     stimInfo(stim).fps_cMass     = sum(freq.*fps)/sum(fps); % center of mass, ie weighted average
%     stimInfo(stim).fps_bwRg_90   = [min(freq(fps>0.1)), max(freq(fps>0.1))]; %bw range
%     stimInfo(stim).fps_bwDist_90 = length(freq(fps>0.1))*fqBin; %bw distribution
%     stimInfo(stim).fps_noiseIdx  = sum(fps)/length(freq); %normalized ((1/length(freq))=tone, 1=noise)
end

%% Compute dynamic bandwidth
for stim = 1:nFiles
    pwrRg = [min(stimInfo(stim).spec(:)), max(stimInfo(stim).spec(:))];
    lvlThr = pwrRg(2)-0.8*diff(pwrRg);
    fqInt = mean(diff(stimInfo(stim).spec_fq));
    
    sigFrq = sum(stimInfo(stim).spec>lvlThr,1);
    sigIdx = sigFrq > (200/fqInt);      % sig time bins must have bw >= 200 Hz
    
    [barY,barX] = hist(sigFrq(sigIdx),1:length(stimInfo(stim).spec_fq));
    bw = (barX*fqInt)/1000;
    bwPct = barY/sum(barY);
    
    stimInfo(stim).dynBW_X       = bw;
    stimInfo(stim).dynBW_Y       = bwPct; % pct song that has a particular bandwidth
    stimInfo(stim).dynBW_cMass   = sum(bw.*bwPct);
    stimInfo(stim).dynBW_dist90  = sum(bwPct>(max(bwPct)*0.1))*mean(diff(bw));
end

%% Compute modulation power spectrum
conPwr = 0.9;
for stim = 1:nFiles
    [mps,mps2,mpsCV,tmf,smf] = makeMPS(stimInfo(stim).wave,stimInfo(stim).fs);
    
    TMF = (-50:50);
    SMF = (0.05:0.05:2); %(MR: Originally (0:0.05:2)
    [X,Y] = meshgrid(tmf,smf);
    [Xq,Yq] = meshgrid(TMF,SMF);
    MPS = interp2(X,Y,mps,Xq,Yq);
    MPS = MPS/sum(MPS(:));
    
    TMC = nansum(MPS,1);
    SMC = nansum(MPS,2);
    
    % get power contour
    XS = sort(MPS(:),'descend');
    cumXS = cumsum(XS)/sum(XS);
    idx = cell(1,length(conPwr));
    for ctr = 1:length(conPwr)
        diffXS = abs(cumXS-conPwr(ctr));
        [~,idx{ctr}] = min(diffXS);
        CON = XS(idx{ctr}(1));
    end
    
%     % log-transform & relativize
%     mpsTemp = log(mps);
%     conTemp = log(mpsContour);
%     
%     mpsLog = (mpsTemp-min(mpsTemp(:)))/range(mpsTemp(:));
%     conLog = (conTemp-min(mpsTemp(:)))/range(mpsTemp(:));
    
    stimInfo(stim).tmf = TMF;
    stimInfo(stim).smf = SMF';
    stimInfo(stim).mps = MPS;
    stimInfo(stim).con = CON;
    stimInfo(stim).tmc = TMC;
    stimInfo(stim).smc = SMC;
end

%% Get acoustic features
% for stim = 1:nFiles
%     wave = stimInfo(stim).wave;
%     if length(wave)<882 % syllable must be >=20msec
%         wave = [wave,zeros(1,882-length(wave))];
%     end
%     amp = wave./max(wave);
%     
%     features = features_10082013(amp,fs);
%     gravity_center = nanmean(features(:,1),1);
%     gravity_var = nanmean(features(:,2),1);
%     entropy = nanmean(features(:,3),1);
%     good_pitch = nanmean(features(:,4),1);
%     
%     stimInfo(stim).gravity = gravity_center;
%     stimInfo(stim).entropy = entropy;
%     stimInfo(stim).good_pitch = good_pitch;
% end

%% Save the stimulus data
save([figDir,'\stimInfo.mat'],'stimInfo','-v7.3');

%% Plot the figure
ampPlot = 'wav'; % {'env','wav'}
for stim = 1:nFiles
    fig = figure;
    set(fig,'PaperPositionMode','auto', ...
        'Units','normalized', ...
        'Position',[0 0 1 1]);
    
    % stimulus envelope/waveform ------------------------------------------
    ax(1) = subplot(3,5,[1,5]);
    hold on;
    
    fs = stimInfo(stim).fs;
    switch ampPlot
        case 'env'
            wave = [0,stimInfo(stim).posEnvY,0,stimInfo(stim).negEnvY,0];
            time = [0,stimInfo(stim).posEnvY,0,stimInfo(stim).negEnvY,0];
            fill(time,wave,'k');
            
        case 'wav'
            wave = stimInfo(stim).wave;
            time = (1:length(wave))/fs;
            plot(ax(1),time,wave,'k-');
            
    end
    xLim = [0,time(end)];
    yLim = [min(wave)-0.05*range(wave), max(wave)+0.05*range(wave)];
    set(ax(1),'XLim',xLim, ...
        'XTick',0:1:floor(time(end)), ...
        'XTickLabel',0:1:floor(time(end)), ...
        'YLim',yLim, ...
        'YTick',[min(wave), 0, max(wave)], ...
        'YTickLabel',[round(min(wave)/0.01)*0.01, 0, round(max(wave)/0.01)*0.01], ...
        'TickLength',[0.005,0.005], ...
        'FontSize',9);
    ylabel('Amplitude','FontSize',9);
    title(stimInfo(stim).name,'Interpreter','none');
    
    % spectrogram ---------------------------------------------------------
    ax(2) = subplot(3,5,[6,10]);
    hold on;
    
    cmap = jet;
    cmap = [0,0,0; 0,0,0.5; cmap];
    colormap(ax(2),cmap);
    
    spec = stimInfo(stim).spec;
    time = (1:1:size(spec,2))/stimInfo(stim).spec_fs;
    freq = stimInfo(stim).spec_fq;
    
    cRg = [min(spec(:)), max(spec(:))];
    cRg = [cRg(1)+range(cRg)*0.01, cRg(2)-range(cRg)*0.01];
    
    imagesc(time,freq,spec,cRg);
    
    set(ax(2),'XLim',[0,time(end)], ...
        'XTick',0:1:floor(time(end)), ...
        'XTickLabel',0:1:floor(time(end)), ...
        'YLim',[0,maxFreq], ...
        'YTick',linspace(minFreq,maxFreq,6), ...
        'YTickLabel',linspace(minFreq,maxFreq,6)/1000, ...
        'TickLength',[0.005 0.005], ...
        'FontSize',9);
    xlabel('Time (sec)','FontSize',9);
    ylabel('Frequency (kHz)','FontSize',9);
    
    linkaxes([ax(1),ax(2)],'x');
    
    % frequency power spectrum --------------------------------------------
    ax(3) = subplot(3,5,11);
    hold on;
    
    frq = [0,stimInfo(stim).frq];
    fps = [0,stimInfo(stim).fps];
    line(frq,fps, ...
        'Color','k');
    xLim = [0,frq(end)];
    yLim = [0,max(fps)];
    axis square;
    set(ax(3),'XLim',xLim, ...
        'XTick',(0:1:frq(end)), ...
        'XTickLabel',(0:1:frq(end)));
    xlabel('Frequency (kHz)');
    ylabel('Prop. Power');
    title('FPS', ...
        'FontWeight','normal');
    
    text(xLim(2)*0.98,yLim(2)*0.95,sprintf('%s%4.2f','ctrMass = ',stimInfo(stim).fps_cMass), ...
        'FontSize',8, ...
        'HorizontalAlignment','Right');
    text(xLim(2)*0.98,yLim(2)*0.85,sprintf('%s%4.2f','noiseIdx = ',stimInfo(stim).fps_noiseIdx), ...
        'FontSize',8, ...
        'HorizontalAlignment','Right');
    
    % dynamic bandwidth ------------------------------------------
    ax(4) = subplot(3,5,12);
    hold on;
    
    frq = [0,stimInfo(stim).dynBW_X,stimInfo(stim).dynBW_X(end)];
    bw  = [0,stimInfo(stim).dynBW_Y,0];
    line(frq,bw, ...
        'Color','k');
    
    xLim = [0,frq(end)];
    yLim = [0,max(bw)];
    axis square;
    set(ax(4),'XLim',xLim, ...
        'XTick',(0:1:frq(end)), ...
        'XTickLabel',(0:1:frq(end)));
    xlabel('Bandwidth (kHz)');
    ylabel('Pct Song');
    title('Dynamic Bandwidth', ...
        'FontWeight','normal');
    
    text(xLim(2)*0.98,yLim(2)*0.95,sprintf('%s%4.2f','ctrMass = ',stimInfo(stim).dynBW_cMass), ...
        'FontSize',8, ...
        'HorizontalAlignment','Right');
    text(xLim(2)*0.98,yLim(2)*0.85,sprintf('%s%4.2f','dist90 = ',stimInfo(stim).dynBW_dist90), ...
        'FontSize',8, ...
        'HorizontalAlignment','Right');
    
    % modulation power spectrum -------------------------------------------
    ax(5) = subplot(3,5,13);
    hold on;
    
    cmap = jet;
    cmap = [0,0,0.5; cmap];
    colormap(ax(5),cmap);
    cMargin = 0.05;
    
    tmf     = stimInfo(stim).tmf;
    smf     = stimInfo(stim).smf;
    mps     = stimInfo(stim).mps;
    mpsLog  = log(mps);
    con     = mps>stimInfo(stim).con;
    
    mpsRel = (mpsLog-min(mpsLog(:)))/range(mpsLog(:));
    imagesc(tmf,smf,mpsRel,[cMargin,1-cMargin]);
    contour(tmf,smf,con,1,'k-', ...
        'LineWidth',2);
    
    axis square;
    set(ax(5),'XLim',[tmf(1),tmf(end)], ...
        'XTick',[-50,0,50], ...
        'XTickLabel',[-50,0,50], ...
        'YLim',[smf(1),smf(end)], ...
        'YTick',[0.05,1,2], ... $(MR: Originally [0,1,2])
        'YTickLabel',[0,1,2]);
    xlabel('{\omega}_{t} (Hz)');
    ylabel('{\Omega}_{s} (cycles/kHz)');
    title('MPS', ...
        'FontWeight','normal');
    
    % temporal modulation curve -------------------------------------------
    ax(6) = subplot(3,5,14);
    hold on;
    
    tmf = stimInfo(stim).tmf;
    tmc = stimInfo(stim).tmc;
    line(tmf,tmc, ...
        'Color','k');
    
    axis square;
    set(ax(6),'XLim',[tmf(1),tmf(end)], ...
        'XTick',[-50,0,50], ...
        'XTickLabel',[-50,0,50], ...
        'YScale','log');
    xlabel('{\omega}_{t} (Hz)');
    ylabel('Rel. Power');
    title('TMC', ...
        'FontWeight','normal');
    
    % spectral modulation curve -------------------------------------------
    ax(7) = subplot(3,5,15);
    hold on;
    
    smf = stimInfo(stim).smf;
    smc = stimInfo(stim).smc;
    line(smc,smf, ...
        'Color','k');
    
    axis square;
    set(ax(7),'XScale','log', ...
        'YLim',[smf(1),smf(end)], ...
        'YTick',(0:0.25:2), ...
        'YTickLabel',(0:0.25:2));
    xlabel('Rel. Power');
    ylabel('{\Omega}_{s} (cycles/kHz)');
    title('SMC', ...
        'FontWeight','normal');
    
    
    
    
    
    
    
    
%     % temporal & spectral modulations -------------------------------------
%     MPS = interp2(X,Y,mps,Xq,Yq,'cubic');
%     
%     upIdx = TMF<0;
%     dnIdx = TMF>0;
%     smIdx = TMF==0;
%     tmIdx = SMF==0;
%     
%     % directionality index
%     dnData = mean(mean(MPS(~tmIdx,dnIdx)));
%     upData = mean(mean(MPS(~tmIdx,upIdx)));
%     dirIdx = (dnData-upData)/(dnData+upData);
%     
%     % temporal modulations
%     tmOffset = 2;
%     tmf     = TMF(smIdx | dnIdx);
%     mTMF    = mean([MPS(:,(smIdx|dnIdx)); fliplr(MPS(:,(upIdx|smIdx)))],1);
%     mTMF    = mTMF/max(mTMF);
%     [pkTMF,tmIdx] = max(mTMF(tmOffset:end));
%     cmTMF   = sum(tmf(tmOffset:end) .* mTMF(tmOffset:end))/sum(mTMF(tmOffset:end)); % center of mass
%     
%     % spectral modulations
%     smOffset = 4;
%     smf    = SMF';
%     mSMF   = mean(MPS,2);
%     mSMF   = mSMF/max(mSMF);
%     [pkSMF,smIdx] = max(mSMF(smOffset:end));
%     cmSMF  = sum(smf(smOffset:end) .* mSMF(smOffset:end))/sum(mSMF(smOffset:end)); % center of mass
%     
%     ax(6) = subplot(3,5,14);
%     hold on;
%     
%     plot(ax(6),repmat(tmf(tmOffset),1,2),[0,1],'k:', ...
%         tmf,mTMF,'k', ...
%         tmIdx,pkTMF+0.1,'r*', ...
%         cmTMF,pkTMF+0.1,'b*');
%     
%     set(ax(6),'XLim',[tmf(1),tmf(end)], ...
%         'XTick',(tmf(1):10:tmf(end)), ...
%         'XTickLabel',(tmf(1):10:tmf(end)), ...
%         'YLim',[0,1], ...
%         'YTick',(0:0.25:1), ...
%         'YTickLabel',(0:0.25:1));
%     xlabel('{\omega}_{t} (Hz)');
%     ylabel('Relative Amplitude');
%     title('TM Power Spectrum');
%     
%     ax(7) = subplot(3,5,15);
%     hold on;
%     
%     plot(ax(7),[0,1],repmat(smf(smOffset),1,2),'k:', ...
%         mSMF,smf,'k', ...
%         pkSMF+0.05,smf(smOffset+smIdx-1),'r*', ...
%         pkSMF+0.05,cmSMF,'b*');
%     
%     set(ax(7),'XLim',[0,1], ...
%         'XTick',(0:0.25:1), ...
%         'XTickLabel',(0:0.25:1), ...
%         'YLim',[smf(1),smf(end)], ...
%         'YTick',(smf(1):0.5:smf(end)), ...
%         'YTickLabel',(smf(1):0.5:smf(end)));
%     xlabel('Relative Amplitude');
%     ylabel('{\Omega}_{t} (cyc/kHz)');
%     title('SM Power Spectrum');
    
     if nFiles > 1
        disp('paused');
%         pause;
        print(fig,'-depsc','-loose','-painters', ... 
            [figDir,'\',stimInfo(stim).name,'.eps']);
        close(fig);
    end
end


function [features] = features_10082013(TS,fs)
%% Feature extraction
handles.siDir = 'C:\Users\Woolley Lab\Documents\MATLAB\Similarity Index';
%input:
% TS - a signal
% fs - sampling frequency
%output: 
% The signal representation in the feature space:
% A matrix with 4 columns: {'Gravity center','Gravity variance','Entropy','Pitch goodness'}

if size(TS,2)>size(TS,1)
    TS = TS';
end

TS = TS-min(TS);
TS = TS/max(TS)*2-1;
TS = TS-mean(TS)+0.005;
if fs==11025
    TS = interp(x,4);
elseif fs==22050
    TS = interp(x,2);
elseif fs==44100
else
    TS = interp1(1:length(TS),TS,1:fs/44100:length(TS));
    TS = TS';
end
fs=44100;

P.fs = fs;
P.pad = 1024;
P.winSize = 409;
P.winStep = 44;
P.spectrum_range = 256;
P.Fmin = 20;
P.Fmax = 200;
P.up_pitch = 13;
P.low_pitch = 55;
P.pitch_HoP = 800;
P.gdn_HoP = 100;
P.up_wiener = -3;
P.cepsNFFT=2^10;
P.minPitchFreq = 400;
P.maxPitchFreq = 6000;
P.totalPowerRange = [300,min(10000,fs/2)];
P.cepsUpperLimit = 6000;
P.PeakDetect = true;

load([handles.siDir,'\tapers.mat']);    %loads Tapers
N = length(TS);

TSM = windowing(TS', P.winSize, P.winStep);
S = 0;
SF = 0;

J1 = (fft(TSM(:,:).*(ones(size(TSM,1),1)*(Tapers(:,1))'),P.pad,2));
J1 = J1(:,1:P.spectrum_range);
J2 = (fft(TSM(:,:).*(ones(size(TSM,1),1)*(Tapers(:,2))'),P.pad,2));
J2 = J2(:,1:P.spectrum_range);

P_Spec = real(J1).^2+real(J2).^2+imag(J1).^2+imag(J1).^2;
f = fs/P.pad*(1:256)';

%==================Gravity center and variance==============
freq_index = [P.Fmin:P.Fmax];
m_power = (real(J1(:,freq_index)).^2+imag(J1(:,freq_index)).^2).*(real(J2(:,freq_index)).^2+imag(J2(:,freq_index)).^2);
[m,n] = size(m_power);
gravityINDEX = (freq_index*m_power')./sum(m_power');
gravity_center = gravityINDEX*fs/P.pad;
gravity_variance = sqrt(sum((repmat(freq_index,m,1)-repmat(gravityINDEX',1,n))'.^2.*m_power')./sum(m_power'))*fs/P.pad;

%===========Wiener entropy===================================
sumlog = sum(log(P_Spec(:,freq_index)+eps),2);
logsum = (sum(P_Spec(:,freq_index),2));

logsum(find(logsum==0)) = length(freq_index);
logsum = log(logsum/length(freq_index)); %divide by the number of frequencies
entropy = (sumlog/length(freq_index))-logsum;
entropy(find(logsum==0)) = 0;

P_Spec = P_Spec';
logp = log(P_Spec);

%========== Pitch goodness===============================
cepstrum1 = abs(fft(logp))';
z = cepstrum1(:,P.up_pitch:P.low_pitch);
[pitchGoodness,m_Pitch_xxx] = sort(z,2);
pitchGoodness = pitchGoodness(:,end);
features = [gravity_center', gravity_variance', entropy, pitchGoodness];

function [windowed_data] = windowing(data,window_size,step_size)
%% from vector to matrix according to window_size. 
%step_size determines how much the windows overlap
windowed_data = [];
size_data = size(data,2);
for start_window = 1:step_size:(size_data-window_size+1);
    windowed_data = [windowed_data; data(start_window:start_window+window_size-1)];
end
clear data;
