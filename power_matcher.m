function equate_power
% Equate the RMS power of all wave files w/in chosen directory
% Use check_power.m to verify the power before and after transformation.
% dms 5/15/12
% updated 4/4/2018 jae - fixed deprecated function calls

%% Set parameters for pass-filtering (Hz)
loBound = 500;
hiBound = 16000;

%% Set the desired SPL:
SPLgoal = 65;

%% DO NOT CHANGE:
% Set the RMS value at a particular intensity (dB SPL)
SPLbase = 60;
RMSbase = 0.01378;

%% Adjust the desired RMS value based on desired SPL
deltaSPL = SPLgoal - SPLbase;
RMSgoal = RMSbase * 10^(deltaSPL/20);

%% Choose a directory
d = uigetdir;
directory = dir(d);
mkdir(d,'MatchedStims');
%savedir = [d,'\MatchedStims'];
savedir = [d,'/MatchedStims']; % mac

%% Loop through each wave file and adjust the RMS
for i = 1:length(directory) % loop through all items in dir
    if endsWith(directory(i).name, '.wav') % make sure it is a .wav file
        filename = directory(i).name;

        % Load the file
        cd(d)
        [y , fs] = audioread(filename);
        
        % Calculate the RMS power
        yRMS = sqrt(mean(y.*y));

        % Scale to the appropriate level
        yN = y*RMSgoal/yRMS;
        
        % Band-pass filter
        %yN = bandpass(yN, [loBound, hiBound], fs);

        % Save the power-matched signal
        cd(savedir)
        filenameSave = filename(1:end-4);  %chop off the .wav
        audiowrite([filenameSave,'_powerMatched.wav'],yN,fs);
    end
end
disp('Done power matching!')