%% Song Training 02/26/2019 by Yow-Tyng (Tim) Yeh
% Modified on 05/08/2019 to replace PsychPortAudio with matlab sound player
% Modified on 05/15/2019 to automatically update speaker number when the computer restarts

% This script is for training young birds to sing (song development). 
% The bird triggers song playback by breaking the IR beam (on the left). The liquid crystal panel becomes transparent and a decoy is revealed at the beginning of the experiment. 
% Total duration of song playback per day is capped at 28 seconds.

% New instances of matlab is required to run additional booths

%% Create a list of random playbacks
training_songs = dir('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\song_training');

fs = 44100;    %audio files sample rate 44.1 kHz
songs_cell = transpose({training_songs(3:end).name});  %get all the filenames and create a cell array

for i=1:(length(training_songs)-2)
    
    path = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\song_training', training_songs(i+2).name);
    info = audioinfo(path);
    songs_cell(i,2) = {info.Duration};    %song duration, column 2 of the array
    songs_cell(i,3) = {audioread(path)};           %audio signal, column 3 of the array
    
end

block_size = int32(60);   %define block size
[nr, nc] = size(songs_cell);  %get the dimension of the array
numRepeat = idivide(block_size, nr, 'ceil');      %find out the number of repetitions to reach block size
repeat_songs = repmat(songs_cell,numRepeat,1);
repeat_random = repeat_songs(randperm(length(repeat_songs)),:);   %shuffle rows
song_list = repeat_random(1:block_size,:);          %take only a defined number of rows to create a block of song list

%% Set up audio parameters and channels
booth = 1;   %set booth number to indicate the booth that is going to run the experiment

switch booth
    case 1
        speaker = 'Playback 1-4 (AudioFire 4) (Windows DirectSound)';
        audioColumn = 1;
    case 2
        speaker = 'Playback 1-4 (AudioFire 4) (Windows DirectSound)';
        audioColumn = 2;
    case 3
        speaker = 'Playback 5-6 (AudioFire 4) (Windows DirectSound)';
        audioColumn = 1;
    case 4
        speaker = 'Playback 5-6 (AudioFire 4) (Windows DirectSound)';
        audioColumn = 2;
    case 5
        speaker = 'Playback 1-4 (2- AudioFire 4) (Windows DirectSound)';
        audioColumn = 1;
    case 6
        speaker = 'Playback 1-4 (2- AudioFire 4) (Windows DirectSound)';
        audioColumn = 2;
    case 7
        speaker = 'Playback 5-6 (2- AudioFire 4) (Windows DirectSound)';
        audioColumn = 1;
    case 8
        speaker = 'Playback 5-6 (2- AudioFire 4) (Windows DirectSound)'; 
        audioColumn = 2;
end

%% Find the correct speaker
% This avoids the problem of reassigning speaker number every time when the computer restarts

audioinfo = audiodevinfo;  % Get audioinfor on speakers
audioInfoMatrix = transpose({audioinfo.output.Name});  %Create a cell array with a list of speakers
audioInfoMatrix(:, 2) = transpose({audioinfo.output.ID});
TF_vector = contains(audioInfoMatrix(:,1), speaker);  %Match the correct speaker according to booth number
matched_row = find(TF_vector);
booth_speaker =  audioInfoMatrix{matched_row, 2}; %Get speaker number

%% Set up arduino
% Connect to the arduino corresponding to the booth number
switch booth
    case 1
        arduino_device = 'COM3';
        disp('booth 1 running');
    case 2
        arduino_device = 'COM4';
        disp('booth 2 running');
    case 3
        arduino_device = 'COM5';
        disp('booth 3 running');
    case 4
        arduino_device = 'COM7';
        disp('booth 4 running');
    case 5
        arduino_device = 'COM8';
        disp('booth 5 running');
    case 6
        arduino_device = 'COM9';
        disp('booth 6 running');
    case 7
        arduino_device = 'COM10';
        disp('booth 7 running');
    case 8
        arduino_device = 'COM6';
        disp('booth 8 running');
        
end

a = arduino(arduino_device);   %connect to arduino

LIGHTPIN = 'D44';          % Light pin

sensor_L = 'D45';          % Sensors
sensor_R = 'D10';

LEDPIN_L = 'D5';           % LED lights
LEDPIN_R = 'D7';

configurePin(a, sensor_L, 'digitalInput');  % Set sensor pin to digital input
configurePin(a, sensor_L, 'pullup');        % Pull up
configurePin(a, sensor_R, 'digitalInput');  % Set sensor pin to digital input
configurePin(a, sensor_R, 'pullup');

configurePin(a, LIGHTPIN, 'digitalOutput');  % Light

configurePin(a, LEDPIN_L, 'PWM'); % LED lights
configurePin(a, LEDPIN_R, 'PWM');  

sensorState = 0;
lastState = 0;

sensorState2 = 0;
lastState2 = 0;

total_playbackTime = 0;
counter = 0;
Done = 0;

%% Actual experiment
while true
    
    now = fix(clock);
    time = now(4);
    
    if time > (7-1) && time < 21 && (Done == 0)  %if time greater than 6 and smaller than 21 and the training session has not started, light is on and the training session starts
        writeDigitalPin(a, LIGHTPIN, 1);
        %disp('light is on');
        
        while total_playbackTime < 28    %after reaching the maximum time duration of 28 sec (Ofer's papers) for song playback per day, the training session stops

            writePWMVoltage(a, LEDPIN_R, 0.2);
            writePWMVoltage(a, LEDPIN_L, 0.2);

            sensorState = readDigitalPin(a, sensor_L);
            sensorState2 = readDigitalPin(a, sensor_R);

            if (sensorState && ~lastState) || (sensorState2 && ~lastState2)  %if True(1) && not false(1)
                %fprintf('Unbroken\n');          %lastState = false (0)

            end

            if (~sensorState && lastState) || (~sensorState2 && lastState2)  %if not false(1) && True(1)
                %fprintf('Broken\n');            %sensorState = false (0) if touched

                audioMatrix = zeros(length(song_list{counter+1, 3}), 2);    %create a zero matrix that is the length of audio signal x 2(L and R channels)
                audioMatrix(:, audioColumn) = song_list{counter+1, 3};      %load the audio signal into the column based on selected booth number
                player = audioplayer(audioMatrix, fs, 24, booth_speaker); %#ok<TNMLP>
                play(player);    %play audio

                pause(song_list{counter+1, 2});   %time-out

                current_time = clock;   %show current time

                disp(['song playback at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                disp(song_list{counter+1, 1});   %display the song that is played

                playbackDuration = song_list{counter+1, 2};
                disp(['song duration: ', num2str(playbackDuration), 'sec']);   %display song duration

                total_playbackTime = total_playbackTime + playbackDuration;
                disp(['total playback time: ', num2str(total_playbackTime), 'sec']);   %show the total duration of playback

                counter = counter+1;
                disp(['playback ', num2str(counter), '-------------------------']); %show the number of times in song playback
                %pause; % For initial training
            end

            lastState = sensorState;
            lastState2 = sensorState2;

        end
        Done = 1;
        
        total_playbackTime = 0;
        counter = 0;
        
        disp('Done.....you have reached the maximum number of song playback!');
        
    elseif time > (7-1) && time < 21 && (Done == 1)   % After the training session is finished during daytime
        writeDigitalPin(a, LIGHTPIN, 1);
        writePWMVoltage(a, LEDPIN_R, 0);
        writePWMVoltage(a, LEDPIN_L, 0);
        
    else
        writeDigitalPin(a, LIGHTPIN, 0);              % After the light goes off, training session resets
        writePWMVoltage(a, LEDPIN_R, 0);
        writePWMVoltage(a, LEDPIN_L, 0);
        Done = 0;
        
        %disp('light is off');
    
    end
end
