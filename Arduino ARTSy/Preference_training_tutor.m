%% Preference 02/26/2019 by Yow-Tyng (Tim) Yeh
% Modified on 05/10/2019 to replace PsychAudio with Matlab sound player

% The bird will elicit song playback by triggering the IR sensor. If the
% playback has positive valence, the bird will continue to trigger the sensor.

% New instances of matlab is required to run additional booths

%% Enter the information about the bird
birdID = 'xx0000';
sex = 'F';
stimulus_L = 'heterospecific songs'; % Specify the type of playback used on the left sensor
stimulus_R = 'conspecific songs';    % Specifiy the type of playback used on the right sensor

%% Set up audio parameters and channels
booth = 1 ;   % Set booth number to indicate the booth that is going to run the experiment

switch booth
    case 1
        speaker = 'Playback 1-4 (AudioFire 4) (Windows DirectSound)';
        tutor = 'b1';
        audioColumn = 1;
    case 2
        speaker = 'Playback 1-4 (AudioFire 4) (Windows DirectSound)';
        tutor = 'b2';
        audioColumn = 2;
    case 3
        speaker = 'Playback 5-6 (AudioFire 4) (Windows DirectSound)';
        tutor = 'b3';
        audioColumn = 1;
    case 4
        speaker = 'Playback 5-6 (AudioFire 4) (Windows DirectSound)';
        tutor = 'b4';
        audioColumn = 2;
    case 5
        speaker = 'Playback 1-4 (2- AudioFire 4) (Windows DirectSound)';
        tutor = 'b5';
        audioColumn = 1;
    case 6
        speaker = 'Playback 1-4 (2- AudioFire 4) (Windows DirectSound)';
        tutor = 'b6';
        audioColumn = 2;
    case 7
        speaker = 'Playback 5-6 (2- AudioFire 4) (Windows DirectSound)';
        tutor = 'b7';
        audioColumn = 1;
    case 8
        speaker = 'Playback 5-6 (2- AudioFire 4) (Windows DirectSound)'; 
        tutor = 'b8';
        audioColumn = 2;
end

%% Create a list of random playbacks
right_songs = dir(['C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\preference_songs\tutor\', tutor]);
left_songs = dir(['C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\preference_songs\tutor\', tutor]);

fs = 44100;    % Audio files sample rate 44.1 kHz
songs_cell_R = transpose({right_songs(3:end).name});  %get all the filenames from right_IR folder and create a cell array
songs_cell_L = transpose({left_songs(3:end).name});  %get all the filenames from left_IR folder and create a cell array

for i=1:(length(right_songs)-2)
    
    path_R = fullfile(['C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\preference_songs\tutor\', tutor], right_songs(i+2).name);
    info_R = audioinfo(path_R);
    songs_cell_R(i,2) = {info_R.Duration};    %song duration, column 2 of the array
    songs_cell_R(i,3) = {audioread(path_R)};           %audio signal, column 3 of the array
    
end

for i=1:(length(left_songs)-2)
    
    path_L = fullfile(['C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\preference_songs\tutor\', tutor], left_songs(i+2).name);
    info_L = audioinfo(path_L);
    songs_cell_L(i,2) = {info_L.Duration};    %song duration, column 2 of the array
    songs_cell_L(i,3) = {audioread(path_L)};           %audio signal, column 3 of the array
    
end

block_size = int32(5000);   %define block size

[nr_R, nc_R] = size(songs_cell_R);  %get the dimension of the array
[nr_L, nc_L] = size(songs_cell_L);

numRepeat_R = idivide(block_size, nr_R, 'ceil');      %find out the number of repetitions to reach block size
numRepeat_L = idivide(block_size, nr_L, 'ceil');

repeat_songs_R = repmat(songs_cell_R,numRepeat_R,1);
repeat_songs_L = repmat(songs_cell_L,numRepeat_L,1);

repeat_random_R = repeat_songs_R(randperm(length(repeat_songs_R)),:);   %shuffle rows
repeat_random_L = repeat_songs_L(randperm(length(repeat_songs_L)),:);

song_list_R = repeat_random_R(1:block_size,:);          %take only a defined number of rows to create a block of song list
song_list_L = repeat_random_L(1:block_size,:); 

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

a = arduino(arduino_device);   % Connect to arduino

SENSORPIN = 'D45';          % Sensor
SENSORPIN2 = 'D10';

LEDPIN_R = 'D7';   % LED lights
LEDPIN_L = 'D5';

% Variables will change:
sensorState1 = 0;
lastState1 = 0;        % Variable for reading the pushbutton status
sensorState2 = 0;
lastState2 = 0;

% Initialize the sensor pin as an input:
configurePin(a, SENSORPIN, 'digitalInput');
configurePin(a, SENSORPIN, 'pullup');
configurePin(a, SENSORPIN2, 'digitalInput');
configurePin(a, SENSORPIN2, 'pullup');

%initialize the LED pin as an output:
configurePin(a, LEDPIN_L, 'PWM');
configurePin(a, LEDPIN_R, 'PWM');

% Counter for right and left sensors
counter_L = 0;
counter_R = 0;
counter = 0; % Trial counter

% Variables to determine whether to inactivate the sensor or not (for
% training purpose); the sensor will be inactivated if L_repeat = 3, 
%L_repeat = 0;
%R_repeat = 0;

% Create data table
data = cell(3, 2);
data(1,:) = {'Left_IR', 'Right_IR'};
data(2,:) = {stimulus_L, stimulus_R};

while counter < block_size
  
          % Read the state of the pushbutton value:
          sensorState1 = readDigitalPin(a, SENSORPIN);
          sensorState2 = readDigitalPin(a, SENSORPIN2);
          
          writePWMVoltage(a, LEDPIN_L, 0.2);
          writePWMVoltage(a, LEDPIN_R, 0.2);
          %if L_repeat < 3
          %writePWMVoltage(a, LEDPIN_L, 0.2);
          %else
          %writePWMVoltage(a, LEDPIN_L, 0);    
          %end 
          
          %if R_repeat < 3
          %writePWMVoltage(a, LEDPIN_R, 0.2);
          %else
          %writePWMVoltage(a, LEDPIN_R, 0);    
          %end 
          
          if (~sensorState1 && lastState1) %&& (L_repeat < 3)  % If not false(1) && True(1)
            disp('Left sensor triggered ---------------------');            %sensorState = false (0) if touched

            audioMatrix = zeros(length(song_list_L{counter_L+1, 3}), 2);    % Create a zero matrix that is the length of audio signal x 2(L and R channels)
            audioMatrix(:, audioColumn) = song_list_L{counter_L+1, 3};      % Load the audio signal into the column based on selected booth number
            player = audioplayer(audioMatrix, fs, 24, booth_speaker); %#ok<TNMLP>
            play(player);    %play audio

            pause(song_list_L{counter_L+1, 2});   %time-out

            counter = counter +1;
            counter_L = counter_L +1;

            %L_repeat = L_repeat +1;  % Left sensor triggered counter +1
            %if R_repeat > 0 % If right sensor triggered counter greater than 0, -1 right repeat counter
            %    R_repeat = R_repeat -1;
            %end

            current_time = clock;   %show current time

            disp(['song playback at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
            '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);

            disp(song_list_L{counter_L+1, 1});   %display the song that is played

            playbackDuration = song_list_L{counter_L+1, 2};
            disp(['song duration: ', num2str(playbackDuration), 'sec']);   %display song duration
            disp(['Left sensor ', num2str(counter_L), ' Right sensor ', num2str(counter_R)]); % Count how many times each sensor got triggered
          end

          if (~sensorState2 && lastState2) %&& (R_repeat < 3)  % If not false(1) && True(1)
            disp('Right sensor triggered ---------------------');            % SensorState = false (0) if touched

            audioMatrix = zeros(length(song_list_R{counter_R+1, 3}), 2);    % Create a zero matrix that is the length of audio signal x 2(L and R channels)
            audioMatrix(:, audioColumn) = song_list_R{counter_R+1, 3};      % Load the audio signal into the column based on selected booth number
            player = audioplayer(audioMatrix, fs, 24, booth_speaker); %#ok<TNMLP>
            play(player);    % Play audio

            pause(song_list_R{counter+1, 2});   % Time-out

            counter = counter +1;
            counter_R = counter_R +1;

            %R_repeat = R_repeat +1;  % Left sensor triggered counter +1
           % if L_repeat > 0 % If right sensor triggered counter greater than 0, -1 right repeat counter
            %    L_repeat = L_repeat -1;
            %end

            current_time = clock;   %show current time

            disp(['song playback at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
            '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);

            disp(song_list_R{counter_R+1, 1});   %display the song that is played

            playbackDuration = song_list_R{counter_R+1, 2};
            disp(['song duration: ', num2str(playbackDuration), 'sec']);   %display song duration

            disp(['Left sensor ', num2str(counter_L), ' Right sensor ', num2str(counter_R)]); % Count how many times each sensor got triggered
          end

          lastState1 = sensorState1;
          lastState2 = sensorState2;


end
 
