%% Go/Nogo 03/04/2019 by Yow-Tyng (Tim) Yeh
% Modified on 05/08/2019 to replace PsychPortAudio with matlab sound player
% Modified on 05/15/2019 to automatically update speaker number when the computer restarts

%The bird first triggers the IR sensor to elicit a song playback that is randomly
%selected. Depending on the playback, the bird has to decide whether to
%touch the sensor again (Go) or not touch the sensor (NoGo). 

%If the bird responds to the Go playback correctly, it gets food (reward) for a period
%of time. If the bird responds to the NoGo playback correctly, it avoids
%the light being shut off (punishment).

%% Enter the information about the bird
birdID = 'xx0000';
sex = 'F';

%% Create a list of random playbacks from go and nogo folders
go_songs = dir('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\Go');
nogo_songs = dir('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\NoGo');

fs = 44100;    %audio files sample rate 44.1 kHz

%get the file names and make them into cell arrays
go_cell = transpose({go_songs(3:end).name});
go_cell(:,2) = {'go'};     %assign the song as go or nogo, column 2 of the go cell array

for i=1:(length(go_songs)-2)
    
    path = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\Go', go_songs(i+2).name);
    info = audioinfo(path);
    go_cell(i,3) = {info.Duration};    %song duration, column 3 of the go cell array
    go_cell(i,4) = {audioread(path)};           %audio signal, column 4 of the go cell array
    
end

nogo_cell = transpose({nogo_songs(3:end).name});
nogo_cell(:,2) = {'nogo'};      %assign the song as go or nogo, column 2 of the nogo cell array

for i=1:(length(nogo_songs)-2)
    
    path = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\NoGo', nogo_songs(i+2).name);
    info = audioinfo(path);
    nogo_cell(i,3) = {info.Duration};    %song duration, column 3 of the nogo cell array
    nogo_cell(i,4) = {audioread(path)};           %audio signal, column 4 of the nogo cell array
    
end

go_nogo = [go_cell; nogo_cell];  %concatenate go and nogo arrays into one

session_size = int32(400);   %define the size of one session
[nr, nc] = size(go_nogo);  %get the dimension of the array
numRepeat = idivide(session_size, nr, 'ceil');      %find out the number of repetitions to reach session size
repeat_songs = repmat(go_nogo,numRepeat,1);
repeat_random = repeat_songs(randperm(length(repeat_songs)),:);   %shuffle rows
song_list = repeat_random(1:session_size,:);          %take only a defined number of rows to create a block of song list

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
%connect to the arduino corresponding to the booth number
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

light_PIN = 'D44';              %light pin
sensor_PIN = 'D45';          %sensor pin
SOLENOID_PIN = 'D46';         %solenoid pin

configurePin(a, sensor_PIN, 'digitalInput');  %set sensor pin to digital input
configurePin(a, sensor_PIN, 'pullup');        %pull up
configurePin(a, SOLENOID_PIN, 'digitalOutput'); %set solenoid pin to digital output
configurePin(a, light_PIN, 'digitalOutput');     %set light pin to digital output

writeDigitalPin(a, light_PIN, 1);   %turn up the light becuase the light turns off automatically when the script initializes

%% Variables for calculations (DO NOT CHANGE)

%Variables for calculating the state of IR sensor
sensorState = 0;
lastState = 0;

%Counter for the number of trials
counter = 0;

%For shifting the focus from constantly pecking the sensor to listening the
%cues
repeat = 0;
repeat_counter = 0;

%Cell array for data plotting
data = cell(session_size, 5);
data(1,:) = {'Hit', 'Miss', 'CR', 'FA', 'go/nogo'};
accuracy_counter = [0, 0, 0, 0];

%% Set up experiment parameters

%Set pre-response duration- the amount of time before IR sensor can be
%triggered again for response
preResponse = 1;

%Set post-stim response duration- after song playback, the amount of time
%for response
%was 4 -> 3 -> 2
postStim = 1;
postStim_loop = postStim/0.015;  %0.015 sec per loop (to avoid complex calculations)

%Set reward duration
%was 2
reward_time = 3;

%Set punishment duration
punishment_time = 20;

%Set null response duration- if the bird does not trigger IR sensor again after song
%playback, null response kicks in and here is the time duration
null_response = 3;

%% Actual experiment
%The experiment will stop after running through the song list
while (counter < session_size)
    
    sensorState = readDigitalPin(a, sensor_PIN);
    
    if (sensorState && ~lastState)  %if True(1) && not false(1)
        %fprintf('Unbroken\n');          %lastState = false (0)
        
    end
    
    %the IR beam is broken and a song playback is triggered by the bird
    if (~sensorState && lastState)  %if not false(1) && True(1)
        %fprintf('Broken\n');            %sensorState = false (0) if touched
        
        if (repeat ==3)    %repeat the wrong trial if repeat == 1 % this setting is only for training
            counter = counter -1;
            %repeat_counter = repeat_counter +1;
        end
        
        audioMatrix = zeros(length(song_list{counter+1, 4}), 2);    %create a zero matrix that is the length of audio signal x 2(L and R channels)
        audioMatrix(:, audioColumn) = song_list{counter+1, 4};      %load the audio signal into the column based on selected booth number
        player = audioplayer(audioMatrix, fs, 24, booth_speaker); %#ok<TNMLP>
        play(player);    %play audio
        
        disp(['trial ', num2str(counter+1),'---------------------']);     %show the number of times in song playback
        
        current_time = clock;   %show current time
        disp(['song playback at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
        '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
        
        disp(song_list{counter+1, 1});   %display the song that is played
        disp(song_list{counter+1, 2});   %display whether the song is go or nogo
        
        playbackDuration = song_list{counter+1, 3};
        disp(['song duration: ', num2str(playbackDuration), 'sec']);   %display song duration
        
        which_one = song_list{counter+1, 2};   %for later string matching in waiting response loop
        
        counter = counter+1;
        waiting_response = 1;
        
        loop = 0;
        
        %waiting for resposne after song playback
        while(waiting_response == 1)
            if loop == 0
                pause(preResponse);   %first time entering the loop results in a pre-response pause
            end
          
            sensorState = readDigitalPin(a, sensor_PIN);
        
            if (sensorState) && (strcmp(which_one, 'go')) && (loop > postStim_loop) 
                current_time = clock;   %show current time
            
                pause(null_response);   %null response duration
                
                accuracy_counter(1, 2) = accuracy_counter(1, 2)+1;
                
                disp(['-> Miss at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                disp(accuracy_counter);
                
                data(counter+1,:) = {0, 1, 0, 0, 'go'};
                repeat = 1;
                waiting_response = 0;  %reset to exit waiting response loop
            
            elseif (~sensorState) && (strcmp(which_one, 'go'))
                current_time = clock;
            
                writeDigitalPin(a, SOLENOID_PIN, 1);  %the feeder is out for a period of time defined by reward time
                pause(reward_time);
            
                writeDigitalPin(a, SOLENOID_PIN, 0);
                
                accuracy_counter(1, 1) = accuracy_counter(1, 1)+1;
                
                disp(['-> Hit at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                disp(accuracy_counter);
            
                data(counter+1,:) = {1, 0, 0, 0, 'go'};
                repeat = 0;
                waiting_response = 0;  %reset
            
            elseif (sensorState) && (strcmp(which_one, 'nogo')) && (loop > postStim_loop) 
                current_time = clock;
            
                pause(null_response);
                
                accuracy_counter(1, 3) = accuracy_counter(1, 3)+1;
                
                disp(['-> Correct response at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                disp(accuracy_counter);
                
                data(counter+1,:) = {0, 0, 1, 0, 'nogo'};
                repeat = 0;
                waiting_response = 0; %reset
            
            elseif (~sensorState) && (strcmp(which_one, 'nogo'))
                current_time = clock;
            
                disp('punishment');
                writeDigitalPin(a, light_PIN, 0);       %turn off the light for a time duration defined by punishment time
                pause(punishment_time);
            
                writeDigitalPin(a, light_PIN, 1);
    
                accuracy_counter(1, 4) = accuracy_counter(1, 4)+1;
                
                disp(['-> False alarm at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                disp(accuracy_counter);
            
                data(counter+1,:) = {0, 0, 0, 1, 'nogo'};
                repeat = 1;
                waiting_response = 0;  %reset
            
            end
           
            loop = loop+1;  %this is for calculating post-stim time
        end
        
    end    
       
    lastState = sensorState;   %to avoid continuous activation of IR sensor when no state change
    
end
disp('end of the experiment!');

%caculate accuracy and error of this session
Hit = double(sum([data{2:end,1}]));
Miss = double(sum([data{2:end,2}]));
CR = double(sum([data{2:end,3}]));
FA = double(sum([data{2:end,4}]));
total_trials = double(session_size);

accuracy = ((Hit + CR)/total_trials)*100;
error = ((Miss + FA)/total_trials)*100;

disp(['Accuracy is ', num2str(accuracy), '%']);
disp(['Error rate is ', num2str(error), '%']);

%% Save experimental data
%save the data as an excel file
data = cell2table(data);
filename = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Data',['go_nogo_data_booth', num2str(booth),'_', birdID, '_', sex, '_', num2str(current_time(1))...
    , '_',num2str(current_time(2)),'_',num2str(current_time(3)), '.xlsx']);
writetable(data, filename);