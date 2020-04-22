%% Song Training 02/26/2019 by Yow-Tyng (Tim) Yeh
% Modified on 05/08/2019 to replace PsychPortAudio with matlab sound player
% Modified on 05/15/2019 to automatically update speaker number when the computer restarts

% This script is for training young birds to sing (song development). 
% The bird triggers song playback by breaking the IR beam (on the left). The liquid crystal panel becomes transparent and a decoy is revealed at the beginning of the experiment. 
% Total duration of song playback per day is capped at 28 seconds.

% New instances of matlab is required to run additional booths

%% Create a list of random playbacks
response_songs = dir('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\Response2sounds');

fs = 44100;    %audio files sample rate 44.1 kHz
songs_cell = transpose({response_songs(3:end).name});  %get all the filenames and create a cell array

for i=1:(length(response_songs)-2)
    
    path = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Testing_songs\Response2sounds', response_songs(i+2).name);
    info = audioinfo(path);
    songs_cell(i,2) = {info.Duration};    %song duration, column 2 of the array
    songs_cell(i,3) = {audioread(path)};           %audio signal, column 3 of the array
    
end

block_size = int32(9);   %define block size
[nr, nc] = size(songs_cell);  %get the dimension of the array
numRepeat = idivide(block_size, nr, 'ceil');      %find out the number of repetitions to reach block size
repeat_songs = repmat(songs_cell,numRepeat,1);
repeat_random = repeat_songs(randperm(length(repeat_songs)),:);   %shuffle rows
song_list = repeat_random(1:block_size,:);          %take only a defined number of rows to create a block of song list

%% Set up audio parameters and channels
booth = 7;   %set booth number to indicate the booth that is going to run the experiment

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

total_playbackTime = 0;
counter = 0;

%% Actual experiment

%writeDigitalPin(a, LQ_PIN, 0);   %make LQ panel transparent to reveal the decoy at the beginning of the experiment

while counter < 9
        if (counter == 0) || (counter == 3) || (counter == 6)
            pause(60*5);
            %pause(5);
        end
        audioMatrix = zeros(length(song_list{counter+1, 3}), 2);    %create a zero matrix that is the length of audio signal x 2(L and R channels)
        audioMatrix(:, audioColumn) = song_list{counter+1, 3};      %load the audio signal into the column based on selected booth number
        player = audioplayer(audioMatrix, fs, 24, booth_speaker); %#ok<TNMLP>
        play(player);    %play audio
        
        pause(song_list{counter+1, 2});   %time-out
        pause(1);
        
        current_time = clock;   %show current time
    
        disp(['song playback at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
        '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm ',num2str(current_time(6)), 's']);
        disp(song_list{counter+1, 1});   %display the song that is played
        
        playbackDuration = song_list{counter+1, 2};
        disp(['song duration: ', num2str(playbackDuration), 'sec']);   %display song duration
        
        total_playbackTime = total_playbackTime + playbackDuration;
        disp(['total playback time: ', num2str(total_playbackTime), 'sec']);   %show the total duration of playback
        
        counter = counter+1;
        disp(['playback ', num2str(counter), '-------------------------']); %show the number of times in song playback
        
end
