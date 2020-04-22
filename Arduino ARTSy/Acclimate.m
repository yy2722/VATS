%% Acclimate 02/26/2019 by Yow-Tyng (Tim) Yeh
%This script is for the acclimation stage prior to conditional procedures that involve a feeder.
%A feeder will present reward (food) for a defined duration of time at random intervals within a range that is set by user.

%New instances of matlab is required to run additional booths
%% Set up arduino

booth = 2;

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

%Initialize arduino
a = arduino(arduino_device);

%Define solenoid pin # and configure the pin to digital output
SOLENOID_PIN = 'D46';
configurePin(a, SOLENOID_PIN, 'digitalOutput');

%Set reward duration
reward_time = 10;

%Set ITI range (min and max)  %if the range is too small, might cause the
%bird to be stressful
ITI_min = 30;
ITI_max = 60;

%% Actual experiment
%The feeder moves after a random time interval and has the reward present for a defined time duration
while true
    writeDigitalPin(a, SOLENOID_PIN, 0);
    interval = randi([ITI_min,ITI_max],1,1);
    disp(interval);
    pause(interval);
    
    
    writeDigitalPin(a, SOLENOID_PIN, 1);
    pause(reward_time);
end