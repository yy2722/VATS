%% Housing 02/26/2019 by Yow-Tyng (Tim) Yeh
%This script is for housing the bird when not running experiment. 
%Arduino will control the onset (7 am)and offset (9 pm) of the housing light.

%New instances of matlab is required to run additional booths

%% Set up arduino
booth =2;

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

a = arduino(arduino_device);

configurePin(a, 'D44', 'digitalOutput');

%% Actual experiment
while true
    
    now = fix(clock);
    time = now(4);
    
    if time > (7-1) && time < 21  %if time greater than 6 and smaller than 21, light is on
        writeDigitalPin(a, 'D44', 1);
        %disp('light is on');
        
    else 
        writeDigitalPin(a, 'D44', 0);
        %disp('light is off');
    
    end
end