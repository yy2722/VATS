%% Shaping-Phase 3 by Yow-Tyng (Tim) Yeh 07/21/2019
% This script is for training bird to touch left sensor for trial
% initiation and then touch right sensor within a defined time period to get food
% reward. This training is prior to audiotor selectivity test.

booth = 1;   %set booth number to indicate the booth that is going to run the experiment

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

sensor_PIN = 'D45';          %sensor pin on right
sensor_PIN2 = 'D10';    

SOLENOID_PIN = 'D46';         %solenoid pin

LED_L_PIN = 'D5';             %LED pin on right
LED_R_PIN = 'D7';  

configurePin(a, sensor_PIN, 'digitalInput');  %set sensor pin to digital input
configurePin(a, sensor_PIN, 'pullup');        %pull up
configurePin(a, sensor_PIN2, 'digitalInput');  %set sensor pin to digital input
configurePin(a, sensor_PIN2, 'pullup');        %pull up

configurePin(a, SOLENOID_PIN, 'digitalOutput'); %set solenoid pin to digital output

configurePin(a, LED_L_PIN, 'PWM');
configurePin(a, LED_R_PIN, 'PWM');

%% Variables for calculations (DO NOT CHANGE)

%Variables for calculating the state of IR sensor
sensorState = 0;
lastState = 0;
sensorState2 = 0;

%Counter for the number of trials
counter = 0;

%% Set up experiment parameters

% Set pre-response duration
preResponse = 1;

% Set post-stim response duration- after left sensor triggered, the amount of time for response
% 1 sec -> 4 sec
postStim = 4;
postStim_loop = postStim/0.015;  %0.015 sec per loop (to avoid complex calculations)

% Set reward duration
reward_time = 3;

% Set null response duration
% -if the bird does not trigger right IR sensor after triggering left sensor, null response kicks in
null_response = 3;

%% Actual experiment
%The experiment will stop after running through the song list
while true
    
    writePWMVoltage(a, LED_L_PIN, 0.2);
    sensorState = readDigitalPin(a, sensor_PIN);
    
    %the IR beam is broken and a song playback is triggered by the bird
    if (~sensorState && lastState)  %if not false(1) && True(1)
        
        disp(['trial ', num2str(counter+1)]);     %show the number of times in song playback
        
        counter = counter+1;
        waiting_response = 1;
        
        loop = 0;
        
        writePWMVoltage(a, LED_L_PIN, 0);
        %waiting for resposne after song playback
        while(waiting_response == 1)
            if loop == 0
                pause(preResponse);   %first time entering the loop results in a pre-response pause
            end
            
            writePWMVoltage(a, LED_R_PIN, 0.2);
            sensorState2 = readDigitalPin(a, sensor_PIN2);
        
            if (sensorState2) && (loop > postStim_loop)
                writePWMVoltage(a, LED_R_PIN, 0);
                
                current_time = clock;   %show current time
            
                pause(null_response);   %null response duration
                
                disp(['- missed at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
                
                waiting_response = 0;  %reset to exit waiting response loop
            
            elseif (~sensorState2)
                writePWMVoltage(a, LED_R_PIN, 0);
                
                current_time = clock;
            
                writeDigitalPin(a, SOLENOID_PIN, 1);  %the feeder is out for a period of time defined by reward time
                pause(reward_time);
            
                writeDigitalPin(a, SOLENOID_PIN, 0);
                
                disp(['- triggered at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
                '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
               
                waiting_response = 0;  %reset
            
            end
           
            loop = loop+1;  %this is for calculating post-stim time
        end
    end    
       
    lastState = sensorState;   %to avoid continuous activation of IR sensor when no state change
    
end
