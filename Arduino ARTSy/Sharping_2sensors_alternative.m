%This script is for the shaping phase. The bird will trigger the IR sensor
%to get reward (food). Once the bird learns the association between
%triggering the sensor and reward, proceed the next phase immediately.

booth = 7;

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
LEDPIN_R = 'D7';
LEDPIN_L = 'D5';

configurePin(a, LEDPIN_L, 'PWM');
configurePin(a, LEDPIN_R, 'PWM');

configurePin(a, SOLENOID_PIN, 'digitalOutput');

sensor_PIN = 'D45';          %sensor
configurePin(a, sensor_PIN, 'digitalInput');  %set sensor pin to digital input
configurePin(a, sensor_PIN, 'pullup');        %pull up

sensor_PIN2 = 'D10';          %sensor
configurePin(a, sensor_PIN2, 'digitalInput');  %set sensor pin to digital input
configurePin(a, sensor_PIN2, 'pullup');        %pull up

sensorState = 0;
lastState =0;
sensorState2 = 0;
lastState2 =0;

%Set reward duration
%reward time = 2 seems to work the best
reward_time = 2 ;

%The feeder moves after the sensor is triggered and has the reward present for a defined time
%duration
counter = 0;
counterL = 0;
counterR = 0;

%Parameters to enable/disable sensors
rightResponse = 0;
leftResponse = 1;

writePWMVoltage(a, LEDPIN_L, 0.2);

while true
    sensorState = readDigitalPin(a, sensor_PIN);
    sensorState2 = readDigitalPin(a, sensor_PIN2);
    
    if (~sensorState && lastState) && (leftResponse == 1)  %if not false(1) && True(1)
        disp('Left IR touched');            %sensorState = false (0) if touched
        writeDigitalPin(a, SOLENOID_PIN, 1);    %solenoid is activated (feeder out)
        pause(reward_time);    %feeder presents for a defined period of time
        writeDigitalPin(a, SOLENOID_PIN, 0);
        
        current_time = clock;   %show current time
    
        disp(['trial ', num2str(counter),'- sensor triggered at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
        '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
 
        counter = counter+1;
        counterL = counterL+1;
        disp(['Left-',num2str(counterL), ' Right-', num2str(counterR)]); 
        
        rightResponse = 1;
        leftResponse = 0;
        
        writePWMVoltage(a, LEDPIN_R, 0.2);
        writePWMVoltage(a, LEDPIN_L, 0);
        
    elseif (~sensorState2 && lastState2) && (rightResponse == 1)
        disp('Right IR touched');            %sensorState = false (0) if touched
        writeDigitalPin(a, SOLENOID_PIN, 1);    %solenoid is activated (feeder out)
        pause(reward_time);    %feeder presents for a defined period of time
        writeDigitalPin(a, SOLENOID_PIN, 0);
        
        current_time = clock;   %show current time
    
        disp(['trial ', num2str(counter),'- sensor triggered at ',num2str(current_time(1)),'-',num2str(current_time(2)), ...
        '-',num2str(current_time(3)),'  ',num2str(current_time(4)),'h ',num2str(current_time(5)), 'm']);
 
        counter = counter+1;
        counterR = counterR+1;
        disp(['Left-',num2str(counterL), ' Right-', num2str(counterR)]);
        
        rightResponse = 0;
        leftResponse = 1;
        
        writePWMVoltage(a, LEDPIN_R, 0);
        writePWMVoltage(a, LEDPIN_L, 0.2);
    end
  
    lastState = sensorState;
    lastState2 = sensorState2;
    
end