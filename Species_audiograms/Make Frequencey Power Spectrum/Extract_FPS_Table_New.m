%%Goes through the stimInfo output from plotStimulus.m and creates a Matrix
%%of all the FPS for each file (aligned vertically for each 0.025kHz
%%frequency interval (columns)

%Creating the relative amplitude matrix
for i = 1:length(stimInfo)
    A = stimInfo(i).FPS; %the data in this matrix at the end of the loop will only apply to the last file (i.e., the i-th stimInfo) 
    B(i,:) = vertcat(A); %rows are each sound file, columns correspond to the FPS for that interval of frequency (0.025kHz to 8kHz)

meanFPS = mean(B,1); %gives you the output meanFPS with the mean FPS from all files at each frequency interval
maxFPS = max(meanFPS(:)); %find max FPS value to be used in relative amplitude calculation below
relAmp = meanFPS/maxFPS; %gives relative amplitude (of power) at each frequency interval across all files

Name = stimInfo(i).name;

%Concatenating the frequency and relative amplitude matrices and saving the
%file

labels = {'Frequency';'Relative Amplitude'}; %labels for FPS table
valFPS = vertcat(stimInfo(i).FRQ,relAmp); %concatenates (vertically) only the values of your FPS table/data
valuesFPS = num2cell(valFPS);
FPS = horzcat(labels,valuesFPS); %appends (concatenates) the labels to your table
xlswrite([Name '_FPS.xlsx'],FPS); %saves the labeled matrix/table as Excel file 

end
