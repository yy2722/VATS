birdID = 'xx0000';   %choose the bird ID for compiling data and data plotting
data_dir = dir(['C:\Users\labadmin\Desktop\Arduino ARTsy\Data\*', birdID,'*.xlsx']);
data_cell = transpose({data_dir.name});

for i=1:(length(data_dir))

    info = fullfile('C:\Users\labadmin\Desktop\Arduino ARTsy\Data', data_dir(i).name);
    data_cell(i) = {(xlsread(info)};
    
end
contat = [data_cell{:}];

%%accuracy = (hit + cr)/trial_number;
%%error = (false_alarm + miss)/trial_number;



%%accuracy_plot = plot(trial_number, accuracy);
%%error_plot = plot(triaL_number, error);
