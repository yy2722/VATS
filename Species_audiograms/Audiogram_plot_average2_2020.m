%% Auditory Selectivity Plotting by Yow-Tyng (Tim) 08/05/2019
% edited by Tim 02/11/20 - includes error bars (SE) and default threshold
% value; averaging threshold values instead of doing cumulative threshold
% values

% This script is for plotting data acquired from auditory selectivity test.
% It will compile all excel data files within the selected folder and
% plot for an audiogram.



path = uigetdir;
dir_folder = dir(path);  % Get directory using user interface

%% Parameters

% These variables are called by functions later
global sound_trials
global num_files
global num_freq

bird_species = 'LF Average Audiogram'; % ID title for the audiogram

num_freq = 9;   % Number of frequencies tested    ** depends on experiment
lvls = 7;        % Number of sound levels tested
sound_trials = num_freq*lvls;         % Total number of trials excluding sham trials
num_files = 5;    % Number of data files per bird (sessions)

num_birds = 6;   % Number of birds analyzed       ** depends on data files being analyzed

correct_threshold = (num_files*num_birds)/2; % Threshold for 50% correct response for average
indiv_cr = (num_files)/2; % CR threshold for each bird

%% Create empty arrays for later calculations

compiled_array = cell(sound_trials, 3);   % Array for storing all 
bird_array = cell(sound_trials, num_birds); % Arrary for storing responses from all the birds
indiv_threshold_cell = cell(num_birds, 1);
all_fq_threshold = zeros(num_freq, num_birds+1);
error_bar = zeros(num_freq, 3);

%% Analyze the data from all the birds

for f = 1:length(dir_folder)
    if endsWith(dir_folder(f).name, '_completed') % only analyze data folders end with _completed
        bird = dir(fullfile(path, dir_folder(f).name));
        bird_array(:, f-2)= sum_trials(bird);        % Get the sum of all trials for each bird
    end
end

bird_array = cell2mat(bird_array);
total_bird_trials = sum(bird_array, 2);

%% Get cumulative minimal dB thresholds

min_threshold_average = get_min_threshold(total_bird_trials, correct_threshold);  % Get the average response thresholds

for b = 1:num_birds
    indiv_threshold_cell(b,:) = {get_min_threshold(bird_array(:,b), indiv_cr)}; % Get the minimal dB thresholds for each bird
end

%% Averaging all thresholds at each frequency

unnested_indiv = vertcat(indiv_threshold_cell{:});
fq_range = [750, 1000:1000:8000]; % the range needs to be adjusted accordingly
input = 0;
average_mat = zeros(numel(fq_range), 4);

for i = fq_range
    input = input +1;
    
    unnested_indiv(unnested_indiv(:,1) == i, 2);
    size = numel(unnested_indiv(unnested_indiv(:,1) == i, 2));
    
    % 02/10/20 - Tim
    % use default threshold 80 dB for birds that cannot detect the frequency to calculataverage
    if size < num_birds
       diff_num = num_birds - size;
       default_add = repelem(80, diff_num)';
       
       mean_threshold = mean([unnested_indiv(unnested_indiv(:,1) == i, 2); default_add]);
       std_threshold = std([unnested_indiv(unnested_indiv(:,1) == i, 2); default_add]);
    
    else
       mean_threshold = mean(unnested_indiv(unnested_indiv(:,1) == i, 2));
       std_threshold = std(unnested_indiv(unnested_indiv(:,1) == i, 2));
    end    
    
    se_threshold = std_threshold/sqrt(num_birds);
    
    average_mat(input, 1) = i;
    average_mat(input, 2) = mean_threshold;
    average_mat(input, 3) = std_threshold;
    average_mat(input, 4) = se_threshold;    
end
%% Error bars

% row = 1;
% 
% 
% for fqRange = [750, 1000:1000:8000] % This needs to be changed to the testing frequency range!
%     
%     all_fq_threshold(row, 1) = fqRange;
%     error_bar(row, 1) = fqRange;
%     
%     for i = 1:num_birds   % Get the response threshold for each bird at each frequency
%         threshold = indiv_threshold_cell{i}(indiv_threshold_cell{i}(:,1) == fqRange,2);
%         
%         if isempty(threshold)  % If the response threshold is above 70 dB (frequency cannot be matched), response threshold equals zero
%             threshold = 0;
%         end
%         
%         all_fq_threshold(row, i+1) = threshold;  % Add the response threshold to each column on a specific frequency (row) 
%     end
%     
%     greater_zero = all_fq_threshold(row, all_fq_threshold(row,:) > 0);  
%     
%     if length(greater_zero) == 1
%         greater_zero_max = 0;
%         greater_zero_min = 0;
%         
%     else
%         greater_zero_max = max(greater_zero(2:end));
%         greater_zero_min = min(greater_zero(2:end));
%         
%     end
%   
%     error_bar(row, 2) = greater_zero_max;  % Highest response threshold among all birds at a specific fq
%     error_bar(row, 3) = greater_zero_min;  % Lowest response threshold among all birds at a specific fq
%     
%     row = row+1;
% end
% 
% % Find rows with overlapping frequencies between the average thresholds and individual thresholds
% [all_fq,thres_fq,error_fq]= intersect(min_threshold_average(:,1), error_bar(:,1), 'rows');
% 
% % Get the differences between nonzero max and min response thresholds and average thresholds
% err_max = abs(nonzeros(error_bar(error_fq,2))' - min_threshold_average(:,2)');
% 
% % If the sum of error bar and average is greater than 70 dB, set max error bar value to 0
% if any((min_threshold_average(:,2)' + err_max) > 70)
%     err_max((min_threshold_average(:,2)' + err_max) > 70) = 0;
% end
% 
% err_min = abs(nonzeros(error_bar(error_fq,3))' - min_threshold_average(:,2)');
% 
% 

%% Plot audiogram

%plot(min_threshold_average(:,1), min_threshold_average(:,2),'r') % Plot frequency vs threshold amplitude
%plot(min_threshold_average(:,1), min_threshold_average(:,2),'b')
%plot(min_threshold_average(:,1), min_threshold_average(:,2),'k')
%hold on

%errorbar(average_mat(:,1), average_mat(:,2), average_mat(:,3), '-or')
plot(average_mat(:,1), average_mat(:,2), '-k');

%errorbar(min_threshold_average(:,1), min_threshold_average(:,2), err_min, err_max, 'o');  % Create error bars
%e.color = 'blue';

%title([bird_species, ' Audiogram'])
title(bird_species)
set(gca, 'FontSize', 20)

xlabel('Frequency', 'FontSize', 20)
xlim([0 10000])

ylabel('dB', 'FontSize', 20)
ylim([0 90])

%legend({'ZF', 'GW', 'LF'},'Location', 'northeast') % Adding legend
%lgd.FontSize = 20;
%% Funtion that gets the frequency for each trial

function y = freq(x)
    split_str = split(x, '_');
    y = split_str(2);
end

%% Function that determines each bird's responses across all the sessions

function each_bird = sum_trials(bird)

global sound_trials
global num_files
global sorted_data

raw_dataSets = zeros(sound_trials, num_files);

    for i = 1:length(bird)
        if endsWith(bird(i).name, '.xlsx')
           filename = bird(i).name;
           data = readtable(filename);

           data_cell = table2cell(data(:,1:2));
           freq_cell = strfind(data_cell(:,1), 'tone');
           freq_index = ~cellfun(@isempty, freq_cell); %#ok<STRCLFH>
           freq_data = data_cell(freq_index,:);
           sorted_data = sortrows(freq_data);

           raw_dataSets(:,i-2) = cell2mat(sorted_data(:, 2));
        end
    end
    
each_bird = num2cell(sum(raw_dataSets,2));
end

%% Function that gets threshold amplitude for each frequency

function final_output = get_min_threshold(trials, cr)

global sorted_data
global num_freq

    freq_str = cellfun(@freq, sorted_data(:,1));  % Get frequency range
    freq_range = cellfun(@str2num, freq_str);
    lvl_range = repmat([10:10:70]', num_freq,1); %#ok<NBRAK> Get dB lvls

    compiled_array(:,1) = num2cell(freq_range);
    compiled_array(:,2) = num2cell(lvl_range);
    compiled_array(:,3) = num2cell(trials);

    over_threshold = compiled_array((trials >= cr),:);

    for i = 1:(length(over_threshold)-1)   % Find minimal threshold dB by covering other dB levels with threshold dB
        if over_threshold{i,1} == over_threshold{i+1,1}
           over_threshold(i+1,:) = over_threshold(i,:);
        end    
    end

    min_threshold_dB = unique(cell2table(over_threshold)); % Get rows with unique values
    min_threshold_dB = table2array(min_threshold_dB);   % Convert data type from table to array

final_output = min_threshold_dB;
end