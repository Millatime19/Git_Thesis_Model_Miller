clear all;
% Load the Obs_Weights.csv file
obs_weights = readtable('Obs_Weights.csv');
Fig_Position = [100, 100, 1300, 600];
% Load the B4_SEN.csv file
b1_sen = readtable('Building_SEN_Run2.csv');
% b1_sen = readtable('Building_SEN1.csv');
% Take the absolute values of the relevant columns in b4_sen
b1_sen.sen_mean = abs(b1_sen.sen_mean);
b1_sen.sen_mean_abs = abs(b1_sen.sen_mean_abs);
b1_sen.sen_std_dev = abs(b1_sen.sen_std_dev);

% Get the list of unique parameters
unique_params = unique(b1_sen.parameter_name);

% Initialize arrays to store the results
final_sen_mean = zeros(length(unique_params), 1);
final_sen_mean_abs = zeros(length(unique_params), 1);
final_sen_std_dev = zeros(length(unique_params), 1);

% Initialize a cell array for final output
result_data = {};

% Loop through each unique parameter
for p = 1:length(unique_params)
    param = unique_params{p};
    
    % Find the rows in b4_sen corresponding to this parameter
    param_rows = strcmp(b1_sen.parameter_name, param);
    
    % Initialize variables to store the weighted sums and counts
    weighted_sen_mean_sum = 0;
    weighted_sen_mean_abs_sum = 0;
    weighted_sen_std_dev_sum = 0;
    total_weight = 0;
    
    % Loop through each row corresponding to the current parameter
    for i = find(param_rows)'
        % Get the observation name for this row
        obs_name = b1_sen{i, 1}; % Column 1 in b4_sen corresponds to observation_name
        
        % Find the corresponding weight in the obs_weights table
        weight_idx = find(strcmp(obs_weights{:, 1}, obs_name));
        
        if ~isempty(weight_idx)
            % Get the weight for this observation
            weight = obs_weights{weight_idx, 3}; % Assuming weight is in column 3
            
            % Get the values from b4_sen for this observation
            sen_mean = b1_sen{i, 4};      % Column 4 is sen_mean
            sen_mean_abs = b1_sen{i, 5};  % Column 5 is sen_mean_abs
            sen_std_dev = b1_sen{i, 6};   % Column 6 is sen_std_dev
            
            % Multiply by the weight and accumulate the sums
            weighted_sen_mean_sum = weighted_sen_mean_sum + weight * sen_mean;
            weighted_sen_mean_abs_sum = weighted_sen_mean_abs_sum + weight * sen_mean_abs;
            weighted_sen_std_dev_sum = weighted_sen_std_dev_sum + weight * sen_std_dev;
            total_weight = total_weight + weight;
        end
    end
    
    % Normalize by the total weight
    if total_weight > 0
        normalized_sen_mean = weighted_sen_mean_sum / total_weight;
        normalized_sen_mean_abs = weighted_sen_mean_abs_sum / total_weight;
        normalized_sen_std_dev = weighted_sen_std_dev_sum / total_weight;
    else
        normalized_sen_mean = 0;
        normalized_sen_mean_abs = 0;
        normalized_sen_std_dev = 0;
    end
    
    % Add the results to the result_data cell array
    n_samples = max(b1_sen.n_samples(param_rows)); % Assuming n_samples is constant for each parameter
    result_data = [result_data; {param, n_samples, normalized_sen_mean, normalized_sen_mean_abs, normalized_sen_std_dev}];
end

% Convert the result_data cell array to a table
result_table = cell2table(result_data, 'VariableNames', {'parameter_name', 'n_samples', 'sen_mean', 'sen_mean_abs', 'sen_std_dev'});
result_table(1, :) = [];
% % Define the new names
% new_names = {'Air Diff', 'FOC', 'Pore-size Dist.', 'SS Bottom', 'SS Top', ...
%              'K Sat.', 'Air Pressure', 'Tot. Porosity', 'Min. Vapor Porosity'};
% new_names = {'Air Diff', 'FOC', 'Pore-size Dist.', 'SS Bottom', 'SS Top', ...
             % 'K Sat.', 'Air Pressure', 'Tot. Porosity', 'Min. Vapor Porosity'};
% Update the first column of result_table with the new names
% result_table.parameter_name = new_names';
% Create a bar chart to visualize the mean absolute sensitivity and standard deviation for each parameter
figure(345);
bar_data = [result_table.sen_mean_abs, result_table.sen_std_dev];
bar(bar_data);
set(gca, 'XTickLabel', result_table.parameter_name, 'XTickLabelRotation', 45);
xlabel('Parameter');
ylabel('Sensitivity');
legend({'Mean Absolute Sensitivity', 'Standard Deviation Sensitivity'}, 'Location', 'northwest');
title('Global Sensitivity Results (Building Model)');
 set(gcf, 'Position', Fig_Position);
% Display the final result table
disp(result_table);

% Assuming you have the data in result_table
sen_mean_abs = result_table.sen_mean_abs;
sen_std_dev = result_table.sen_std_dev;
param_names = result_table.parameter_name;

% Create the scatter plot
figure(341);
scatter(sen_mean_abs, sen_std_dev, 50, [0.6, 0, 0], '^', 'filled'); % Red triangles

% Add text annotations for each point (use param_names)
for i = 1:length(sen_mean_abs)
    text(sen_mean_abs(i), sen_std_dev(i), param_names{i}, 'FontSize', 10, 'Color', 'black');
end

% Set axis limits based on the min/max of the data
mx = max([max(sen_mean_abs), max(sen_std_dev)]);
mn = min([min(sen_mean_abs), min(sen_std_dev)]);

% Plot diagonal line
hold on;
plot([mn, mx], [mn, mx], 'k--'); % Diagonal dashed line

% Set the axis limits
xlim([mn, mx]);
ylim([mn, mx]);

% Add labels and grid
xlabel('Mean Absolute Sensitivity'); % Use LaTeX interpreter for the special character
ylabel('Standard Deviation Sensitivity'); % Use LaTeX interpreter for the special character
title('Standard Deviation vs Absolute Mean Sensitivity (Building Model)');
grid on;
 set(gcf, 'Position', Fig_Position);
% Adjust layout
set(gca, 'Box', 'on');

hold off


% Assuming you have the data in result_table
sen_mean = result_table.sen_mean;  % Mean sensitivity for each parameter
param_names = result_table.parameter_name;  % Parameter names

% Calculate the total mean sensitivity
total_sensitivity = sum(sen_mean);

% Calculate the percent contribution of each parameter
percent_sensitivity = (sen_mean / total_sensitivity) * 100;

% Create a bar chart
figure(689);
bar(percent_sensitivity, 'FaceColor', [0.2, 0.2, 0.5]); % Dark blue bars

% Add labels and grid
set(gca, 'XTickLabel', param_names, 'XTickLabelRotation', 45); % Rotate x-axis labels for readability
xlabel('Parameter');
ylabel('Percent Contribution to Total Sensitivity');
title('Percent Total Mean Sensitivity by Parameter (Building Model)');
grid on;

% Adjust layout
set(gca, 'Box', 'on');
set(gca, 'YScale', 'log');
 set(gcf, 'Position', Fig_Position);
