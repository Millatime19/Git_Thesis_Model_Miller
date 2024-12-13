% Define the summary file names for Building models
summary_files = {'B0_summary.csv', 'B1_summary.csv', 'B2_summary.csv', 'B3_summary.csv', 'B4_summary.csv'};
% Read the Building_PEST_Phi data for plotting
% Define number of parameters
num_parameters = [0 1 2 4 5]; % Corresponds to 1 parameter, 2 parameters, ..., 5 parameters
phi_building_data = [PEST_Phi_B1;PEST_Phi_B2;PEST_Phi_B3;PEST_Phi_B4] 
building_phi_values = phi_building_data{:, 2}; % Extract the Phi PEST values
building_phi_parameters = num_parameters; % Define parameter set for Building Phi data
% Initialize an empty cell array to store data from all summary files
combined_data = [];

% Define the headers for each test
headers = {'TCE', 'PCE', 'Cis-1,2-DCE', 'Toluene', 'VC', 'Xylenes', 'Ethylbenzene', 'Total Sum'};

% Iterate over each summary file and extract the data
for i = 1:length(summary_files)
    % Read the summary data from the current file
    data = readtable(summary_files{i}, 'ReadVariableNames', false, 'ReadRowNames', true);
    
    % Extract the relevant data (all rows, columns 2 to end)
    test_data = table2array(data(:, 1:end));
    
    % Append the test data to the combined array
    combined_data = [combined_data, test_data];
end

% Define the row names for the combined data
row_names = {'Residuals', 'Weighted', 'Weighted Average', 'RMSE', 'Weighted RMSE', 'Phi Weighted', 'Phi'};


% Create the headers for the combined table with test identifiers
combined_headers = {};
for i = 1:length(summary_files)
    for j = 1:length(headers)
        combined_headers = [combined_headers, strcat(headers{j}, '_B', num2str(i-1))];
    end
end

% Convert the combined data to a table for easier viewing and export
combined_table = array2table(combined_data, 'VariableNames', combined_headers, 'RowNames', row_names);

% Display the combined table
disp(combined_table);

% Export the combined table to a CSV file
writetable(combined_table, 'Combined_Building_Summary_Horizontal.csv', 'WriteRowNames', true);
% Load the combined data table
combined_table = readtable('Combined_Building_Summary_Horizontal.csv', 'ReadVariableNames', true, 'ReadRowNames', true);

% Extract the data for each metric (each row in the combined table)
residuals = table2array(combined_table('Residuals', :));
weighted = table2array(combined_table('Weighted', :));
weighted_avg = table2array(combined_table('Weighted Average', :));
rmse = table2array(combined_table('RMSE', :));
weighted_rmse = table2array(combined_table('Weighted RMSE', :));
phi_weight = table2array(combined_table('Phi Weighted', :));
phi = table2array(combined_table('Phi', :));
% Define number of parameters
num_parameters = [0 1 2 4 5]; % Corresponds to 1 parameter, 2 parameters, ..., 5 parameters

% Define contaminant names for plotting and legend purposes
contaminant_names = headers;


% Reshape the data into 5 columns, each representing a parameter set
residuals = reshape(residuals, [], length(num_parameters));
weighted = reshape(weighted, [], length(num_parameters));
weighted_avg = reshape(weighted_avg, [], length(num_parameters));
rmse = reshape(rmse, [], length(num_parameters));
weighted_rmse = reshape(weighted_rmse, [], length(num_parameters));
phi_weight = reshape(phi_weight, [], length(num_parameters));
phi = reshape(phi, [], length(num_parameters));

% Create the plots for each metric
figure(239);
% Define number of parameters
num_parameters = [0 1 2 4 5]; % Corresponds to 1 parameter, 2 parameters, ..., 5 parameters
% Plot Residuals
subplot(4, 2, 1);
hold on;
plot_handles = plot(num_parameters, residuals, '-o');
set(plot_handles(end), 'Color', 'k'); % Set the last plot ("Sum") to black
title('Residuals');
xlabel('Number of Parameters');
ylabel('Residuals');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted Residuals
subplot(4, 2, 2);
hold on;
plot_handles = plot(num_parameters, weighted, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted Residuals');
xlabel('Number of Parameters');
ylabel('Weighted Residuals');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted Average
subplot(4, 2, 3);
hold on;
plot_handles = plot(num_parameters, weighted_avg, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted Average');
xlabel('Number of Parameters');
ylabel('Weighted Average');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot RMSE
subplot(4, 2, 4);
hold on;
plot_handles = plot(num_parameters, rmse, '-o');
set(plot_handles(end), 'Color', 'k');
title('RMSE');
xlabel('Number of Parameters');
ylabel('RMSE');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted RMSE
subplot(4, 2, 5);
hold on;
plot_handles = plot(num_parameters, weighted_rmse, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted RMSE');
xlabel('Number of Parameters');
ylabel('Weighted RMSE');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Phi Weighted
subplot(4, 2, 6);
hold on;
plot_handles = plot(num_parameters, phi_weight, '-o');
set(plot_handles(end), 'Color', 'k');
title('Phi Weighted (Sum of Weighted Squared Errors)');
xlabel('Number of Parameters');
ylabel('Phi Weighted');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Phi (regular) and Building Phi data
subplot(4, 2, 7);
hold on;
plot_handles = plot(num_parameters, phi, '-o'); % Plot phi normally
set(plot_handles(end), 'Color', 'k'); 
plot(building_phi_parameters, building_phi_values, '-*k', 'MarkerSize', 8); % Black line with asterisk markers
title('Phi (Sum of Squared Errors)');
xlabel('Number of Parameters');
ylabel('Phi');
set(gca, 'YScale', 'log'); % Set y-axis to log scale
% Update the legend to include Building Phi data
legend([contaminant_names, 'Building PEST Phi'], 'Location', 'best'); 
hold off;

% Adjust the overall layout to accommodate the plots
sgtitle('PEST Results and Testing Parameters Comparison - Building');

% Define the summary file names for Outside models
outside_summary_files = {'O0_summary.csv', 'O1_summary.csv', 'O2_summary.csv', 'O3_summary.csv', 'O4_summary.csv'};

% Initialize an empty array to store combined data for the outside runs
outside_combined_data = [];

% Define headers for each test (common for all runs)
headers = {'TCE', 'PCE', 'Cis-1,2-DCE', 'Toluene', 'VC', 'Xylenes', 'Ethylbenzene', 'Total Sum'};

% Iterate over each outside summary file and extract data
for i = 1:length(outside_summary_files)
    % Read the summary data from the current outside file
    data = readtable(outside_summary_files{i}, 'ReadVariableNames', false, 'ReadRowNames', true);
    
    % Extract relevant data (all rows, columns 2 to end)
    test_data = table2array(data(:, 1:end));
    
    % Append the test data to the combined array
    outside_combined_data = [outside_combined_data, test_data];
end

% Define the row names for the combined data
row_names = {'Residuals', 'Weighted', 'Weighted Average', 'RMSE', 'Weighted RMSE', 'Phi Weighted', 'Phi'};

% Create headers for the combined table with test identifiers
outside_combined_headers = {};
for i = 1:length(outside_summary_files)
    for j = 1:length(headers)
        outside_combined_headers = [outside_combined_headers, strcat(headers{j}, '_O', num2str(i-1))];
    end
end

% Convert combined data to a table for easier viewing and export
outside_combined_table = array2table(outside_combined_data, 'VariableNames', outside_combined_headers, 'RowNames', row_names);

% Display the combined table
disp(outside_combined_table);

% Export the combined table to a CSV file
writetable(outside_combined_table, 'Combined_Outside_Summary_Horizontal.csv', 'WriteRowNames', true);

% Load the combined data table for plotting
outside_combined_table = readtable('Combined_Outside_Summary_Horizontal.csv', 'ReadVariableNames', true, 'ReadRowNames', true);

% Extract the data for each metric (each row in the combined table)
residuals_out = table2array(outside_combined_table('Residuals', :));
weighted_out = table2array(outside_combined_table('Weighted', :));
weighted_avg_out = table2array(outside_combined_table('Weighted Average', :));
rmse_out = table2array(outside_combined_table('RMSE', :));
weighted_rmse_out = table2array(outside_combined_table('Weighted RMSE', :));
phi_weight_out = table2array(outside_combined_table('Phi Weighted', :));
phi_out = table2array(outside_combined_table('Phi', :));

% Define number of parameters
num_parameters = [0 1 2 4 5]; % Corresponds to 1 parameter, 2 parameters, ..., 5 parameters
phi_outside_data = [PEST_Phi_O1;PEST_Phi_O2;PEST_Phi_O3;PEST_Phi_O4] 
building_phi_values = phi_building_data{:, 2}; % Extract the Phi PEST values
building_phi_parameters = num_parameters; % Define parameter set for Building Phi data
% Define contaminant names for plotting and legend purposes
contaminant_names = headers;

% Reshape data into 5 columns, each representing a parameter set
residuals_out = reshape(residuals_out, [], length(num_parameters));
weighted_out = reshape(weighted_out, [], length(num_parameters));
weighted_avg_out = reshape(weighted_avg_out, [], length(num_parameters));
rmse_out = reshape(rmse_out, [], length(num_parameters));
weighted_rmse_out = reshape(weighted_rmse_out, [], length(num_parameters));
phi_weight_out = reshape(phi_weight_out, [], length(num_parameters));
phi_out = reshape(phi_out, [], length(num_parameters));

% Create plots for each metric
figure(240);

% Plot Residuals
subplot(4, 2, 1);
hold on;
plot_handles = plot(num_parameters, residuals_out, '-o');
set(plot_handles(end), 'Color', 'k'); % Set the last plot ("Sum") to black
title('Residuals - Outside');
xlabel('Number of Parameters');
ylabel('Residuals');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted Residuals
subplot(4, 2, 2);
hold on;
plot_handles = plot(num_parameters, weighted_out, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted Residuals - Outside');
xlabel('Number of Parameters');
ylabel('Weighted Residuals');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted Average
subplot(4, 2, 3);
hold on;
plot_handles = plot(num_parameters, weighted_avg_out, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted Average - Outside');
xlabel('Number of Parameters');
ylabel('Weighted Average');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot RMSE
subplot(4, 2, 4);
hold on;
plot_handles = plot(num_parameters, rmse_out, '-o');
set(plot_handles(end), 'Color', 'k');
title('RMSE - Outside');
xlabel('Number of Parameters');
ylabel('RMSE');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Weighted RMSE
subplot(4, 2, 5);
hold on;
plot_handles = plot(num_parameters, weighted_rmse_out, '-o');
set(plot_handles(end), 'Color', 'k');
title('Weighted RMSE - Outside');
xlabel('Number of Parameters');
ylabel('Weighted RMSE');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Phi Weighted
subplot(4, 2, 6);
hold on;
plot_handles = plot(num_parameters, phi_weight_out, '-o');
set(plot_handles(end), 'Color', 'k');
title('Phi Weighted (Sum of Weighted Squared Errors) - Outside');
xlabel('Number of Parameters');
ylabel('Phi Weighted');
legend(contaminant_names, 'Location', 'best');
hold off;

% Plot Phi (regular)
subplot(4, 2, 7);
hold on;
plot_handles = plot(num_parameters, phi_out, '-o'); % Plot phi normally
set(plot_handles(end), 'Color', 'k');
title('Phi (Sum of Squared Errors) - Outside');
xlabel('Number of Parameters');
ylabel('Phi');
set(gca, 'YScale', 'log'); % Set y-axis to log scale
legend(contaminant_names, 'Location', 'best');
hold off;

% Adjust the overall layout to accommodate the plots
sgtitle('PEST Results and Testing Parameters Comparison - Outside');
