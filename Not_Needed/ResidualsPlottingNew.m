% GLM2 Analysis Script
% Load and analyze GLM2 data with time series

% Read the data
data = readtable('GLM2.csv');

% Extract the time vector
time = data.Time;

% Extract modeled and measured concentrations
measured = data.Measured;
modeled = data.Modelled;
weights = data.Weight;

% Get unique observation groups (assuming format mod_obsX_Y where Y is 1-7)
names = data.Name;
contaminant_indices = cellfun(@(x) str2double(x(end)), names);
unique_contaminants = unique(contaminant_indices);
num_contaminants = length(unique_contaminants);

% Initialize arrays for each contaminant
measured_by_contaminant = zeros(length(unique(time)), num_contaminants);
modeled_by_contaminant = zeros(length(unique(time)), num_contaminants);
residuals_by_contaminant = zeros(length(unique(time)), num_contaminants);
weighted_residuals_by_contaminant = zeros(length(unique(time)), num_contaminants);
unique_times = unique(time);

% [Previous code until data reorganization remains the same]

% Reorganize data by contaminant and filter out 999999 values
for i = 1:num_contaminants
    contam_mask = contaminant_indices == i;
    
    % Store data for each time point
    for t = 1:length(unique_times)
        time_mask = time == unique_times(t);
        mask = find(contam_mask & time_mask, 1); % Get the index of the first matching row
        
        if ~isempty(mask)
            % Only store if measured value is not 999999
            if measured(mask) == 999999
                measured_by_contaminant(t,i) = NaN;
                modeled_by_contaminant(t,i) = NaN;
                residuals_by_contaminant(t,i) = NaN;
                weighted_residuals_by_contaminant(t,i) = NaN;
            else
                measured_by_contaminant(t,i) = measured(mask);
                modeled_by_contaminant(t,i) = modeled(mask);
                residuals_by_contaminant(t,i) = modeled(mask) - measured(mask);
                weighted_residuals_by_contaminant(t,i) = (modeled(mask) - measured(mask)) * weights(mask);
            end
        else
            measured_by_contaminant(t,i) = NaN;
            modeled_by_contaminant(t,i) = NaN;
            residuals_by_contaminant(t,i) = NaN;
            weighted_residuals_by_contaminant(t,i) = NaN;
        end
    end
end

% Define selected contaminants
selected_contaminants = {'TCE', 'PCE', 'Cis-1,2-DCE', 'VC'};
selected_indices = [1, 2, 3, 5];  % Corresponding indices in the data

% Figure 1: Time series of measured vs modeled concentrations
figure(1);
for i = 1:length(selected_indices)
    subplot(2,2,i);
    valid_data = ~isnan(measured_by_contaminant(:,selected_indices(i)));
    plot(unique_times(valid_data), measured_by_contaminant(valid_data,selected_indices(i)), 'bo-', 'DisplayName', 'Measured');
    hold on;
    plot(unique_times(valid_data), modeled_by_contaminant(valid_data,selected_indices(i)), 'r.-', 'DisplayName', 'Modeled');
    title(selected_contaminants{i});
    xlabel('Time');
    ylabel('Concentration');
    legend('show');
    grid on;
    hold off;
end
sgtitle('Measured vs Modeled Concentrations');

% Figure 2: Residuals over time
figure(2);
hold on;
plot_handles = [];
for i = 1:length(selected_indices)
    valid_data = ~isnan(residuals_by_contaminant(:,selected_indices(i)));
    h = plot(unique_times(valid_data), residuals_by_contaminant(valid_data,selected_indices(i)), '-o');
    plot_handles = [plot_handles, h];
end
xlabel('Time');
ylabel('Residuals (Modeled - Measured)');
title('Residuals over Time');
legend(selected_contaminants, 'Location', 'best');
grid on;
hold off;

% Figure 3: Weighted residuals over time
figure(3);
hold on;
plot_handles = [];
for i = 1:length(selected_indices)
    valid_data = ~isnan(weighted_residuals_by_contaminant(:,selected_indices(i)));
    h = plot(unique_times(valid_data), weighted_residuals_by_contaminant(valid_data,selected_indices(i)), '-o');
    plot_handles = [plot_handles, h];
end
xlabel('Time');
ylabel('Weighted Residuals');
title('Weighted Residuals over Time');
legend(selected_contaminants, 'Location', 'best');
grid on;
hold off;
% [Previous code remains the same until before the summary statistics section]

% Calculate summary statistics for each contaminant
sum_residuals = zeros(length(selected_indices), 1);
sum_weighted_residuals = zeros(length(selected_indices), 1);
RMSE = zeros(length(selected_indices), 1);
weighted_RMSE = zeros(length(selected_indices), 1);
phi = zeros(length(selected_indices), 1);
phi_weighted = zeros(length(selected_indices), 1);

for i = 1:length(selected_indices)
    valid_data = ~isnan(residuals_by_contaminant(:,selected_indices(i)));
    
    % Calculate basic statistics
    sum_residuals(i) = sum(residuals_by_contaminant(valid_data,selected_indices(i)));
    sum_weighted_residuals(i) = sum(weighted_residuals_by_contaminant(valid_data,selected_indices(i)));
    
    % Calculate RMSE
    RMSE(i) = sqrt(mean(residuals_by_contaminant(valid_data,selected_indices(i)).^2));
    weighted_RMSE(i) = sqrt(mean(weighted_residuals_by_contaminant(valid_data,selected_indices(i)).^2));
    
    % Calculate Phi (sum of squared residuals)
    phi(i) = sum(residuals_by_contaminant(valid_data,selected_indices(i)).^2);
    phi_weighted(i) = sum(weighted_residuals_by_contaminant(valid_data,selected_indices(i)).^2);
end

% Calculate total statistics
total_sum_residuals = sum(sum_residuals);
total_sum_weighted_residuals = sum(sum_weighted_residuals);
total_RMSE = sqrt(mean([residuals_by_contaminant(:,selected_indices)].^2, 'omitnan'));
total_weighted_RMSE = sqrt(mean([weighted_residuals_by_contaminant(:,selected_indices)].^2, 'omitnan'));
total_phi = sum(phi);
total_phi_weighted = sum(phi_weighted);

% Display summary statistics
disp('Summary Statistics for Each Contaminant:');
for i = 1:length(selected_indices)
    fprintf('\n%s:\n', selected_contaminants{i});
    fprintf('Residual Sum: %.10f\n', sum_residuals(i));
    fprintf('Weighted Sum: %.10f\n', sum_weighted_residuals(i));
    fprintf('RMSE: %.10f\n', RMSE(i));
    fprintf('Weighted RMSE: %.10f\n', weighted_RMSE(i));
    fprintf('Phi: %.10f\n', phi(i));
    fprintf('Weighted Phi: %.10f\n', phi_weighted(i));
end

% Display totals
fprintf('\nTotal Statistics:\n');
fprintf('Total sum of residuals: %.10f\n', total_sum_residuals);
fprintf('Total sum of weighted residuals: %.10f\n', total_sum_weighted_residuals);
fprintf('Total RMSE: %.10f\n', total_RMSE);
fprintf('Total Weighted RMSE: %.10f\n', total_weighted_RMSE);
fprintf('Total Phi: %.10f\n', total_phi);
fprintf('Total Weighted Phi: %.10f\n', total_phi_weighted);