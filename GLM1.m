% Read the CSV file
data = readtable('GLM2.csv');

% Extract data
measured = data.Measured;
modelled = data.Modelled;
weights = data.Weight;
names = data.Name;

% Create contaminant categories based on the last digit of the name
contaminantTypes = zeros(size(names));
for i = 1:length(names)
    % Extract the last character of each name
    nameParts = split(names{i}, '_');
    lastDigit = str2double(nameParts{end});
    contaminantTypes(i) = lastDigit;
end

% Remove entries where Measured is 999999 (placeholder values)
valid_indices = measured ~= 999999;
measured = measured(valid_indices);
modelled = modelled(valid_indices);
contaminantTypes = contaminantTypes(valid_indices);

% Create figure
figure('Position', [100, 100, 800, 600], 'Name', 'GLM2');

% Define colors and names for each contaminant type
colors = {[0.8500, 0.3250, 0.0980],  % TCE (red-orange)
          [0, 0.4470, 0.7410],       % PCE (blue)
          [0.4660, 0.6740, 0.1880],  % CIS (green)
          [0.4940, 0.1840, 0.5560],  % Toluene (purple)
          [0.9290, 0.6940, 0.1250],  % VC (yellow)
          [0.6350, 0.0780, 0.1840],  % Xylenes (dark red)
          [0.3010, 0.7450, 0.9330]}; % Ethylbenzene (light blue)

contaminantNames = {'TCE', 'PCE', 'CIS', 'Toluene', 'VC', 'Xylenes', 'Ethylbenzene'};

% Plot each contaminant type separately
hold on;
legendEntries = [];
for i = 1:7
    idx = contaminantTypes == i;
    if any(idx)
        h = scatter(measured(idx), modelled(idx), 50, colors{i}, 'filled', 'MarkerFaceAlpha', 0.6);
        legendEntries = [legendEntries h];
    end
end

% Set axes to logarithmic scale
set(gca, 'XScale', 'log', 'YScale', 'log');

% Get current axis limits after plotting the data
xlims = xlim;
ylims = ylim;
min_val = min(xlims(1), ylims(1));
max_val = max(xlims(2), ylims(2));

% Add 1:1 line using axis limits
h_line = plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 1.5);
legendEntries = [legendEntries h_line];

% Add labels and title
ylabel('Modeled VIMS Concentrations (\mug/m^3)');
xlabel('Measured Surface Vapor Concentrations (\mug/m^3)');
title('Modeled vs Measured Observations PESTPP-GLM Optimized Parameter Set');

% Add grid
grid on;

% Adjust axes to be equal and square
axis square;

% Add legend with contaminant names
legend([legendEntries], [contaminantNames, '1:1 Line'], 'Location', 'northwest');

% Calculate and add R-squared value
valid_data = ~isnan(measured) & ~isnan(modelled);
R = corrcoef(measured(valid_data), modelled(valid_data));
R2 = R(1,2)^2;
text(min_val + (max_val-min_val)*0.1, max_val*0.9, sprintf('R^2 = %.3f', R2), 'FontSize', 10);

% Set x-axis limit
xlim([0,1000000]);

hold off;