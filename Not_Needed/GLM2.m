% Read the CSV file
data = readtable('GLM1.csv');

% Extract Measured and Modelled data
measured = data.Measured;
modelled = data.Modelled;
weights = data.Weight;

% Remove entries where Measured is 999999 (placeholder values)
valid_indices = measured ~= 999999;
measured = measured(valid_indices);
modelled = modelled(valid_indices);

% Create the scatter plot
figure('Position', [100, 100, 800, 600], 'Name', 'GLM2');
scatter(measured, modelled, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;

% Set axes to logarithmic scale
set(gca, 'XScale', 'log', 'YScale', 'log');

% Get current axis limits after plotting the data
xlims = xlim;
ylims = ylim;
min_val = min(xlims(1), ylims(1));
max_val = max(xlims(2), ylims(2));

% Add 1:1 line using axis limits
plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 1.5);

% Add labels and title
xlabel('Modeled VIMS Concentrations (\mug/m^3)');
ylabel('Measured Surface Vapor Concentrations (\mug/m^3)');
title('Modeled vs Measured Observations PESTPP-GLM Optimized Parameter Set ');

% Add grid
grid on;

% Adjust axes to be equal and square
axis square;

% Add legend
legend('Data Points', '1:1 Line', 'Location', 'northwest');

% Calculate and add R-squared value
valid_data = ~isnan(measured) & ~isnan(modelled);
R = corrcoef(measured(valid_data), modelled(valid_data));
R2 = R(1,2)^2;
text(min_val + (max_val-min_val)*0.1, max_val*0.9, sprintf('R^2 = %.3f', R2), 'FontSize', 10);
xlim([0,1000000])
hold off;