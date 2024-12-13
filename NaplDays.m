% Data
days = [6000 6500 7000 7500 8000 9000 10000];
NAPL = [0 51 212 402 592 972 1352];

% Create the plot
figure;
plot(NAPL, days,  'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');

% Add labels and title
ylabel('Shallow Soil Total Mass Soil Concentration (mg/kg)');
xlabel('Days NAPL Present');
title(' Total Mass Concentration of Soil vs NAPL Persistence');

% Add grid
grid on;

% Make it look nicer
set(gca, 'FontSize', 8);
box on;

% If you want to adjust axis limits
% xlim([4500 10500]);
% ylim([-100 1400]);

%%

% % % % % % % Percentage of Total Soil Mass Removed % % % % % % 

% Data
runtime = [0, 3300, 10000, 20000, 30000, 40000, 50000, 100000];
percent_removed = [0 , 30, 47, 60, 67, 71, 75, 85];

% Create the figure
figure('Position', [100, 100, 800, 500]);

% Plot with a line and markers
plot(runtime, percent_removed, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');

% Add labels and title
xlabel('VIMS Runtime (Days)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Percent Removed (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('VIMS Runtime vs. Percent Soil Source Removed', 'FontSize', 14, 'FontWeight', 'bold');

% Add grid
grid on;

% Format x-axis to use comma separator for thousands
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.0f');

% Adjust axis limits to add some padding
xlim([0, max(runtime)*1.05]);
ylim([0, max(percent_removed)*1.05]);

% Add data labels
for i = 1:length(runtime)
    text(runtime(i), percent_removed(i)+1, [num2str(percent_removed(i)) '%'], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end