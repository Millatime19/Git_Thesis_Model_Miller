% Read the data
data = readtable('Extracted_PCETCE.csv');

% Convert dates to datetime format
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM');

% Calculate days since first measurement
days_elapsed = days(dates - dates(1));

% Convert ug/sec to ug/day (multiply by seconds in a day)
PCE_daily = data.PCE_Discharge * 86400;  % 86400 seconds in a day
TCE_daily = data.TCE_Discharge * 86400;

% Calculate time differences in days
days_between = days(diff([dates; dates(end) + calmonths(6)]));  % Convert to numeric days

% Calculate cumulative mass (in ug)
PCE_cumulative = cumsum(PCE_daily .* days_between);
TCE_cumulative = cumsum(TCE_daily .* days_between);

% Define darker colors
PCE_color = [0, 0, 0.6];  % Darker blue
TCE_color = [0.6, 0, 0];  % Darker red

% Create figure 1: Daily extraction rates
figure('Position', [100, 100, 800, 500], 'Name', 'Daily Extraction Rates');
semilogy(days_elapsed, PCE_daily, '-o', 'Color', PCE_color, 'LineWidth', 1.5, 'MarkerFaceColor', PCE_color);
hold on;
semilogy(days_elapsed, TCE_daily, '-s', 'Color', TCE_color, 'LineWidth', 1.5, 'MarkerFaceColor', TCE_color);
grid on;
xlabel('Days Since Start');
ylabel('Extraction Rate (µg/day)');
title('PCE and TCE Daily Extraction Rates');
legend('PCE', 'TCE', 'Location', 'best');

% Create figure 2: Cumulative mass extracted
figure('Position', [100, 100, 800, 500], 'Name', 'Cumulative Mass Extracted');
plot(days_elapsed, PCE_cumulative/1e6, '-o', 'Color', PCE_color, 'LineWidth', 1.5, 'MarkerFaceColor', PCE_color);
hold on;
plot(days_elapsed, TCE_cumulative/1e6, '-s', 'Color', TCE_color, 'LineWidth', 1.5, 'MarkerFaceColor', TCE_color);
grid on;
xlabel('Days Since Start');
ylabel('Cumulative Mass Extracted (g)');
title('PCE and TCE Cumulative Mass Extracted');
legend('PCE', 'TCE', 'Location', 'best');

% Print total mass extracted
fprintf('Total PCE mass extracted: %.2f g\n', PCE_cumulative(end)/1e6);
fprintf('Total TCE mass extracted: %.2f g\n', TCE_cumulative(end)/1e6);

% Print total days of operation
fprintf('Total days of operation: %.0f days\n', days_elapsed(end));