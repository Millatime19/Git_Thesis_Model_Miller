% Define q values and concentrations (using small value instead of 0)
q_values = [1e-8, 8e-5, 5e-4, 1e-3, 3e-3];  % Changed 0 to 1e-8 for plotting

% PCE data
pce_thesis = [52875, 12000, 7140, 4810, 2090];
pce_thesis_source = [13400, 13700, 13400, 13400, 13000];  % Your model GW table results
pce_je_subslab = 460;
pce_je_source = 12200;

% TCE data
tce_thesis = [444, 387, 230, 155, 67.5];
tce_thesis_source = [441, 441, 431, 431, 420];  % Your model GW table results
tce_je_subslab = 26;
tce_je_source = 505;

% DCE data
dce_thesis = [950000, 955000, 498000, 336000, 146000];
dce_thesis_source = [953000, 955000, 932000, 932000, 909000];  % Your model GW table results
dce_je_subslab = 63000;
dce_je_source = 947000;

% VC data
vc_thesis = [235, 189, 99.1, 66.0, 29.0];
vc_thesis_source = [189, 189, 185, 185, 180];  % Your model GW table results
vc_je_subslab = 20000;
vc_je_source = 265000;

% Define darker colors
thesis_color = [0 0 0.7];  % Darker blue
thesis_source_color = [0.7 0 0.7];  % Purple for thesis source
subslab_color = [0 0.5 0];  % Darker green
source_color = [0.7 0 0];   % Darker red

% Create figure with 4 subplots
figure('Position', [100, 100, 1200, 800])

% PCE subplot
subplot(2,2,1)
semilogx(q_values, pce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2)
hold on
semilogx(q_values, pce_thesis_source, '.-', 'Color', thesis_source_color, 'MarkerSize', 20, 'LineWidth', 2)
yline(pce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2)
yline(pce_je_source, '--', 'Color', source_color, 'LineWidth', 2)
hold off
grid on
title('PCE Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([1e2 1e5])
legend('Thesis Subslab', 'Thesis Source', 'JE Subslab', 'JE Source', 'Location', 'best')
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')

% TCE subplot
subplot(2,2,2)
semilogx(q_values, tce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2)
hold on
semilogx(q_values, tce_thesis_source, '.-', 'Color', thesis_source_color, 'MarkerSize', 20, 'LineWidth', 2)
yline(tce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2)
yline(tce_je_source, '--', 'Color', source_color, 'LineWidth', 2)
hold off
grid on
title('TCE Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([1e1 1e4])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')

% DCE subplot
subplot(2,2,3)
semilogx(q_values, dce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2)
hold on
semilogx(q_values, dce_thesis_source, '.-', 'Color', thesis_source_color, 'MarkerSize', 20, 'LineWidth', 2)
yline(dce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2)
yline(dce_je_source, '--', 'Color', source_color, 'LineWidth', 2)
hold off
grid on
title('DCE Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([1e4 1e7])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')

% VC subplot
subplot(2,2,4)
semilogx(q_values, vc_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2)
hold on
semilogx(q_values, vc_thesis_source, '.-', 'Color', thesis_source_color, 'MarkerSize', 20, 'LineWidth', 2)
yline(vc_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2)
yline(vc_je_source, '--', 'Color', source_color, 'LineWidth', 2)
hold off
grid on
title('VC Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([1e1 1e6])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')

% Adjust subplot spacing
set(gcf, 'Position', [100, 100, 1200, 800])