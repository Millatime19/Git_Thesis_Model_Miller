clear all; 
% Define q values and concentrations (using small value instead of 0)
q_values = [1e-8, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1];

% PCE data
pce_thesis = [52875,11626,	4812,	701,	146,	74];
pce_je_subslab = 460;
pce_je_source = 12200;
pce_thesis_source = 13797; 
% TCE data
tce_thesis = [444,	376, 155, 23,	5,	2];
tce_je_subslab = 26;
tce_je_source = 505;
tce_thesis_source = 445; 

% DCE data
dce_thesis = [950000,	811807,	336030,	48978,	10211,	5143];
dce_je_subslab = 63000;
dce_je_source = 947000;
dce_thesis_source = 963363; 
% VC data
vc_thesis = [281892,	237545,	98326,	14331,	2988,	1502];
vc_je_subslab = 20000;
vc_je_source = 265000;
vc_thesis_source = 281829; 

% Define colors and markers
thesis_color = [0 0 0.7];    % Darker blue
thesis_source_color = [0 0 0];  % Black
subslab_color = [0 0.5 0];    % Darker green
source_color = [0.7 0 0];     % Darker red

% Best q 
pce_best = [0.01 701];
tce_best = [0.01 23];
dce_best = [0.01 48978];
vc_best = [0.01 14331];

figure(462)
% PCE subplot
subplot(2,2,1)
h1 = semilogx(q_values, pce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2);
hold on
h2 = semilogx(q_values(1), pce_thesis_source, 'kd', 'MarkerSize', 5, 'LineWidth', 2);
h3 = yline(pce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2);
h4 = yline(pce_je_source, '--', 'Color', source_color, 'LineWidth', 2);
h5 = plot(pce_best(1,1), pce_best(1,2), 'k*', 'MarkerSize', 15);
legend([h1 h2 h3 h4 h5], {'Thesis Subslab', 'Thesis Source', 'JE Subslab', 'JE Source', 'Best'}, 'Location', 'best')
grid on
title('PCE Vapor Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([0 1e5])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')
hold off

% TCE subplot
subplot(2,2,2)
h1 = semilogx(q_values, tce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2);
hold on
h2 = semilogx(q_values(1), tce_thesis_source, 'kd', 'MarkerSize', 5, 'LineWidth', 2);
h3 = yline(tce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2);
h4 = yline(tce_je_source, '--', 'Color', source_color, 'LineWidth', 2);
h5 = plot(tce_best(1,1), tce_best(1,2), 'k*', 'MarkerSize', 15);
legend([h1 h2 h3 h4 h5], {'Thesis Subslab', 'Thesis Source', 'JE Subslab', 'JE Source', 'Best'}, 'Location', 'best')
grid on
title('TCE Vapor Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([0 1e4])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')
hold off

% DCE subplot
subplot(2,2,3)
h1 = semilogx(q_values, dce_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2);
hold on
h2 = semilogx(q_values(1), dce_thesis_source, 'kd', 'MarkerSize', 5, 'LineWidth', 2);
h3 = yline(dce_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2);
h4 = yline(dce_je_source, '--', 'Color', source_color, 'LineWidth', 2);
h5 = plot(dce_best(1,1), dce_best(1,2), 'k*', 'MarkerSize', 15);
legend([h1 h2 h3 h4 h5], {'Thesis Subslab', 'Thesis Source', 'JE Subslab', 'JE Source', 'Best'}, 'Location', 'best')
grid on
title('DCE Vapor Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([0 1e7])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')
hold off

% VC subplot
subplot(2,2,4)
h1 = semilogx(q_values, vc_thesis, '.-', 'Color', thesis_color, 'MarkerSize', 20, 'LineWidth', 2);
hold on
h2 = semilogx(q_values(1), vc_thesis_source, 'kd', 'MarkerSize', 5, 'LineWidth', 2);
h3 = yline(vc_je_subslab, '--', 'Color', subslab_color, 'LineWidth', 2);
h4 = yline(vc_je_source, '--', 'Color', source_color, 'LineWidth', 2);
h5 = plot(vc_best(1,1), vc_best(1,2), 'k*', 'MarkerSize', 15);
legend([h1 h2 h3 h4 h5], {'Thesis Subslab', 'Thesis Source', 'JE Subslab', 'JE Source', 'Best'}, 'Location', 'best')
grid on
title('VC Vapor Concentrations', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('q (cm/s)', 'FontSize', 11)
ylabel('Concentration (µg/m³)', 'FontSize', 11)
set(gca, 'YScale', 'log')
ylim([0 1e6])
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off')
hold off

% Adjust subplot spacing
set(gcf, 'Position', [100, 100, 1200, 800])