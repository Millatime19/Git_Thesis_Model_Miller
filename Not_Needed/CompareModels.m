% Given values for Henry's constants and groundwater concentrations (ug/L)
H_TCE = 2.56E-01;     
GW_TCE = 1.9;         

H_PCE = 5.07E-01;     
GW_PCE = 26.9;        

H_CIS = 1.14E-01;     
GW_CIS = 8211.7;      

H_VC = 9.57E-01;      
GW_VC = 291.1;          

% Calculated soil gas concentrations (ug/m^3)
SoilGas_TCE = GW_TCE * H_TCE * 1000;    
SoilGas_PCE = GW_PCE * H_PCE * 1000;    
SoilGas_CIS = GW_CIS * H_CIS * 1000;    
SoilGas_VC = GW_VC * H_VC * 1000;       

% Thesis data (modeled site sub-slab vapor concentrations)
thesis_data =  [5.66, 34, 16033, 5317];

% JE Model Surface Vapor Concentration (indoor air) and Sub Slab Vapor Concentration
JE_IA_conc = [1.4, 7.80e-02, 1.90e+02, 6.00E+01];  % in ug/m^3
JE_sub_slab_conc = [4.60e+02, 2.60e+01, 6.30e+04, 9.95E+02]; % in ug/m^3

% Sub Slab Vapor Sampling Concentrations prior to VIMS Operations
Measured_SS_conc = [2.15e+05, 6.23e+04, 5.24e+04, 0.00e+00];

% Contaminant names and concentrations
contaminants = {'PCE', 'TCE', 'Cis-1,2-DCE', 'Vinyl Chloride'};
calculated_concentrations = [SoilGas_PCE, SoilGas_TCE, SoilGas_CIS, SoilGas_VC];

% Plot
figure(423);
hold on;

% Plot JE model data with dots
plot(1:4, JE_IA_conc, 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'JE IA Data');
% Plot JE model data
plot(1:4, JE_sub_slab_conc, 'p', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'JE SS Data');

% Plot JE model data with dots
plot(1:4, Measured_SS_conc, 'x', 'MarkerSize', 8, 'Color','black','LineWidth', 1.5, 'DisplayName', 'Site SS Data');

% Plot calculated data with squares
plot(1:4, calculated_concentrations, 's', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Calculated Data');

% Plot thesis data with stars
plot(1:4, thesis_data, '*', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Thesis Data');

% Customize plot
set(gca, 'XTick', 1:4, 'XTickLabel', contaminants);
set(gca, 'YScale', 'log');
xlabel('Contaminant');
ylabel('Soil Gas Concentration (\mug/m^3)');
title('Comparison of Soil Gas Concentrations');
legend('Location', 'northeast');
grid on;
hold off;
