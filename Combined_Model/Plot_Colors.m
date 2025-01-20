function [current_run,Figure_Count,colors,darker_colors, Screening_PCE,...
    Screening_TCE,colors_solvents,colors_hydrocarbons, Sim_ID] = Plot_Colors()

Sim_ID        =    1;          % Simulation ID
Sim_ID; current_run=1; 
Figure_Count=1; 
plot_markers = {'o', 'd', 's', 'o', '^', 'd', 's'};
colors = [
    0.0000 0.4470 0.7410;  % Blue
    0.2500 0.2500 0.2500   % Gray
    0.8500 0.3250 0.0980;  % Red
    0.0000 0.4470 0.7410;  % Blue
    0.4660 0.6740 0.1880;  % Green
    0.2500 0.2500 0.2500   % Gray
    0.8500 0.3250 0.0980;  % Red
    
];

darker_colors = [
 0.0000 0.3576 0.5928;  % Darker Blue
 0.2000 0.2000 0.2000;  % Darker Gray
 0.6800 0.2600 0.0784;  % Darker Red
  0.0000 0.3576 0.5928; % Darker Blue
 0.3728 0.5392 0.1504;  % Darker Green
  0.2000 0.2000 0.2000; % Darker Gray
 0.5080 0.0624 0.1472;  % Darker Dark Red
];
Screening_PCE = [-4000 6000; 7000 6000]; % Superfund Site PCE Screening Level for Vapor ug/m3
Screening_TCE = [-4000 293; 7000 293]; % Superfund Site PCE Screening Level for Vapor ug/m3
colors_solvents = darker_colors([1,2,3,5], :);
colors_hydrocarbons = darker_colors([4, 6, 7], :); 