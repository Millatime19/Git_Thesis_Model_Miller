clear all; 
data = load('plot.mat'); % Load the results from the plot.mat file
PCE_VISL = 6000; % PCE VISL ug/L
TCE_VISL = 293; % TCE VISL ug/L
CIS_VISL = 5840; % PCE VISL ug/L
VC_VISL = 929; % TCE VISL ug/L
PCE_VISL_20 = 0.2 * PCE_VISL;  % 20% of 6000 ug/L
TCE_VISL_20 = 0.2 * TCE_VISL;   % 20% of 293 ug/L
CIS_VISL_20 = 0.2 * CIS_VISL;  % 20% of 5840 ug/L
VC_VISL_20 = 0.2 * VC_VISL;   % 20% of 929 ug/L

% Extract all fields into the current workspace
fields = fieldnames(data);
for i = 1:length(fields)
    assignin('base', fields{i}, data.(fields{i}));
end
data = load('Model_Export.mat');  % Load the results from the plot.mat file
% Extract all fields into the current workspace
fields = fieldnames(data); 
for i = 1:length(fields)
    assignin('base', fields{i}, data.(fields{i}));
end
nc=size(MC,2);
Fig_Position = [100, 100, 1300, 600];
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
MW1 = load("MW_1.csv"); % Load monitor well contaminant data (ug/L Dissolved Phase) 
Model_Start = min(MW1(:,1)); MW1(:,1) = MW1(:,1) - Model_Start; % Adjust the starting dates
dates = unique(MW1(:,1));contaminants = unique(MW1(:,2)); % Get unique dates and contaminants
MW_1_transformed = zeros(length(dates), length(contaminants) + 1); % Initialize new matrix for plotting
MW_1_transformed(:,1) = dates; % Fill in the dates column

% % % % % MW-1 GW Concentrations Plotting (Chlorinated Solvents % % % % %
MW1_startDate = datetime(2004, 7, 1);  % Convert days to datetime values  Define the starting date% July 1, 2004
MW_1_x_dates = MW1_startDate + days(MW_1_transformed(:,1)); % Define the starting date
MW1_startDate = datetime(2004, 7, 1);  % Convert days to datetime values  Define the starting date% July 1, 2004
MW_1_x_dates = MW1_startDate + days(MW_1_transformed(:,1)); % Define the starting date
% Setup VIMS Construction Date and Sampling Dates 
SS_Data = load("SS8_Conc_Data.csv"); % VIMs Sub-System Vapor Exhaust Concentration Data
SSV_Data = load("SubSlabVapor.csv");
% Setup VIMS Construction Date and Sampling Dates 
MW1_Start = min(MW1(:,1));MW1_dates = MW1(:, 1); MW1_dates = MW1_dates - MW1_Start;
VIMS_Install = 0; VIMS_Sample = min(SS_Data(:,1))-MW1_Start;
Start_GWSample = min(MW1_dates); VIMS_Sample_Start = VIMS_Install-MW1_Start;
Solvents_MW = [2 3 4 6]; % Used for plotting 
Hydrocarbons_MW = [5 7 8]; 
Fig_Position = [100, 100, 1300, 600];
L           =   24.5;       % Vertical domain length units in ft
L           =   L*12*2.54;  % Vertical domain length units in cm 
nin         =   15;         % Specify number of internal spatial nodes
dx          =   L/nin;      % Specify dx based on L and nin. 
for i = 1:length(dates)% Loop through each date and contaminant
    for j = 1:length(contaminants)
        
        idx = find(MW1(:,1) == dates(i) & MW1(:,2) == contaminants(j));% Find the corresponding concentration
            if ~isempty(idx)
            MW_1_transformed(i, j+1) = MW1(idx, 3);
        else
            MW_1_transformed(i, j+1) = NaN; % or 0, depending on how you want to handle missing data
        end
    end
end
MW_1_transformed = sortrows(MW_1_transformed, 1); % Sort the matrix by date

MW2 = load("MW_2.csv"); % Load monitor well contaminant data (ug/L Dissolved Phase) 
Model_Start = min(MW2(:,1)); MW2(:,1) = MW2(:,1) - Model_Start; % Adjust the starting dates
dates = unique(MW2(:,1));contaminants = unique(MW2(:,2)); % Get unique dates and contaminants
MW_2_transformed = zeros(length(dates), length(contaminants) + 1); % Initialize new matrix for plotting
MW_2_transformed(:,1) = dates; % Fill in the dates column
Solvents_MW = [2 3 4 6]; % Used for plotting 
MW_2_transformed = zeros(length(dates), length(contaminants) + 1); % Initialize new matrix for plotting
MW_2_transformed(:,1) = dates; % Fill in the dates column
for i = 1:length(dates)% Loop through each date and contaminant
    for j = 1:length(contaminants)
        
        idx = find(MW2(:,1) == dates(i) & MW2(:,2) == contaminants(j));% Find the corresponding concentration
            if ~isempty(idx)
            MW_2_transformed(i, j+1) = MW2(idx, 3);
        else
            MW_2_transformed(i, j+1) = NaN; % or 0, depending on how you want to handle missing data
        end
    end
end
MW_2_transformed = sortrows(MW_2_transformed, 1); 
SubSlabVapor_Data = load("SubSlabVapor.csv"); %  % Load monitor well contaminant data (ug/L Dissolved Phase) 
Model_Start = min(SubSlabVapor_Data(:,1)); SubSlabVapor_Data(:,1) = SubSlabVapor_Data(:,1) - Model_Start; % Adjust the starting dates
dates = unique(SubSlabVapor_Data(:,1));contaminants = unique(SubSlabVapor_Data(:,2)); % Get unique dates and contaminants
SubSlabVapor_Data_transformed = zeros(length(dates), length(contaminants) + 1); % Initialize new matrix for plotting
SubSlabVapor_Data_transformed(:,1) = dates+50; % Fill in the dates column

for i = 1:length(dates)% Loop through each date and contaminant
    for j = 1:length(contaminants)
        
        idx = find(SubSlabVapor_Data(:,1) == dates(i) & SubSlabVapor_Data(:,2) == contaminants(j));% Find the corresponding concentration
            if ~isempty(idx)
            SubSlabVapor_Data_transformed(i, j+1) = SubSlabVapor_Data(idx, 3);
        else
            SubSlabVapor_Data_transformed(i, j+1) = NaN; % or 0, depending on how you want to handle missing data
        end
    end
end
SubSlabVapor_Data_transformed = sortrows(SubSlabVapor_Data_transformed, 1); % Sort the matrix by date

Vapor_All = [SubSlabVapor_Data_transformed;VIMS_Data];


% % % % % % % %     Plot Surface Node Vapor Phase Concentrations     % % % % % % % %
% (Chlorinated Solvents) %%%%%%
    figure(199)
    
    Solvents_SSD = [1 2 3 5]; 
    plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,[2 3 4 6]),'LineWidth', 2)
    hold on
   for i = [2 3 4 6]
        % Plot each contaminant with a different marker  
        plot(Vapor_All(:,1), Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured VIMS vapor data
   end  
yline(PCE_VISL, '--', 'PCE VISL,                               ', 'LineWidth', 1.5, 'Color', [0 0 0]) % Darkest black
yline(TCE_VISL, '--', ' TCE VISL', 'LineWidth', 1.5, 'Color', [0 0 0.5]) % Darker blue
yline(CIS_VISL, '--', 'Cis-1,2-DCE VISL', 'LineWidth', 1.5, 'Color',[0.5 0 0] ) % Darker red
yline(VC_VISL, '--', 'VC VISL', 'LineWidth', 1.5, 'Color', [0 0.4 0]) % Darker green

% yline(PCE_VISL_20, '--k', '20% PCE VISL', 'Color', [0 0 0], 'LineWidth', 1.5)  % Black dashed line for PCE 20 % VISL
% yline(TCE_VISL_20, '--b', '20% TCE VISL', 'Color', [0 0 0.5], 'LineWidth', 1.5)  % Blue dashed line for TCE 20 % VISL
% yline(CIS_VISL_20, '--k', '20% CIS VISL', 'Color', [0.5 0 0], 'LineWidth', 1.5)  % Black dashed line for CIS 20 % VISL
% yline(VC_VISL_20,  '--b', '20% TC VISL',  'Color', [0 0.4 0], 'LineWidth', 1.5)  % Blue dashed line for VC 20 % VISL
hold off
    ylabel('Vapor Conc. (ug/m^3)'); xlabel('Days');
    set(gca, 'ColorOrder', colors_solvents); 
    x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    title('Vapor Phase Concentrations at Ground Surface (Combined Model)','GLM Optimized Parameter Set - Chlorinated Solvents');
    legend('Modeled TCE','Modeled PCE','Modeled cis-1-2,DCE','Modeled Vinyl Chloride', ...
    'Measured TCE','Measured PCE','Measured CIS','Measured VC',...
    'Location', 'northeast');
     set(gca, 'YScale', 'log'); grid on; ylim([1e-1 2e8]); xlim([BC(1,2), 3300]);
     % set(gca, 'YScale', 'log'); grid on; ylim([1E4 1E6]); xlim([40 60]);
     set(gcf, 'Position', Fig_Position);
     % gtext('\bf PCE Screening Level = 6,000 \mug/m^3'); gtext('\bf TCE Screening Level = 293 \mug/m^3');





%%
     % % % % % % % %   Plot of NAPL at surface and node and lowest node during simulation      % % % % % % % % 
   
    figure (118) % 
    subplot(3,1,1) 
    % Top Node NAPL History (At surface; Node 1)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,2),'o','Color',[0, 0, 0.5])
    xlim([BC(1,2),BC(3,3)]);
    hold on
    ylim([3,4])
    title('NAPL History - At Surface');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    subplot(3,1,2) % Bottom Node NAPL History (Just above water table; Node 16)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,16),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - Above Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
     subplot(3,1,3) % Bottom Node NAPL History (At water table; Node 17)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,17),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - At Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    % save_figure_with_sim_no(Sim_No, 'Phase_History');


%%

% % % % % % %      Soil Moisture Plot     % % % % % % %
figure(110)
    plot(theta(:,2),(dx/100)*(0.5+(0:(nin))),'-o')
    set(gca,'Ydir','reverse')
    axis([0 max(por) 0 (L+dx/2)/100]);
    title(['Soil Moisture - Time = ',num2str(elapsed/86400),' days'])
    xlabel('Vol. Soil Moisture'); ylabel('Depth (m)');
    % save_figure_with_sim_no(Sim_No, 'Soil Moisture');

    %%

     % % % % % % % %     Plot above WT (node 15) Vapor Phase Concentrations     % % % % % % % %
% (Chlorinated Solvents) %%%%%%
    figure(100)
    
    Solvents_SSD = [1 2 3 5]; 
    plot(G_WT_ug_m3(:,1)/86400/365,G_WT_ug_m3(:,[2 3 4 6]),'LineWidth', 2)
    hold on
   for i = [2 3 4 6]
        % Plot each contaminant with a different marker  
        plot(Vapor_All(:,1)/365, Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured VIMS vapor data
   end  
   
    %  for i = [2 3 4 6]
    %     % Plot each contaminant with a different marker
    %     plot(SubSlabVapor_Data_transformed(:,1), SubSlabVapor_Data_transformed(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured Subslab vapor data
    % end  
    plot(Screening_PCE(:,1), Screening_PCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site PCE Screening Level
    plot(Screening_TCE(:,1), Screening_TCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site TCE Screening Level
    ylabel('Vapor Conc. (ug/m^3)'); xlabel('Years');
    set(gca, 'ColorOrder', colors_solvents); 
    x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    title('Baseline  Vapor Phase Concentrations just above Water Table (Building Model)','Chlorinated Solvents, q = 2.78E-04 cm/sec');
    legend('Modeled TCE','Modeled PCE','Modeled cis-1-2,DCE','Modeled Vinyl Chloride', ...
    'Measured TCE','Measured PCE','Measured CIS','Measured VC',...
    'Location', 'northeast');
     set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e8]); xlim([BC(1,2),30000/365]);
     set(gcf, 'Position', Fig_Position);
     gtext('\bf PCE Screening Level = 6,000 \mug/m^3'); gtext('\bf TCE Screening Level = 293 \mug/m^3');
     hold off;
     %% 
% % % % % % % %     Plot Surface Node Vapor Phase Concentrations     % % % % % % % %
% (Chlorinated Solvents) %%%%%%
    figure(200)
    
    Solvents_SSD = [1 2]; 
    plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,[2 3]),'LineWidth', 2); 
    hold on
    for i = [2 3]
        % Plot each contaminant with a different marker  
        plot(Vapor_All(:,1), Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured VIMS vapor data
    end  
    plot(Screening_PCE(:,1), Screening_PCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site PCE Screening Level
    plot(Screening_TCE(:,1), Screening_TCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site TCE Screening Level
    ylabel('Vapor Conc. (ug/m^3)'); xlabel('Days');
    set(gca, 'ColorOrder', colors_solvents([1,2],:)); 
    x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    title('Baseline  Vapor Phase Concentrations at Ground Surface (Building Model)','Chlorinated Solvents, Soil Source Only');
    legend('Modeled TCE','Modeled PCE','Measured TCE','Measured PCE','Location', 'northeast');
    set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e8]); xlim([BC(1,2),3255]);
    set(gcf, 'Position', Fig_Position);
    gtext('\bf PCE Screening Level = 6,000 \mug/m^3'); gtext('\bf TCE Screening Level = 293 \mug/m^3');
    hold off;

%% 


% % % % % % % %     Plot Surface Node Vapor Phase Concentrations     % % % % % % % %
% (Chlorinated Solvents) %%%%%%
    figure(200)
    
    Solvents_SSD = [1 2 3 5]; 
    plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,[2 3 4 6]),'LineWidth', 2)
    hold on
   for i = [2 3 4 6]
        % Plot each contaminant with a different marker  
        plot(Vapor_All(:,1), Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured VIMS vapor data
    end  
    plot(Screening_PCE(:,1), Screening_PCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site PCE Screening Level
    plot(Screening_TCE(:,1), Screening_TCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site TCE Screening Level
    ylabel('Vapor Conc. (ug/m^3)'); xlabel('Days');
    set(gca, 'ColorOrder', colors_solvents); 
    x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    title('Baseline  Vapor Phase Concentrations at Ground Surface (Outside Model)','Chlorinated Solvents, Post VIMS Deactivation');
    legend('Modeled TCE','Modeled PCE','Modeled cis-1-2,DCE','Modeled Vinyl Chloride', ...
    'Measured TCE','Measured PCE','Measured CIS','Measured VC',...
    'Location', 'northeast');
     set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e8]); xlim([BC(1,2),2000]);
     set(gcf, 'Position', Fig_Position);
     gtext('\bf PCE Screening Level = 6,000 \mug/m^3'); gtext('\bf TCE Screening Level = 293 \mug/m^3');
     hold off;

     %%
% % % % % % % %     Plot Surface Node Vapor Phase Concentrations     % % % % % % % %
% (Petroleum Hydrocarbons) %%%%%%        
figure(203)
        plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,[4 6 7]),'LineWidth', 2)
        hold on
        Hydrocarbons_VIMS = [5 7 8]; 
          for i = Hydrocarbons_VIMS
        % Plot each contaminant with a different marker
        plot(SubSlabVapor_Data_transformed(:,1)+200, SubSlabVapor_Data_transformed(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured Subslab vapor data
    end  
    for i = Hydrocarbons_VIMS
        % Plot each contaminant with a different marker  
        plot(VIMS_Data(:,1), VIMS_Data(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',10,'MarkerFaceColor',darker_colors(i-1,:)+.1,'MarkerEdgeColor','black'); % Plot measured VIMS vapor data
    end  
     
        % axis([min(C_Top(:,1))/86400 max(C_Top(:,1))/86400 0 max(max(VIMS_Data(:,3)))]);set(gca, 'YScale', 'log'); 
        xlabel('Days'); ylabel('Vapor Conc. (ug/m^3)');grid on;
        title('Baseline Vapor Phase Concentrations at Ground Surface (Combined Model)','PHCs, Post VIMS Deactivation');
        set(gca, 'ColorOrder', colors_hydrocarbons);
        set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e8]); xlim([BC(1,2),3251]);
         legend('Modeled Toluene','Modeled Total Xylenes','Modeled Ethylbenzene', ...
        'Measured Toluene','Measured Xylenes','Measured Ethylbenzene','Location', 'northeast');
         set(gcf, 'Position', Fig_Position);
         hold off
%% 

% % % % % % % % Plot Groundwater Node Dissolved Phase Concentrations
% (Chorinated Solvents) %%%%%%
figure (121)
     
    plot((1/86400)*C_Top(:,1),W_GW_ug_L(:,[2 3 4 6]),'LineWidth', 2)
    hold on
    for i = [2 3 4 6]
        % Plot each contaminant with a different marker
        plot(MW_1_transformed(:,1)-4042, MW_1_transformed(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',8);
    end

    set(gca, 'ColorOrder', colors_solvents);
    set(gcf, 'Position', Fig_Position);
    axis([min(C_Top(:,1))/86400 3275 0 max(max(MW1(:,3)))]);set(gca, 'YScale', 'log'); 
    xlabel('Days'); ylabel('Dissolved Conc. (ug/L)');grid on;
    title('Dissolved Chlorinated Solvents Concentration at Water Table');
    legend('Modeled TCE','Modeled PCE','Modeled cis-1-2,DCE','Modeled Vinyl Chloride', ...
    'Measured TCE','Measured PCE','Measured CIS','Measured VC',...
    'Location', 'eastoutside');
    figure_name = 'GW_ModMeas_CS'; % Change this manually for each figure
    filename = sprintf('%s%s%d', figure_name);
    saveas(gcf,  filename);
    filename = sprintf('%s.svg', figure_name);
    saveas(gcf,  filename);
    hold off

    % % % % % % % % Plot Groundwater Node Dissolved Phase Concentrations
% (Petroleum Hydrocarbons) %%%%%%
    figure (122)
    Hydrocarbons_MW = [5 7 8]; 
    plot((1/86400)*C_Top(:,1),W_GW_ug_L(:,[5 7 8]),'LineWidth', 2)
    hold on
    for i = Hydrocarbons_MW
        % Plot each contaminant with a different marker
        plot(MW_1_transformed(:,1), MW_1_transformed(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',8);
    end
    set(gca, 'ColorOrder', colors_hydrocarbons);
    set(gcf, 'Position', Fig_Position);
    axis([min(C_Top(:,1))/86400 3275 0 max(max(MW_1_transformed(:,3)))]);set(gca, 'YScale', 'log'); 
    xlabel('Days'); ylabel('Dissolved Conc. (ug/L)');grid on;
    title('Dissolved Petroleum Hydrocarbons Concentrations at Water Table');
    legend('Modeled Toluene','Modeled Total Xylenes','Modeled Ethylbenzene', ...
    'Measured Toluene','Measured Xylenes','Measured Ethylbenzene','Location', 'northeast');
    hold off
  %% 
  

  % % % % % % % %     Plot Surface Node Vapor Phase Concentrations     % % % % % % % %
% (Chlorinated Solvents PCE and TCE Only) %%%%%%
    figure(2040)
    
    % Solvents_SSD = [1 2]; 
    Solvents_MW = [2 3];
    plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,[2 3]),'LineWidth', 2)
    hold on
    Vapor_All = [SubSlabVapor_Data_transformed;VIMS_Data]
    for i = [2 3];
        % Plot each contaminant with a different marker  
        plot(Vapor_All(:,1), Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',8); % Plot measured VIMS vapor data
    end  
    plot(Screening_PCE(:,1), Screening_PCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site PCE Screening Level
    plot(Screening_TCE(:,1), Screening_TCE(:,2),'-.','Color', [0.3, 0.3, 0.3],'LineWidth',1.5); % Plot Superfund Site TCE Screening Level
    ylabel('Vapor Concentration (ug/m^3)'); xlabel('Days');
    x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    title('Chlorinated Solvents Vapor Phase Concentrations at Surface');
    legend('Modeled TCE','Modeled PCE','Measured TCE','Measured PCE',...
    'Location', 'northeast');
     set(gca, 'YScale', 'log');grid on; ylim([1e1 2e8]); xlim([BC(1,2),BC(3,3)-4]);
     set(gcf, 'Position', Fig_Position);
     gtext('PCE Screening Level = 6,000 \mug/m^3'); gtext('TCE Screening Level = 293 \mug/m^3');
    set(gca, 'ColorOrder', colors_solvents(1:2,:));
     hold off;


       %% 
       
% % % % % % % %     Plot All Surface Node Vapor Phase Concentrations     % % % % % % % %
% (All Contaminants) %%%%%%     
figure(204)
colors = colors([1,2,3,4,5,6,7], :);
        plot(G_Top_ug_m3(:,1)/86400,G_Top_ug_m3(:,2:end),'LineWidth', 2)
        hold on
        % plot(VIMS_Data(:,1), VIMS_Data(:,2:end),'--o','LineWidth', 1, 'MarkerSize',8);
         for i = [2 3 4 5 6 7 8]
        % Plot each contaminant with a different marker  
            plot(Vapor_All(:,1), Vapor_All(:,i), [plot_markers{i-1} '--'],'LineWidth', 1, 'MarkerSize',8); % Plot measured VIMS vapor data
        end  
        set(gca, 'YScale', 'log');  xlabel('Days'); ylabel('Vapor Conc. (ug/m^3)');grid on;
        title('Petroleum Hydrocarbons Vapor Phase Concentrations at Surface');
        set(gca, 'ColorOrder', colors);
        set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e6]); xlim([BC(1,2),BC(3,3)-4]);
            legend('Modeled Toluene','Modeled Total Xylenes','Modeled Ethylbenzene', ...
        'Measured Toluene','Measured Xylenes','Measured Ethylbenzene','Location', 'eastoutside');
        set(gcf, 'Position', Fig_Position);
            hold off
    
    %% 
    


     % % % % % % % %   Plot of NAPL at surface and node and lowest node during simulation      % % % % % % % % 
   
    figure (118) % 
    subplot(3,1,1) 
    % Top Node NAPL History (At surface; Node 1)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,2),'o','Color',[0, 0, 0.5])
    xlim([BC(1,2),BC(3,3)]);
    hold on
    ylim([3,4])
    title('NAPL History - At Surface');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    subplot(3,1,2) % Bottom Node NAPL History (Just above water table; Node 16)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,16),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - Above Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
     subplot(3,1,3) % Bottom Node NAPL History (At water table; Node 17)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,17),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - At Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    % save_figure_with_sim_no(Sim_No, 'Phase_History');

%% 

% % % % % Plot the model output and measured data for PEST vs time -
% % % % %  DIRECTLY FROM PEST REC FILE
% Load data from CSV file copy paste from rec file. 
data = readtable('PEST_Results_recfile.csv'); % Make sure the file path is correct

% Extract relevant columns
Names = data.Name; % Observation names
Measured = data.Measured; % Measured values
Modeled = data.Modelled; % modeled values

% Identify NaN values in the Measured data and update modeled values
nanIdx = Measured == 999999; % Assuming 999999 represents NaN values in Measured
Modeled(nanIdx) = NaN; % Set corresponding modeled values to NaN
figure (453)
hold on; grid on;
% Determine the range for the 1:1 line and plot 
minValue = min([Measured; Modeled], [], 'omitnan'); % Minimum value of both arrays, ignoring NaNs
maxValue = max([Measured; Modeled], [], 'omitnan'); % Maximum value of both arrays, ignoring NaNs
plot([minValue, maxValue], [minValue, maxValue], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Line');
xlabel('Measured'); ylabel('modeled');
set(gca, 'XScale', 'log', 'YScale', 'log');
scatter(Measured(:,1), Modeled(:,1))
legend('TCE', 'PCE', 'CIS', 'Toluene', 'VC', 'Xylenes', 'Ethylbenzene', 'Location', 'best');
title('Comparison of Modeled and Measured Observations with Time (PEST .rec File)');


% Create a time vector assuming sequential time points (you can adjust as needed)
Time = (1:length(Measured))'; % Assuming time is sequential based on the data length
% Define colors for each contaminant type (1-7)
colors = lines(7); % Generates 7 distinct colors
% Identify NaN values in the Measured data and update modeled values
nanIdx = Measured == 999999; % Assuming 999999 represents NaN values in Measured
Measured(nanIdx) = NaN; 
Modeled(nanIdx) = NaN; % Set corresponding modeled values to NaN


% Prepare figure for Measured vs. Time and Modeled vs. Time using .rec file
% data 
figure (4934); %
hold on;
scatter(Time, Measured, 60, 'b', 'filled', 'DisplayName', 'Measured');% Plot Measured vs. Time
scatter(Time, Modeled, 60, 'r', 'filled', 'DisplayName', 'Modeled');% Plot Modeled vs. Time
grid on; xlabel('Time'); ylabel('Concentration');set(gca, 'YScale', 'log');
title('Measured and Modeled vs. Time');
legend('Measured', 'Modeled', 'Location', 'best');% Add legend to the second plot
title('Comparison of Modeled and Measured Observations with Time (PEST .rec File)');

hold off;
    % figure_name = 'PlotObsTimeModelTest'; % Change this manually for each figure
    % filename = sprintf('%s%s%d', figure_name, current_run);
    % saveas(gcf,  filename);
    % hold off;
%% 
% % %  % Plot the model output and measured data for PEST vs time
figure (201)
model_data = Model_Data_Export(:, 2:end); % Extract the analytical data from columns 2 to the end and combine arrays
% for plotting
vims_data = VIMS_Data(:, 2:end);  new_vims_data = vims_data; new_model_data = model_data;
nan_positions = isnan(vims_data); % Find the positions of NaNs in vims_data
new_model_data(nan_positions) = NaN; % Replace corresponding positions in model_data with NaN
plot(Model_Data_Export(:, 1),new_model_data(:,:),'-o','DisplayName', 'Modeled');
hold on
plot(VIMS_Data(:, 1), VIMS_Data(:, 2:end), '*',  'DisplayName', 'Measured');
set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e4]); xlim([BC(1,2),BC(3,3)]);
title([current_run, ' Comparison of Modeled and Measured Observations with Time (Model Test)']);
legend('TCE', 'PCE', 'cis-1,2-DCE', 'Toluene', 'Vinyl Chloride', 'Xylenes', 'Ethylbenzene', 'Location', 'best');
hold off
    figure_name = 'PlotObsTimeModelTest'; % Change this manually for each figure
    filename = sprintf('%s%s%d', figure_name, current_run);
    saveas(gcf,  filename);
    hold off;

% % %  % Plot the model output and measured data for PEST 
figure(202); % Plot each column pair on a single plot
num_columns = size(new_model_data, 2); % Get the number of columns
hold on; % Keep all data on the same plot
Obs_VIMS = []; % Store x-values (VIMS data)
Obs_Model = []; % Store y-values (Model data)

% Loop through each column to plot
for i = 1:num_columns
    % Plot each pair of data
    scatter(new_vims_data(:, i), new_model_data(:, i), 'filled', 'DisplayName', ['Column ', num2str(i)]);
    
    % Save the plotted data
    Obs_VIMS = [Obs_VIMS; new_vims_data(:, i)];
    Obs_Model = [Obs_Model; new_model_data(:, i)];
end
minValue = min([Obs_VIMS; Obs_Model], [], 'omitnan');
maxValue = max([Obs_VIMS; Obs_Model], [], 'omitnan');
plot([minValue, maxValue], [minValue, maxValue], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Line'); % Plot the 1:1 line
set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');% Set log scale for both axes
xlabel('Model Data');ylabel('VIMS Data');
title([current_run, ' Comparison of Model Data vs VIMS Data (Model Test)']);
legend('TCE', 'PCE', 'cis-1,2-DCE', 'Toluene', 'Vinyl Chloride', 'Xylenes', 'Ethylbenzene', 'Location', 'best');
grid on;hold off;
    figure_name = 'CompareObsModelTest'; % Change this manually for each figure
    filename = sprintf('%s%s%d', figure_name, current_run);
    saveas(gcf,  filename);
    hold off;

%% 
% Plot the measured times vs model ouputs. 
Obs_VIMS_Model = [Obs_VIMS; Obs_Model];
% Get the number of columns in Reshape_Model_PEST_Export
[~, num_columns] = size(Obs_VIMS_Model);

% Create a new figure
figure (7654);

% Hold on to plot multiple series
hold on;

% Create a cell array to store legend entries
legendEntries = cell(1, num_columns);

% Loop through each column of Reshape_Model_PEST_Export
for i = 1:num_columns
    % Plot column 1 of measured_times against column i of Reshape_Model_PEST_Export
    plot(measured_times(:,1), Obs_VIMS_Model(:,i), 'o-');
    
    % Create legend entry for this series
    legendEntries{i} = ['Model Output ', num2str(i)];
end
% minValue = min([Measured; Modeled], [], 'omitnan'); % Minimum value of both arrays, ignoring NaNs
% maxValue = max([Measured; Modeled], [], 'omitnan'); % Maximum value of both arrays, ignoring NaNs
% plot([minValue, maxValue], [minValue, maxValue], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Line');
% Add labels and title
xlabel('Measured Times');
ylabel('Model Outputs');
title('Measured Times vs Model Outputs');

% Add legend
legend(legendEntries, 'Location', 'best');

% Add grid
grid on;

% Release the hold
hold off;

% Adjust the figure size
set(gcf, 'Position', [100, 100, 800, 600]);
plot(VIMS_Data_New(:, 1), VIMS_Data_New(:, 2:end), '*-', 'LineWidth', 2, 'DisplayName', 'Measured');
hold off; 








    %% 

%  % Degradation Plot 
% % figure (43)
% %     plot(W_Deg_plot(1:n,1)/86400,W_Deg_plot(1:n,2:end),'-o')
% %     hold on
% %     xlabel('Time (days)');
% %     ylabel('Change in Dissolved Phase Concentration (ug/L)');
% %     legend('PCE', 'TCE', 'CIS', 'VC');
% %     title('Change in COC Concentration Over Time - Semi-Implicit Degradation Method');
% %     grid on;
% %     hold off
% %     W_deg_sum = sum(W_Deg_plot(1:n,2:end))




     % % % % % % % %      Vapor Concentration Plot    % % % % % % % %

  figure(111)
    Cplotnew=[Gnew(:,1:(nc-1)); G(nin+1,1:(nc-1))]; % Soil Gas Concentration Plot
    semilogx(1e12*Cplotnew.*MW',(dx/100)*(0.5+(0:(nin))),'-o')
    hold on 
    set(gca,'Ydir','reverse')
    xlabel('Soil Gas Conc. [mg/L]'); ylabel('Depth (m)')
    axis([1e-2 1e9 0 (L+dx/2)/100])
    legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
     title(['Vapor Con. - Time = ',num2str(elapsed/86400),' days'])
     line([x_pce_vapor x_pce_vapor], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'GW Goal PCE')
    line([x_tce_vapor x_tce_vapor], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal TCE')
    line([x_vc_vapor x_vc_vapor], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'Vapor Goal VC')
  % line([x_toluene_vapor x_toluene_vapor], y_limits1, 'Color', color_Toluene, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Toluene')
  % line([x_xylenes_vapor x_xylenes_vapor], y_limits1, 'Color', color_Xylenes, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Xylenes')
  % line([x_ethylbenzene_vapor x_ethylbenzene_vapor], y_limits1, 'Color', color_Ethylbenzene, 'LineStyle', '--', 'DisplayName', 'EPA Vapor Goal Ethylbenzene')
    hold off 

%   % save_figure_with_sim_no(Sim_No, 'Vapor');
% 

%% 

    % % % % % % % %     Total Mass Concentration Plot    % % % % % % % %
    figure(1133)
    semilogx((1e6/SDENS)*MC(:,1:ncomp).*MW',(dx/100)*(0.5+(0:(nin))),'-+')
    set(gca,'Ydir','reverse')
    axis([1e-5 1000 0 (L+dx/2)/100])
    xlabel('Total Mass Conc. [mg/kg]'); ylabel('Depth (m)')
    legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
    title(['Tot. Mass - Time = ',num2str(elapsed/86400),' days'])  
    hold on
    % line([x_pce_soil x_pce_soil], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'Max Soil PCE')
    % line([x_tce_soil x_tce_soil], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'Max Soil TCE')
    % line([x_vc_soil x_vc_soil], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'Max Soil VC')
    % line([x_toluene_soil x_toluene_soil], y_limits, 'Color', color_Toluene, 'LineStyle', '--', 'DisplayName', 'Max Soil Toluene')
    % line([x_xylenes_soil x_xylenes_soil], y_limits, 'Color', color_Xylenes, 'LineStyle', '--', 'DisplayName', 'Max Soil Xylenes')
    % line([x_ethylbenzene_soil x_ethylbenzene_soil], y_limits, 'Color', color_Ethylbenzene, 'LineStyle', '--', 'DisplayName', 'Max Soil Ethylbenzene')
    hold off
    % save_figure_with_sim_no(Sim_No, 'TotalSoilConc');
    %% 
figure(113)
% Create 2x2 subplot layout
darkOrange = [0.8500 0.3250 0.0980]; % Dark orange color
darkBlue = [0 0.4470 0.7410];        % Dark blue color

% PCE subplot (2)
subplot(2,2,1)
semilogx((1e6/SDENS)*MC(:,2).*MW(2),(dx/100)*(0.5+(0:(nin))),'-+','Color',darkBlue)
set(gca,'Ydir','reverse')
axis([1e-5 1000 0 (L+dx/2)/100])
xlabel('Total Mass Conc. [mg/kg]')
ylabel('Depth (m)')
legend('PCE')
title(['PCE - Time = ',num2str(elapsed/86400),' days'])

% TCE subplot (1)
subplot(2,2,2)
semilogx((1e6/SDENS)*MC(:,1).*MW(1),(dx/100)*(0.5+(0:(nin))),'-+','Color',darkBlue)
set(gca,'Ydir','reverse')
axis([1e-5 1000 0 (L+dx/2)/100])
xlabel('Total Mass Conc. [mg/kg]')
ylabel('Depth (m)')
legend('TCE')
title(['TCE - Time = ',num2str(elapsed/86400),' days'])

% CIS subplot (3)
subplot(2,2,3)
semilogx((1e6/SDENS)*MC(:,3).*MW(3),(dx/100)*(0.5+(0:(nin))),'-+','Color',darkBlue)
set(gca,'Ydir','reverse')
axis([1e-5 1000 0 (L+dx/2)/100])
xlabel('Total Mass Conc. [mg/kg]')
ylabel('Depth (m)')
legend('CIS')
title(['CIS - Time = ',num2str(elapsed/86400),' days'])

% VC subplot (5)
subplot(2,2,4)
semilogx((1e6/SDENS)*MC(:,5).*MW(5),(dx/100)*(0.5+(0:(nin))),'-+','Color',darkBlue)
set(gca,'Ydir','reverse')
axis([1e-5 1000 0 (L+dx/2)/100])
xlabel('Total Mass Conc. [mg/kg]')
ylabel('Depth (m)')
legend('VC')
title(['VC - Time = ',num2str(elapsed/86400),' days'])

% Overall title
sgtitle({['Domain Total Mass Concentrations'], ...
    ['Time = ',num2str(elapsed/86400),' days']}, ...
    'FontWeight', 'bold')
%% 

    %   % % % % % % % %    Dissolved Concentration Plot    % % % % % % % %
figure(112)
    semilogx(W_ug_L,(dx/100)*(0.5+(0:(nin))),'-+')
    set(gca,'Ydir','reverse')
    axis([1e-5 2000 0 (L+dx/2)/100])
    xlabel('Dissolved Conc. [mg/L]'); ylabel('Depth (m)')
    title(['Diss. Conc. - Time = ',num2str(elapsed/86400),' days'])
    hold on
    legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
    % line([x_pce_diss x_pce_diss], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal PCE')
    % line([x_tce_diss x_tce_diss], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'GW Goal TCE')
    % line([x_cis_diss x_cis_diss], y_limits, 'Color', color_Cis, 'LineStyle', '--', 'DisplayName', 'GW Goal CIS')
    % line([x_vc_diss x_vc_diss], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'GW Goal VC')
    hold off 


%% 
% % % % % % % %   Plot of Sum of Dissolved Conc. at node above water table during simulation vs Infiltration      % % % % % % % %    


    % subplot(2,1,2)
    % 
    % semilogy((1/86400)*C_Top(1:n,1),C_Top(:,2))
    % hold on
    % xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    % xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    % xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start 
    % text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    % text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    % text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    % 
    % plot(I(:,1)/86400,86400*I(:,2))
    % xlabel('Time (d)'); ylabel('Rain rate (cm/d)');
    % xlim([0 ttot/86400])
    % title('Infiltration')
    % save_figure_with_sim_no(Sim_No, 'Diss_bottom');    
%% 

     % % % % % % % %   Plot of NAPL at surface and node and lowest node during simulation      % % % % % % % % 
   
    figure (118) % 
    subplot(3,1,1) 
    % Top Node NAPL History (At surface; Node 1)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,2),'o','Color',[0, 0, 0.5])
    xlim([BC(1,2),BC(3,3)]);
    hold on
    ylim([3,4])
    title('NAPL History - At Surface');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    subplot(3,1,2) % Bottom Node NAPL History (Just above water table; Node 16)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,16),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - Above Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
     subplot(3,1,3) % Bottom Node NAPL History (At water table; Node 17)
    plot((1/86400)*PHASE_History(:,1), PHASE_History(:,17),'o','Color',[0, 0.5, 0])  
    hold on
    xlim([BC(1,2),BC(3,3)]);
    ylim([3,4])
    title('NAPL History - At Water Table');
    xlabel('Days');
    ylabel('Phases Present');
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off
    % save_figure_with_sim_no(Sim_No, 'Phase_History');
 
    CtotCurrent_All = sum(1e6/SDENS)*MC(:,1:ncomp).*MW';  %  Ctot for each node and contaminant at the end of the simulation
    for i = 1:nin+1
    CtotCurrent_Sum(i) = sum(CtotCurrent_All(i,:)); % Sum of Ctot at end of simulation for each node
    end 
    CtotCurrent_Sum = CtotCurrent_Sum';
%% 

 % % % % % % %      Soil Moisture Plot     % % % % % % %
figure(110)
    plot(theta(:,2),(dx/100)*(0.5+(0:(nin))),'-o')
    set(gca,'Ydir','reverse')
    axis([0 max(por) 0 (L+dx/2)/100]);
    title(['Soil Moisture - Time = ',num2str(elapsed/86400),' days'])
    xlabel('Vol. Soil Moisture'); ylabel('Depth (m)');
    % save_figure_with_sim_no(Sim_No, 'Soil Moisture');

      % % % % % % % %     Mass Extraction Rate during Simulation Plot   % % % % % % % %
%     figure(116)
%     title(['Mass Extraction Rate - Time = ',num2str(elapsed/86400),' days (Inside)'])
%     subplot(2,1,1)
%     semilogy((1/86400)*massremoved(1:n-1,1),(MW(1)/dx)*diff(massremoved(1:n,2:3))./diff(massremoved(1:n,1)))  % PCE
%     hold on
%     xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
%     xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
%     xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
%     text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%     xlabel('Time (d)'); ylabel('Mass extraction rate (gm/cc/d)')
%     legend('TCE', 'PCE')
%     title('Mass Extraction Rate vs Time')
%     hold off
%     subplot(2,1,2)
%     plot(I(:,1)/86400,86400*I(:,2))
%     xlabel('Time (d)'); ylabel('Rain rate (cm/d)');
%     xlim([0 ttot/86400])
%      title('Infiltration vs Time')
%     % save_figure_with_sim_no(Sim_No, 'MassExtractInfiltration');

 % % % % % % % %     Mass Removed during Simulation Plot   % % % % % % % % 
 figure(117)

    semilogy((1/86400)*massremoved(1:n,1),(MW(1)/dx)*massremoved(1:n,1:7))
    hold on
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    xlabel('Time (d)'); ylabel('Cumulative Mass extraction (gm/cc soil)')
    title(['Mass Extracted - Time = ',num2str(elapsed/86400),' days (Inside)'])
    legend('TCE', 'PCE','CIS','Toluene', 'VC', ' Tot Xy' ,'Ethylbenzene')
    hold off
    % save_figure_with_sim_no(Sim_No, 'MassExtract');
  %% 

     
%     % % % % % % % %   Plot of Ctot at current time (end of simulation)    % % % % % % % % 
%     figure(119)
%     semilogx(CtotCurrent_Sum,(dx/100)*(0.5+(0:(nin))),'-+')
%     set(gca,'Ydir','reverse')
%     axis([(min(CtotCurrent_Sum(:,1))) 1000 0 (L+dx/2)/100])
%     xlabel('Total Soil Conc. [mg/kg]'); ylabel('Depth (m)')
%          legend('Total Mass Concentration');
%          title(['Modeled Tot. Mass = ',num2str(elapsed/86400),' days (Inside-VIMS)'])
%     CtotCurrent_Top = sum(1e6*MW(1:ncomp)'.*G_Top(n,2:ncomp+1)) % 
%     save_figure_with_sim_no(Sim_No, 'C_tot_Nodes');
% 
% 
%  % % % % % % % %   Plot of Ctot for all contaminants top during simulation      % % % % % % % % 
% figure (120)
%     subplot (2,1,1)
%     semilogy((1/86400)*C_Top(1:n,1),C_Top_All(1:n,1:ncomp),'-')
%     hold on
%     ylabel('Total Mass Conc. (mg/kg)'); xlabel('Time (Days)');
%     xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
%     xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
%     xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start 
%     text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%     legend('Model TCE','Model PCE','Model cis-1-2,DCE','Model Toluene','Model Vinyl Chloride','Model Xylenes')
%     title(['Tot. Mass Top Node - Time = ',num2str(elapsed/86400),' days (Inside - Baseline)'])
%     axis([min(C_Top(:,1))/86400 max(C_Top(:,1))/86400 1E-5 max(max(C_Top_All(:,:)))])
%     hold off
%     % xlim([min(C_Top(:,1))/86400 max(C_Top(:,1))/86400]);ylim([1E-5 max(max(C_Top_All(:,:)))])
%     subplot (2,1,2)
%     semilogy((1/86400)*C_Top(1:n,1),C_Bottom_All(1:n,1:ncomp),'-')
%     hold on 
%     xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
%     xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
%     xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start 
%     text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%     text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%     ylabel('Total Mass Conc. (mg/kg)'); xlabel( 'Time (Days)');
%     legend('Model TCE','Model PCE','Model cis-1-2,DCE','Model Toluene','Model Vinyl Chloride','Model Xylenes')
%     axis([min(C_Top(:,1))/86400 max(C_Top(:,1))/86400 1E-5 max(max(C_Bottom_All(:,:)))]);
%     title(['Tot. Mass Bottom Node - Time = ',num2str(elapsed/86400),' days (Inside - Baseline)'])
%     hold off
%     % save_figure_with_sim_no(Sim_No, 'C_tot_TopBottom');
