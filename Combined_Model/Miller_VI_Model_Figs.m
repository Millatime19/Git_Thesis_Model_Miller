clear all; 
data = load('plot.mat'); % Load the results from the plot.mat file
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
% % % % % % % %     Plot Modeled Surface Vapor & SS8 Exhaust Vapor Data    % % % % % % % %
    figure(3549)
    plot((1/86400)*G_Top(1:n,1),1e12*MW(1:ncomp)'.*G_Top(1:n,2:ncomp+1),'-') % Top Node Soil Vapor Conc. vs Time convert from mol/cc to ug/L (conversion = 1e-9)
    hold on
    SS8_Data = load("SS8_Conc_Data.csv");
    dates = SS8_Data(:, 1);Model_Start = min(SS8_Data(:,1));dates = dates - Model_Start+VIMS_Exhaust + 90; %+BC(1,3); % Extract the columns
    contaminantIDs = SS8_Data(:, 2);concentrations = SS8_Data(:, 3); % Setup IDs for loop
    uniqueContaminants = unique(contaminantIDs);colors = lines(length(uniqueContaminants));% Unique contaminant ID
        for i = 1:length(uniqueContaminants)% Loop through each unique contaminant and plot its data
            contaminant = uniqueContaminants(i);
            idx = contaminantIDs == contaminant;% Extract data for the current contaminant
            contaminantDates = dates(idx); contaminantConcentrations = concentrations(idx);
            plot(contaminantDates, contaminantConcentrations, 'o-', 'DisplayName', ['Contaminant ' num2str(contaminant)], 'Color', colors(i, :));% Plot data
        end
     SubSlabVapor_Data = load("SubSlabVapor.csv"); % 
    dates = SubSlabVapor_Data(:, 1);Model_Start = min(SubSlabVapor_Data(:,1));dates = dates - Model_Start+VIMS_Exhaust;%+BC(1,3); % Extract the columns
    contaminantIDs = SubSlabVapor_Data(:, 2);concentrations = SubSlabVapor_Data(:, 3); % Setup IDs for loop
    uniqueContaminants = unique(contaminantIDs);colors = lines(length(uniqueContaminants));% Unique contaminant ID
        for i = 1:length(uniqueContaminants)% Loop through each unique contaminant and plot its data
            contaminant = uniqueContaminants(i);
            idx = contaminantIDs == contaminant;% Extract data for the current contaminant
            contaminantDates = dates(idx); contaminantConcentrations = concentrations(idx);
            plot(contaminantDates, contaminantConcentrations, 'o-', 'DisplayName', ['Contaminant ' num2str(contaminant)], 'Color', colors(i, :));% Plot data
        end
        x_limits = get(gca, 'XLim'); % Get the current x-axis limits
    line(x_limits, [x_pce_vapor x_pce_vapor], 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal PCE')
    line(x_limits, [x_tce_vapor x_tce_vapor], 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal TCE')
    line(x_limits, [x_vc_vapor x_vc_vapor], 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'Vapor Goal VC')
    % line(x_limits, [x_toluene_vapor x_toluene_vapor], 'Color', color_Toluene, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Toluene')
    % line(x_limits, [x_xylenes_vapor x_xylenes_vapor], 'Color', color_Xylenes, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Xylenes')
    % line(x_limits, [x_ethylbenzene_vapor x_ethylbenzene_vapor], 'Color', color_Ethylbenzene, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Ethylbenzene')
    xline(x1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS ON
    xline(x2, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);   % Line for VIMS OFF
    xline(x3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Line for Surface Soil Contamination Start
    text(x1, 10, 'VIMS ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x2, 10, 'VIMS OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(x3, 10, 'Soil Contam.', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    ylabel('Concentration (mg/L)');
    if Building_Model==1
        title('Surface Vapor Concentration of Contaminants Over Time (Building Model)');
    else 
        title('Surface Vapor Concentration of Contaminants Over Time (Outside Model)');
    end
    legend('Model TCE','Model PCE','Model cis-1-2,DCE','Model Toluene','Model Vinyl Chloride','Model Xylenes',...
        'Model Ethylbenzene','Exhaust TCE','Exhaust PCE','Exhaust CIS','Exhaust Toluene','Exhaust VC','Exhaust Xylenes','Exhaust Ethylbenzene');
    set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e6]); xlim([BC(1,2),BC(3,3)]);
    % save_figure_with_sim_no(Sim_No, 'SurfVaporCompare');
    hold off;

    % Test the model output and measured data for PEST 
% figure (233)
% plot(model_time,model_data)
% hold on
% 
% plot (measured_times,new_SS8_Data(:,2:ncomp),'*-')
% set(gca, 'YScale', 'log');grid on; ylim([1e-1 2e4]); xlim([BC(1,2),BC(3,3)]);
% hold off
%     saveas(gcf, 'TestFig.svg')  
%  % 
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


% 
  % 
  %    % % % % % % % %      Vapor Concentration Plot    % % % % % % % %
  % 
  % figure(111)
  %   Cplotnew=[Gnew(:,1:(nc-1)); G(nin+1,1:(nc-1))]; % Soil Gas Concentration Plot
  %   semilogx(1e6*Cplotnew.*MW',(dx/100)*(0.5+(0:(nin))),'-o')
  %   hold on 
  %   set(gca,'Ydir','reverse')
  %   xlabel('Soil Gas Conc. [mg/L]'); ylabel('Depth (m)')
  %   axis([1e-5 1000 0 (L+dx/2)/100])
  %   legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
  %    title(['Vapor Con. - Time = ',num2str(elapsed/86400),' days'])
  %    line([x_pce_vapor x_pce_vapor], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'GW Goal PCE')
  %   line([x_tce_vapor x_tce_vapor], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal TCE')
  %   line([x_vc_vapor x_vc_vapor], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'Vapor Goal VC')
  % % line([x_toluene_vapor x_toluene_vapor], y_limits1, 'Color', color_Toluene, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Toluene')
  % % line([x_xylenes_vapor x_xylenes_vapor], y_limits1, 'Color', color_Xylenes, 'LineStyle', '--', 'DisplayName', 'Vapor Goal Xylenes')
  % % line([x_ethylbenzene_vapor x_ethylbenzene_vapor], y_limits1, 'Color', color_Ethylbenzene, 'LineStyle', '--', 'DisplayName', 'EPA Vapor Goal Ethylbenzene')
    % hold off 

%   % save_figure_with_sim_no(Sim_No, 'Vapor');
% 
% %   % % % % % % % %    Dissolved Concentration Plot    % % % % % % % %
% figure(112)
%     semilogx(W_mg_L,(dx/100)*(0.5+(0:(nin))),'-+')
%     set(gca,'Ydir','reverse')
%     axis([1e-5 1000 0 (L+dx/2)/100])
%     xlabel('Dissolved Conc. [mg/L]'); ylabel('Depth (m)')
%     title(['Diss. Conc. - Time = ',num2str(elapsed/86400),' days'])
%     hold on
%     legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
%     line([x_pce_diss x_pce_diss], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'Vapor Goal PCE')
%     line([x_tce_diss x_tce_diss], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'GW Goal TCE')
%     line([x_cis_diss x_cis_diss], y_limits, 'Color', color_Cis, 'LineStyle', '--', 'DisplayName', 'GW Goal CIS')
%     line([x_vc_diss x_vc_diss], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'GW Goal VC')
%     hold off 
%     % save_figure_with_sim_no(Sim_No, 'Diss');
% 
%     % % % % % % % %     Total Mass Concentration Plot    % % % % % % % %
%     figure(113)
%     semilogx((1e6/SDENS)*MC(:,1:ncomp).*MW',(dx/100)*(0.5+(0:(nin))),'-+')
%     set(gca,'Ydir','reverse')
%     axis([1e-5 1000 0 (L+dx/2)/100])
%     xlabel('Total Mass Conc. [mg/kg]'); ylabel('Depth (m)')
%     legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
%     title(['Tot. Mass - Time = ',num2str(elapsed/86400),' days'])  
%     hold on
%     line([x_pce_soil x_pce_soil], y_limits, 'Color', color_PCE, 'LineStyle', '--', 'DisplayName', 'Max Soil PCE')
%     line([x_tce_soil x_tce_soil], y_limits, 'Color', color_TCE, 'LineStyle', '--', 'DisplayName', 'Max Soil TCE')
%     % line([x_vc_soil x_vc_soil], y_limits, 'Color', color_VC, 'LineStyle', '--', 'DisplayName', 'Max Soil VC')
%     % line([x_toluene_soil x_toluene_soil], y_limits, 'Color', color_Toluene, 'LineStyle', '--', 'DisplayName', 'Max Soil Toluene')
%     % line([x_xylenes_soil x_xylenes_soil], y_limits, 'Color', color_Xylenes, 'LineStyle', '--', 'DisplayName', 'Max Soil Xylenes')
%     % line([x_ethylbenzene_soil x_ethylbenzene_soil], y_limits, 'Color', color_Ethylbenzene, 'LineStyle', '--', 'DisplayName', 'Max Soil Ethylbenzene')
%     hold off
%     % save_figure_with_sim_no(Sim_No, 'TotalSoilConc');