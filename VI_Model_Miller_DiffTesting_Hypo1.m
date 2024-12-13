clear all;
[current_run,Figure_Count,colors,darker_colors, Screening_PCE,...
    Screening_TCE,colors_solvents,colors_hydrocarbons, Sim_ID] = Plot_Colors(); 

% Boundary Conditions and Sampling Dates
L           =     24.5;       % Vertical domain length units in ft
L           =     L*30.48;    % Convert vertical domain length units to cm 
nin         =     15;         % Specify number of internal spatial nodes
dx          =     L/nin;      % Specify dx based on L and nin. 
n           =     0;          % Start at t = 0. Loop begins with "n = n + 1". 
Building_Model =  1;          % Building Model = 1, Outside Model = 0
PEST_Optimized =  0;          % Set 1 if optimized PEST parameters available. 
SoilhasExecuted = 1;          % Initialize as 1. Sets contaminated soil to only trigger once, then changes to 0. 
[Contam_Soil ,VIMS_Present, I_Present, Bare_ground,VIMS_Sample,Bio_Deg] = Set_Model(Building_Model); % Load Model Type Conditions
% Read initial parameters from 'Parameters.csv'
param_data = readtable('Parameters.csv', 'ReadVariableNames', true);
for i = 1:height(param_data)
    assignin('base', param_data.Parameter{i}, param_data.Value(i));
end

% Update parameters with PEST values if file exists
if isfile('Model_Input_All.txt')
    values = dlmread('Model_Input_All.txt');
    PEST_Params = {'por', 'por_v_min','FOC', 'airdiff', 'K_sat', 'BC_lambda', 'P_air', ...
                   'Conc_Soil_1', 'Conc_Soil_2', 'ctot_16','theta_min','q_air_Slab','q_air_VIMS','q_air_Open','q_var'};
    for i = 1:min(length(values), length(PEST_Params))
        assignin('base', PEST_Params{i}, values(i));
    end

end

if PEST_Optimized ==1 
values = dlmread('Optimization_B0.txt');                   % Read values from Model_Input.txt
PEST_Params_Opt = {'por','FOC',  'airdiff', 'K_sat', 'BC_lambda', 'P_air','Conc_Soil_1', 'Conc_Soil_2','ctot_16','Finalphi'};  % Values for optimized parameters for Outside PEST - All Params
for i = 1:length(values)
     assignin('base', PEST_Params_Opt{i}, values(i)); 
end
end 
%  Assign Ctot to the base workspace based on parameter inputs for each node. 
Ctot = [ctot_1 ctot_2 ctot_3 ctot_4 ctot_5 ctot_6 ctot_7 ctot_8 ctot_9 ctot_10 ctot_11 ctot_12 ctot_13 ctot_14 ctot_15 ctot_16]';

% % Make Ctot only GW Source 
% Ctot(1:nin,1) = 0;  
% Load Input Parameters
filename = 'Screening_Info_2.csv'; param_data = readtable(filename, 'ReadVariableNames', true); % Load the CSV file
for i = 1:height(param_data)
    param_name = param_data.Parameter{i}; param_value = param_data.Value(i); assignin('base', param_name, param_value); 
end

% Initialize parameter variables 
por = por*ones(nin+1,1);                                                    % Specify total porosity at nodes
theta_max = por(1) - por_v_min; theta=zeros(nin+1,2);    % Inital moisture content max, current. 
theta(:,2) =theta_initial*ones(nin+1,1);                                    % Inital water content throughout domain 
theta(:,1) = dx*(0.5:1:(nin+0.5))';                                         % Vertical distance for each node
theta(nin+1,2)=theta_max;                                                   % Max soil moisture at water table
D_nodes=zeros(nin+1,1); D_half=zeros(nin,1);                                % Setup Effective Diffusion array for finite difference approx. 
lwaves=zeros(1,3); twaves=zeros(1,3);                                       % Initialize for the kinematic wave infiltration function

% Setup Infiltration Array if Infiltration present
if I_Present ==1                                                % If I present then load in infiltration. 
I = load('Site_precip_cm.csv'); I = flipud(I);                  % Load infiltration array (Col 1 = days, Col 2 = cm/day). 
I(:,1)=86400*I(:,1);I(:,2)=(1/86400)*I(:,2); I = I(:,1:2);      % Convert days to seconds and cm/day to cm/sec
end 

% Put in hydrostatic soil moisture according to Brooks-Corey:
apple=find(theta(:,1)<(L-P_air));                               % Nodes above air-entry pressure
theta(apple,2)=theta_min+(theta_max-theta_min)*(P_air./(L-theta(apple,1))).^BC_lambda;
theta(nin+1,2)=theta_max;                                       % Water table

% Function to initialize chemical data and arrays
MASSF_soil_PEST = [Soil_Frac_1,Soil_Frac_2,Soil_Frac_3,Soil_Frac_4,Soil_Frac_5,Soil_Frac_6,Soil_Frac_7]';
MASSF_GW_PEST = [GW_Frac_1,GW_Frac_2,GW_Frac_3,GW_Frac_4,GW_Frac_5,GW_Frac_6,GW_Frac_7]';
[MC,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,MASSF_GW,nodes_contam,nc,ncomp,G,G_new,W,...
         G_Top,G_WT,W_WT,C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
         CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,elapsed]...
         = Initialize_Model(MASSF_soil_PEST, MASSF_GW_PEST, por,theta,dx,nin,DENSLIQ,SDENS,...
         TEMP,FOC,Ctot);% Initalize arrarys related to phase concentrations using inital parameters.

% Function to setup stress periods
[BC, dt, dt_max, ttot, VIMS_Exhaust, n_Soil_Contam,ntsteps,W_Deg] = Stress_Periods(Conc_Soil_1,Conc_Soil_2,Contam_Soil,...
    VIMS_Present, VIMS_Sample,Bare_ground,nc);                      % Function to setup time stress periods and contamination based on given criteria. 

if Building_Model == 1 
    if n==0 % Slab condition 
       q_air = q_air_Slab; 
    end 
    if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % Slab condition 
       q_air = q_air_Slab; 
    end
    if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
     q_air = q_air_VIMS;
    end
end   

if Building_Model == 0 
     q_air = q_air_Open; 
end



% tic
% % % % % % % % % % % % % % % % % % Main loop % % % % % % % % %  % % % % % % % % % 
while elapsed < ttot
  n=n+1;
    % Figure out water contents, with infiltration
    if I_Present == 1 % If I_present = 1, induce infiltration 
    [theta,lwaves,twaves] = Kinwave_Theta(theta_max,theta_min,K_sat,BC_lambda,P_air,theta, ...
                                          lwaves,twaves,elapsed,dt,I,L);
    end
    theta(nin+1,2)=theta_max; % Force saturation at water table
   
    if Bio_Deg == 1
    % % % Calculate Biodegradation based on first order rate constants for
    % chlorianted solvents. % % 
    [W,W_Deg]=Biodeg_Solvents(W,lambda_years,elapsed,nin); 

    end 
    % % % Calculate the multi-phase equilibrium at each node. % % 
    % Molar concentrations (mol/cm^3) for: G is gas conc (mol/cc), W is dissolved (mol/cc),
    % MC is total (mol/cc).  Then calculate void space and Deff at all nodes.
    [Retardation,totmass,G,W,MC,D_nodes,void,PHASE,calc_H,theory_H] = Equilibrium_Function(G,W,MC,por,por_v_min,MW,VP,ACT,KD,SOL,...
    airdiff,DENSLIQ,SDENS,TEMP);
    G_old=G;                                % Update gas diffusion with previous timestep +

    % % % Diffuse vapor using finite differences.   % 
    % D_nodes=D_nodes*0+0.4;                        % Force effective diffusion. Set for variable diffusion 
    % D_nodes(nin+1) = 0.02; 
    % D_nodes=logspace(2,-2,nin+1)';                % Set for variable diffusion
    q_send=q_air*(1+2*q_var*(-0.5+rand()));         % This sends Q+-20% 
    [G_new] = Vapor_Diffusion(G,D_nodes,dt,dx,q_send); 
    G(1:nin,:)=G_new(:,:);                          % Update G with G_New
    flux=(.5/dx)*(D_nodes(1:end-1)+D_nodes(2:end)).*diff(G);
    massremoved(n,1) = elapsed; 
    massremoved(n+1,2:end) = massremoved(n,2:end)+dx*q_send*G_new(1,1:ncomp); 
       
    % Update Q if VIMS is on or off. 
    if Building_Model == 1                          
      if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % Slab condition 
       q_air = q_air_Slab; 
       end
        if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
        q_air = q_air_VIMS; 
        end
        if elapsed/86400 > BC(2,3)  % Turn off VIMS  
        q_air = q_air_Slab; 
        end
    end   

    % q_air = 2.73e-4;%1e-3; % Remove - used for testing q. 1e-3
    % Advect the soil water concentrations if Infiltration is present 
    % if I_Present ==1
    % [W_new] = Water_Advect(theta,W,K_sat,theta_max,theta_min,BC_lambda,I,elapsed,dt,dx);
    % % W_Save=(W_new(1:nin,1:ncomp)-W(1:nin,1:ncomp));
    % % MC(1:nin,1:ncomp)=MC(1:nin,1:ncomp)+W_Save;
    % end 
    % Figure out change in gas molar conc and remove from total molar conc.
    G_save=(G_new-G_old(1:nin,:)).*void(1:nin);
    MC(1:nin,:)=MC(1:nin,:)+G_save;
    MC(:,end)=(1/18)*theta(:,2);  % Keep track of water moles per cm^3 soil.
   
    % Update Soil Boundary Condition 
    if Contam_Soil == 1 && elapsed > 48*86400 && SoilhasExecuted ==1
    Surf_Contam = zeros(nodes_contam,1);  
    Ctot_Surf(1:length(Surf_Contam)) = linspace(BC(2,4),BC(2,5),nodes_contam); % Avg soil contamination
    elapsed_Soil_Contam = elapsed/86400;
    % Labels based on VIMS operation and soil contamination dates. 
    x3 = elapsed_Soil_Contam; % Used for plotting label. 
        for i = 1:nodes_contam
        MC(i,1:ncomp) = Ctot_Surf(i)*(SDENS/1e6)*MASSF_soil(1:ncomp)./MW(1:ncomp); % Update molar concentrations using total mass conc. from contaminated surface soil. (Nodes 1-4)
        end
    SoilhasExecuted = 0; % Set the flag to 0 so this block doesn't run again
    end
    
    % % % % % % Save Phase concentrations for plotting and update time % % % % % %
   [G_Top, G_Top_ug_m3, G_ug_m3, G_WT, G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_WT,W_ug_L, W_GW_ug_L,  CtotCurrent_All, CtotCurrent_Sum, dt, elapsed] = ...
    Save_Model_Data(G_Top, G_WT, W_WT,C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    elapsed, ncomp,G_new,n,SDENS,MC,MW,PHASE,W,dt,dt_max,nin);
    
    if Bio_Deg == 1 
    W_Deg_plot(n,1) = elapsed; 
    W_Deg_plot(n,2:end) = W_Deg*1e6.*MW([2 1 3 5],:)';    % Store the change in concentration after degradation occurs. 
    end
        if min(por-theta(:,2))+1e-6 < por_v_min

        disp('Pore Vapor Min Exceeded End')

        end 


 
     
     % %     % % % % % % %      Kin Wave Plot     % % % % % % %
     % figure(321)
     % plot(twaves(:,2),L-twaves(:,1),'bd')
     % hold on
     % plot(lwaves(:,2),L-lwaves(:,1),'ro')
     % plot(theta(:,2),L-theta(:,1),'+-k')
     % axis([0 theta_max min(L-theta(:,1)) L])
     % ylabel('Height above WT (cm)'); xlabel('Vol. Soil Moisture')
     % title(['Lead and Trailing Waves; Time = ',num2str(elapsed/86400),' days'])
     % legend('Trailing wave','Leading wave')
     % drawnow
     % hold off

end 


[Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta, Conc_Soil_1, Conc_Soil_2, ctot_16, K_sat, BC_lambda,...
    P_air, por, airdiff,FOC, n, ncomp, L, G_Top, G_Top_ug_m3,G_WT,G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All,...
    PHASE_History, W_ug_L, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
    MW, VIMS_Exhaust, G_new, G_ug_m3, nin, BC, Building_Model, Sim_ID, dx, SDENS, MC); 

    % %% 
    % % figure(333) Plot G and Flux
    % figure(333)
    % % plot(G(:,1:(nc-1)),L-dx*(0:(nin)),'-o'); hold on
    % % plot (flux,L-dx*(0:(nin-1)),'-x')
    % plot (flux(:,2),L-dx*(0:(nin-1)),'-x'); hold on
    % plot(G_ug_m3(:,2),L-dx*(0:(nin-1)),'-o');
    % legend('Conc.','flux')
    % % set(gca, 'XScale', 'log');
    
    % % axis([0 5*max(G(nin+1,1:ncomp)) 0 (nin)*dx])
    % legend( 'Flux PCE','Conc. PCE')
    % title(['G and flux - Time = ',num2str(elapsed/86400),' days'])
    % xlabel('G (\mug/m^3), flux (moles/cm^2/sec)'); ylabel('Depth above water table (cm)')
    % drawnow
    % hold off
    % figure(3333)

    %% 

% % Create 2x2 subplot layout
darkOrange = [0.8500 0.3250 0.0980]; % Dark orange color
darkBlue = [0 0.4470 0.7410];        % Dark blue color
figure (232)
subplot(2,2,2) % TCE (column 1)
plot(flux(:,1),dx*(0:(nin-1)),'-x','Color',darkOrange); hold on  % Changed L-dx to just dx
plot(G_ug_m3(:,1),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Flux TCE','Conc. TCE')
title(['TCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3), flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')  % Add this line
hold off

subplot(2,2,1) % PCE (column 2)
plot(flux(:,2),dx*(0:(nin-1)),'-x','Color',darkOrange); hold on
plot(G_ug_m3(:,2),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Flux PCE','Conc. PCE')
title(['PCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3), flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')
hold off

subplot(2,2,3) % CIS (column 3)
plot(flux(:,3),dx*(0:(nin-1)),'-x','Color',darkOrange); hold on
plot(G_ug_m3(:,3),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Flux CIS','Conc. CIS')
title(['CIS - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3), flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')
hold off

subplot(2,2,4) % VC (column 5)
plot(flux(:,5),dx*(0:(nin-1)),'-x','Color',darkOrange); hold on
plot(G_ug_m3(:,5),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Flux VC','Conc. VC')
title(['VC - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3), flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')
hold off
% Adjust spacing between subplots
sgtitle({['Domain Flux & Vapor Conc.'], ...
         ['Time = ',num2str(elapsed/86400),' days']}, ...
         'FontWeight', 'bold')




% 
% %% 
% 
% % Define colors
% darkOrange = [0.8500 0.3250 0.0980]; % Dark orange color
% darkBlue = [0 0.4470 0.7410];        % Dark blue color
% 
% % First figure for Domain Flux
% figure(232)
% subplot(2,2,2) % TCE
% plot(flux(:,1),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
% legend('Flux TCE')
% title('TCE')
% xlabel('Flux (moles/cm^2/sec)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,1) % PCE
% plot(flux(:,2),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
% legend('Flux PCE')
% title('PCE')
% xlabel('Flux (moles/cm^2/sec)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,3) % CIS
% plot(flux(:,3),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
% legend('Flux CIS')
% title('Cis-1,2-DCE')
% xlabel('Flux (moles/cm^2/sec)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,4) % VC
% plot(flux(:,5),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
% legend('Flux VC')
% title('Vinyl Chloride')
% xlabel('Flux (moles/cm^2/sec)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% sgtitle({['Domain Flux'], ...
%     ['Time = ',num2str(elapsed/86400),' days,'],...
%     ['q = 0.001 cm/sec'] }, ...
%     'FontWeight', 'bold')
% % Second figure for Concentration Gradients
% figure(233)
% subplot(2,2,2) % TCE
% plot(G_ug_m3(:,1),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
% legend('Conc. TCE')
% title('TCE')
% xlabel('Concentration (\mug/m^3)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,1) % PCE
% plot(G_ug_m3(:,2),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
% legend('Conc. PCE')
% title('PCE')
% xlabel('Concentration (\mug/m^3)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,3) % CIS
% plot(G_ug_m3(:,3),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
% legend('Conc. CIS')
% title('Cis-1,2-DCE')
% xlabel('Concentration (\mug/m^3)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% subplot(2,2,4) % VC
% plot(G_ug_m3(:,5),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
% legend('Conc. VC')
% title('Vinyl Chloride')
% xlabel('Concentration (\mug/m^3)')
% ylabel('Depth bgs (cm)')
% set(gca,'YDir','reverse')
% 
% sgtitle({['Vapor Concentration Gradient'], ...
%     ['Time = ',num2str(elapsed/86400),' days,'], ...
%     ['q = 0.001 cm/sec']}, ...
%     'FontWeight', 'bold')
% 
% %% 
% figure(34)
% plot((1/86400)*massremoved(1:n,1),massremoved(1:n,2)*MW(1),'-o')
% hold on
% plot((1/86400)*massremoved(1:n,1),massremoved(1:n,3)*MW(2),'-o')
%     legend('Mass Removed')
%     xlabel('Days'); ylabel('Mass Removed (grams)')
%     legend('PCE Removed', 'TCE Removed')
% hold off
% 
% % Molecular weights in g/mol
% MW_PCE = 165.83;
% MW_TCE = 131.39;
% 
% % First plot the model data (Fig 34)
% % Converting moles to g: moles * g/mol * 1e6 µg/g / 1e6 µg/g = grams
% figure('Position', [100, 100, 800, 500], 'Name', 'Mass Removal Comparison');
% plot((1/86400)*massremoved(1:n,1), massremoved(1:n,2)*MW(1), '-o', 'Color', [0, 0, 0.6], 'LineWidth', 1.5, 'DisplayName', 'PCE Model');
% hold on
% plot((1/86400)*massremoved(1:n,1), massremoved(1:n,3)*MW(2), '-s', 'Color', [0.6, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'TCE Model');
% 
% % Now add the site data
% data = readtable('Extracted_PCETCE.csv');
% dates = datetime(data.Date, 'InputFormat', 'yyyy-MM');
% days_elapsed = days(dates - dates(1)) + 150; % Add 150 days offset
% 
% % Calculate cumulative mass
% days_between = days(diff([dates; dates(end) + calmonths(6)]));
% PCE_cumulative = cumsum(data.PCE_Discharge .* 86400 .* days_between)/1e6; % µg to g
% TCE_cumulative = cumsum(data.TCE_Discharge .* 86400 .* days_between)/1e6; % µg to g
% 
% % Plot site data
% plot(days_elapsed, PCE_cumulative, '--o', 'Color', [0, 0, 0.8], 'LineWidth', 1.5, 'DisplayName', 'PCE Site');
% plot(days_elapsed, TCE_cumulative, '--s', 'Color', [0.8, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'TCE Site');
% 
% % Formatting
% xlabel('Days');
% ylabel('Mass Removed (g)');
% title('Comparison of Modeled and Measured Mass Removal');
% legend('Location', 'best');
% grid on;
% 
% % Print final values
% fprintf('Final PCE mass removed (Site): %.2f g\n', PCE_cumulative(end));
% fprintf('Final TCE mass removed (Site): %.2f g\n', TCE_cumulative(end));
% fprintf('Final PCE mass removed (Model): %.2f g\n', massremoved(n,2)*MW_PCE*1e6/1e6);
% fprintf('Final TCE mass removed (Model): %.2f g\n', massremoved(n,3)*MW_TCE*1e6/1e6);
% 
% % Convert model data from moles/sec to ug/sec
% % moles/sec * g/mol * 1e6 ug/g = ug/sec
% PCE_model_rate = diff(massremoved(1:n,2))./diff(massremoved(1:n,1)) * MW_PCE * 1e6;
% TCE_model_rate = diff(massremoved(1:n,3))./diff(massremoved(1:n,1)) * MW_TCE * 1e6;
% time_model = (1/86400)*massremoved(2:n,1); % Skip first point due to diff
% 
% % Site data is already in ug/sec
% data = readtable('Extracted_PCETCE.csv');
% dates = datetime(data.Date, 'InputFormat', 'yyyy-MM');
% days_elapsed = days(dates - dates(1)) + 150; % Add 150 days offset
% 
% % Create figure for rates
% figure('Position', [100, 100, 800, 500], 'Name', 'Mass Removal Rates');
% semilogy(time_model, PCE_model_rate, '-o', 'Color', [0, 0, 0.6], 'LineWidth', 1.5, 'DisplayName', 'PCE Model');
% hold on
% semilogy(time_model, TCE_model_rate, '-s', 'Color', [0.6, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'TCE Model');
% 
% % Plot site data
% semilogy(days_elapsed, data.PCE_Discharge, '--o', 'Color', [0, 0, 0.8], 'LineWidth', 1.5, 'DisplayName', 'PCE Site');
% semilogy(days_elapsed, data.TCE_Discharge, '--s', 'Color', [0.8, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'TCE Site');
% 
% % Formatting
% xlabel('Days');
% ylabel('Mass Removal Rate (µg/sec)');
% title('Comparison of Modeled and Measured Mass Removal Rates');
% legend('Location', 'best');
% grid on;
% 
% % Print some statistics
% fprintf('Model PCE rate range: %.2e to %.2e µg/sec\n', min(PCE_model_rate), max(PCE_model_rate));
% fprintf('Model TCE rate range: %.2e to %.2e µg/sec\n', min(TCE_model_rate), max(TCE_model_rate));
% fprintf('Site PCE rate range: %.2e to %.2e µg/sec\n', min(data.PCE_Discharge), max(data.PCE_Discharge));
% fprintf('Site TCE rate range: %.2e to %.2e µg/sec\n', min(data.TCE_Discharge), max(data.TCE_Discharge));
% % 
% % 
% %     figure(111)
% %     Cplotnew=[G_new(:,1:(nc-1)); G(nin+1,1:(nc-1))];
% %     G_new_Solvents = G_new(:,[2 1 3 5]);
% %     Cplotnew_Solvents=[G_new_Solvents; G(nin+1,[2 1 3 5])]; 
% %     semilogx(1e12 * Cplotnew_Solvents(:,1) .* MW(2,:)', (dx/100) * (0.5 + (0:nin)), '-o', 'MarkerFaceColor', darker_colors(2, :)+.1,'MarkerEdgeColor','black','MarkerSize', 8,'LineWidth',1.5);%, ...
% %     hold on
% %     semilogx(1e12 * Cplotnew_Solvents(:,2) .* MW(1,:)', (dx/100) * (0.5 + (0:nin)), '-o',  'MarkerFaceColor', darker_colors(1, :)+.1,'MarkerEdgeColor','black','MarkerSize', 8,'LineWidth',1.5);%,...
% %     semilogx(1e12 * Cplotnew_Solvents(:,3) .* MW(3,:)', (dx/100) * (0.5 + (0:nin)), '-o',  'MarkerFaceColor', darker_colors(3, :)+.1,'MarkerEdgeColor','black','MarkerSize', 8,'LineWidth',1.5);%, ...
% %     semilogx(1e12 * Cplotnew_Solvents(:,4) .* MW(5,:)', (dx/100) * (0.5 + (0:nin)), '-o',  'MarkerFaceColor', darker_colors(5, :)+.1,'MarkerEdgeColor','black','MarkerSize', 8,'LineWidth',1.5);%, ...
% % 
% %     set(gca,'Ydir','reverse')
% %     xlabel('Soil Gas Conc. [ug/m^3]'); ylabel('Depth (m)')
% %     axis([1e-1 1e10 0 (L+dx/2)/100])
% %     set(gca, 'ColorOrder', colors_solvents([2 ,1, 3, 4],:)); 
% %     legend('PCE','TCE','cis-1-2,DCE','Vinyl Chloride');
% %     % legend('TCE','PCE','cis-1-2,DCE','Toluene','Vinyl Chloride','Xylenes','Ethylbenzene');
% %     title(['Vapor Con. - Time = ',num2str(elapsed/86400),' days'])
% %     line([6000 6000], [-10 10], 'Color', colors_solvents(2,:), 'LineStyle', '--','LineWidth', 2, 'DisplayName', 'PCE VISL')
% %     line([298 298], [-10 10], 'Color', colors_solvents(1,:), 'LineStyle', '--', 'LineWidth', 2,'DisplayName', 'TCE VISL')
% %     hold off
% %     grid on;
% %% 
%     % figure(333)
%     % % plot(G(:,1:(nc-1)),L-dx*(0:(nin)),'-o'); hold on
%     % % plot (flux,L-dx*(0:(nin-1)),'-x')
%     % plot(G_new(:,1),L-dx*(0:(nin-1)),'-o'); hold on
%     % plot (flux(:,1),L-dx*(0:(nin-1)),'-x')
%     % legend('Conc.','flux')
%     % % set(gca, 'XScale', 'log');
%     % % axis([0 5*max(G(nin+1,1:ncomp)) 0 (nin)*dx])
%     % legend('Conc. TCE', 'Flux TCE')
%     % title(['G and flux - Time = ',num2str(elapsed/86400),' days'])
%     % xlabel('G, flux'); ylabel('Depth above water table (cm)')
%     % drawnow
%     % hold off
