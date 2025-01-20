clear all;

% Set Model Colors and IDs
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

% Step 1: Read initial parameters from 'Parameters.csv'
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
q_air_Slab = 10^(q_air_Slab);  % PEST uses log(q_air_Slab), take out of log
q_air_VIMS = 10^(q_air_VIMS);  % PEST uses log(q_air_VIMS), take out of log 


if PEST_Optimized ==1 
values = dlmread('Optimization_B0.txt');                   % Read values from Model_Input.txt
PEST_Params_Opt = {'por','FOC',  'airdiff', 'K_sat', 'BC_lambda', 'P_air','Conc_Soil_1', 'Conc_Soil_2','ctot_16','Finalphi'};  % Values for optimized parameters for Outside PEST - All Params
for i = 1:length(values)
     assignin('base', PEST_Params_Opt{i}, values(i)); 
end
end 
%  Assign Ctot to the base workspace based on parameter inputs for each node. 
Ctot = [ctot_1 ctot_2 ctot_3 ctot_4 ctot_5 ctot_6 ctot_7 ctot_8 ctot_9 ctot_10 ctot_11 ctot_12 ctot_13 ctot_14 ctot_15 ctot_16]';
%  Assign Ctot to the base workspace based on parameter inputs for each node. 
% ctot_16=0;
Ctot = [ctot_1 ctot_2 ctot_3 ctot_4 ctot_5 ctot_6 ctot_7 ctot_8 ctot_9 ctot_10 ctot_11 ctot_12 ctot_13 ctot_14 ctot_15 ctot_16]';
% Conc_Soil_1 = 10000;
% Conc_Soil_2 = 10000; 
% % % Make Ctot only GW Source 
 % Ctot(1:nin,1) = 0;  
% Load Input Parameters
filename = 'Screening_Info_2.csv'; param_data = readtable(filename, 'ReadVariableNames', true); % Load the CSV file
for i = 1:height(param_data)
    param_name = param_data.Parameter{i}; param_value = param_data.Value(i); assignin('base', param_name, param_value); 
end

% Initialize parameter variables 
por = por*ones(nin+1,1);                                        % Specify total porosity at nodes
theta_max = por(1) - por_v_min; theta_slab=zeros(nin+1,2);      % Inital moisture content max, current. 
theta_slab(:,2) =theta_initial*ones(nin+1,1);                   % Inital water content throughout domain 
theta_slab(:,1) = dx*(0.5:1:(nin+0.5))';                        % Vertical distance for each node
theta_slab(nin+1,2)=theta_max;                                  % Max soil moisture at water table
D_nodes_slab=zeros(nin+1,1); D_half_slab=zeros(nin,1);          % Setup Effective Diffusion array for finite difference approx. 
lwaves=zeros(1,3); twaves=zeros(1,3);                           % Initialize for the kinematic wave infiltration function

% Put in hydrostatic soil moisture according to Brooks-Corey:
apple=find(theta_slab(:,1)<(L-P_air));                               % Nodes above air-entry pressure
theta_slab(apple,2)=theta_min+(theta_max-theta_min)*(P_air./(L-theta_slab(apple,1))).^BC_lambda;
theta_slab(nin+1,2)=theta_max;                                       % Water table

% Copy stuff to bare-ground model
theta_bare=theta_slab;  D_nodes_bare=D_nodes_slab;  D_half_bare=D_half_slab;

% Function to initialize chemical data and arrays
MASSF_soil_PEST = [Soil_Frac_1,Soil_Frac_2,Soil_Frac_3,Soil_Frac_4,Soil_Frac_5,Soil_Frac_6,Soil_Frac_7]';
MASSF_GW_PEST = [GW_Frac_1,GW_Frac_2,GW_Frac_3,GW_Frac_4,GW_Frac_5,GW_Frac_6,GW_Frac_7]';

% Initalize arrarys related to phase concentrations using inital parameters.
[MC_bare,MC_Before,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G_bare,W_bare,...
 G_Top,G_WT,W_WT,C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
         CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,massremoved_dt,massremoved_g,massremoved_g_tot, elapsed]...
 = Initialize_Model(MASSF_soil_PEST, MASSF_GW_PEST, por,theta_bare,dx,nin,DENSLIQ,SDENS,...
         TEMP,FOC,Ctot);

[MC_slab,MC_Before,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G_slab,W_slab,...
 G_Top,G_WT,W_WT,C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
         CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,massremoved_dt,massremoved_g,massremoved_g_tot, elapsed]...
 = Initialize_Model(MASSF_soil_PEST,MASSF_GW_PEST,por,theta_slab,dx,nin,DENSLIQ,SDENS,...
 TEMP,FOC,Ctot);


% Setup Infiltration Array if Infiltration present
if I_Present ==1                                                % If I present then load in infiltration. 
I = load('Site_precip_cm.csv'); I = flipud(I);                  % Load infiltration array (Col 1 = days, Col 2 = cm/day). 
I(:,1)=86400*I(:,1);I(:,2)=(1/86400)*I(:,2); I = I(:,1:2);      % Convert days to seconds and cm/day to cm/sec
end 

% Function to setup stress periods
[BC, dt, dt_max, ttot, VIMS_Exhaust, n_Soil_Contam,ntsteps,W_Deg] = Stress_Periods(Conc_Soil_1,Conc_Soil_2,Contam_Soil,...
    VIMS_Present, VIMS_Sample,Bare_ground,nc);                      % Function to setup time stress periods and contamination based on given criteria. 
elapsed = 0;

if Building_Model == 1 
    if n==0 % Slab condition 
       q_air = q_air_Slab; 
    end 
    if VIMS_Present==1
    if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % Slab condition 
       q_air = q_air_Slab; 
    end
    if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
     q_air = q_air_VIMS;
    end
    end

end   

if Building_Model == 0 
     q_air = q_air_Open; 
end
elapsed = 0;
% % % % % Main loop  % % % % % 
while elapsed < ttot 
  n=n+1;
    % Figure out water contents, with infiltration
    if I_Present == 1 % If I_present = 1, induce infiltration 
    [theta_bare,lwaves,twaves] = Kinwave_Theta(theta_max,theta_min,K_sat,BC_lambda,P_air,theta_bare, ...
                                          lwaves,twaves,elapsed,dt,I,L);
    end
    theta_bare(nin+1,2)=theta_max; % Force saturation at water table
    
    if Bio_Deg == 1% Calculate Biodegradation of dissolved phase based on first order rate constants for chlorianted solvents. 
    [W_bare,W_Deg]=Biodeg_Solvents(W_bare,lambda_years,elapsed,nin); 
    [W_slab,W_Deg]=Biodeg_Solvents(W_slab,lambda_years,elapsed,nin); 
    end 
    % Calculate the multi-phase equilibrium at each node.
    % Molar concentrations (mol/cm^3) for: G is gas conc (mol/cc), W is dissolved (mol/cc),
    % MC is total (mol/cc).  Then calculate void space and Deff at all nodes.
    [G_slab,W_slab,MC_slab,D_nodes_slab,void_slab,PHASE,calc_H,theory_H,Retardation,totmass] = ...
        Equilibrium_Function(G_slab,W_slab,MC_slab,por,por_v_min,MW,VP,ACT,KD,SOL,...
    airdiff,DENSLIQ,SDENS,TEMP);
    [G_bare,W_bare,MC_bare,D_nodes_bare,void_bare,PHASE,calc_H,theory_H,Retardation,totmass] = ...
        Equilibrium_Function(G_bare,W_bare,MC_bare,por,por_v_min,MW,VP,ACT,KD,SOL,...
    airdiff,DENSLIQ,SDENS,TEMP);
    Gold_slab=G_slab;                                % Update gas diffusion with previous timestep 
    Gold_bare=G_bare;                                % Update gas diffusion with previous timestep 
    
    % Diffuse vapor using finite differences. See Vapor Diffuse Function. 
    q_send_slab=q_air; % Add q_vapor for sublab and VIMS 
    q_send_bare = q_air_Open; % Add q_vapor for outside 
    [Gnew_slab]=Vapor_Diffuse(G_slab,D_nodes_slab,dt,dx,q_send_slab); 
    [Gnew_bare]=Vapor_Diffuse(G_bare,D_nodes_bare,dt,dx,q_send_bare); 
    G_slab(1:nin,:)=Gnew_slab(:,:); % Update G with G_New
    G_bare(1:nin,:)=Gnew_bare(:,:); % Update G with G_New

    
    flux=(.5/dx)*(D_nodes_slab(1:end-1)+D_nodes_slab(2:end)).*diff(G_slab);
    massremoved(n,1) = elapsed; 
    massremoved(n+1,2:end) = massremoved(n,2:end)+dx*q_send_slab*Gnew_slab(1,1:ncomp); 

        % Update q if VIMS is on or off. 
    if Building_Model == 1                          
      if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % Slab condition 
       q_air = q_air_Slab; 
      end
      if VIMS_Present ==1
            if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
                q_air = q_air_VIMS; 
            end
            if elapsed/86400 > BC(2,3)  % Turn off VIMS  
                q_air = q_air_Slab; 
            end 
        end
    end   
    % Advect the soil water concentrations if Infiltration is present 
    if I_Present ==1
    [W_new_bare] = water_advect(theta_bare,W_bare,K_sat,theta_max,theta_min,BC_lambda,I,elapsed,dt,dx);
    W_Save_bare=(W_new_bare(1:nin,1:ncomp)-W_bare(1:nin,1:ncomp));
    MC_bare(1:nin,1:ncomp)=MC_bare(1:nin,1:ncomp)+W_Save_bare;
    end 

    % Figure out change in gas molar conc and remove from total molar conc.
    G_save_slab=(Gnew_slab-Gold_slab(1:nin,:)).*void_slab(1:nin);
    MC_slab(1:nin,:)=MC_slab(1:nin,:)+G_save_slab;
    MC_slab(:,end)=(1/18)*theta_slab(:,2);  % Keep track of water moles per cm^3 soil.
    G_save_bare=(Gnew_bare-Gold_bare(1:nin,:)).*void_bare(1:nin);
    MC_bare(1:nin,:)=MC_bare(1:nin,:)+G_save_bare;
    MC_bare(:,end)=(1/18)*theta_bare(:,2);  % Keep track of water moles per cm^3 soil.

 % Update Soil Boundary Condition 
    if Contam_Soil == 1 && elapsed > 48*86400 && SoilhasExecuted ==1
    Surf_Contam = zeros(nodes_contam,1);  
    Ctot_Surf(1:length(Surf_Contam)) = linspace(BC(2,4),BC(2,5),nodes_contam); % Avg soil contamination
    elapsed_Soil_Contam = elapsed/86400;
    % Labels based on VIMS operation and soil contamination dates. 
    x3 = elapsed_Soil_Contam; % Used for plotting label. 
        for i = 1:nodes_contam
        MC_slab(i,1:ncomp) = Ctot_Surf(i)*(SDENS/1e6)*MASSF_soil(1:ncomp)./MW(1:ncomp); % Update molar concentrations using total mass conc. from contaminated surface soil. (Nodes 1-4)
        MC_Before = MC_slab; % Save starting molar conentration. 
        
        end
    SoilhasExecuted = 0; % Set the flag to 0 so this block doesn't run again
    end

   
    % % % % % % Save Phase concentrations for plotting and update time % % % % % %
    % Right now this is just from the slab portion ...

    % [G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    % W_ug_L, W_GW, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed] = ...
    % Save_Model_Data(G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    % W_ug_L,W_GW, W_Deg_plot, W_Deg, elapsed, ncomp,Gnew_slab,n,SDENS,MC_slab,MW,PHASE,W_slab,dt,dt_max,nin);
  
       [G_Top, G_Top_ug_m3, G_ug_m3, G_WT, G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_WT,W_ug_L, W_GW_ug_L,  CtotCurrent_All, CtotCurrent_Sum, dt, elapsed] = ...
    Save_Model_Data(G_Top, G_WT, W_WT,C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    elapsed, ncomp,Gnew_slab,n,SDENS,MC_slab,MW,PHASE,W_slab,dt,dt_max,nin);
    
end 
% [Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta_bare, Conc_Soil_1,...
%     Conc_Soil_2, ctot_16, K_sat, BC_lambda, P_air, por, airdiff,FOC, n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
%     PHASE_History, W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
%     MW, VIMS_Exhaust, Gnew_slab, nin, BC, Building_Model, Sim_No, dx, SDENS, MC_slab); 
% [Model_Data_Export] = Export__Model_Data(n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
%     PHASE_History, W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
%     MW, VIMS_Exhaust, Gnew_slab, nin, BC, Building_Model, Sim_No, dx, SDENS, MC_slab);
% [Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta_slab,theta_bare, Conc_Soil_1, Conc_Soil_2, ctot_16, K_sat, BC_lambda,...
%     P_air, por, airdiff,FOC, n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
%     PHASE_History, W_ug_L, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
%     MW, VIMS_Exhaust, nin, BC, Building_Model, Sim_ID, dx, SDENS, MC_slab); 


[Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta_slab,theta_bare, Conc_Soil_1, Conc_Soil_2, ctot_16, K_sat, BC_lambda,...
    P_air, por, airdiff,FOC, n, ncomp, L, G_Top, G_Top_ug_m3,G_WT,G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All,...
    PHASE_History, W_ug_L, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
    MW, VIMS_Exhaust, Gnew_slab, G_ug_m3, nin, BC, Building_Model, Sim_ID, dx, SDENS, MC_slab,MC_Before); 

% toc
% model_time = toc
%%
   % %     % % % % % % %      Kin Wave Plot     % % % % % % %
     figure(321)
     plot(twaves(:,2),L-twaves(:,1),'bd')
     hold on
     plot(lwaves(:,2),L-lwaves(:,1),'ro')
     plot(theta_bare(:,2),L-theta_bare(:,1),'+-k')
     axis([0 theta_max min(L-theta_bare(:,1)) L])
     ylabel('Height above WT (cm)'); xlabel('Vol. Soil Moisture')
     title(['Lead and Trailing Waves; Time = ',num2str(elapsed/86400),' days'])
     legend('Trailing wave','Leading wave')
     drawnow
     hold off
%%
% Create 2x2 subplot layout
darkOrange = [0.8500 0.3250 0.0980]; % Dark orange color
darkBlue = [0 0.4470 0.7410];        % Dark blue color

% Figure 1 - Domain Concentration Flux
figure(2245)
subplot(2,2,2) % TCE
plot(flux(:,1),dx*(0:(nin-1)),'-x','Color',darkOrange);
legend('Flux TCE')
title(['TCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,1) % PCE
plot(flux(:,2),dx*(0:(nin-1)),'-x','Color',darkOrange);
legend('Flux PCE')
title(['PCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,3) % CIS
plot(flux(:,3),dx*(0:(nin-1)),'-x','Color',darkOrange);
legend('Flux CIS')
title(['CIS - Time = ',num2str(elapsed/86400),' days'])
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,4) % VC
plot(flux(:,5),dx*(0:(nin-1)),'-x','Color',darkOrange);
legend('Flux VC')
title(['VC - Time = ',num2str(elapsed/86400),' days'])
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

sgtitle({['Domain Concentration Flux'], ...
 ['Time = ',num2str(elapsed/86400),' days']}, ...
'FontWeight', 'bold')

% Figure 2 - Domain Vapor Concentrations
figure(2246)
subplot(2,2,2) % TCE
plot(G_ug_m3(:,1),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Conc. TCE')
title(['TCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,1) % PCE
plot(G_ug_m3(:,2),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Conc. PCE')
title(['PCE - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,3) % CIS
plot(G_ug_m3(:,3),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Conc. CIS')
title(['CIS - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,4) % VC
plot(G_ug_m3(:,5),dx*(0:(nin-1)),'-o','Color',darkBlue);
legend('Conc. VC')
title(['VC - Time = ',num2str(elapsed/86400),' days'])
xlabel('G (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

sgtitle({['Domain Vapor Conc. Concentrations'], ...
 ['Time = ',num2str(elapsed/86400),' days']}, ...
'FontWeight', 'bold')



%% 

% Define colors
darkOrange = [0.8500 0.3250 0.0980]; % Dark orange color
darkBlue = [0 0.4470 0.7410];        % Dark blue color

% First figure for Domain Flux
figure(232)
subplot(2,2,2) % TCE
plot(flux(:,1),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
legend('Flux TCE')
title('TCE')
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,1) % PCE
plot(flux(:,2),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
legend('Flux PCE')
title('PCE')
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,3) % CIS
plot(flux(:,3),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
legend('Flux CIS')
title('Cis-1,2-DCE')
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,4) % VC
plot(flux(:,5),dx*(0:(nin-1)),'-x','Color',darkOrange,'LineWidth',1.5)
legend('Flux VC')
title('Vinyl Chloride')
xlabel('Flux (moles/cm^2/sec)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

sgtitle({['Domain Flux'], ...
    ['Time = ',num2str(elapsed/86400),' days,'],...
    ['q = 0.001 cm/sec'] }, ...
    'FontWeight', 'bold')
% Second figure for Concentration Gradients
figure(233)
subplot(2,2,2) % TCE
plot(G_ug_m3(:,1),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
legend('Conc. TCE')
title('TCE')
xlabel('Concentration (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,1) % PCE
plot(G_ug_m3(:,2),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
legend('Conc. PCE')
title('PCE')
xlabel('Concentration (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,3) % CIS
plot(G_ug_m3(:,3),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
legend('Conc. CIS')
title('Cis-1,2-DCE')
xlabel('Concentration (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

subplot(2,2,4) % VC
plot(G_ug_m3(:,5),dx*(0:(nin-1)),'-o','Color',darkBlue,'LineWidth',1.5)
legend('Conc. VC')
title('Vinyl Chloride')
xlabel('Concentration (\mug/m^3)')
ylabel('Depth bgs (cm)')
set(gca,'YDir','reverse')

sgtitle({['Vapor Concentration Gradient'], ...
    ['Time = ',num2str(elapsed/86400),' days,'], ...
    ['q = 0.001 cm/sec']}, ...
    'FontWeight', 'bold')
