clear all;
% tic
Figure_Count = 145; 
current_run = 'CombinedOutside';
Sim_ID = 1; 
% Boundary Conditions and Sampling Dates
L           =   609;        % Length units in cm (20 ft at site = 20*12*2.54 609 cm) 
nin         =   15;         % Specify number of internal spatial nodes
dx          =   L/nin;      % Specify dx based on L and nin. 
n           =   0;          % Start at t = 0. Loop begins with "n = n + 1". 
Building_Model = 0;         % Building Model = 1, Outside Model = 0
PEST_Optimized = 0;         % Set 1 if optimized PEST parameters available. 
[Contam_Soil ,VIMS_Present, I_Present, Bare_ground,VIMS_Sample,Spin_Up,Bio_Deg] = Set_Model(Building_Model); % Load Model Type Conditions

% Step 1: Read initial parameters from 'Parameters.csv'
param_data = readtable('Parameters.csv', 'ReadVariableNames', true);
for i = 1:height(param_data)
    assignin('base', param_data.Parameter{i}, param_data.Value(i));
end

% Update parameters with PEST values if file exists
if isfile('Model_Input_All.txt')
    values = dlmread('Model_Input_All.txt');
    PEST_Params = {'por', 'por_v_min','FOC', 'airdiff', 'K_sat', 'BC_lambda', 'P_air', ...
                   'Conc_Soil_1', 'Conc_Soil_2', 'ctot_16'};
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

% Initalize arrarys related to phase concentrations using inital parameters.
[MC_bare,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G_bare,W_bare,...
 G_Top, C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
 CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,elapsed]...
 = Initialize_Model(MASSF_soil_PEST,por,theta_bare,dx,nin,DENSLIQ,SDENS,...
 TEMP,FOC,Spin_Up,Ctot);

[MC_slab,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G_slab,W_slab,...
 G_Top, C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
 CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,elapsed]...
 = Initialize_Model(MASSF_soil_PEST,por,theta_slab,dx,nin,DENSLIQ,SDENS,...
 TEMP,FOC,Spin_Up,Ctot);


% Setup Infiltration Array if Infiltration present
if I_Present ==1                                                % If I present then load in infiltration. 
I = load('Site_precip_cm.csv'); I = flipud(I);                  % Load infiltration array (Col 1 = days, Col 2 = cm/day). 
I(:,1)=86400*I(:,1);I(:,2)=(1/86400)*I(:,2); I = I(:,1:2);      % Convert days to seconds and cm/day to cm/sec
end 

% Function to setup stress periods
[BC, dt, dt_max, ttot, VIMS_Exhaust, n_Soil_Contam,ntsteps] = Stress_Periods(Conc_Soil_1,Conc_Soil_2,Contam_Soil,...
    VIMS_Present, VIMS_Sample,Bare_ground,nc);                        % Function to setup time stress periods and contamination based on given criteria. 
% tic

% % % % % Main loop  % % % % % 
while elapsed < ttot 
  n=n+1;
    % Figure out water contents, with infiltration
    if I_Present == 1 % If I_present = 1, induce infiltration 
    [theta_bare,lwaves,twaves] = Kinwave_Theta(theta_max,theta_min,K_sat,BC_lambda,P_air,theta_bare, ...
                                          lwaves,twaves,elapsed,dt,I,L);
    end
    theta_bare(nin+1,2)=theta_max; % Force saturation at water table
    
    % Calculate Biodegradation based on first order rate constants for
    % chlorianted solvents. 
    [W_bare,W_Deg]=Biodeg_Solvents(W_bare,lambda_years,elapsed,nin); 
    [W_slab,W_Deg]=Biodeg_Solvents(W_slab,lambda_years,elapsed,nin); 

    % Calculate the multi-phase equilibrium at each node.
    % Molar concentrations (mol/cm^3) for: G is gas conc (mol/cc), W is dissolved (mol/cc),
    % MC is total (mol/cc).  Then calculate void space and Deff at all nodes.
    [G_slab,W_slab,MC_slab,D_nodes_slab,void_slab,PHASE] = ...
        Equilibrium_Function(G_slab,W_slab,MC_slab,por,MW,VP,ACT,KD,SOL,airdiff,DENSLIQ,SDENS,TEMP);
    [G_bare,W_bare,MC_bare,D_nodes_bare,void_bare,PHASE] = ...
        Equilibrium_Function(G_bare,W_bare,MC_bare,por,MW,VP,ACT,KD,SOL,airdiff,DENSLIQ,SDENS,TEMP);
    Gold_slab=G_slab;                                % Update gas diffusion with previous timestep 
    Gold_bare=G_bare;                                % Update gas diffusion with previous timestep 
    
    % Diffuse vapor using finite differences. See Vapor Diffuse Function. 
    [Gnew_slab]=Vapor_Diffuse(G_slab,D_nodes_slab,dt,dx); 
    [Gnew_bare]=Vapor_Diffuse(G_bare,D_nodes_bare,dt,dx); 

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

    % Mass Removal and udpate total mass concentrations for VIMS on/OFF
    [massremoved, MC_slab,MC_bare,Gnew_slab,Gnew_bare,x1,x2,x3] = Boundary_Conditions(BC,n,dx,MC_slab,MC_bare,...
           elapsed,massremoved,Contam_Soil,n_Soil_Contam,...
           nodes_contam,SDENS,ncomp,MASSF_soil,por,theta_slab,theta_bare,Gnew_slab,Gnew_bare,MW,x1,x2,x3);
%pause

%  Old version    
%    [massremoved, MC_slab,MC_bare,x1,x2,x3] = Boundary_Conditions(BC,n,dx,MC,...
%    elapsed,massremoved,Contam_Soil,n_Soil_Contam,...
%    nodes_contam,SDENS,ncomp,MASSF_soil,por,theta,G_new,void,MW,...
%    x1,x2,x3); % Function to update molar concentration and labels depending on VIMS operation status and soil contam. 
    
    % % % % % % Save Phase concentrations for plotting and update time % % % % % %
    % Right now this is just from the slab portion ...

    [G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_ug_L, W_GW, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed] = ...
    Save_Model_Data(G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_ug_L,W_GW, W_Deg_plot, W_Deg, elapsed, ncomp,Gnew_slab,n,SDENS,MC_slab,MW,PHASE,W_slab,dt,dt_max,nin);
  
end 
% [Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta_bare, Conc_Soil_1,...
%     Conc_Soil_2, ctot_16, K_sat, BC_lambda, P_air, por, airdiff,FOC, n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
%     PHASE_History, W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
%     MW, VIMS_Exhaust, Gnew_slab, nin, BC, Building_Model, Sim_No, dx, SDENS, MC_slab); 
% [Model_Data_Export] = Export__Model_Data(n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
%     PHASE_History, W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
%     MW, VIMS_Exhaust, Gnew_slab, nin, BC, Building_Model, Sim_No, dx, SDENS, MC_slab);
[Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta_slab,theta_bare, Conc_Soil_1, Conc_Soil_2, ctot_16, K_sat, BC_lambda,...
    P_air, por, airdiff,FOC, n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
    PHASE_History, W_ug_L, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
    MW, VIMS_Exhaust, nin, BC, Building_Model, Sim_ID, dx, SDENS, MC_slab); 

% toc
% model_time = toc
