clear all;
% tic

% Boundary Conditions and Sampling Dates
L           =   609;        % Length units in cm (20 ft at site = 20*12*2.54 609 cm) 
nin         =   15;         % Specify number of internal spatial nodes
dx          =   L/nin;      % Specify dx based on L and nin. 
n           =   0;          % Start at t = 0. Loop begins with "n = n + 1". 
Building_Model = 0;         % Building Model = 1, Outside Model = 0
PEST_Optimized = 0;         % Set 1 if optimized PEST parameters available. 
[Contam_Soil ,VIMS_Present, I_Present, Bare_ground,VIMS_Sample,Spin_Up,Bio_Deg] = Set_Model(Building_Model); % Load Model Type Conditions

% Read in inital parameters and associated values in 'Parameters.csv' then
% update the parameters with PEST values (Model_Input_New.txt'). 
param_data = readtable('Parameters.csv', 'ReadVariableNames', true); % Load parameter names and values
for i = 1:height(param_data), assignin('base', param_data.Parameter{i}, param_data.Value(i)); end
        values = dlmread('Model_Input_2.txt');                       % Read values from Model_Input.txt
% values = dlmread('Model_Input_New.txt'); % Read values from Model_Input.txt
% varNames_All = {'DENSLIQ', 'SDENS', 'TEMP', 'FOC', 'airdiff', 'por', 'theta_initial', ...
            % 'BC_lambda', 'ctot_1', 'ctot_2', 'ctot_3', 'ctot_4', 'ctot_5', 'ctot_6', ...
            % 'ctot_7', 'ctot_8', 'ctot_9', 'ctot_10', 'ctot_11', 'ctot_12', 'ctot_13', ...
            % 'ctot_14', 'ctot_15', 'ctot_16', 'Conc_Soil_1', 'Conc_Soil_2', 'Soil_Frac_1', ...
            % 'Soil_Frac_2', 'Soil_Frac_3', 'Soil_Frac_4', 'Soil_Frac_5', 'Soil_Frac_6', 'Soil_Frac_7'};
% varNames_1 = { 'por','FOC', 'airdiff', 'Conc_Soil_1', 'Conc_Soil_2','ctot_16'}; % Param Run 1. 
varNames_2 = { 'por','FOC', 'airdiff'}; % PEST parameter Set 2. 
varNames_3 = { 'por','FOC', 'airdiff'}; % Pest parameter Set 3. 
% Ctot = param_data(9:24); % Directly assign ctot values 

for i = 1:length(values)
    assignin('base', varNames_2{i}, values(i)); 
end

if PEST_Optimized ==1 
    for i = 1:height(param_data), assignin('base', param_data.Parameter{i}, param_data.Value(i)); 
    end
values = dlmread('Model_Input_PEST_Opt.txt');                   % Read values from Model_Input.txt
varNames_Opt = { 'por','FOC', 'airdiff', 'Conc_Soil_1', 'Conc_Soil_2','ctot_16'}; % Values for optimized parameters for PEST run #1. 
for i = 1:length(values)
    assignin('base', varNames_Opt{i}, values(i)); 
end
end

% Assign Ctot to the base workspace based on parameter inputs for each node. 
Ctot = [ctot_1 ctot_2 ctot_3 ctot_4 ctot_5 ctot_6 ctot_7 ctot_8 ctot_9 ctot_10 ctot_11 ctot_12 ctot_13 ctot_14 ctot_15 ctot_16]';

% Load Input Parameters
filename = 'Screening_Info_2.csv'; param_data = readtable(filename, 'ReadVariableNames', true); % Load the CSV file
for i = 1:height(param_data)
    param_name = param_data.Parameter{i}; param_value = param_data.Value(i); assignin('base', param_name, param_value); 
end

% Initialize parameter variables 
por = por*ones(nin+1,1);                                        % Specify total porosity at nodes
theta_max = por(1) - por_v_min; theta=zeros(nin+1,2);           % Inital moisture content max, current. 
theta(:,2) =theta_initial*ones(nin+1,1);                        % Inital water content throughout domain 
theta(:,1) = dx*(0.5:1:(nin+0.5))';                             % Vertical distance for each node
theta(nin+1,2)=theta_max;                                       % Max soil moisture at water table
D_nodes=zeros(nin+1,1); D_half=zeros(nin,1);                    % Setup Effective Diffusion array for finite difference approx. 
lwaves=zeros(1,3); twaves=zeros(1,3);                           % Initialize for the kinematic wave infiltration function

% Put in hydrostatic soil moisture according to Brooks-Corey:
apple=find(theta(:,1)<(L-P_air));                               % Nodes above air-entry pressure
theta(apple,2)=theta_min+(theta_max-theta_min)*(P_air./(L-theta(apple,1))).^BC_lambda;
theta(nin+1,2)=theta_max;                                       % Water table

% Function to initialize chemical data and arrays
MASSF_soil_PEST = [Soil_Frac_1,Soil_Frac_2,Soil_Frac_3,Soil_Frac_4,Soil_Frac_5,Soil_Frac_6,Soil_Frac_7]';

% Initalize arrarys related to phase concentrations using inital parameters.
[MC,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G,W,...
 G_Top, C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
 CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,elapsed]...
 = Initialize_Model(MASSF_soil_PEST,por,theta,dx,nin,DENSLIQ,SDENS,...
 TEMP,FOC,Spin_Up,Ctot);

% Setup Infiltration Array if Infiltration present
if I_Present ==1                                                % If I present then load in infiltration. 
I = load('Site_precip_cm.csv'); I = flipud(I);                  % Load infiltration array (Col 1 = days, Col 2 = cm/day). 
I(:,1)=86400*I(:,1);I(:,2)=(1/86400)*I(:,2); I = I(:,1:2);      % Convert days to seconds and cm/day to cm/sec
end 

% Function to setup stress periods
[BC, dt, dt_max, ttot, VIMS_Exhaust, n_Soil_Contam,ntsteps] = Stress_Periods(Conc_Soil_1,Conc_Soil_2,Contam_Soil,...
    VIMS_Present, VIMS_Sample,Bare_ground);                      % Function to setup time stress periods and contamination based on given criteria. 
% tic
% % % % % Main loop  % % % % % 
while elapsed < ttot 
  n=n+1;
    % Figure out water contents, with infiltration
    if I_Present == 1 % If I_present = 1, induce infiltration 
    [theta,lwaves,twaves] = Kinwave_Theta(theta_max,theta_min,K_sat,BC_lambda,P_air,theta, ...
                                          lwaves,twaves,elapsed,dt,I,L);
    end
    theta(nin+1,2)=theta_max; % Force saturation at water table
    
    % Calculate Biodegradation based on first order rate constants for
    % chlorianted solvents. 
    [W,W_Deg]=Biodeg_Solvents(W,lambda_years,elapsed,nin); 

    % Calculate the multi-phase equilibrium at each node.
    % Molar concentrations (mol/cm^3) for: G is gas conc (mol/cc), W is dissolved (mol/cc),
    % MC is total (mol/cc).  Then calculate void space and Deff at all nodes.
    [G,W,MC,D_nodes,void,PHASE] = Equilibrium_Function(G,W,MC,por,MW,VP,ACT,KD,SOL,airdiff,DENSLIQ,SDENS,TEMP);
    Gold=G;                                % Update gas diffusion with previous timestep 
    
    % Diffuse vapor using finite differences. See Vapor Diffuse Function. 
    [G_new]=Vapor_Diffuse(G,D_nodes,dt,dx); 

    % Advect the soil water concentrations if Infiltration is present 
    if I_Present ==1
    [W_new] = water_advect(theta,W,K_sat,theta_max,theta_min,BC_lambda,I,elapsed,dt,dx);
    W_Save=(W_new(1:nin,1:ncomp)-W(1:nin,1:ncomp));
    MC(1:nin,1:ncomp)=MC(1:nin,1:ncomp)+W_Save;
    end 
    % Figure out change in gas molar conc and remove from total molar conc.
    G_save=(G_new-Gold(1:nin,:)).*void(1:nin);
    MC(1:nin,:)=MC(1:nin,:)+G_save;
    MC(:,end)=(1/18)*theta(:,2);  % Keep track of water moles per cm^3 soil.
    % Mass Removal and udpate total mass concentrations for VIMS on/OFF
    [massremoved, MC,x1,x2,x3] = Boundary_Conditions(BC,n,dx,MC,...
    elapsed,massremoved,Contam_Soil,n_Soil_Contam,...
    nodes_contam,SDENS,ncomp,MASSF_soil,por,theta,G_new,void,MW,...
    x1,x2,x3); % Function to update molar concentration and labels depending on VIMS operation status and soil contam. 
    
    % % % % % % Save Phase concentrations for plotting and update time % % % % % %
    [G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum,dt, elapsed] = ...
    Save_Model_Data(G_Top, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_ug_L, W_GW, W_Deg_plot, W_Deg, elapsed, ncomp,G_new,n,SDENS,MC,MW,PHASE,W,dt,dt_max,nin);
  
end 

[Model_Data_Export] = Export__Model_Data(n, ncomp, L, G_Top, C_Top, C_Top_All, C_Bottom_All,...
    PHASE_History, W_ug_L, W_GW, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
    MW, VIMS_Exhaust, G_new, nin, BC, Building_Model, Sim_No, dx, SDENS, MC);

% toc
% model_time = toc

