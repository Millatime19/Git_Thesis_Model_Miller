function [Model_Data_Export] = Export__Model_Data(current_run,Figure_Count,theta, Conc_Soil_1,...
    Conc_Soil_2, ctot_16, K_sat, BC_lambda, P_air, por, airdiff,FOC, n, ncomp, L, G_Top,G_Top_ug_m3,G_WT,G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All,...
    PHASE_History, W_ug_L, W_GW_ug_L, W_Deg_plot, CtotCurrent_All, CtotCurrent_Sum, dt, elapsed,...
    MW, VIMS_Exhaust, G_new, G_ug_m3, nin, BC, Building_Model, Sim_ID, dx, SDENS, MC)
   
    % Populate Model_Export structure with all required variables
    Model_Export.current_run = current_run; 
    Model_Export.Figure_Count = Figure_Count;
    Model_Export.G_Top = G_Top;
    Model_Export.G_WT = G_WT;
    Model_Export.C_Top = C_Top;
    Model_Export.C_Top_All = C_Top_All;
    Model_Export.C_Bottom_All = C_Bottom_All;
    Model_Export.PHASE_History = PHASE_History;
    Model_Export.W_ug_L = W_ug_L;
    Model_Export.W_GW_ug_L = W_GW_ug_L;
    Model_Export.W_Deg_plot = W_Deg_plot;
    Model_Export.CtotCurrent_All = CtotCurrent_All;
    Model_Export.CtotCurrent_Sum = CtotCurrent_Sum;
    Model_Export.dt = dt;
    Model_Export.elapsed = elapsed;
    Model_Export.MW = MW;
    Model_Export.VIMS_Exhaust = VIMS_Exhaust;
    Model_Export.G_Top_ug_m3 = G_Top_ug_m3;
    Model_Export.G_WT_ug_m3 = G_WT_ug_m3;
    Model_Export.Gnew = G_new;
    Model_Export.G_ug_m3 = G_ug_m3;
    Model_Export.n = n;
    Model_Export.nin = nin;
    Model_Export.ncomp = ncomp;
    Model_Export.BC = BC;
    Model_Export.Building_Model = Building_Model;
    Model_Export.Sim_ID = Sim_ID;
    Model_Export.dx = dx;
    Model_Export.SDENS = SDENS;
    Model_Export.MC = MC;
    Model_Export.airdiff = airdiff;
    Model_Export.K_sat = K_sat;
    Model_Export.BC_lambda = BC_lambda;
    Model_Export.P_air = P_air;
    Model_Export.BC = BC;
    Model_Export.por = por;
    Model_Export.theta = theta;
    Model_Export.Conc_Soil_1 = Conc_Soil_1;
    Model_Export.Conc_Soil_2= Conc_Soil_2;
    Model_Export.ctot_16 = ctot_16; 
    % Model_Export.I = I; 
    Model_Export.L = L; 
    % Load the SS8 data from CSV file and adjust dates
    VIMS_Data = load("SS8_Conc_Data.csv");
   
    dates = VIMS_Data(:, 1);
    Model_Start = min(VIMS_Data(:, 1));
    dates = dates - Model_Start + VIMS_Exhaust + 90; % Adjust dates based on VIMS_Exhaust
    contaminantIDs = VIMS_Data(:, 2);
    concentrations = VIMS_Data(:, 3);

    % Find unique dates and contaminants
    unique_dates = unique(dates);
    unique_contaminants = unique(contaminantIDs);

    % Initialize a matrix to store the reshaped data
    VIMS_Data = NaN(length(unique_dates), length(unique_contaminants) + 1);
    VIMS_Data(:, 1) = unique_dates;

    % Loop through each unique contaminant and fill the corresponding column
    for i = 1:length(unique_contaminants)
        contaminant = unique_contaminants(i);
        for j = 1:length(unique_dates)
            date = unique_dates(j);
            idx = (dates == date) & (contaminantIDs == contaminant);
            if any(idx)
                VIMS_Data(j, i + 1) = concentrations(idx);
            end
        end
    end

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
    
    Model_Export.VIMS_Data = VIMS_Data;
    Vapor_All = [SubSlabVapor_Data_transformed;VIMS_Data];

 % % % % %  Create Export File of modeled observations for PEST % % % % %  
    model_time = G_Top(:, 1)*(1/86400);                                  % Extract model time (dates) from the first column of G_Top_ug_m3
    model_data= G_Top_ug_m3;                                             % Extract modeled concentrations from columns 2 to 9 of G_Top_ug_m3
    measured_times = Vapor_All(:,1);                                     % Extract the measured observation times from new_SS8_Data % Ensure this is a column vector with dates/times
    % % % % %  Create Export File of modeled observations for PEST % % % % %  
    model_time = G_Top(:, 1)*(1/86400);                                  % Extract model time (dates) from the first column of G_Top_ug_m3
    model_data= G_Top_ug_m3;                                             % Extract modeled concentrations from columns 2 to 9 of G_Top_ug_m3
    measured_times = Vapor_All(:,1);                                     % Extract the measured observation times from new_SS8_Data % Ensure this is a column vector with dates/times
    Model_PEST_Out = interp1(model_time, model_data, measured_times);    % Interpolate the modeled concentrations at the measured times
    Model_PEST_Out = Model_PEST_Out(:,2:end);                            % Grab only the analytical data (no time) 
    Model_Obs_New = reshape(Model_PEST_Out', [], 1);                     % Reshape the 17x7 array into a single column
    filename = 'Modeled_Obs_New.txt';                                    % Save the reshaped array to a text file 
    writematrix(Model_Obs_New, filename, 'Delimiter', 'tab');            % Write the data to the text file

    % % % % Save other data - need to clean up 
    Model_Data_Export = [measured_times, Model_PEST_Out]; % Create Data_Export array with measured times and interpolated data
    Model_Export.Model_Data_Export = Model_Data_Export; 
    Model_Export.FOC = FOC; PEST_Export.FOC = FOC;
    Model_Export.model_time = model_time;  PEST_Export.model_time = model_time; 
    Model_Export.model_data = model_data; PEST_Export.model_data = model_data; 
    Model_Export.measured_times = measured_times; PEST_Export.measured_times = measured_times;  
    Model_Export.Model_Obs_New = Model_Obs_New; 

    % Save Model_Export structure to a .mat file
    save('Model_Export.mat', '-struct', 'Model_Export');
    % save('Miller_Model_Export.txt', '-struct', 'Model_Export', '-ascii');
    plot.nc_colors = ncomp;

    % Define colors and other stuff for plotting. 
    plot.color_TCE = [0 0.4470 0.7410];          % Blue
    plot.color_PCE = [0.8500 0.3250 0.0980];     % Red
    plot.color_Cis = [0.9290 0.6940 0.1250];     % Orange
    plot.color_Toluene = [0.4940 0.1840 0.5560]; % Purple
    plot.color_VC = [0.4660 0.6740 0.1880];      % Green
    plot.color_Xylenes = [0.5294 0.8078 0.9216]; % Sky Blue
    plot.color_Ethylbenzene = [0.5 0 0];         % Maroon

    % Y-axis limits
    plot.y_limits_1 = [0, L * 2];

    % Remedial Goal for vapor concentration: ug/m^3
    plot.x_tce_vapor = 293;
    plot.x_pce_vapor = 6000;
    plot.x_vc_vapor = 4.00 ;
    plot.x_toluene_vapor = 2.20E-01;
    plot.x_xylenes_vapor = 4.40E-03;
    plot.x_ethylbenzene_vapor = 4.90E-05;

    % Dissolved concentrations
    plot.x_tce_diss = 7.40E-06;
    plot.x_pce_diss = 6.50E-05;
    plot.x_cis_diss = 7.00E-02;
    plot.x_vc_diss = 2.00E-03;

    % Soil concentrations
    plot.x_tce_soil = 2.82E-04; % Max Soil Sample
    plot.x_pce_soil = 1.49E-02; % Max Soil Sample   
    plot.x_vc_soil = 0.023; % EPA Residential 
    plot.x_toluene_soil = 8.80E-6; % Max Soil Sample
    plot.x_xylenes_soil = 100; % EPA Residential 
    plot.x_ethylbenzene_soil = 7.8; % EPA Residential 

    % Other values
    plot.x1 = 1;
    plot.x2 = 1;
    plot.x3 = 1;
    plot.x_limits = [0 800];
    plot.y_limits = [0, L * 2];
    
    % Save variables to a .mat file
    save('plot.mat', '-struct', 'plot');

end