function [BC, dt, dt_max, ttot, VIMS_Exhaust, n_Soil_Contam,ntsteps,W_Deg] = Stress_Periods(Conc_Soil_1,Conc_Soil_2,Contam_Soil,...
    VIMS_Present, VIMS_Sample,Bare_ground,nc)
    % Setup Stress Period Conditions 
    Max_Days = 3300; % Total Model Runtime (Days) 3300
    VIMS_Days = 3300; % Total VIMS Runtime (Days) 3300
    BC_1_Start = 1; BC_1_End = BC_1_Start+1; BC_1_Soil1 = 0; BC_1_Soil4 = 0; BC_1_GW = 0;
    BC_2_Soil1 = Contam_Soil * Conc_Soil_1; % Avg Soil Conc. 20.44; % Max Soil Conc. (mg/kg)
    BC_2_Soil4 = Contam_Soil * Conc_Soil_2; % Avg Soil Conc. 8.63; % Max Soil Conc. (mg/kg)
    BC_2_GW = 0; BC_3_Soil1 = 0; BC_3_Soil4 = 0; BC_3_GW = 0;
    
    if VIMS_Present == 0 
        BC_2_Start = Max_Days; BC_2_End = BC_2_Start; % Start and End of VIMS
        BC_3_Start = BC_2_End ; BC_3_End = BC_3_Start ;
        ttot = BC_3_End * 86400; % Runtime for no VIMS. 
        VIMS_Exhaust = VIMS_Sample; % Used for plotting 
    else
        BC_2_Start = VIMS_Sample; % Start of VIMS. (01/01/2016 - 806)
        BC_2_End = VIMS_Sample + VIMS_Days; % End of VIMS/Just after last Sampling Event (1/3/2024 - 2638)
        BC_3_Start = BC_2_End + 1; BC_3_End = Max_Days; 
        ttot = BC_3_End * 86400; % Total runtime
        VIMS_Exhaust = VIMS_Sample; % Used for plotting 
    end 
    if Bare_ground ==1 
        BC_2_Start = BC_1_End+1; % Start of VIMS. (01/01/2016 - 806)
        BC_2_End = Max_Days; % End of VIMS/Just after last Sampling Event (1/3/2024 - 3711)
        BC_3_Start = BC_2_End + 1; BC_3_End = BC_3_Start+1; 
        ttot = BC_3_End * 86400; % Total runtime
        VIMS_Exhaust = VIMS_Sample; % Used for plotting 
    end
    % Boundary Condition Setup: [Stress Period #, Start Date, End Date, Node 1 Mass Contam, Node 2 Mass Contam, GW Mass Contam]
    BC = [1 BC_1_Start BC_1_End BC_1_Soil1 BC_1_Soil4 BC_1_GW; % Spin Up (VIMS OFF) 
          2 BC_2_Start BC_2_End BC_2_Soil1 BC_2_Soil4 BC_2_GW; % VIMS ON 
          3 BC_3_Start BC_3_End BC_3_Soil1 BC_3_Soil4 BC_3_GW]; % VIMS OFF 

    dt = 0.5 * 86400;  dt_max = 5 * 86400;   % Timestep min and max seconds

    if Contam_Soil == 0 % Setup introduction of contaminated soil
        n_Soil_Contam = BC(3,3) + 200;  % Set date for soil contamination. 
    else 
        % n_Soil_Contam = BC(1,2) + 1;  % Start soil contamination early
        n_Soil_Contam = BC(2,2);  % Start soil contamination just before VIMS
    end
    ttot = BC_3_End * 86400; % Total model run
    ntsteps=ceil(ttot/dt); % Number of timesteps required (not including 0.4 dt) 
    W_Deg = zeros(nc); % Setp plot data foor degraded solvents 
end
