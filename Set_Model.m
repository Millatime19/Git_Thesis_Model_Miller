function [Contam_Soil ,VIMS_Present, I_Present, Bare_ground,VIMS_Sample,Bio_Deg] = Set_Model(Building_Model)
VIMS_Sample = 251;% First VIMS Exhaust Sample (2000 06/29/2016) 250

if Building_Model == 1 % 1 if Building Model, 0 if Outside Model. 
Contam_Soil  =   1;    % 1 if surface soil contamination is present, 0 if no contaminated soil 
VIMS_Present =   1;    % 1 if VIMS is be present, 0 if no VIMS 
I_Present    =   0;    % 1 if infiltration is present, 0 if no infiltration 
Bare_ground  =   0;    % 1 if bare ground which acts as infinite sink
Bio_Deg      =   1;    % 1 if Degradation on, 2 if off. 
else 
Contam_Soil  =   0;    % 1 if surface soil contamination is present, 0 if no contaminated soil 
VIMS_Present =   0;    % 1 if VIMS is be present, 0 if no VIMS 
I_Present    =   1;    % 1 if infiltration is present, 0 if no infiltration 
Bare_ground  =   1;    % 1 if bare ground which acts as a sink
Bio_Deg      =   1;    % 1 if Degradation on, 2 if off. 
end


