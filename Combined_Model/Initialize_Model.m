function[MC,MC_Before,KD,ACT,VP,SOL,MW,SDENS,MASSF_soil,nodes_contam,nc,ncomp,G,W,...
G_Top,G_WT,W_WT,C_Top, C_Top_All, C_Bottom_All, W_ug_L, W_GW, CtotCurrent_All,...
         CtotCurrent_Sum, PHASE_History, W_ChlSol, W_Deg_plot,lambda_years,massremoved,massremoved_dt,massremoved_g,massremoved_g_tot, elapsed]...
         = Initialize_Model(MASSF_soil_PEST, MASSF_GW_PEST, por,theta,dx,nin,DENSLIQ,SDENS,...
         TEMP,FOC,Ctot)
% G,W,T,N are concentrations in gas, water, total, and NAPL phases, solubility (mg/L), Kow, and a flag for printing results.

% Initalize model domain for function. 
elapsed=0;
soil_contam_depth=150;  % Depth interval for contaminated Soil Contamination Depth (cm)
n_solvents = 4;         % number of solvents used for biodegradation 
nodes_contam = round(soil_contam_depth/dx+0.5); 
% Load chem data for soil:
% Columns are :mass fraction, molec. weight (g/mol), boiling point (C), 
% Vapor pressure (atm) at STP, solubility (mg/L), Kow, flag for printing results.

% % 25 Degrees C % % % % 
% Remember to adjust the soil temp in parameters.csv  
chem_data = ... 
[0.01  131.4	87	    0.091	1280	407.38   1; ... % 1 TCE
 0.600 165.8	121.1	0.024	206	    2511.89  1; ... % 2 PCE
 0.364 96.95	60.2	0.263	6410	72.44    1; ... % 3 Cis-1,2-DCE
 0.001 92.14	110.6	0.037	526	    537.03   1; ... % 4 Toluene
 0.021 62.5	    -13.4	3.329	2763	23.99    1; ... % 5 Vinyl Chloride
 0.003 106.17   140	    0.011	106	    1318.26  1; ... % 6 Total Xylenes
 0.001 106.168	136.2	0.013	177	    1412.54  1];    % 7 Ethylbenzene

NCOMP=size(chem_data,1);        % Number of chemical components (i.e. PCE, TCE, 1,2-DCE, etc.)
% MASSF_soil=chem_data(:,1);      % Masss Fraction for soil contamination 
MASSF_soil = MASSF_soil_PEST;   % Allow PEST to change Mass Fraction for Soil. 

% Separate mass fractions for GW contamination:
% MASSF_GW = [0.006 0.003 0.922 0.003 0.0056 0.008 0.002]';    % Mass fraction average from MW-1 GW samples collected between 2016-2024)
MASSF_GW = MASSF_GW_PEST;
MASSF_GW=MASSF_GW./sum(MASSF_GW);                           % Mass fraction for groundwater contamination 
MW=chem_data(:,2); BP=chem_data(:,3);
VP=chem_data(:,4); SOL=chem_data(:,5); KD=chem_data(:,6); 
ICOMPFL=logical(chem_data(:,7));
ACT=zeros(size(SOL));
TEMP=TEMP+273;                                              % Convert temperature from Celclius to Kelvin 

% Make a vector of initial total soil concentrations (mg/kg)
nnodes=length(por);

% Setup array for molar concentration stored at each node
MC=zeros(nnodes-1,NCOMP+1);

%    Numerical constants etc.
eps=1e-10; tiny=1.0d-20;
TMOLE=0.0;

REF=273+25;
mwoil=0.0;

for I=1:NCOMP
    BPK=BP(I)+273.0;
    if(abs(BPK-REF)>eps)
        VP(I)=VP(I)*exp((BPK*REF/(BPK-REF))* ...
                          (1.0/TEMP-1.0/REF)*log(VP(I)));
    end
end 
    MASSF_soil=MASSF_soil./sum(MASSF_soil); % Correct the mass fractions so they sum to one, then calculate the average molecular weight of the oil.
	mwoil=sum(MASSF_soil.*MW);              % Only do this for shallow soil DNAPL

for	I=1:NCOMP
       if(SOL(I)>0.0)
	    ACT(I)=55.55*MW(I)/(SOL(I)/1000.0);
       else 
		ACT(I)=-SOL(I)*mwoil/(18.0*DENSLIQ);
       end
	 if(VP(I)>1.0) 
       ACT(I)=ACT(I)/VP(I);
     end
   KD(I)=0.63*KD(I)*FOC;
end
VPWAT=0.023*exp((373*293/(373-293))*(1/TEMP-1/293)* ...
               log(0.023));

% Calculate the molar concentration at each node. 
for i=1:nodes_contam  % 
    MC(i,1:NCOMP)=(Ctot(i)*SDENS/1.0d06)*MASSF_GW(1:NCOMP)./MW(1:NCOMP);
end
MC(1:nnodes-1,NCOMP+1)=(1/18).*theta(1:nnodes-1,2);   % Put in the molar conc of water

% Calculate  the molar concentration for lowest (water table) node
for i=nodes_contam+1:nnodes
        MC(i,1:NCOMP)=(Ctot(i)*SDENS/1.0d06)*MASSF_GW(1:NCOMP)./MW(1:NCOMP);
end

% Initalize arrays used in model functions
MC(nnodes,NCOMP+1)=(1/18).*theta(nnodes,2);   % Put in the molar conc of water
nc = size(MC,2);  ncomp=nc-1;               % nc = quantitiy of compounds incl. H2O. ncomp = quantitiy of VOCs. 
G=zeros(nin+1,nc); W=zeros(nin+1,nc);       % Save vapor conc.(mol/cc) and water conc.(mol/cc) during simulation. 
% G(:,nc)=VPWAT/(82.1*TEMP);                % Tested including H2O vapor concentration, but effect on model output is minimual. Not included in model.
G_Top=zeros(100,1+ncomp);                   % Save vapor conc. at top node (surface) throughout simulation. 
G_WT=zeros(100,1+ncomp); 
W_WT=zeros(100,1+ncomp); 
C_Top=zeros(2,1+ncomp);                     % Save sum of total mass conc. (sum of all contaminants) at top node (surface) throughout simulation. 
C_Top_All=zeros(100,ncomp);                 % Save total mass conc. for each of the contaminants at top node (surface) throughout simulation. 
C_Bottom_All=zeros(100,ncomp);              % Save total mass conc. for each of the contaminants at bottom node (water table) throughout simulation. 
W_ug_L=zeros(100,ncomp);                    % Save dissolved phase concentrarion (ug/L).
W_GW=zeros(100,ncomp);                      % Save dissolved phase concentration (ug/L).
CtotCurrent_All=zeros(100,ncomp);           % Current total mass conc. for each contaminant at all nodes.
CtotCurrent_Sum=zeros(2,1+ncomp);           % Current sum of total mass conc. of all contaminants for each node. 
PHASE_History = zeros(100,nin+2);           % Save phase history for all nodes throughout the simulation.\
W_ChlSol =  W(:,[2 1 3 5]);                 % W_ChlSol = Initalize array of chlorinated solvent dissolved phase concentrations. PCE, TCE, cis-1,2-DCE, VC
W_Deg_plot=zeros(1,5); W_Deg_plot(1,2:end)=W_ChlSol(round(nin/2),:); % Initalize array for saving the concentration of chlorinated solvent degraded with each timestep.  
W_Deg_Sum_plot=zeros(1,5); 
lambda_years = [1.1 1.2 1.2 0.0];            % Biodegradation lambda values from article (1/years)
massremoved=zeros(100,1+ncomp);              % Used to save time elapsed (seconds) and VIMS mass removal throughout simulation
massremoved_dt= zeros(size(massremoved)); 
massremoved_g= zeros(size(massremoved)); 
massremoved_g_tot = zeros(size(massremoved)); 
G_new=G(1:nin); 
MC_Before = MC; 
end % end function

