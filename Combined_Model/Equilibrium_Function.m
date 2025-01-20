function [G,W,MC,D_nodes,void,PHASE,calc_H,theory_H,Retardation,totmass] = Equilibrium_Function(G,W,MC,por,por_v_min,MW,VP,ACT,KD,SOL,...
    airdiff,DENSLIQ,SDENS,TEMP)

% To start out, I'll make this function work for 1-D domains passed here
% G,W,T,N are concentrations in gas, water, total, and NAPL phases

% % Some constants particular to this contamination event
% DENSLIQ     =   1.6;     % density of NAPL in grams/cc
% SDENS       =   1.5;     % density of dry soil in gm/cc
% TEMP        =   15;      % temperature of soil (deg C)
% airdiff     =   0.084;   % free-air diffusion coeff. (cm^2/s)

TEMP=TEMP+273;   % Convert to Kelvin

%    Numerical constants etc.
eps=1e-10; tiny=1.0d-20;
nc=size(MC,2);      % total number of compounds incl. water
NCOMP=nc-1;        % Number of non-water species
nnodes=size(MC,1);
TMOLE=0.0;
G=0*G; % Why?? 
VPWAT=0.023*exp((373*293/(373-293))*(1/TEMP-1/293)* ...
               log(0.023));

void=zeros(nnodes,1);
DELWAT=zeros(nnodes,1);
PHASE=3*ones(nnodes,1);
%EVAC=VOID;
%EVAC(EVAC<0)=eps;
CONCMASS=zeros(nnodes,1);
CONCMOLE=zeros(nnodes,1);
totmass=zeros(nnodes,1);

%blah2=MC(:,1:NCOMP)
% Any previous error that produces negative molar concentration zeroed out
MC(MC<0)=0;
for i=1:nnodes
   totmass(i)=(1e6/SDENS)*sum(MC(i,1:NCOMP).*MW(1:NCOMP)');
end
void=por-totmass.*(SDENS/(1.0d06*DENSLIQ)) - MC(:,nc)*18.0;  %***
void(void<0)=1e-6;

%pordiff=airdiff./(por.*por);

%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%   CHECK FOR FREE PHASE IN ANY OF THE COMPOUNDS.
%   IF ALL ARE 3-PHASE, CALCULATE VAPOR CONC. OF EACH
%   COMPOUND AND RETURN. IF 4-PHASE, ITERATE TO THE PHASE DISTRIBUTION 
%	Only check those compounds that are immiscible

DELWAT(MC(:,nc)>0)=1;  % Check to see if there is soil water 
C1=zeros(size(por)); C2=zeros(size(por)); C4=zeros(size(por));
for N=1:NCOMP
%	if(sol(n).gt.0.0)then    % Only for immiscible compounds
    CONCMASS=CONCMASS+MC(:,N)*MW(N);
    CONCMOLE=CONCMOLE+MC(:,N);
    C1=VP(N)*void/(82.1*TEMP);
    C2=MC(:,nc)/ACT(N);
%	    if(sol(n).gt.0.0)C2=MC(I,J,K,nc)/ACT(N)
    C3=KD(N)*SDENS*DELWAT/(18.0*ACT(N));
    C4=MC(:,N)./(C1+C2+C3);
%  if the compound is miscible, don't count it - doesn't force a NAPL
    if SOL(N)>0
        PHASE(C4>1)=4;
        G(:,N)=C4*(VP(N)/(82.1*TEMP));  % get gas conc for 3-phase
%    if (MC(I,J,K,N).GT.0.0) 
%        RETARD(N)=MC(I,J,K,N)/(CI(I,J,K,N)*VOID)
%    end
    end
end
% PHASE
%
%CCCCCCCCCCCCCCCCCCCCCCC
%
%     IF THERE IS A SEPARATE PHASE, MAKE A FIRST GUESS THAT
%     THE SEPARATE PHASE MOLE FRACTION (X(N)) EQUALS THE TOTAL 
%     MOLES OF N DIVIDED BY THE TOTAL CONTAMINANT MOLES, OR
%     CI*RT/VP(N) = MI/CONCMOLE.
%
%       IF(PHASE(I,J,K).gt.3.5) THEN
%	 DO 20 N=1,NCOMP
%	    CI(I,J,K,N)=MC(I,J,K,N)*VP(N)/(CONCMOLE*82.1*TEMP)
%  20    CONTINUE
apple=find(PHASE>3.5);
for i=1:length(apple)
    ii=apple(i);
    G(ii,1:NCOMP)=MC(ii,1:NCOMP).*VP'./(CONCMOLE(ii)*82.1*TEMP);
%
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%
%    THIS IS THE ITERATIVE LOOP WHERE, ESSENTIALLY, 
%    MI(N)=(CONST+MCHC)CI(N). SO SOLVE FOR MCHC, THEN CI(N), ETC. 
%    UNTIL SUM(CI(N)*R*T/VP(N)), WHICH = SUM(X(N))=1.000
%
	 faze3=0;
	 for NN=1:NCOMP
	   c1=void(ii);
	   c2=MC(ii,nc)*TEMP*82.1/(ACT(NN)*VP(NN));
	   c3=KD(NN)*SDENS*DELWAT(ii)*82.1*TEMP/(ACT(NN)*VP(NN)*18.0);
	   faze3=faze3+(G(ii,NN)*(c1+c2+c3));
     end
	 MCHC(ii)=CONCMOLE(ii)-faze3;

% ---- WITH NEW MCHC, RECALCULATE CI(N), AND CHECK THE SUM OF X(N)

	 NUMITER=0; sumx=0.0;
     while(abs(sumx-1.0)>0.01)
        NUMITER=NUMITER+1;
	    for N=1:NCOMP
	      k1=82.1*TEMP/VP(N);
	      k2=void(ii);
   	      k3=MCHC(ii)*k1;
	      k4=MC(ii,nc)*k1/ACT(N);
	      k5=KD(N)*SDENS*DELWAT(ii)*k1/(ACT(N)*18.0);
	      G(ii,N)=MC(ii,N)/(k2+k3+k4+k5);
        end

%    CHECK THAT THE SUM OF FREE PRODUCT MOLE FRACTIONS = 1.0
%    IF NOT, ADJUST THE G(N) VALUES AND REITERATE.
	   
       sumx = 82.1*TEMP*sum( (G(ii,1:NCOMP)./(VP(1:NCOMP))' ) );
	   MCHC(ii)=MCHC(ii)*sumx;
     end

end      % Done going through the 4-phase nodes 

%MC(:,NCOMP+1)=(1/18)*theta(:);   % Put in the molar conc of water
%G(:,NCOMP+1)=VPWAT/(82.1*TEMP);      % H2O gas concentration


for i=1:nnodes   % calculate the dissolved conc. in moles/cc
%   W(i,1:NCOMP)=1e6*TEMP*82.1*G(i,1:NCOMP).*MW(1:NCOMP)'./ ...
%                      (18.0*ACT(1:NCOMP).*VP(1:NCOMP))';  % Dissolved conc. in mg/L
   W(i,1:NCOMP)=TEMP*82.1*G(i,1:NCOMP)./ ...
                       (18.0*ACT(1:NCOMP).*VP(1:NCOMP))';  % Dissolved conc. in mole/cc


pordiff=airdiff./(por.*por);  
void=por-totmass.*(SDENS/(1.0d06*DENSLIQ)) - MC(:,nc)*18.0;  %***
void(void<0)=1e-6;
D_nodes= pordiff.*void.^3.333; % Calculate effective diffusion using Millington Quirk 
theory_H=18.0*ACT.*VP/(TEMP*82.1);
calc_H=G(nnodes,:)./W(nnodes,:);
Retardation = MC(:,1:NCOMP)./(G(:,1:NCOMP).*void);
Retardation_Mean = mean(Retardation);
end