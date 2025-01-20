function [G_new] = Vapor_Diffuse(G,D_nodes,dt,dx,q)

% As written, specifies a reflecting boundary at top (i.e., pavement)
% and a specified concentration at the bottom (like water table).
% The lower specified conc. are held in G(nin+1,:). 
% ALso note that I commented out 2 lines of the diffusion function.  
% The q in the code and the write-up is the volumetric VIMS rate Q divided
% by its effective area.  Units need to agree in the model (minutes? seconds?
% I can’t rememer), but it would be, for example, 10 cfm over a 10 x 10 foot
% area = 10 ft^3/min/100 ft^2 = 0.1 ft/minute.  The flux out of the model can
% be calculated (like I did in the example), which goes back into QC by
% multiplying by that area again. 

 nin = size(G,1)-1; nc= size(G,2); G_new=zeros(nin,nc); 
 % Get internode D: D(1) is D(1 1/2) etc. 
 D_half(1:nin) = 0.5*(D_nodes(1:nin)+D_nodes(2:nin+1));
 %D_half(1:nin) = sqrt(D_nodes(1:nin).*D_nodes(2:nin+1));  % Geom. mean
 r=dt/dx/dx;
 A=diag(-r*D_half(1:end-1),-1)+diag(-r*D_half(1:end-1),+1);
 maindiag=1+r*D_half;   
 maindiag(2:end)=maindiag(2:end)+r*D_half(1:end-1);
 A=A+diag(maindiag);

% Adjust for BCs - without q in equation
 %A(1,1)=A(1,1)+r*D_half(1);     % reflecting
 %A(1,2)=2*A(1,2);               % reflecting

 A(1,1)=A(1,1)+q*dt/dx;           % source/sink at node 1

 G(nin,:)=G(nin,:)+r*D_half(nin).*G(nin+1,:);  % constant at bottom
    
 % solve:
 G_new=A\G(1:nin,:);
 
end  % gas diffuse function