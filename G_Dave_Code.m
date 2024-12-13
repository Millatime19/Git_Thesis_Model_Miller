clear all; 
L           =   24.5;       % Vertical domain length units in ft
L           =   L*12*2.54;  % Vertical domain length units in cm  
L           =   15;
nin         =   15;         % Specify number of internal spatial nodes
dz          =   L/nin;      % Specify dz based on L and nin. 
lower_bdy_val = 5;          % 
dt=1; % Seconds
ntsteps = 500;
% dt_max = 100; 

% Setup flowrate - q. Q over area (L/T) (cm/s)
q_air_VIMS = 0.01;       % VIMS condition  
q_air_Open = 10000;      % Open Air condition
q_air_Slab = 0.00001;    % Slab  condtion
Surface = 1;             % 1 for slab, 2 for VIMS, 3 for Surface. 
%q_air=dz/dt
q_var=0.0; %


% Setup Effective Diffusion Coefficient
% D_nodes=logspace(2,-2,nin+1)'; % Set for variable diffusion 
D_nodes=(1e-1*ones(nin+1,1)); % Set for constant diffusion 

% G=lower_bdy_val+zeros(nin+1,1);  % Start at steady-state
G=0+zeros(nin+1,1);                % Start clean
G(nin+1)=lower_bdy_val;
Gnew=G(1:nin);       

mass_removed=zeros(1,ntsteps+1);

% Calculate analytical solution at specific distances
C0 = lower_bdy_val;    % Initial concentration at water table
D_avg = mean(D_nodes);     % Average diffusion coefficient (or use a specific value)
t_max = dt * ntsteps;  % Maximum time seconds (total simulation time)
z_values = dz * (0:nin);           % Depths to calculate
C_analytical = C0 * erfc(z_values ./ (2 * sqrt(D_avg * t_max)));
C_analytical= (flip(C_analytical)');

for n=1:ntsteps
    if Surface == 1
        q_air = q_air_Slab; 
    end
    if Surface == 2
        q_air = q_air_VIMS; 
    end
    if Surface == 3
        q_air = q_air_Open; 
    end
    q_send=q_air*(1+2*q_var*(-0.5+rand()));  % This sends Q+-20% 
    Gnew=gas_diffuse(G,D_nodes,dt,dz,q_send); 
    G(1:nin)=Gnew;
    G(1)
 %   G(1:2)=0;     % force zero top boundary
    flux=(.5/dz)*(D_nodes(1:end-1)+D_nodes(2:end)).*diff(G); % I need to understand this. 
    mass_removed(n+1)=mass_removed(n)+flux(1)*dt; % Mass out of the model  
   
% Plot analytical solution
% figure(34)
% plot (flux,L-dz*(0:(nin-1)),'-x')
% hold on 
% plot(C_analytical, L-dz*(0:(nin)), '-r', 'LineWidth', 2); hold on
% plot(G, L-dz*(0:(nin)), '-o');       % Plot numerical solution
% legend('Vapor Flux', 'Analytical Solution', 'Numerical Solution')
% xlabel('Concentration')
% ylabel('Depth (cm)')
% % set(gca, 'XScale', 'log');
% hold off
% drawnow

     % dt = min(1.01*dt, dt_max); 
end
figure(33)
    plot(G,L-dz*(0:(nin)),'-o'); hold on
    plot (flux,L-dz*(0:(nin-1)),'-x')
    legend('Conc.','flux')
    % set(gca, 'XScale', 'log');
    %axis([ 0 5*lower_bdy_val 0 (nin)*dz])
    drawnow
    hold off

    % Plot G and analytical solution
figure(34)
plot (flux,L-dz*(0:(nin-1)),'-x')
hold on 
plot(C_analytical, L-dz*(0:(nin)), '-r', 'LineWidth', 2); hold on
plot(G, L-dz*(0:(nin)), '-o');       % Plot numerical solution
legend('Vapor Flux', 'Analytical Solution', 'Numerical Solution')
xlabel('Concentration')
ylabel('Depth (cm)')
% set(gca, 'XScale', 'log');
hold off
drawnow

figure(35)
plot(dt*(0:ntsteps),mass_removed,'-o')


function [Gnew] = gas_diffuse(G,D_nodes,dt,dx,q)

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

 nin = size(G,1)-1; nc= size(G,2); Gnew=zeros(nin,nc); 
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
 Gnew=A\G(1:nin,:);

end  % gas diffuse function


% function [Gnew] = gas_diffuse(G,D_nodes,dt,dx)
% 
% % As written, specifies a reflecting boundary at top (i.e., pavement)
% % and a specified concentration at the bottom (like water table).
% % The specified concentrations are held in G(nin+1,:).
% % So only calculate Gnew(1:nin).
% 
%  nin = size(G,1)-1; nc= size(G,2);     % Number of internal nodes and compounds
%  Gnew=zeros(nin,nc);                   % Only solving on internal nodes
%  %D_half(1:nin) = 0.5*(D_nodes(1:nin)+D_nodes(2:nin+1)); % arith. mean
%  D_half(1:nin) = sqrt(D_nodes(1:nin).*D_nodes(2:nin+1));  % Geom. mean
%  r=dt/dx/dx;
%  A=diag(-r*D_half(1:end-1),-1)+diag(-r*D_half(1:end-1),1);
%  maindiag=zeros(1,nin);   
%  maindiag(2:end)=1+r*D_half(1:end-1)+r*D_half(2:end);
%  maindiag(1)=1+(2*dt/dx/dx)*D_half(1);    % treat D(-1/2) as = D(+1/2) 
%  %maindiag(1)=maindiag(1) + 0;    % treat D(-1/2) as = 0 
% 
%  A=A+diag(maindiag);
% 
%     % Adjust for BCs
%  %A(1,1)=A(1,1)+r*D_half(1); 
%  %A(1,2)=2*A(1,2);                % reflecting boundary at top
%  A(1,2)=-(2*dt/dx/dx)*D_half(1);                % reflecting boundary at top
% 
% % sum(A,2)
%  G(nin,:)=G(nin,:)+r*D_half(nin).*G(nin+1,:);   % making the known vector with BC info
% 
%  Gnew=A\G(1:nin,:);
%  Gnew(nin+1)=G(nin+1);
%  %Gnew=inv(A)*G(1:nin,:);
% 
% end  % gas diffuse function
