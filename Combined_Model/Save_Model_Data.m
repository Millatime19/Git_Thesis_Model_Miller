function[G_Top, G_Top_ug_m3, G_ug_m3, G_WT, G_WT_ug_m3, C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    W_WT, W_ug_L, W_GW_ug_L,  CtotCurrent_All, CtotCurrent_Sum, dt, elapsed] = ...
    Save_Model_Data(G_Top, G_WT, W_WT,C_Top, C_Top_All, C_Bottom_All, PHASE_History,...
    elapsed, ncomp,G_new,n,SDENS,MC,MW,PHASE,W,dt,dt_max,nin)
% % % % % % Save Phase concentrations for plotting and update time % % % % % %
    G_Top(n,1) = elapsed;                                 % Save top node vapor phase concentration during simulation (mol/cc)
    G_Top(n,2:end) = G_new(1,1:ncomp);                    
    G_WT(n,1) = elapsed;                                 % Save WT node vapor phase concentration during simulation (mol/cc)
    G_WT(n,2:end) = G_new(15,1:ncomp); 
    W_WT(n,1) = elapsed;                                     % Save top node vapor phase concentration during simulation (mol/cc)
    W_WT(n,2:end) = W(16,1:ncomp);   
    C_Top(n,1) = elapsed;                                 % Save total mass (mg/kg)
    C_Top(n,2) = sum((1e6/SDENS)*MC(1,1:ncomp).*MW');     % Total mass (mg/kg) (top node)
    C_Top_All(n,:) = (1e6/SDENS)*MC(1,1:ncomp).*MW';      % Save separate masses for each contaminant (top node)
    C_Bottom_All(n,:) = (1e6/SDENS)*MC(nin,1:ncomp).*MW'; % Save bottom node total mass conc.                  % Save separate masses for each contaminant (second to last node)
    PHASE_History(n,1) = elapsed; 
    PHASE_History(n,2:end) = PHASE';                                        % Save NAPL History Info
    G_Top_ug_m3 = zeros(length(G_Top),1+ncomp);
    G_Top_ug_m3(:,1) = G_Top(:,1); 
    G_Top_ug_m3(:,2:end) = G_Top(:, 2:end) .* (MW' * 1e12); % Conversion factor: for mol/cc to µg/m³  
    
    G_WT_ug_m3 = zeros(length(G_WT),1+ncomp);
    G_WT_ug_m3(:,1) = G_WT(:,1); 
    G_WT_ug_m3(:,2:end) = G_WT(:, 2:end) .* (MW' * 1e12); % Conversion factor: 1e9 for mol/cc to µg/m³  
    G_ug_m3 = G_new(:,1:ncomp).*(MW' * 1e12);
    
    
    W_ug_L=W(:,1:ncomp); W_ug_L=1e9*W_ug_L.*MW';                            % Save dissolved phase (all nodes)
    W_GW_ug_L = zeros(length(W_WT),1+ncomp);
    W_GW_ug_L(:,1) = W_WT(:,1); 
    W_GW_ug_L(:,2:end) = W_WT(:, 2:end).* (MW' * 1e9);
    
    % W_GW_ug_L(n,:) = W_ug_L(16,1:ncomp);                                         % Save dissolved phase (second to last node)
    dt = min(1.01*dt, dt_max);                                              % Update dt 
    elapsed = elapsed+dt;                                                   % Keep track of elapsed time
    CtotCurrent_All = sum(1e6/SDENS)*MC(:,1:ncomp).*MW';                    % Ctot for each node and contaminant at the end of the simulation
    for i = 1:nin+1
    CtotCurrent_Sum(i) = sum(CtotCurrent_All(i,:));                        % Sum of Ctot at end of simulation for each node
    end 
    CtotCurrent_Sum = CtotCurrent_Sum';

   end 