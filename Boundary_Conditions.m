function  [MC,x1,x2,x3] = Boundary_Conditions(BC,n,...
           elapsed,Contam_Soil,n_Soil_Contam,...
           nodes_contam,SDENS,ncomp,MASSF_soil,MW)
% % This function updates and keeps track of the mass removed, molar concentration and labels
% for "VIMs On", "VIMS Off" and "Soil Contamination" depending on the user
% set boundary conditions. 

% G(1:nin)=G_new(:,1:ncomp);
% flux=(.5/dx)*(D_nodes(1:end-1)+D_nodes(2:end)).*diff(G);
% massremoved(n+1)=massremoved(n)+flux(1)*dt;
% q_send=q_air*(1+2*q_var*(-0.5+rand()));  % This sends Q+-20% put this into boudnary conditions? 
% MC(1,1:ncomp)=MC(1,1:ncomp)-Gnew(1,1:ncomp)*void(1); % Remove node 1 mass from total molar concentration. 
% Update vapor discharge vector conditions for Building Model and Outside
% Models based on elapsed time. 
% if Building_Model == 1 
%     % if n==BC(1,2) % Slab condition 
%     %    q_air = q_air_Slab; 
%     % end 
%     if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % Slab condition 
%        q_air = q_air_slab; 
%     end
%     if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
%      x1 = BC(2,2);x2 = BC(2,3); % Used for plotting label. 
%      q_air = q_air_VIMS; 
%     end
% end   

% if n==BC(1,2) % Pre VIMS Segment 
%         massremoved(n,2:end)= 0;
% end 
% if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % No mass removed if VIMS is off
%        massremoved(n,2:end) = massremoved(n-1,2:end);
% end

if Contam_Soil == 1  
    if n == 19 %  
        Surf_Contam = zeros(nodes_contam,1);  
        n
        Ctot_Surf(1:length(Surf_Contam)) = linspace(BC(2,4),BC(2,5),nodes_contam); % Avg soil contamination
        elapsed_Soil_Contam = elapsed/86400;
        % Labels based on VIMS operation and soil contamination dates. 
        x3 = elapsed_Soil_Contam; % Used for plotting label. 
        for i = 1:nodes_contam
        MC(i,1:ncomp) = Ctot_Surf(i)*(SDENS/1e6)*MASSF_soil(1:ncomp)./MW(1:ncomp); % Update molar concentrations using total mass conc. from contaminated suface soil. (Nodes 1-4)
        end
    end
end


end % End function