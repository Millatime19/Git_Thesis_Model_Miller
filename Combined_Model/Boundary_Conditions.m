function  [massremoved, MC_slab,MC_bare,Gnew_slab,Gnew_bare,x1,x2,x3] = Boundary_Conditions(BC,n,dx,MC_slab,MC_bare,...
           elapsed,massremoved,Contam_Soil,n_Soil_Contam,...
           nodes_contam,SDENS,ncomp,MASSF_soil,por,theta_slab,theta_bare,Gnew_slab,Gnew_bare,MW,x1,x2,x3)
% % This function updates and keeps track of the mass removed, molar concentration and labels
% for "VIMs On", "VIMS Off" and "Soil Contamination" depending on the user
% set boundary conditions. 

massremoved(n,1)=elapsed; % Update mass removed time column. 

rel_perm_slab=((por-theta_slab(:,2))./por).^3.33;   % Millington-Quirk relative permeability factor, unchanging under slab
rel_perm_bare=((por-theta_bare(:,2))./por).^3.33;   % Millington-Quirk relative permeability factor

% For now just use rel_perm to make flow-weighted averages
%flow_slab=sum(rel_perm_slab(1:1)/sum(rel_perm_slab(1:1)+rel_perm_bare(1:1));
%flow_bare=sum(rel_perm_bare(1:1)/sum(rel_perm_slab(1:1)+rel_perm_bare(1:1));

if n==BC(1,2) % Pre VIMS Segment 
        massremoved(n,2:end)= 0;
end 
if elapsed/86400 >= BC(1,2) && elapsed/86400 <= BC(2,2) % No mass removed if VIMS is off
       massremoved(n,2:end) = massremoved(n-1,2:end);
end

if Contam_Soil == 1  
    if n == n_Soil_Contam %  Incorporate surface contamination just before VIMS is on
        Surf_Contam = zeros(nodes_contam,1);  
        Ctot_Surf(1:length(Surf_Contam)) = linspace(BC(2,4),BC(2,5),nodes_contam); % Avg soil contamination
        elapsed_Soil_Contam = elapsed/86400;
        % Labels based on VIMS operation and soil contamination dates. 
        x3 = elapsed_Soil_Contam; % Used for plotting label. 
        for i = 1:nodes_contam
            % Update molar concentrations using total mass conc. from contaminated suface soil. (Nodes 1-4)
            MC_slab(i,1:ncomp) = Ctot_Surf(i)*(SDENS/1e6)*MASSF_soil(1:ncomp)./MW(1:ncomp); 
        end
    end
end

if elapsed/86400 > BC(2,2) && elapsed/86400 < BC(2,3) % Turn on VIMS  
     x1 = BC(2,2);x2 = BC(2,3); % Used for plotting label. 
%     massremoved(n,2:end)=massremoved(n-1,2:end)+dx*(por(1)-theta(1,2))*Gnew(1,1:ncomp); 
% Could make the flow a function of permeabilities, but for now just the
% ratio of concentrations Qbot_slab is from node 2 to 1; Qbot_bare is bare
% node 2 to 1; Qtop_bare is ATM to node 1 bare ground, all relative

     Qbot_slab=0.5*(rel_perm_slab(1)+rel_perm_slab(2));
     Qbot_bare=0.5*(rel_perm_bare(1)+rel_perm_bare(2));
     Qtop_bare=rel_perm_bare(1);
     denom=Qbot_slab+Qbot_bare+Qtop_bare;
     Qbot_slab=Qbot_slab/denom;
     Qbot_bare=Qbot_bare/denom;
     Qtop_bare=Qtop_bare/denom;  % Q's sum to one

     G1_bare = Qbot_bare*Gnew_bare(2,:)./(Qbot_bare+Qtop_bare);

     G1_slab = (Gnew_bare(1,:)*(Qtop_bare+Qbot_bare) + Qbot_slab*Gnew_slab(2,:))/ ...
               (Qtop_bare+Qbot_bare+Qbot_slab);

     massremoved(n,2:end)=massremoved(n-1,2:end)+...
                          dx * (por(1)-theta_slab(1,2))*G1_slab(1:ncomp);
 
     MC_slab(1,1:ncomp)=MC_slab(1,1:ncomp)-(Gnew_slab(1,1:ncomp)-G1_slab(1:ncomp))*(por(1)-theta_slab(1,2)); % Remove node 1 mass from MC if VIMS on
     MC_bare(1,1:ncomp)=MC_bare(1,1:ncomp)-(Gnew_bare(1,1:ncomp)-G1_bare(1:ncomp))*(por(1)-theta_bare(1,2)); % Remove node 1 mass from MC if VIMS on

         % Matt - need to put the newly updated top node G back into
         % Gnew_slab and Gnew_bare and pass back to main program

     Gnew_bare(1,:)=G1_bare;
     Gnew_slab(1,:)=G1_slab;

end
 
  if elapsed/86400 > BC(3,2)   % VIMS OFF  
     massremoved(n,2:end) = massremoved(n-1,2:end); % Keep total mass removed the same once VIMS off. 
  end
end
