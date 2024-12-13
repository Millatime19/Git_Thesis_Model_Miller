function [W] = Water_Advect(theta,W,K_sat,theta_max,theta_min,BC_lambda,I,time,dt,dx)

theta(:,2)=min(theta_max, theta(:,2));   % catch any errors here
Sy=theta_max-theta_min;                  % notational simplicity

BC_eps=(2+3*BC_lambda)/BC_lambda;

I_current=interp1(I(:,1),I(:,2),time);   % Current infilt. rate
I_current=min([K_sat I_current]);

q_top=zeros(size(W,1),1);
ncomp=size(W,2)-1;
q_top(1)=I_current;                      % only needed if contam water entering

q_top(2:end)=K_sat*((0.5*(theta(1:end-1,2)+theta(2:end,2))-theta_min)/Sy).^BC_eps;

% Courant=max((dt/dx)*q_up./theta(:,2));
% Simple upwind-weighting:
% Flux in from tops:
W(2:end,1:ncomp)  = W(2:end,1:ncomp)+(dt/dx)*q_top(2:end).*(W(1:end-1,1:ncomp)./theta(2:end,2));
% Flux out from bottoms:
W(1:end-1,1:ncomp)= W(1:end-1,1:ncomp)-(dt/dx)*q_top(2:end).*(W(1:end-1,1:ncomp)./theta(1:end-1,2));
% Flux in=flux out at WT (bottom)
W(end,1:ncomp)=W(end,1:ncomp)-(dt/dx)*q_top(end)*W(end-1,1:ncomp)./theta(end); 
end
