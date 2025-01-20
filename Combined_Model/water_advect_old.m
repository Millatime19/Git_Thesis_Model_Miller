function [W] = water_advect(theta,W,K_sat,theta_max,theta_min,BC_lambda,I,elapsed,dt,dx)
% This function moves kinematic waves (waves) according to infiltration times
% series (I), and maps it to the finite difference grid (theta).  The
% "grid" for the waves is probably finer than the finite-difference grid.
% Make sure to pass air entry pressure 

theta(:,2)=min(theta_max, theta(:,2));     % catch any errors here

BC_eps=(2+3*BC_lambda)/BC_lambda;
Sy=theta_max-theta_min;                 % notational simplicity

I_current=interp1(I(:,1),I(:,2),elapsed);  % Current infilt. rate
I_current=min([K_sat I_current]);

q_up=zeros(size(W,1),1);
ncomp=size(W,2)-1;
q_up(1)=I_current;                       % only needed if contam water entering
q_up(2:end)=K_sat*((0.5*(theta(1:end-1,2)+theta(2:end,2))-theta_min)/Sy).^BC_eps;

%Courant=max((dt/dx)*q_up./theta(:,2));

% Simple upwind-weighting:
W(2:end,1:ncomp)  = W(2:end,1:ncomp)+(dt/dx)*q_up(2:end).*W(1:end-1,1:ncomp);      % Flux in from tops
W(1:end-1,1:ncomp)= W(1:end-1,1:ncomp)-(dt/dx)*q_up(2:end).*W(1:end-1,1:ncomp);    % Flux out from bottoms
W(end,1:ncomp)=W(end,1:ncomp)-(dt/dx)*q_up(end)*W(end,1:ncomp);                    % Flux in=flux out at WT (bottom)
%Wnew(:,1:ncomp)=W(:,1:ncomp)./theta(:,2)

end
