function[W,W_Deg]=Biodeg_Solvents(W,lambda_years,elapsed,nin)

% This function takes a row vector of N concentrations and calculates
% the chain reaction C1 -> C2 -> ... CN, with a corresponding
% row vector of rates k1, k2, ... kN, where C1 -> C2 has rate k1, etc.
C =  W(:,[2 1 3 5]);

% This uses mixed implicit/explicit for 2nd order accuracy.
lambda_sec = lambda_years/31536000; % [PCE, TCE, CIS, VC] (sec^-1)
j = 1;
nc = length(C(1,:));
C_deg = zeros(nc);
C_0 = C(round(nin/2),:);
while j < nin
C_node = C(j,:);
% Make matrix
H = 1 + 0.5*elapsed*lambda_sec;
J = -0.5*elapsed*lambda_sec;
A = diag(H)+diag(J(1:end-1),-1);

% Make know vector
E = 1 - 0.5*elapsed*lambda_sec;
F = 0.5*elapsed*lambda_sec;
known = C_node.*E;
known(2:end) = known(2:end) + F(1:end-1).*C_node(1:end-1);
known=known';

C_node=A\known;
C(j,:)=C_node';
j = j+1;
end
C_1 = C(round(nin/2),:); 
C_deg =C_1-C_0; 

% % Check for 99% degradation of CIS and save the time only once - Doesn't
% work, need to update. 
%     if C(1, 3) <= C_99 && ~recorded_99
%         Time_days_99 = elapsed/86400; % Save the time (days) when 99% of CIS is degraded, convert to hours
%         recorded_99 = true; % Set the flag to true so it doesn't record again
%         fprintf('99%% degradation of CIS reached at %.2f days\n', Time_days_99); % Export a note
%     end
W(:,[2 1 3 5]) = C; % Save new dissolved phase concentrations for dissolved solvents after first order deg.  
W_Deg = C_deg;  % Save change dissolved phase concentrations for dissolved solvents after first order deg. 
end 