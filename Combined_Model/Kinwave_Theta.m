function [theta,lwaves,twaves] = Kinwave_Theta(theta_max,theta_min,K_sat,BC_lambda,P_air,theta, ...
                                          lwaves,twaves,elapsed,dt,I,L)
% This function moves kinematic waves (waves) according to infiltration times
% series (I), and maps it to the finite difference grid (theta).  The
% "grid" for the waves is probably finer than the finite-difference grid.
% Make sure to pass air entry pressure 
BC_eps=(2+3*BC_lambda)/BC_lambda;
n_trail=10;
n_sit=50;  % DAB new stationary (hydrostatic) "waves"
small=1e-10;
Sy=theta_max-theta_min;     % notational simplicity

if(elapsed<small)              % Initialize profile with trailing waves
    % Setup theta inital with a 1x100 grid from surface (0) to the water table. 
    % Inital setup yeilds column 2 being equal to the inital theta value (0.2)

    theta_initial=zeros(100,2);
    theta_initial(:,1)=linspace(1e-6,L,100);
    theta_initial(:,2)=interp1(theta(:,1),theta(:,2),theta_initial(:,1),'linear','extrap');
    % Setup trailing wave matrix with 3 columns. 
    % Column 1 = Depth 
    % Column 2 = Theta 
    % Column 3 = Velocity  
    twaves=zeros(100,3);
    twaves(:,1:2)=theta_initial;
    lwaves(1,2)=theta_min; lwaves(1,1)=1e-3; 
    lwaves(1,:)=[]; % Deletes the first row of l-wave at the end.
    
    I_current=interp1(I(:,1),I(:,2),elapsed,'previous');  % Current infilt. rate Columns =sec, column 2 = cm/sec. 
    I_last=interp1(I(:,1),I(:,2),elapsed); % Last time step infiltration (previous infilt. rate)
    theta_current=theta_min+Sy*(I_current/K_sat)^(1/BC_eps); % Update theta_current based on infiltration. 
    new_wave=zeros(n_trail,3); 
    new_wave(:,2) = linspace(theta_current,theta_initial(1,2),n_trail);
    new_wave(:,1) = linspace(0,1e-6-small,n_trail);
    twaves=[new_wave; twaves]; % Reset twaves with new wave pasted the top of twaves. 
else
    % If time is after t=0, then update I current and I last, both equal to 0.5. 
    % Check if I_current < K_sat, Check if I_last < K_sat. 
    % If no change then update theta last with theta max. This keeps checking until there 
    % is a change in infilration and theta
    % current - theta last is no longer zero. 
    I_current=interp1(I(:,1),I(:,2),elapsed,'linear','extrap');  % Current infilt. rate
    I_last=interp1(I(:,1),I(:,2),elapsed-dt,'linear','extrap'); % last time step infilt. rate
end
I_current=max(I_current,0);
I_last=max(I_last,0);

if I_current<K_sat 
   theta_current=theta_min+Sy*(I_current/K_sat)^(1/BC_eps);
   theta_current=max([theta_current theta_min]);
   theta_current=min([theta_max theta_current]);
else
   theta_current=theta_max;
end
if I_last<K_sat 
   theta_last=theta_min+Sy*(I_last/K_sat)^(1/BC_eps);
   theta_last=max([theta_last theta_min]);
   theta_last=min([theta_max theta_last]);
else
   theta_last=theta_max;
end
% If the water content is greater than the prevoius timestep, add infilatration with leading wave. 
% Add on a new row to lwaves for 0 ft (start depth at top of domain), current water content, and 0 velocity. 
if theta_current-theta_last > small   % Increase in infiltration
    % Throw a leading wave [depth, theta, v]
    %disp('making a lead wave')
    lwaves=[0, theta_current 0; lwaves];
    
    
end

% Decrease infiltration if the previous timestep water content is greater than the current water contnet. 
% Then add new wave using theta_last. 30 rows long. Column 1 is depth. 
% Column 2 theta linspace going from current theta down to theta min. 
% Basically an array of the thetas as long as the drying front needs to be. 
% Then tack on the new wave to the top of the twaves array. 
if theta_last-theta_current > small  % Decrease in infiltration
%    disp('Throwing a big trailing wave');
    % throw a trailing wave with n_trail segments, each row [depth, theta, v]
    new_wave=zeros(n_trail,3);
    new_wave(:,2) = linspace(theta_current,theta_last,n_trail);
    new_wave(:,1) = linspace(0,1e-6,n_trail);
%    new_wave(:,1) = theta_min+Sy*(P_air./new_wave(:,1)).^BC_lambda;
    twaves=[new_wave; twaves];
end
%  If no infiltration then add a trailing wave to the top of the domain
%  (new_wave on top of twaves). 
if abs(theta_last-theta_current)<=small & theta_current>theta_min;  % No change in infiltration
    % throw a single trailing wave [depth, theta, v]
    % disp('should be here?');
    new_wave=[0 theta_current 0];
    twaves=[new_wave; twaves];
end

% Calculate velocity of trailing wave
    twaves(:,3) = (BC_eps*K_sat/Sy)*((twaves(:,2)-theta_min)/Sy).^(BC_eps-1);

% Calculate velocities of all lead waves. First locate lead waves below the
% current lead wave. 
for i=1:size(lwaves,1) % Search through all the lead waves
    depth_now=lwaves(i,1); % Keep track of current depth of lead wave thats being looped
    blah=[lwaves; twaves]; % Combine leadwaves with trailing waves to compare their velocity. 
    apple=find(blah(:,1)>(depth_now+1e-6)); % Find all the points that are below the lead wave and recall theta. 
    [x, idx]=min(blah(apple,1)-depth_now); % Indexing used to find the closest wave to current. 
     % Find next trailing wave segment ahead of this lead wave
    if(~isempty(apple)) %  % Now loop through all the thetas atpoints below the lead wave (these are thetas) 
        theta_dn=blah(apple(idx),2); % Grab theta from the lowest point. 
    else
        theta_dn=theta_min;   % Catch an error if this lead wave is deepest of all
    end
    % Use newly calculated thetas to calculate K_up and K_dn
    K_up=K_sat*((lwaves(i,2)-theta_min)/Sy)^BC_eps; 
    K_dn=K_sat*((theta_dn-theta_min)/Sy)^BC_eps;
%    catch any zero/zero
    lwaves(i,3)=(K_up-K_dn)/(lwaves(i,2)-theta_dn); % Then use the Kup Kdn and Thetaup and Thetadn to find the velocity
    if isnan(lwaves(i,3))
        lwaves(i,3)=0.0;
    end
end

% Now mark the depths for each wave at the new timestep. With conv. of momentum. X = VT 
% The size of the matrix is 1 x many; with the size of the lead waves and trailing waves profile.  
% These locations will be used to determine if waves take over other waves. 
new_depths_lead  = lwaves(:,1)+lwaves(:,3)*dt;

new_depths_trail = twaves(:,1)+twaves(:,3)*dt;

% see if any lead waves overtake trailing wave segments
for i=1:size(lwaves,1)
    old_depth=lwaves(i,1);
    new_depth=new_depths_lead(i,1);
    kill=find(twaves(:,1)>=old_depth & new_depths_trail(:,1)<=new_depth); % All waves that are greater than the old depth, but less than the new depth. 
    twaves(kill,:)=[]; % Removes the waves marked "kill" 
    new_depths_trail(kill,:)=[]; % Now remove the depths from new depths trailing waves, since those are deleted. 
end

% see if any lead waves overtake lead waves
% Run loop for a variable amount of time, equal to the size of lwaves. Note that lwaves will be cut down during the loop, 
% as lead waves are overtaken by the current lead wave, and removed from the lead wave matrix.  The length of the loop
% is adjusted as you go though each wave that is overtaken. Example: If the current wave is 5, and it overtakes waves 6 and 7,
% the next in line would have been wave 8, but since wave 8 is now wave 6, since the previous wave 6 and 7 are deleted in the array. 
% This is why the number of loops required is equal to the size of l_waves. 

length_loop=size(lwaves,1);
i=1;
while i<length_loop
    new_depth=new_depths_lead(i,1);
    kill=find(new_depth>new_depths_lead(:,1)); % Locatioins where new depth is greater than new depths lead lead wave.  
    kill(kill<=i)=[];
    lwaves(kill,:)=[];
    new_depths_lead(kill,:)=[];
    length_loop=length_loop-length(kill);
    i=i+1;
end

% See if any trailing waves will overtake the updated lead waves
% Use the spatially highest to set new theta of lead wave and delete the
% trailers. this occurs when a shock takes place, and the wetting front is
% moving into dry soil with lower K, thus the wetting front is moving
% slower than what is behind it. 
for i=1:size(lwaves,1)
    candidates=find(twaves(:,1)<lwaves(i,1) & new_depths_trail(:,1)>new_depths_lead(i,1));
    best=min(candidates);
    if(size(best,1)>0)
        lwaves(i,2)=twaves(best,2);
        twaves(candidates,:)=[];
        new_depths_trail(candidates,:)=[];
    end
end

% Go ahead and finalize the lead wave positions now:
lwaves(:,1)=new_depths_lead(:,1);

% Finalize trailing wave segments positions
twaves(:,1) = new_depths_trail(:,1);

% DAB this is new: calculate the moisture content of every wave relative
% to the hydrostatic moisture profile.  Kill all waves that are dryer than 
% (are absorbed by) the upward capillary profile

sitwaves=zeros(n_sit,3);
sitwaves(:,1)=linspace(0,L,n_sit);
sitwaves(:,2)=theta_max;
apple=find(sitwaves(:,1)<(L-P_air));
sitwaves(apple,2)= theta_min+(theta_max-theta_min)*(P_air./(L-sitwaves(apple,1))).^BC_lambda;

hydro_stat=interp1(sitwaves(:,1),sitwaves(:,2),lwaves(:,1));
lwaves(:,2)=max(lwaves(:,2),hydro_stat);

hydro_stat=interp1(sitwaves(:,1),sitwaves(:,2),twaves(:,1));
twaves(:,2)=max(twaves(:,2),hydro_stat);


% Kill allwaves deeper than the water table
lwaves(lwaves(:,1)>L,:)=[];
lwaves(lwaves(:,1)==L,:)=[];
twaves(twaves(:,1)>L,:)=[];
twaves(twaves(:,1)==L,:)=[];

% Finally, map the fine-scale water contents to the coarser grid used for
% transport etc.
%all=[lwaves; twaves; sitwaves];

all=[lwaves; twaves];

all=[-1e-6 all(1,2) all(1,3); all];
[blah,ia,ic]=unique(all(:,1));
%all=[all; WT+1e-6 all(end,2) all(end,3)];
all=all(ia,:);


%theta(:,2) = interp1(all(:,1),all(:,2),theta(:,1),'nearest');
theta_wave = interp1(all(:,1),all(:,2),theta(:,1),'linear','extrap');
theta_sit  = interp1(sitwaves(:,1),sitwaves(:,2),theta(:,1),'linear','extrap');
theta(:,2) = max(theta_wave,theta_sit);
theta(:,2) = min(theta_max,theta(:,2));
theta(:,2) = max(theta_min,theta(:,2));

 % figure(300)
 % plot(theta(:,2),L-theta(:,1),'+-k')
 % hold on
 % plot(twaves(:,2),L-twaves(:,1),'bd')
 % plot(lwaves(:,2),L-lwaves(:,1),'ro')
 %  axis([0 theta_max min(L-theta(:,1)) L])
 % ylabel('Height above WT (cm)'); xlabel('Vol. Soil Moisture')
 % title(['Lead and Trailing Waves; Time = ',num2str(elapsed/86400),' days'])
 % legend('Soil Moisture','Trailing wave','Leading wave')
 % hold off
 % drawnow

end
