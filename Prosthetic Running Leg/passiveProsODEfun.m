function dtheta_dt = passiveProsODEfun(t,theta,b,k,theta_sp,L,GRFvals,GRFtime)

% First, we need to interpolate within our GRF data to get an estimate of
% GRF at the current time point.
F = [interp1(GRFtime,GRFvals,t)]; % **Fill this in with a function that interpolates GRFvals 
            %   at time t. Hint: look up "interp1".**


dtheta_dt = [F/b*L*cosd(theta)-k/b*(theta-theta_sp)]; %**Fill in your differential equation from part a. Be sure
                %  to use only the variables we've defined in this
                %  function: theta, b, k, theta_sp, L and F.