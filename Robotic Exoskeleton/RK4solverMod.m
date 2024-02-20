function [tout, yout] = RK4solver(ODEfunIn, solveTime, y0, stepSize)
% This is an implementation of the fourth-order Runge_Kutta method. Inputs are as follows:
% - ODEfunIn is a function handle to a differential equation
% - solveTime is the range of interest of the independent variable
% - y0 is the initial condition
% - stepSize is step size for the solver, in the same units as solveTime

% Define the values of t at which we will evaluate our function
tVals = solveTime(1):stepSize:solveTime(2);

% Allocate an array
n = length(y0);
yout = zeros(n, length(tVals));

% Start us off at our initial conditions
yout(:,1) = y0;

% Amount of steps taken
nSteps = width(yout)-1;
             
for stepNum = 1:nSteps
    
    % Grab the last value of y and current value of t.
    y_n = yout(:,stepNum);
    t_n = tVals(stepNum);
    
    % Calculate the dydt's using the ODE function passed in
    kn1=ODEfunIn(t_n, y_n);
    kn2=ODEfunIn(t_n+1/2*stepSize, y_n+1/2*stepSize*kn1);
    kn3=ODEfunIn(t_n+1/2*stepSize, y_n+1/2*stepSize*kn2);
    kn4=ODEfunIn(t_n+stepSize, y_n+stepSize*kn3);
    
    % Calculate y_nplus1 using the Runge-Kutta method, and store it in the output array
    y_nplus1 = y_n+stepSize*(kn1+2*kn2+2*kn3+kn4)/6;
    yout(:, stepNum+1) = y_nplus1; 
end

% Rename this to be more intuitive as output
tout = tVals;