%% Problem 3b
% In this problem, we are going to use ODE45 to explore design prospects
% for a robotic tail.

% As before, define some things to make our plots look nice
fs = 14; % Font size
fn = 'Arial'; %Font name
lw = 2; % Linewidth
primColor = [39 116 174]/255; %UCLA Blue
secColor = [255 209 0]/255; %UCLA Gold

% Our first objective is to solve our system numerically, to see how our 
% tail handles to different external torques. Let's start by setting up 
% our system.

% Set the time range and parameter values
t_range = [0 30];
m = 0.4; % **Replace brackets with the correct value.**
L = 0.5; % **Replace brackets with the correct value.**
k = 0.2; % **Replace brackets with the correct value.**
b = 0.3; % **Replace brackets with the correct value.**
J = m*L^2; % **Calculate the rotational inertia J here from your other variables.**

% I'll help with initial conditions, because this one is a bit tricky.
theta_0 = -1; % **Replace brackets with initial conditions on theta.**
thetaD_0 = 0; % **Replace brackets with initial conditions on theta dot.**
y_0 = [theta_0;
       thetaD_0];

% Here, we need to create the external torque trajectories. I've done the 
% time vector and the homogeneous case for you.
torTime = t_range(1):0.01:t_range(2);
T1 = zeros(size(torTime));
T2 = 30.*torTime.*exp(1).^(-10.*torTime); % **Replace brackets with the nonhomogeneous external torque.**

% Use ODE45 to simulate the system response. Make sure you put your
% differential equation into robotailODEfun.m first! Hint: you'll need to 
% call ODE45 twice (once for each torque input). 

% **ADD CODE HERE TO FIND NUMERICAL SOLUTIONS FOR THE HOMOGENEOUS AND
% NONHOMOGENEOUS CASES**

[t_H, y_H] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T1,k,b,J), torTime, y_0);
[t_NH, y_NH] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2,k,b,J), torTime, y_0);

% Plot your results on top of each other. Both results should be plotted on
% the same axes in UCLA blue. The homogeneous case should be a dashed
% line, and the non-homogeneous a solid line. Include a legend in the 
% top-right corner, and use our standard formatting throughout. Be sure to 
% label your axes, including units.

figure(1)
% **ADD CODE HERE TO PLOT YOUR RESULTS**
plot(t_H, y_H(:,1), 'color', primColor, 'linestyle', '- -','LineWidth', lw);
hold on;
plot(t_NH, y_NH(:,1), 'Color', primColor, 'LineWidth', lw);
legend('Homogeneous Case', 'Non-Homogeneous Case');
ylabel('Angle (rad)');
xlabel('Time (s)');
grid on;
set(gca, 'FontSize', fs, 'FontName', fn, 'linewidth', lw, 'box', 'off');

%% Problem 3c

% In this problem, our goal is to simulate system behavior as we change one
% parameter at a time. To get you started, I've done the first one for you.

% Let's test k values that are half the original, the same as the original,
% and twice the original.
testKs = [0.5 1 2]*k;

% Set up some colors that I'll use in my plots.
cols = [primColor;
        0 0 0;
        secColor];
    
% Set up the figure, create a subplot. Note: if subplot is new, try
% entering "help subplot" into the command window.
figure(2); subplot(2,2,1)
hold on    

% Iterate through values of k, simulating and plotting each time
for kNum = 1:length(testKs)
    k_cur = testKs(kNum);
    [t_cur, y_cur] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2,...
                            k_cur, b, J),t_range, y_0);
    plot(t_cur,y_cur(:,1),'color',cols(kNum,:),'linewidth',lw)
end
hold off

% Apply some formatting and add a legend
legend({['k = ' num2str(testKs(1))],['k = ' num2str(testKs(2))],...
    ['k = ' num2str(testKs(3))]},'box','off','location','northeast')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off')

% **ADD CODE HERE TO DO THE SAME THING AS ABOVE FOR THE REMAINING
% PARAMETERES b, m, AND L. Don't forget to switch the subplot number!**
%Parameter b
testBs=[0.5 1 2]*b;
subplot(2,2,2);
hold on;
for bNum = 1:length(testBs)
    b_cur = testBs(bNum);
    [t_cur, y_cur] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2, k, b_cur, J),t_range, y_0);
    plot(t_cur,y_cur(:,1),'color',cols(bNum,:),'linewidth',lw)
end
legend({['b = ' num2str(testBs(1))],['b = ' num2str(testBs(2))],['b = ' num2str(testBs(3))]},'box','off','location','northeast')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');
hold off;

%Parameter m
testMs=[0.5 1 2]*m;
subplot(2,2,3);
hold on;
for mNum = 1:length(testMs)
    m_cur = testMs(mNum);
    [t_cur, y_cur] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2, k, b, m_cur*L^2),t_range, y_0);
    plot(t_cur,y_cur(:,1),'color',cols(mNum,:),'linewidth',lw)
end
legend({['m = ' num2str(testMs(1))],['m = ' num2str(testMs(2))],['m = ' num2str(testMs(3))]},'box','off','location','northeast')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');
hold off;

%Parameter L
testLs=[0.5 1 2]*L;
subplot(2,2,4);
hold on;
for lNum = 1:length(testLs)
    l_cur = testLs(lNum);
    [t_cur, y_cur] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2, k, b, m*l_cur.^2),t_range, y_0);
    plot(t_cur,y_cur(:,1),'color',cols(lNum,:),'linewidth',lw)
end
legend({['L = ' num2str(testLs(1))],['L = ' num2str(testLs(2))],['L = ' num2str(testLs(3))]},'box','off','location','northeast')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');
hold off;

% Now, we want to optimize our b value to reduce RMSE. I've added the time
% vector for you, and set up the range of bVals to test. You do the rest.
% Don't forget to interpolate each result to the testTime vector before
% calculating the RMSE.
testTime = 0:0.01:30;
bVals = 0.01:0.02:1;
RMSEs = zeros(size(bVals));

% **ADD CODE HERE TO ITERATE THROUGH THE bVals, CALCULATE RMSE FOR
% EACH, AND STORE THOSE VALUES**

for bVal = 1:length(bVals)
    b_cur = bVals(bVal);
    [time_cur, theta_cur] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2, k, b_cur, J),t_range, y_0);
    thetaAdj_cur = interp1(time_cur, theta_cur, testTime,'pchip','extrap');
    RMSE_cur = [sqrt(mean(thetaAdj_cur(:,1).^2));];
    RMSEs(bVal) = RMSE_cur;
end

% Now that we've optimized, take a look at the landscape. Plot your RMSEs
% versus the bVals that produced them. Use our standard formatting.
figure(3)
% **ADD CODE HERE TO PLOT YOUR RESULTS.**
plot(bVals, RMSEs, 'LineWidth', lw);
xlabel('Damping Constant, b (Nms/rad)');
ylabel('RMSE (rad)');
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off')

% Let's look at the optimized results. Find the minimum value of RMSE, as
% well as the bVal that produced it. Re-run your simulation using this
% bVal, and plot the resultant theta (UCLA Gold) on the same plot with of 
% the nonhomogeneous solution (UCLA Blue) from 3b. Use our standard
% formatting, and label your axes. Include a legend in the top-right
% corner.
%**ADD CODE HERE TO FIND OPTIMIZED RMSE, AND THE bVal THAT PRODUCED IT.**
[minRSMEval, minRMSEindex] = (min(RMSEs))
optbVal = bVals(minRMSEindex)
figure(4)
%**ADD CODE HERE TO PLOT YOUR RESULTS.**

[t_Hc, y_Hc] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T1,k,optbVal,J), torTime, y_0);
[t_NHc, y_NHc] = ode45(@(t,y) roboTailODEfun(t,y,torTime,T2,k,optbVal,J), torTime, y_0);
plot(t_Hc, y_Hc(:,1), 'color', secColor, 'LineWidth', lw);
hold on;
plot(t_NHc, y_NHc(:,1), 'Color', primColor, 'LineWidth', lw);
legend('Homogeneous Case', 'Non-Homogeneous Case');
ylabel('Angle (rad)');
xlabel('Time (s)');
grid on;
set(gca, 'FontSize', fs, 'FontName', fn, 'linewidth', lw, 'box', 'off');