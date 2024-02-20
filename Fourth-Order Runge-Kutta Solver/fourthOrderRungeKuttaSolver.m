%% Problem 4b
% In this problem, we are going to compare the performance of the RK4
% method to Euler's method.

%Beginning of code
clear;
clc;
close;

%Variables that make our plots look nice
fs = 14; % Font size
fn = 'Arial'; %Font name
lw = 2; % Linewidth
primColor = [39 116 174]/255; %UCLA Blue
secColor = [255 209 0]/255; %UCLA Gold

% First, we want to numerically solve our implicit solution to the differential equation.
impSolFun = @(a,t) 0.05*(log(abs(a))+3*a)+3.4/pi*cos(pi/2*t)+t-0.05*(log(0.01)+.03)-3.4/pi;

% For this part, we're only interested in looking at a limited time window
tImplicit = 0.75:0.001:2;

% Now we are going to interatively use fzero() to find a numerical solution 
% to the implicit equation
aImplicit = zeros(size(tImplicit));

for tNum = 1:length(tImplicit)    

    t_cur = tImplicit(tNum);
    % Here, we define an anonymous function "afun" that take only the real
    % part of our implicit solution function, to make it easy for fzero to
    % solve.

    afun = @(a) real(impSolFun(a,t_cur));
   
    % We provide an initial guess in the ballpark
    a_guess = 0.01;
    
    % fzero() does the magic. Google fzero() to see how this works.
    a_cur = fzero(afun,a_guess);

    % Store the current value of a
    aImplicit(tNum) = a_cur;
    
end

%Vector of varying step sizes to be used for Euler's and the RK4 method.
stepSizeVals = [0.1, 0.075, 0.05, 0.01, 0.005];

a_0 = 0.01;
t_range = [0,2]; 

%Euler's Method Graphs
figure;
subplot(1,2,1);
plot(tImplicit, aImplicit, 'linewidth', lw, 'Color', 'k');
legendIDs = cell(1,length(stepSizeVals)+1);
legendIDs{1} = 'Implicit Solution';
hold on
ctr = 1;
for stepSizeNum = 1:length(stepSizeVals)
    curStepSize = stepSizeVals(stepSizeNum);
    [t_cur, a_cur] = eulerSolver(@(t,a) muscleActODEfun(t,a), t_range, a_0, curStepSize);
    plot(t_cur, a_cur, 'linewidth',lw-0.5)
    ctr = ctr + 1;
    legendIDs{ctr} = ['Step Size = ' num2str(curStepSize) ' s'];
end
hold off;
title("Euler's Method");
xlabel('Time (s)');
ylabel('Muscle Activation Level');
xlim([.75, 2]);
ylim([0, 1.2]);
grid on;
legend(legendIDs,'location','northwest','box','on');
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');

%RK4 Method Graphs
subplot(1,2,2);
plot(tImplicit, aImplicit, 'linewidth', lw, 'Color', 'k');
legendIDs = cell(1,length(stepSizeVals)+1);
legendIDs{1} = 'Implicit Solution';
hold on;
ctr = 1;
for stepSizeNum = 1:length(stepSizeVals)
    curStepSize = stepSizeVals(stepSizeNum);
    [t_cur, a_cur] = RK4solver(@(t,a) muscleActODEfun(t,a), t_range, a_0, curStepSize);
    plot(t_cur, a_cur, 'linewidth',lw-0.5)
    ctr = ctr + 1;
    legendIDs{ctr} = ['Step Size = ' num2str(curStepSize) ' s'];
end
hold off;
title("RK4 Method");
xlabel('Time (s)');
ylabel('Muscle Activation Level');
xlim([.75, 2]);
ylim([0, 1.2]);
grid on;
legend(legendIDs,'location','northwest','box','on');
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');

%Muscle Activation Differential Equation Function
function da_dt = muscleActODEfun(t,a)
    da_dt = (1.7*a*sin(pi/2*t)-a)/(.1*(.5+1.5*a));
end



