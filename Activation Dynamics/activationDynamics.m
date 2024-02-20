%% Problem 2b
% In this problem, we are going to use several numerical tools to find a
% solution to a complex initial value problem.

% As before, define some things to make our plots look nice
fs = 14; % Font size
fn = 'Arial'; %Font name
lw = 2; % Linewidth
primColor = [39 116 174]/255; %UCLA Blue
secColor = [255 209 0]/255; %UCLA Gold


% First, we want to numerically solve our implicit solution to the
% differential equation. Hint: the ln function in matlab is "log()"
impSolFun = @(a,t) 0.05*(log(abs(a))+3*a)+3.4/pi*cos(pi/2*t)+t-0.05*(log(0.01)+.03)-3.4/pi;
                        % **put F_imp in place of the brackets, and delete 
                        % the brackets. Do not include the "=0" part of
                        % your eqn - only enter the left hand side. The
                        % only variables in this equation should be a and
                        %t**
% For this part, we're only interested in looking at a limited time window
tImplicit = 0.75:0.001:2;

% Now we are going to interatively use fzero() to find a numerical solution 
% to the implicit equation

aImplicit = zeros(size(tImplicit));
uImplicit = zeros(size(tImplicit));
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
    
    % Add some code here to calculate the current value of u(t), and store
    % it in the uImplicit array at the correct index.
    u_cur = 1.7*a_cur*sin(pi/2*t_cur);
    uImplicit(tNum) = u_cur;
end

% Plot your results. uImplicit should be plotted first as a function of
% time, in UCLA Blue. aImplicit should then be plotted against t on the
% same plot, in black. Use our standard formatting from the last problem
% set, with a legend in the top-left corner that identifies the two
% curves.
figure(1)
% ** PLOTTING CODE HERE **
plot(tImplicit, uImplicit, 'linewidth', lw, 'Color', primColor);
hold on;
plot(tImplicit, aImplicit, 'linewidth', lw, 'Color', 'k');
xlim([.75, 2]);
xlabel('Time (s)');
ylabel('Excitation/Activation Level');
title('Neural Excitation and Muscle Activation vs. Time');
grid on
legend({'Neural Excitation','Muscle Activation'}, 'box','off')
set(gca, 'FontSize', fs, 'FontName', fn, 'linewidth', lw, 'box', 'off')

%% Problem 2c
% Now, define initial conditions
%**Fill in initial conditions here**
a_0 = 0.01;

% Tell the solver how long to run
%**Fill in the start and end times here**
t_range = [0,2]; % should have the format [start_time, end_time]

% Set the step size for your euler solver
stepSize = 0.005;

% In this line, we run the solver. First, we'll need to edit 
% muscleActODEfun.m, and build our solver in eulerSolver.m.
[t_euler, a_euler] = eulerSolver(@(t,a) muscleActODEfun(t,a),...
                        t_range, a_0, stepSize);

                    
% Let's plot some results
figure(2)
plot(tImplicit,aImplicit,'k--','linewidth',lw)
hold on
plot(t_euler,a_euler,'linewidth',lw-0.5)
hold off
legend({'implicit solution','euler solution'},'box','off','location','northwest');
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off');


% What happens if we change the step size?

% Let's choose some step sizes to try.
stepSizeVals = [0.01, 0.005, 0.001, 0.0001];

% Create a figure on which to put everything
figure(3)
plot(tImplicit,aImplicit,'k','linewidth',lw);
legendIDs = cell(1,length(stepSizeVals)+1);
legendIDs{1} = 'implicit solution';
hold on
ctr = 1;

% Now we'll try them all out, and see what we get.
for stepSizeNum = 1:length(stepSizeVals)
    curStepSize = stepSizeVals(stepSizeNum);
    % Put code here to use your euler solver to solve your equation
    % again, using the current step size. Call your outputs t_cur and
    % a_cur
    % ** YOUR CODE HERE **
    [t_cur, a_cur] = eulerSolver(@(t,a) muscleActODEfun(t,a),...
                        t_range, a_0, curStepSize);
    
    % Add to the plot
    plot(t_cur, a_cur, 'linewidth',lw-0.5)
    ctr = ctr + 1;
    legendIDs{ctr} = ['step size = ' num2str(curStepSize) ' s'];
end
hold off
legend(legendIDs,'location','northwest','box','off')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off')

%% Problem 2d

% Put code here to solve the equation with ODE45. Call your output
% variables tODE45 and aODE45
% ** YOUR CODE HERE**
[tODE45, aODE45] = ode45(@(t,a) muscleActODEfun(t,a),t_range, a_0); 

% Plot aODE45 against tODE45, as well as your implicit solution.
figure(4)
plot(tImplicit,aImplicit,'k--','linewidth',lw)
hold on
plot(tODE45,aODE45,'linewidth',lw-1)                
                    
% Find ODE45's average step size                    
averageStepODE = mean(diff(tODE45)); % ** Replace the brackets with code to calculate the
                     % average step size from your ODE45 results. Hint:
                     % look at the functions mean() and diff().

% Now, we're going to run the Euler again, with the average step size from
% ODE45
[t_bigStep, a_bigStep] = eulerSolver(@(t,a) muscleActODEfun(t,a), t_range, a_0, averageStepODE); 
                             % ** Fill in the brackets with code to solve
                             % your equation using your eulerSolver
                             % function, with a stepSize equal to the
                             % average step size you just calculated for
                             % ODE45.

% Add this to your plot
figure(4)
hold on
plot(t_bigStep,a_bigStep,'linewidth',lw-1)
hold off
legend({'implicit solution','ode45','euler'},'location','northwest','box','off')
set(gca,'FontName',fn,'FontSize',fs,'linewidth',lw,'box','off')

%}