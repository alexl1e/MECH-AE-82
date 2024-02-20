%% Problem 5a
%Build a mathematical model and analytically solve the resulting system of equations.

%Beginning of code
clc
clear
close all

%Declaring variables
R=1.03; %Units: Ohms
L=2.04e-4; %Units: H
Kt=4.48e-2; %Units: Nm/A
Kv=7.1*pi; %Units: rad/Vs
Jm=1.01e-5; %Units: kg/m^2
N=1; %Unitless
bm=0.001; %Units: Nms/rad

%Creating matrix
a=-R/L;
b=-N/Kv/L;
c=Kt/Jm/N;
d=-bm/Jm;
A = [a, b; c, d];
[V, D] = eig(A);
B=[0;10];
c=V\B;

%Displaying Output
disp("Eigenvalues:");
disp("r1: " + num2str(D(1,1)));
disp("r2: " + num2str(D(2,2)));
fprintf("\n");
disp("Eigenvectors:");
disp("v1:");
disp(V(:,1));
disp("v2:");
disp(V(:,2));
disp("Arbitrary Constants:")
disp("c1: " + num2str(c(1)));
disp("c2: " + num2str(c(2)));
%% Problem 5b
%Update your RK4 to handle systems of differential equations, and use it to solve the homogeneous equation.

%Variables that make our plots look nice
fs = 14; % Font size
fn = 'Arial'; %Font name
lw = 2; % Linewidth
primColor = [39 116 174]/255; %UCLA Blue
secColor = [255 209 0]/255; %UCLA Gold

%Phase Portrait with Step Size of 0.001s
t_range1 = [0, 1];
ICs1  = [0; 10];
stepSize1 = 0.001; %Units: s
[t1, vals1] = RK4solverMod(@(t, iwout) exoBootsODEfun(t, iwout), t_range1, ICs1, stepSize1);
i1 = vals1(1,:);
wout1 = vals1(2, :);
figure;
plot(wout1, i1, 'linewidth', lw, 'Color', primColor);
title("Phase Portrait with Step Size of "+num2str(stepSize1)+"s");
xlabel('Angular Velocity (rad/s)');
ylabel('Current (A)');
grid on;

%Phase Portrait with Step Size of 0.0001s
stepSize2 = 0.0001; %Units: s
[t2, vals2] = RK4solverMod(@(t, iwout) exoBootsODEfun(t, iwout), t_range1, ICs1, stepSize2);
i2 = vals2(1, :);
wout2 = vals2(2,:);
figure;
plot(wout2, i2, 'linewidth', lw, 'Color', primColor);
title("Phase Portrait with Step Size of "+num2str(stepSize2)+"s");
xlabel('Angular Velocity (rad/s)');
ylabel('Current (A)');
grid on;

%Phase Portrait with Oscillation
Rmod=1.03/10; %Units: Ohms
[t3, vals3] = RK4solverMod(@(t, iwout) exoBootsODEfun(t, iwout, Rmod), t_range1, ICs1, stepSize2);
i3 = vals3(1,:);
wout3 = vals3(2,:);
figure;
plot(wout3, i3, 'linewidth', lw, 'Color', primColor);
title("Phase Portrait with Oscillation");
xlabel('Angular Velocity (rad/s)');
ylabel('Current (A)');
grid on;

%Plotting Angluar Velocity vs. Time and Current vs. Time
figure;
subplot(2,1,1);
plot(t2, wout2, 'linewidth', lw, 'Color', primColor);
title("Angular Velocity vs. Time");
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
grid on;
subplot(2,1,2)
plot(t2, i2, 'linewidth', lw, 'Color', primColor);
title("Current vs. Time");
xlabel('Time (s)');
ylabel('Current (A)');
grid on;

%% Problem 5c
%Solve the non-homogeneous equation while tracking an output trajectory.

%Importing data
exoBootOutput = load('exoBootOutput.mat');
exoAngle = exoBootOutput.exoAngle;
exoTorque = exoBootOutput.exoTorque;
exoTime = exoBootOutput.exoTime;

%Transmission Ratio = 1
%Using RK4 on the Nonhomogeneous ODE
N4=1;
t_range4 = [0, 0.75];
ICs4 = [0,0,0];
stepSize4 = 0.0001; %Units: s
[t4, vals4] = RK4solverMod(@(t, iwtheta) exoBootNHODEfun(t, iwtheta, N4, exoTorque, exoAngle, exoTime), t_range4, ICs4, stepSize4);
thetaout4=vals4(3,:);
wout4 = vals4(2,:);
i4 = vals4(1,:);

%Plotting Output Angle vs. Time, Angluar Velocity vs. Time, and Current vs. Time
%Output Angle vs. Time
figure;
subplot(3,1,1);
plot(t4, thetaout4, 'linewidth', lw, 'Color', primColor);
title("Output Angle vs. Time");
xlabel('Time (s)');
ylabel('Output Angle (rad)');
xlim([0, .75]);
grid on;
hold on;
theta_des4 = interp1(exoTime,exoAngle,t4);
plot(t4, theta_des4, 'linewidth', lw, 'Color', secColor);
%Angular Velocity vs. Time
subplot(3,1,2);
plot(t4, wout4, 'linewidth', lw, 'Color', primColor);
title("Angular Velocity vs. Time");
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
xlim([0, .75]);
grid on;
%Current vs. Time
subplot(3,1,3);
plot(t4, i4, 'linewidth', lw, 'Color', primColor);
title("Current vs. Time");
xlabel('Time (s)');
ylabel('Current (A)');
xlim([0, .75]);
grid on;

%Transmission Ratio = 5
%Using RK4 on the Nonhomogeneous ODE
N5=5;
[t5, vals5] = RK4solverMod(@(t, iwtheta) exoBootNHODEfun(t, iwtheta, N5, exoTorque, exoAngle, exoTime), t_range4, ICs4, stepSize4);
thetaout5=vals5(3,:);
wout5 = vals5(2,:);
i5 = vals5(1,:);

%Plotting Output Angle vs. Time, Angluar Velocity vs. Time, and Current vs. Time
%Output Angle vs. Time
figure;
subplot(3,1,1);
plot(t5, thetaout5, 'linewidth', lw, 'Color', primColor);
title("Output Angle vs. Time");
xlabel('Time (s)');
ylabel('Output Angle (rad)');
xlim([0, .75]);
grid on;
hold on;
theta_des5 = interp1(exoTime,exoAngle,t5);
plot(t5, theta_des5, 'linewidth', lw, 'Color', secColor);
%Angular Velocity vs. Time
subplot(3,1,2);
plot(t5, wout5, 'linewidth', lw, 'Color', primColor);
title("Angular Velocity vs. Time");
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
xlim([0, .75]);
grid on;
%Current vs. Time
subplot(3,1,3);
plot(t5, i5, 'linewidth', lw, 'Color', primColor);
title("Current vs. Time");
xlabel('Time (s)');
ylabel('Current (A)');
xlim([0, .75]);
grid on;

%Transmission Ratio = 100
%Using RK4 on the Nonhomogeneous ODE
N6=100;
[t6, vals6] = RK4solverMod(@(t, iwtheta) exoBootNHODEfun(t, iwtheta, N6, exoTorque, exoAngle, exoTime), t_range4, ICs4, stepSize4);
thetaout6=vals6(3,:);
wout6 = vals6(2,:);
i6 = vals6(1,:);

%Plotting Output Angle vs. Time, Angluar Velocity vs. Time, and Current vs. Time
%Output Angle vs. Time
figure;
subplot(3,1,1);
plot(t6, thetaout6, 'linewidth', lw, 'Color', primColor);
title("Output Angle vs. Time");
xlabel('Time (s)');
ylabel('Output Angle (rad)');
xlim([0, .75]);
grid on;
hold on;
theta_des6 = interp1(exoTime,exoAngle,t6);
plot(t6, theta_des6, 'linewidth', lw, 'Color', secColor);
%Angular Velocity vs. Time
subplot(3,1,2);
plot(t6, wout6, 'linewidth', lw, 'Color', primColor);
title("Angular Velocity vs. Time");
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
xlim([0, .75]);
grid on;
%Current vs. Time
subplot(3,1,3);
plot(t6, i6, 'linewidth', lw, 'Color', primColor);
title("Current vs. Time");
xlabel('Time (s)');
ylabel('Current (A)');
xlim([0, .75]);
grid on;

%% Problem 5d
%Optimize transmission ratio.

%Graphing the winding losses vs the transmission ratio
N7 = 2:2:24;
t_range7=[0,0.75];
ICs7 = [0,0,0];
stepSize7 = .0001; %Units: s
Pw7=zeros(size(N7));
R=1.03; %Units: Ohms
count=0;
for n = N7
    [t7, vals7] = ode45(@(t, iwtheta) exoBootNHODEfun(t, iwtheta, n, exoTorque, exoAngle, exoTime), t_range7, ICs7);
    i7 = vals7(:, 1);
    PwVals=i7.^2*R;
    avgPw = mean(PwVals);
    count = count+1;
    Pw7(count) = avgPw;
end
figure;
plot(N7, Pw7, 'linewidth', lw, 'color', primColor);
title("Winding Losses vs. Transmission Ratio");
xlabel('Transmission Ratio');
ylabel('Winding Losses');
grid on;

%Transmission Ratio = NOpt
%Using RK4 on the Nonhomogeneous ODE
[minPw, index]=min(Pw7);
NOpt = N7(index);
disp("Optimal Transmission Ratio = " + num2str(NOpt));
[tOpt, valsOpt] = RK4solverMod(@(t, iwtheta) exoBootNHODEfun(t, iwtheta, NOpt, exoTorque, exoAngle, exoTime), t_range7, ICs7, stepSize7);
thetaoutOpt=valsOpt(3,:);
woutOpt = valsOpt(2,:);
iOpt = valsOpt(1,:);

%Plotting Output Angle vs. Time, Angluar Velocity vs. Time, and Current vs. Time
%Output Angle vs. Time
figure;
subplot(3,1,1);
plot(tOpt, thetaoutOpt, 'linewidth', lw, 'Color', primColor);
title("Output Angle vs. Time");
xlabel('Time (s)');
ylabel('Output Angle (rad)');
xlim([0, .75]);
grid on;
hold on;
theta_desOpt = interp1(exoTime,exoAngle,tOpt);
plot(tOpt, theta_desOpt, 'linewidth', lw, 'Color', secColor);
%Angular Velocity vs. Time
subplot(3,1,2);
plot(tOpt, woutOpt, 'linewidth', lw, 'Color', primColor);
title("Angular Velocity vs. Time");
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
xlim([0, .75]);
grid on;
%Current vs. Time
subplot(3,1,3);
plot(tOpt, iOpt, 'linewidth', lw, 'Color', primColor);
title("Current vs. Time");
xlabel('Time (s)');
ylabel('Current (A)');
xlim([0, .75]);
grid on;

%Calculating peak current
disp("Peak Current = " +num2str(max(iOpt))+ " amps");