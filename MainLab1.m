% Contributors:
% Course Number: ASEN 3801
% File Name: MainLab1
% Last Updated: 8/26/25

clc
clear
close all

%% Nonzero Initial Conditions
w = 1;
x = 1;
y = 1;
z = 1;

%% Derivative Equations
wdot = -9*w + y;
xdot = 4*w*x*y - x^2;
ydot = 2*w - x - 2*z;
zdot = x*y - y^2 - 3*(z^3);

%% tspan and tolerance
tspan = [0 20];

%% Calling ODE45
[tw,w] = ode45(@(tw,w) wdot,tspan,w);
[tx,x] = ode45(@(tx,x) xdot,tspan,x);
[ty,y] = ode45(@(ty,y) ydot,tspan,y);
[tz,z] = ode45(@(tz,z) ydot,tspan,z);

%% Plotting ODE45
plot1(tw,w); %plotting w vs t
plot2(tx,x); %plotting x vs t
plot3(ty,y); %plotting y vs t
plot4(tz,z); %plotting z vs t