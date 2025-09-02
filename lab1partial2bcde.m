% Contributors:
% Course Number: ASEN 3801
% File Name: MainLab1
% Last Updated: 8/26/25

clc; clear; close all;

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
% tspan = [0 20];

%% Calling ODE45
% [tw,w] = ode45(@(tw,w) wdot,tspan,w);
% [tx,x] = ode45(@(tx,x) xdot,tspan,x);
% [ty,y] = ode45(@(ty,y) ydot,tspan,y);
% [tz,z] = ode45(@(tz,z) zdot,tspan,z);

%% Plotting ODE45
% figure(); hold on;
% plot(tw,w); %plotting w vs t
% plot(tx,x); %plotting x vs t
% plot(ty,y); %plotting y vs t
% plot(tz,z); %plotting z vs t
% legend();

%% Q2 Parts B,C,D,E
% Part A equation to base future parts off of

% function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind_vel) function from 2a

%% Part B, calling atmos data for boulder 1655m
% Initial conditions need for object EOM
rho = stdatmo(1655); % Returns density in kg/m^3

Cd = 0.6; % Coefficient of Drag
diameter = 2.0 / 100; % Diameter of sphere [cm --> m]
A = pi * (diameter/2)^2; % Cross sectional area of sphere
m = 50 / 1000; % Mass [g --> kg]
g = 9.81; % Gravity [m/s^2]
wind_vel = 0; % Wind velocity [m/s]

%% Part C, given initial conditions, solve 2a equation
% Wind 0 in this part

% statevector: [x,y,z,vx,vy,vz]
tspan = [0 20]; % Time span 0 to 20 time units
initialcond = [0;0;0;0;20;-20]; % At origin (m), moving 20 m/s east and upwards
tol = 1e-8; % Tolerance for ode45 call

% Calling ode to find solution, set event to end ode when object reaches ground again.
function [value, isterminal, direction] = groundhit(t,statevector)
    % Event function to detect when the height of the object returns to 0
    value = statevector(3); % When object returns to z = 0 after the first time
    isterminal = 1; % To indicate to ode to stop running simulation
    direction = 0;
end

options = odeset('RelTol',tol,'AbsTol',tol,'Events',@groundhit); % Setting tolerance options and termination event
[t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind_vel),tspan,initialcond,options); % Call of ode function

% 3D Plot of object trajectory
figure();
plot3(statevector(:,1),statevector(:,2),statevector(:,3));
xlabel('X/North (m)');
ylabel('Y/East (m);');
zlabel('Z/down (m)');
title('Question 2c, 3D Trajectory');
set(gca, 'ZDir', 'reverse'); % Flipping z-axis

%% Starter code for 2D
% Same initial conditions as 2c, varying wind speed, plot all on 1 figure

% Creating vector of wind speed from 0 to 20 m/s
windspeedvec = linspace(0,20,21); % [m/s], o to ~45 mph winds
names = ["zero","one","two","three","four","five","six","seven","eight","nine","ten","eleven","twelve","thirteen","fourteen","fifteen","sixteen","seventeen","eighteen","nineteen","twenty"]; % List of names corresponding to each wind speed

% Looping through windspeed vector storing each call in question 2d
% structure and assigning name to each sub speed of wind vector with x,y,z
% information inside.
for i = 1:length(windspeedvec)
    [t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,windspeedvec(i)),tspan,initialcond,options);
    q2d.(names(i)).x = statevector(:,1); % X component
    q2d.(names(i)).y = statevector(:,2); % Y component
    q2d.(names(i)).z = statevector(:,3); % Z component
end

% Creating a plot with all wind speed trajectories on one
figure(); hold on;
% Loop through all wind speed vectors and call information to plot from
% structure
for i = 1:length(windspeedvec)
    plot3(q2d.(names(i)).x,q2d.(names(i)).y,q2d.(names(i)).z);
end
xlabel('X/North (m)');
ylabel('Y/East (m)');
zlabel('Z/down (m)');
title("Question 2d, Variable Wind Speeds on trajectories");
set(gca,'ZDir','reverse'); % Flipping z-axis
legend(names,location="eastoutside"); % Adds legend using name list corresponding to wind speeds
view(45,45);

% Part d1, creating plot of horizontal displacement vs wind speed
% In meters of deflection per m/s wind speed
figure(); hold on;
for i = 1:length(names)
    % scatter(q2d.(names(i)).x(end),windspeedvec(i)); % Scatter plot for naming points
    tempx(i) = q2d.(names(i)).x(end);
end
plot(tempx,windspeedvec); % Connecting Dots
xlabel("Meters of Deflection (m)");
ylabel("Wind Speed (m/s)");
% legend(names,location="best");
title("Question 2dp1, X Displacement vs Wind Speed");
xlim([-1 8]);
ylim([-1 21]);

% Part d2, creating plot of total distance
% Total distance from origin dependsing on wind speed in meters of
% deflection per m/s of wind speed
figure(); hold on;
for i = 1:length(names)
    q2d.(names(i)).dtot(i) = sqrt(q2d.(names(i)).x(end)^2 + q2d.(names(i)).y(end)^2);
    % scatter(q2d.(names(i)).dtot(end),windspeedvec(i));
    tempx2(i)= q2d.(names(i)).dtot(end);
end
plot(tempx2,windspeedvec);
xlabel("Meters of Deflection (m)");
ylabel("Wind Speed (m/s)");
title("Total ")
ylim([-1 21]);
% lgd = legend(names,location="best");
% title(lgd,"Wind Speed in m/s");

%% Question 2e
geopotential_altitude = [0 2000 4000 6000 8000 10000]; % [m], Vector of altitudes
altitudeNames = ["zero","two","four","six","eight","ten"]; % list of names corresponding to each altitude, x1000
for i = 1:length(altitudeNames)
    rhos(i) = stdatmo(geopotential_altitude(i));
end

% Creating figure for part 1 of e, 1 plot, multiple curves, each curve is 1
% geopotential altiude comparing total distance vs windspeed
for i = 1:length(altitudeNames) % Altitude loop
    for ii = 1:length(names) % Wind speed loop
        [t,statevector] = ode45(@(t,x) objectEOM(t,x,rhos(i),Cd,A,m,g,windspeedvec(ii)),tspan,initialcond,options);
        q2e.(altitudeNames(i)).(names(ii)).x = statevector(:,1); % X
        q2e.(altitudeNames(i)).(names(ii)).y = statevector(:,2); % Y
        q2e.(altitudeNames(i)).(names(ii)).z = statevector(:,3); % Z
        % Calculating Total Distances
        xend = q2e.(altitudeNames(i)).(names(ii)).x(end);
        yend = q2e.(altitudeNames(i)).(names(ii)).y(end);
        q2e.(altitudeNames(i)).dtot(ii) = sqrt(xend^2 + yend^2);
    end
end

figure(); hold on;
for i = 1:length(altitudeNames)
    plot(q2e.(altitudeNames(i)).dtot, windspeedvec);
end
xlabel("Distance (m)");
ylabel("Wind Speed (m/s)");
title("Distance vs Wind Speed of Various Geopotential Altitudes");
lgd = legend(altitudeNames,location="best");
title(lgd,"Altitude in 1000s of Meters")

% Plot 2 for part e, minimum distance between origin and landing point as
% function of geopotential altitude
figure(); hold on;
for i = 1:length(names)
    for ii = 1:length(altitudeNames)
        % scatter(q2e.(altitudeNames(i)).dtot,geopotential_altitude(i));
        tempdmin(ii) = q2e.(altitudeNames(ii)).dtot(i);
    end
    plot(tempdmin,geopotential_altitude)
end

xlabel("Minimum Distance (m)");
ylabel("Geopotential Altitude (m)");
title("Minimum Distance Between Origin and Landing Point as Function of Altitude");
ylim([-1000 11000]);
lgd = legend(names,location="eastoutside");
title(lgd,"Wind Speed (m/s)");

%% EOM Function
function xDot = objectEOM(t, x, rho, Cd, A, m, g, vWind)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    v = x(4:6);
    vRel = v-vWind;
    airspeed = norm(vRel);
    D = 0.5*rho*(airspeed^2)*A*Cd;
    fDrag = -D*(vRel/airspeed);
    fGrav = [0;0;m*g];
    acc = (fDrag+fGrav)/m;
    xDot = [v; acc];
end