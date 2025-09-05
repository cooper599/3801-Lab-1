% Contributors: Cooper, Ridley, Alessandro, Gunnar
% Course Number: ASEN 3801
% File Name: MainLab1
% Last Updated: 8/26/25

clc; clear; close all;
%% Q1 Parts A & B
% For Part A calculate the solution to the dynamical system and plot the solution
ic = [1; 1; 1; 1]; % Nonzero Initial Conditions

% Generating the tspan and tolerance vectors
tspan = [0 20];
tol = [1e-2; 1e-4; 1e-6; 1e-8; 1e-10; 1e-12];

% Vectors for saving the t=20 value for each tolerance
finalWVec = [];
finalXVec = [];
finalYVec = [];
finalZVec = [];

% Calling ODE45, runs through each tolerance to get the needed values
for i = 1:length(tol) % Runs through each tolerance and generates the solution 
    options = odeset(RelTol = tol(i), AbsTol = tol(i));
    [t, statevec] = ode45(@(t, statevec) odefnc(t, statevec), tspan, ic, options);
    finalWVec = [finalWVec, statevec(end,1)];
    finalXVec = [finalXVec, statevec(end,2)];
    finalYVec = [finalYVec, statevec(end,3)];
    finalZVec = [finalZVec, statevec(end,4)];

    if i == 4 % Grabs the values for tol=1*10^-8 to plot
        w = statevec(:, 1);
        x = statevec(:, 2);
        y = statevec(:, 3);
        z = statevec(:, 4);

        figure(1) % Subplot for the graphs
        subplot(4, 1, 1)
        plot(t, w) %plotting w vs t
        title('w vs t')
        xlabel('t (n.d.)')
        ylabel('w (n.d.)')
        subplot(4, 1, 2)
        plot(t, x); %plotting x vs t
        title('x vs t')
        xlabel('t (n.d.)')
        ylabel('x (n.d.)')
        subplot(4, 1, 3)
        plot(t, y); %plotting y vs t
        title('y vs t')
        xlabel('t (n.d.)')
        ylabel('y (n.d.)')
        subplot(4, 1, 4)
        plot(t, z); %plotting z vs t
        title('z vs t')
        xlabel('t (n.d.)')
        ylabel('z (n.d.)')
        sgtitle('{\bf Question 1a, Dynamical System vs Time, Tolerances: 1x10^{-8}}', 'FontSize', 10)
        print("Lab1Q1a", "-dpng");
    end
end

%% Q1 Part B, finding differences between the tolerance and reference for the final value of the solution
% Creating each field of the table 
differenceTaken = ["w-w_r"; "x-x_r"; "y-y_r"; "z-z_r"];
tol1 = [abs(finalWVec(1) - finalWVec(6)); abs(finalXVec(1) - finalXVec(6));...
    abs(finalYVec(1) - finalYVec(6)); abs(finalZVec(1) - finalZVec(6))];%10^-2
tol2 = [abs(finalWVec(2) - finalWVec(6)); abs(finalXVec(2) - finalXVec(6));...
    abs(finalYVec(2) - finalYVec(6)); abs(finalZVec(2) - finalZVec(6))];%10^-4
tol3 = [abs(finalWVec(3) - finalWVec(6)); abs(finalXVec(3) - finalXVec(6));...
    abs(finalYVec(3) - finalYVec(6)); abs(finalZVec(3) - finalZVec(6))];%10^-6
tol4 = [abs(finalWVec(4) - finalWVec(6)); abs(finalXVec(4) - finalXVec(6));...
    abs(finalYVec(4) - finalYVec(6)); abs(finalZVec(4) - finalZVec(6))];%10^-8
tol5 = [abs(finalWVec(5) - finalWVec(6)); abs(finalXVec(1) - finalXVec(6));...
    abs(finalYVec(5) - finalYVec(6)); abs(finalZVec(5) - finalZVec(6))];%10^-10

% Generating the table with appropriate rows and columns
differences = table(differenceTaken, tol1, tol2, tol3, tol4, tol5);


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
lgd = legend(names,location="eastoutside"); % Adds legend using name list corresponding to wind speeds
title(lgd,"Wind Speed (m/s)");
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
    tempx2(i)= q2d.(names(i)).dtot(end);
end
plot(tempx2 - min(tempx2(1)),windspeedvec);
xlabel("Meters of Deflection (m)");
ylabel("Wind Speed (m/s)");
title("Total Distance from Origin to Landing Location vs Wind Speed in m per m/s")
xlim([ -5 2]);
ylim([-1 21]);
legend(names,location="best");

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
% Equations of Motion of Object through air accounting for drag, gravity,
% and initial conditions.
%  Takes in properties from object Cd, A, m, and starting conditions. Takes
%  in air density, windspeed, and gravity. And returns the velocity and
%  acceleration of the object passing it into the original statevector of
%  position and velocity.
    v = x(4:6);
    vRel = v-vWind;
    airspeed = norm(vRel);
    D = 0.5*rho*(airspeed^2)*A*Cd;
    fDrag = -D*(vRel/airspeed);
    fGrav = [0;0;m*g];
    acc = (fDrag+fGrav)/m;
    xDot = [v; acc];
end

%% Part F
T1 = 0.5*m*(sqrt(20^2+20^2));
m = 50e-3:10e-2:1;
V = sqrt(T1./(0.5*m));
spd = 1:length(m);
vel = 1:length(m);
figure();
hold on;


for i=1:length(m)
    spd(i) = sqrt(2*T1/m(i));
    vel(i) = sqrt(0.5*(spd(i)^2));
    sVec = [0;0;0;0;vel(i);-vel(i)];
    tspan = [0 20];
    wind_vel = 0; % Wind velocity [m/s]
    tol = 1e-8; % Tolerance for ode45 call
    options = odeset('RelTol',tol,'AbsTol',tol,'Events',@groundhit); % Setting tolerance options and termination event
    [t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m(i),g,wind_vel),tspan,sVec,options); % Call of ode function
    plot3(statevector(:,1),statevector(:,2),statevector(:,3));
    set(gca, 'ZDir', 'reverse');
    view(45,45);
    title('Question 2f, varying mass vs distance (constant KE), windspeed 0m/s')
    xlabel('X/North (m)');
    ylabel('Y/East (m)');
    zlabel('Z/down (m)');
end

legendEntries = cell(size(m));
for i = 1:length(m)
    legendEntries{i} = ['m = ', num2str(m(i)), 'kg'];
end
legend(legendEntries)

figure()
hold on;

for i=1:length(m)
    spd(i) = sqrt(2*T1/m(i));
    vel(i) = sqrt(0.5*(spd(i)^2));
    sVec = [0;0;0;0;vel(i);-vel(i)];
    tspan = [0 20];
    wind_vel = 10; % Wind velocity [m/s]
    tol = 1e-8; % Tolerance for ode45 call
    options = odeset('RelTol',tol,'AbsTol',tol,'Events',@groundhit); % Setting tolerance options and termination event
    [t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m(i),g,wind_vel),tspan,sVec,options); % Call of ode function
    plot3(statevector(:,1),statevector(:,2),statevector(:,3));
    set(gca, 'ZDir', 'reverse');
    view(45,45);
    title('Question 2f, varying mass vs distance (constant KE), windspeed 10m/s')
    xlabel('X/North (m)');
    ylabel('Y/East (m)');
    zlabel('Z/down (m)');
end

legendEntries = cell(size(m));
for i = 1:length(m)
    legendEntries{i} = ['m = ', num2str(m(i)), 'kg'];
end
legend(legendEntries)

figure()
hold on;

for i=1:length(m)
    spd(i) = sqrt(2*T1/m(i));
    vel(i) = sqrt(0.5*(spd(i)^2));
    sVec = [0;0;0;0;vel(i);-vel(i)];
    tspan = [0 20];
    wind_vel = 20; % Wind velocity [m/s]
    tol = 1e-8; % Tolerance for ode45 call
    options = odeset('RelTol',tol,'AbsTol',tol,'Events',@groundhit); % Setting tolerance options and termination event
    [t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m(i),g,wind_vel),tspan,sVec,options); % Call of ode function
    plot3(statevector(:,1),statevector(:,2),statevector(:,3));
    set(gca, 'ZDir', 'reverse');
    view(45,45);
    title('Question 2f, varying mass vs distance (constant KE), windspeed 20m/s')
    xlabel('X/North (m)');
    ylabel('Y/East (m)');
    zlabel('Z/down (m)');
end

legendEntries = cell(size(m));
for i = 1:length(m)
    legendEntries{i} = ['m = ', num2str(m(i)), 'kg'];
end
legend(legendEntries)

figure()
hold on;

for i=1:length(m)
    spd(i) = sqrt(2*T1/m(i));
    vel(i) = sqrt(0.5*(spd(i)^2));
    sVec = [0;0;0;0;vel(i);-vel(i)];
    tspan = [0 20];
    wind_vel = 40; % Wind velocity [m/s]
    tol = 1e-8; % Tolerance for ode45 call
    options = odeset('RelTol',tol,'AbsTol',tol,'Events',@groundhit); % Setting tolerance options and termination event
    [t,statevector] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m(i),g,wind_vel),tspan,sVec,options); % Call of ode function
    plot3(statevector(:,1),statevector(:,2),statevector(:,3));
    set(gca, 'ZDir', 'reverse');
    view(45,45);
    title('Question 2f, varying mass vs distance (constant KE), windspeed 40m/s')
    xlabel('X/North (m)');
    ylabel('Y/East (m)');
    zlabel('Z/down (m)');
end

legendEntries = cell(size(m));
for i = 1:length(m)
    legendEntries{i} = ['m = ', num2str(m(i)), 'kg'];
end
legend(legendEntries)

%% Functions
% Function for problem 1 dynamical system
function dstatevec_dt = odefnc(t, statevec)
% Houses the derivative functions that we are trying to solve. Takes the
% value of the statevector at that instance and solves the following
% derivative values
    w = statevec(1);
    x = statevec(2);
    y = statevec(3);
    z = statevec(4);

    dstatevec_dt = [-9*w + y; 
                    4*w*x*y - x^2;
                    2*w - x - 2*z; 
                    x*y - y^2 - 3*(z^3)]; % dw/dt, dx/dt, dy/dt, dz/dt
end
