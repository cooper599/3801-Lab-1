% Contributors: Ridley 
% Course Number: ASEN 3801
% File Name: MainLab1
% Last Updated: 09/02/25

clc
clear
close all

%% Nonzero Initial Conditions
ic = [1; 1; 1; 1];

%% tspan and tolerance
tspan = [0 20];
tol = [1e-2; 1e-4; 1e-6; 1e-8; 1e-10; 1e-12];

% Vectors for saving the t=20 value for each tolerance
finalWVec = [];
finalXVec = [];
finalYVec = [];
finalZVec = [];

%% Calling ODE45
% Runs through each tolerance to get the needed values
for i = 1:length(tol)
    options = odeset(RelTol = tol(i), AbsTol = tol(i));
    [t, statevec] = ode45(@(t, statevec) odefnc(t, statevec), tspan, ic, options);
    finalWVec = [finalWVec, statevec(end,1)];
    finalXVec = [finalXVec, statevec(end,2)];
    finalYVec = [finalYVec, statevec(end,3)];
    finalZVec = [finalZVec, statevec(end,4)];

    if i == 4
        w = statevec(:, 1);
        x = statevec(:, 2);
        y = statevec(:, 3);
        z = statevec(:, 4);

        figure(1)
        subplot(4, 1, 1)
        plot(t, w) %plotting w vs t
        title('w vs t')
        xlabel('t')
        ylabel('w')
        
        subplot(4, 1, 2)
        plot(t, x); %plotting x vs t
        title('x vs t')
        xlabel('t')
        ylabel('x')
        
        subplot(4, 1, 3)
        plot(t, y); %plotting y vs t
        title('y vs t')
        xlabel('t')
        ylabel('y')
        
        subplot(4, 1, 4)
        plot(t, z); %plotting z vs t
        title('z vs t')
        xlabel('t')
        ylabel('z')

        sgtitle('Dynamic System Solutions: Tolerance 1x10^{-8}')
        print("1a_Plot", "-dpng");
    end
end
 
%% Creating a table of differences for each tolerance
% Creating each field of the table 
differenceTaken = ["w-w_r"; "x-x_r"; "y-y_r"; "z-z_r"];
tol1 = [abs(finalWVec(1) - finalWVec(6)); abs(finalXVec(1) - finalXVec(6));...
    abs(finalYVec(1) - finalYVec(6)); abs(finalZVec(1) - finalZVec(6))];
tol2 = [abs(finalWVec(2) - finalWVec(6)); abs(finalXVec(2) - finalXVec(6));...
    abs(finalYVec(2) - finalYVec(6)); abs(finalZVec(2) - finalZVec(6))];
tol3 = [abs(finalWVec(3) - finalWVec(6)); abs(finalXVec(3) - finalXVec(6));...
    abs(finalYVec(3) - finalYVec(6)); abs(finalZVec(3) - finalZVec(6))];
tol4 = [abs(finalWVec(4) - finalWVec(6)); abs(finalXVec(4) - finalXVec(6));...
    abs(finalYVec(4) - finalYVec(6)); abs(finalZVec(4) - finalZVec(6))];
tol5 = [abs(finalWVec(5) - finalWVec(6)); abs(finalXVec(1) - finalXVec(6));...
    abs(finalYVec(5) - finalYVec(6)); abs(finalZVec(5) - finalZVec(6))];

differences = table(differenceTaken, tol1, tol2, tol3, tol4, tol5);

%% Functions
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