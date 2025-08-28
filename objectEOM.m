function [xDot] = objectEOM(t, x, rho, Cd, A, m, g, vWind)
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