function [attitude321] = EulerAngles321(DCM)
alpha = atan2(-DCM(3,1), sqrt((DCM(1, 1))^2)+ DCM(2, 1)^2);
beta = atan2(DCM(2,1), DCM(1,1));
gamma = atan2(DCM(3,2), DCM(3,3));
attitude321 = [alpha, beta, gamma]';