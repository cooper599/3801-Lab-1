% Contributors: Cooper, Ridley, Alessandro, Gunnar
% Course Number: ASEN 3801
% File Name: MainLab2
% Last Updated: 8/26/25

clc; close all; clear all;

%% Functions
% Function 1
function [t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData(filename)
    %{
    Inputs:
        filename - csv file with aspen position data
    Outputs:
        t_vec - 1xn vector in seconds
        av_pos_inert - 3xn matrix of position vectors of the aerospace vehicle in
        frame E, meters
        av_att - 3xn matrix of vectors of 3-2-1 Euler angles describing the
        attitude of the aerospace vehicle relative to frame E
        (roll,pitch,yaw), radians
        tar_pos_inert - 3xn position vector of the target in frame E, m
        tar_att - 3xn matrix of vectors of the 3-2-1 Euler angles describing
        the attitude of the target relative to frame E (roll,pitch,yaw),
        radians
    %}
    rawdata = readmatrix(filename);
    % Double check if should be divided by 100
    pos_av_aspen = (rawdata(4:end,12:14)./1000)'; % 3xn matrix of position vectors of aerospace vehicle in frame N
    att_av_aspen = rawdata(4:end,9:11)'; % 3xn matrix of vectors of helical angles describing the attitude of the aerospace vehicle in frame N
    pos_tar_aspen = (rawdata(4:end,6:8)./1000)'; % 3xn matrix of position vectors of target in frame N
    att_tar_aspen = rawdata(4:end,3:5)'; % 3xn matrix of vector of helical angles describing attitude of the target in frame N
    [av_pos_inert, av_att, tar_pos_inert, tar_att] = ConvertASPENData(pos_av_aspen,att_av_aspen,pos_tar_aspen,att_tar_aspen);
    t_vec = zeros(1,length(rawdata)-6); % Preallocating Array
    for i = 1:length(rawdata) - 3
        t_vec(i) = (1/100) * i - (1/100); % Multiply frames by (1/100 Hz) subtracting first to normalize to 0 seconds initially
    end
end

% Function 2
function [DCM] = RotationMatrix321(attitude321)
    mat1 = [1 0 0; 0 cos(attitude321(1)) sin(attitude321(1)); 0 -sin(attitude321(1)) cos(attitude321(1))];
    mat2 = [cos(attitude321(2)) 0 -sin(attitude321(2)); 0 1 0; sin(attitude321(2)) 0 cos(attitude321(2))];
    mat3 = [cos(attitude321(3)) sin(attitude321(3)) 0; -sin(attitude321(3)) cos(attitude321(3)) 0; 0 0 1];
    DCM = mat1 * mat2 * mat3;
end

% Function 3
function DCM = RotationMatrix313(attitude313)
    %{
    Inputs:
        attitude313 - 3x1 vector with the 3-1-3 Euler Angles in the form
        attitude313 = [alpha,beta,gamma]^T
    Outputs:
        DCM - rotation matrix calculated from the Euler Angles
    %}
    alpha = attitude313(1); beta = attitude313(2); gamma = attitude313(3);
    R3alpha = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
    R1beta = [1 0 0; 0 cos(beta) sin(beta); 0 -sin(beta) cos(beta)];
    R3gamma = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
    DCM = R3alpha * R1beta * R3gamma;
end

% Function 4
function [attitude321] = EulerAngles321(DCM)
    alpha = atan2(-DCM(3,1), sqrt((DCM(1, 1))^2)+ DCM(2, 1)^2);
    beta = atan2(DCM(2,1), DCM(1,1));
    gamma = atan2(DCM(3,2), DCM(3,3));
    attitude321 = [alpha, beta, gamma]';
end

% Function 5
function attitude313 = EulerAngles313(DCM)
    %{
    Inputs: 
        DCM - rotation matrix
    Outputs:
        attitude313 - 3x1 vector with 3-1-3 Euler angles in form attitude
        313 = [alpha,beta,gamma]^T
    %}
    alpha = atan2(DCM(1,3)/DCM(2,3));
    beta = acos(DCM(3,3));
    gamma = atan2(DCM(3,1)/-DCM(3,2));
    attitude313 = [alpha beta gamma]';
end

%% Plotting
% Function 1 Call
[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData("3801_Sec001_Test1.csv");