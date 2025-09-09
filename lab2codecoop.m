% Contributors: Cooper, Ridley, Alessandro, Gunnar
% Course Number: ASEN 3801
% File Name: MainLab2
% Last Updated: 8/26/25

clc; close all; clear all;

%% Data Cleaning
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

[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData("3801_Sec001_Test1.csv");