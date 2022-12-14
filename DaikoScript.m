%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/kierancurrie-cathey/Documents/University/Matlab/MATLAB/APD/DAIKO_trajectory_v2.csv
%
% Auto-generated by MATLAB on 24-Nov-2022 12:40:59

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["TimeElapseds", "Latdeg", "Londeg", "Geoaltitudem", "ImpactPressurePa", "MeasStaticTempK"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
DAIKOtrajectoryv2 = readtable("/Users/kierancurrie-cathey/Documents/University/Matlab/MATLAB/APD/DAIKO_trajectory_v2.csv", opts);


%% Clear temporary variables
clear opts