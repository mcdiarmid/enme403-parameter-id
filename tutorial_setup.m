% 
% Clean up
clc
close all
clear

% Data
T_data=[25.34,24.28,24.26,25.24,24.82,30.01,28.36,27.00,25.89,25.18, ...
  25.60,25.46,24.78,37.61,32.00,28.41,27.68,27.02,25.03,26.53,25.22, ...
  33.69,30.66,28.08,27.51,26.90,25.97,25.45,24.42,24.79,25.11,24.06, ...
  24.51,25.10,24.40,24.62,25.16,25.10,25.06,26.22,25.34];
fs = 4;
dt_data = 1/fs;
t_start = 0; t_end = 10;
t_meas = t_start:dt_data:t_end;
t_impulse = [1.05, 3.05, 5.05];

% Other constants
parent = [2, 200e3, 400e3, 300e3];
T_amb = 25;
V = 5; % Litres
Cp = 4181; % Specific heat capacity of water
[~, rows] = size(T_data);
cols = 4;

% Matrix initialization
A = zeros(rows, cols);
b = zeros(rows, 1);
