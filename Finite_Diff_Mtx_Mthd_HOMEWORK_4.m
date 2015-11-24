clc; close all; %clears the command window and closes the plot window when the code is rerun
% This is a program to solve for the potential in a square cylinder
%%
% Essential coefficients and constants
Vo = 100.; % Voltage potential
a = 2.; % 1/2 Length in the x-direction
b = 3.; % Total length in the y-direction
nxd = 20; % number of data points in the x-direction
nyd = 15; % number of data points in the y-direction
nx=nxd-2; % computed data points - without boundary points
ny=nyd-2;
n = nx*ny; % number of node points
ah = 2*a/(nxd-1); % distance between node points in the x-direction
bk = b/(nyd-1); % distance between node points in the y-direction
alpa = ah/bk; % alpha matrix variable
alpa2 = alpa^2;
bta2 =(1.+alpa2)/2.; % beta matrix variable
V = zeros(n,1); %zero the excitation vector
F = zeros(n,n); % zero the matrix
for i = 1:n
    for j = 1:n
        if i == j
            F(i,j) = 4.*bta2;
        elseif j == i-1 || j == i+1
            F(i,j) = -alpa2;
        elseif j == i-ny || j == i+ny
            F(i,j) = -1.;
        end
    end
end