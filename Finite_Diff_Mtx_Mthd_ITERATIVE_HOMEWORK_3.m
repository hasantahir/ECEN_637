clc; close all; %clears the command window and closes the plot window when the code is rerun
% This is a program to solve for the potential in a square cylinder
%%
% Essential coefficients and constants
Vo = 100.; % Voltage potential
a = 2.; % 1/2 Length in the x-direction
b = 3.; % Total length in the y-direction
nxd = 10; % number of data points in the x-direction
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
Phi = zeros(nxd,nyd);
R = zeros(nxd,nyd);
Error = zeros(nxd,nyd);
%%% Set Parameters for the iteration
k_max = 100;
err = 1e-4;
E = 1;
n_iter = 0;
w = 1.5;

%% ***************
% Initial Guess of the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nxd
    for j = 1 : nyd
        Phi(j,i) = Vo/2;
    end
end
%% 
%Fill the F matrix and V vector
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
% null out some off diagonal elements that have been set = alpa2 above
for k = ny:ny:n-ny 
    F(k,k+1) = 0.;
    F(k+1,k) = 0.;
end
% Fill the V vector
for k = ny:ny:n
    V(k,1) = alpa2*Vo;
end
%% 
while E >= err
    for k = 1 : k_max
        for i = 2 : nxd-1
            for j = 2 : nyd-1
                R(j,i) = 1/(4*bta2)*(Phi(j-1,i) + Phi(j+1,i)...
                    + alpa2*Phi(j,i-1) + alpa2*Phi(j,i+1)) - Phi(j,i);
                Phi(j,i) = Phi(j,i) + w*R(j,i);
                Err(j,i) = abs(R(j,i))/abs(Phi(j,i));
            end
        end
        n_iter = n_iter + 1;
        E = max(Err(:));
    end
end

           