% Hasan Tahir Abbas
% ECEN 637 Homework 4: Modes of a waveguide
% 10/01/2015
%
clear;
close all; %clears the command window and closes the plot window when the code is rerun
 colorbar
%
% This is a program computes and plots first seven TM modes of a
% rectangular wavguide. The mode propagation constants obtained numerically
% are compared with the analytical values.
%%
% Essential coefficients and constants
a = 2.; % 1/2 Length in the x-direction
b = 3.; % Total length in the y-direction
nxd = 40; % number of data points in the x-direction
nyd = 35; % number of data points in the y-direction
nx=nxd-2; % computed data points - without boundary points
ny=nyd-2;
n = nx*ny; % number of node points
ah = 2*a/(nxd-1); % distance between node points in the x-direction
bk = b/(nyd-1); % distance between node points in the y-direction
alpa = 1; % alpha matrix variable
alpa2 = alpa^2;
bta2 =(1.+alpa2)/2.; % beta matrix variable
A = zeros(n,n); % zero the matrix
mod_num = 7; % number of modes to calculate
% tick_array = zeros(mod_num,1);
%
%
%% Analytical Calculation of the modes
a_x = 4; % Actual Dimension in x
b_y = 3; % Actual Dimension in y
M = 3; % Calculate first seven modes
N = 3; % Calculate first seven modes
KK = zeros(M,N); % Mode matrix
row = zeros(M*N,1); % Row order of modes
col = zeros(M*N,1); % Column order of modes
for i = 1 : M
    for j = 1 : N
        KK(i,j) = sqrt( (i*pi/a_x)^2 + (j*pi/b_y)^2 ); % Definition of TM modes
    end
end
% Sort the modes in ascending order and make a vector out of it
%
sorted = reshape(sort(KK(:)),[],1);
%
%
for i = 1 : length(sorted)
    [row(i),col(i)] = find(KK==sorted(i)); % Find the indices of modes
end
Mode_order = [row,col]; % TM mode order
%%
%Fill the A matrix
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = 4.*bta2;
        elseif j == i-1 || j == i+1
            A(i,j) = -alpa2;
        elseif j == i-ny || j == i+ny
            A(i,j) = -1.;
        end
    end
end
% null out some off diagonal elements that have been set = alpa2 above
for k = ny:ny:n-ny
    A(k,k+1) = 0.;
    A(k+1,k) = 0.;
end
%% Solve the eigen-value problem to obtain modes and propagation constants
[Modes, lambda] = eig(A,'vector');
kc = sqrt(lambda)./ah;
for it= 1 : mod_num
    phi_1 = Modes(:,it); % Each column is a mode
    Phi = zeros(nyd,nxd);
    
    %
    %% Reshape and place the mode in to the waveguide
    %
    jt = 0;
    for i = 1:nxd
        for j = 1:nyd
            if i >1 && i<nxd  && j == nyd
                Phi(j,i)=0;
            elseif i>1 && i<nxd && j>1 && j<nyd
                jt=jt+1;
                Phi(j,i) = phi_1(jt); % Place values inside the rectangle
            end
        end
    end
    %
    %% Plotting routine
    figure(it)
    axis tight
    yd = 0:bk:b;
    xd = -a:ah:a;
    [xdg, ydg] = meshgrid(xd, yd);
    contour(xdg,ydg,Phi,'LineWidth',1.5)
    hold on
    [U,V] = gradient(Phi);
    colormap hsv
    %
    % Quiver Plot
    q = quiver(xdg,ydg,U,V,'LineWidth',1.1,'AutoScaleFactor',1);
    c = q.Color;
    q.Color = 'black';
    
    matFileName = sprintf('Contour plot of TM_{%d%d} mode', row(it) , col(it)); % Create Title with mode numbers
    matFileName1 = sprintf('Contour_plot_of_TM_%d%d_mode', row(it) , col(it)); % Create Title with mode numbers
    title(matFileName,...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',10,...
        'Interpreter','tex');
    tick_array{it} = strcat(num2str(row(it)),num2str(col(it))); % Tick marks for the error plot x-axis
    xlabel('x (meters)','FontWeight','bold') % x-axis label
    ylabel('y (meters)','FontWeight','bold') % y-axis label
     legend('H-field', 'E-field','Location','southoutside','Orientation','horizontal'); % legend
     colorbar
    saveas(gcf,[matFileName1,'.eps'],'epsc')
end
%
% Plot of Error in Analytical and Computed Propagation Constants
figure()
Error = (abs(kc(1:mod_num) - sorted(1:mod_num)))./abs(sorted(1:mod_num));
bar(Error)
ax = gca;
xlabel('Mode Number','HorizontalAlignment','center','FontWeight','bold') % x-axis label
ax.XTickLabel = tick_array;
ylabel('Error',... % y-axis label
    'HorizontalAlignment','center','FontWeight','bold',...
    'Interpreter','latex')
title('Error between Analytical and Computed values of Modes',... % title label
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','tex');
Error_plot = sprintf('Error_plot'); % Create Title with mode numbers
saveas(gcf,[Error_plot,'.eps'],'epsc')