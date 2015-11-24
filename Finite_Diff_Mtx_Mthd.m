clear
clc; close all;%clears the command window and closes the plot window when the code is rerun
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
%% ***************
% **************** Analytical/Exact Part of calculation
% For the analytical solution
x = linspace(-a,a,nxd);
y = linspace(0,b,nyd);
temp = zeros(nyd,nxd); % used for storing summation
tol_ana = 1e-6; % Summation error criterion
%
for i = 1 : nxd
    for j = 1 : nyd
        sum = 0;
        q = exp(-pi*(b-y(j))/(2*a));
        if q == 1 %Modify to avoid singularities later
            q =1-1e-5;
        end
        p = pi/(2*a)*(x(i)+a);
        psi = atan(2*q*sin(p)/(1-q^2));
        for m = 1 : 100
            tem = sum;
            N = 2*m-1; % Sum over only odd terms
            alpha_n = N*pi/(2*a);
            sum = sum + 1/N*((sinh(alpha_n*y(j))/sinh(alpha_n*b) - q^N).*sin(N*p));
            if ( abs(sum-tem)/abs(tem) <= tol_ana )      % Check for summation convergence
                break
            end
        end
        temp(j,i) = sum +psi/2;
    end
end

Phi_analytical = 4*Vo/pi*temp; % Exact Solution
% Check if any element of error matrix is undefined, replace it by zero
for i = 1 : nxd
    for j = 1 : nyd
        if isnan(Phi_analytical(j,i))
            Phi_analytical(j,i) = 0;
        end
    end
end
%% ***************
% Iterative Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 1.5;  % Start with omega = 1.5
max_it = 300; % Maximum number of iterations
tol_it = 1e-4; % Iterative Method Convergence Criteria
Phi_iter = zeros(n,1); % Initial Guess

Vnrm2 = norm (V);
Res = V - F * Phi_iter;  % Residual 
error_norm = norm ( Res ) / Vnrm2; % Error 
%
%  Matrix Splitting for Successive Over Relaxation
% ref. Numerical Methods and Computing, Sixth ed. 
%      Kincaid and Cheney pg. 331
D = diag(diag(F));
U = -triu(F,1);
L = -tril(F,-1);
M =(D - w*L);
N = w*U + (1-w)*D;
%
for iter = 1 : max_it
    
    x_1 = Phi_iter; 
    Phi_iter = M \ (N*Phi_iter + w*V);  %  Update the approximation.
    error_norm = norm(Phi_iter - x_1)/norm(Phi_iter); %Compute the error.
    
    if ( error_norm <= tol_it ) %Check for convergence
        break
    end
    
end

% w = 1: .05: 1.9; len_w = length(w); len_iter = zeros(length(w)); for g =
% 1 : len_w [ Phi_iter, error_norm, len_iter(g), flag ] = sorr ( F,
% zeros(n,1), V, w(g), 500, 1e-4 ); end
%% ***************
% Direct Method Part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=F\V;
%% Incorporate the node position vectors and potentials to include the
%  boundary potentials in a 2-D matrix.

nn = nxd*nyd;% extra dimensions for the boundary points
Phi = zeros(nyd,nxd); %Solution matrix - In this case the first index
%is associated with the row or Y-coordinate and the
%second index with the column or x-coordinate.
Phi_direct = zeros(nyd,nxd);
Phi_iterative = zeros(nyd,nxd);
jt = 0;
for i = 1:nxd
    for j = 1:nyd
        if i >1 && i<nxd  && j == nyd
            Phi(j,i)=Vo;
            Phi_direct(j,i)=Vo;
            Phi_iterative(j,i)=Vo;
        elseif i>1 && i<nxd && j>1 && j<nyd
            jt=jt+1;
            Phi(j,i) = Phi_iter(jt);
            Phi_direct(j,i) = Z(jt);
            Phi_iterative(j,i)=Phi_iter(jt);
        end
    end
end

%% ***************
% ****************
% Calculating Error

% Zero the Phi_analytical matrix elements if <1e-10
% for i = 1 : nxd
%     for j = 1 : nyd
%         if Phi_analytical(j,i) < 1e-10
%             Phi_analytical(j,i) = 0;
%         end
%     end
% end
Res = Phi_analytical - Phi_iterative; %Define Residual
Error = abs(Res(2:nyd-1,2:nxd-1))./abs(Phi_analytical(2:nyd-1,2:nxd-1)); % Exclude Boundaries
% Check if any element of error matrix is undefined, replace 
% it by zero
% for i = 1 : nxd
%     for j = 1 : nyd
%         if isnan(Error(j,i))
%             Error(j,i) = 0;
%         end
%     end
% end
%% ***************
% ****************
% Bonus Part
% Difference between Direct and Iterative Method
Diff = Phi_direct - Phi_iterative;
Diff = abs(Diff(2:nyd-1,2:nxd-1)); % Excluding Boundaries

%% ***************
% ****************
% Plotting Routines
%Plot the functions and the contours
yd = 0:bk:b;
xd = -a:ah:a;
[xdg, ydg] = meshgrid(xd, yd);
%
%
figure('Name','Surface Plot of Iterative Solution')
S1 = surf(xdg,ydg,Phi_iterative,'EdgeColor','interp','FaceLighting','gouraud');
shading interp
grid on
% title('Numerical Solution of Trough Potential')
xlabel('x (meters)'); % X-axis label
h=get(gca,'xlabel');
set(h,'rotation',15)
ylabel('y (meters)'); % y-axis label
h=get(gca,'ylabel');
set(h,'rotation',-25)
zlabel('Amplitude (Volts)') %z-axis label
set(gcf,'Color','white'); %Set background white      
matFileName = sprintf('ECEN6637_hw3_Surf_numerical'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
%
%
figure('Name','Contour Plot of Iterative Solution')
contour(xdg,ydg,Phi_iterative, 'ShowText','on','LineWidth',1.5)
% title('Trough Equipotential Lines')
xlabel('x (meters)') % x-axis label
ylabel('y (meters)') % y-axis label
set(gcf,'Color','white'); %Set background white      
matFileName = sprintf('ECEN6637_hw3_Contour_iterative'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
%
%
figure('Name','Contour Plot of Numerical Solution')
contour(xdg,ydg,Phi_analytical, 'ShowText','on','LineWidth',1.5)
% title('Trough Equipotential Lines')
xlabel('x (meters)') % x-axis label
ylabel('y (meters)') % y-axis label
set(gcf,'Color','white'); %Set background white      
matFileName = sprintf('ECEN6637_hw3_Contour_analytical'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
%
%
figure('Name','Surface Plot of Analytical Solution')
surf(xdg,ydg,Phi_analytical,'EdgeColor','interp','FaceLighting','gouraud')
shading interp
grid on
% title('Analytical Solution of Trough Potential')
xlabel('x (meters)'); % X-axis label
h=get(gca,'xlabel');
set(h,'rotation',15)
ylabel('y (meters)'); % y-axis label
h=get(gca,'ylabel');
set(h,'rotation',-25)
zlabel('Amplitude (Volts)') %z-axis label  
set(gcf,'Color','white'); %Set background white 
matFileName = sprintf('ECEN6637_hw3_Surf_analytical'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
%
%
figure('Name','Surface Plot of Iterative Solution')
surf(xdg,ydg,Phi_iterative,'EdgeColor','interp','FaceLighting','gouraud')
shading interp
grid on
% title('Analytical Solution of Trough Potential')
xlabel('x (meters)'); % X-axis label
h=get(gca,'xlabel');
set(h,'rotation',15)
ylabel('y (meters)'); % y-axis label
h=get(gca,'ylabel');
set(h,'rotation',-25)
zlabel('Amplitude (Volts)') %z-axis label  
set(gcf,'Color','white'); %Set background white 
matFileName = sprintf('ECEN6637_hw3_Surf_iterative'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
%
%
figure('Name','Error Plot')
[~, ~] = meshgrid(xd-2, yd-2);
surf(Error,'EdgeColor','interp','FaceLighting','gouraud')
% title('Error between Analytical and Numerical Solutions discarding boundaries')
xlabel('x (meters)'); % X-axis label
h=get(gca,'xlabel');
set(h,'rotation',15)
ylabel('y (meters)'); % y-axis label
h=get(gca,'ylabel');
set(h,'rotation',-25)
zlabel('Amplitude (Volts)') %z-axis label  
set(gcf,'Color','white'); %Set background white 
shading interp
matFileName = sprintf('ECEN6637_hw3_Error'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
e = max(Error(:)); % Maximum Local Error Value
X = sprintf('Maximum value of local error is: %d',e);
disp(X) % Print on screen error value
%
%
figure('Name','Difference between Direct and Iterative Method Solutions')
[xdg, ydg] = meshgrid(xd-2, yd-2);
surf(abs(Phi_direct - Phi_iterative),'EdgeColor','interp','FaceLighting','gouraud')
% title('Error between Analytical and Numerical Solutions discarding boundaries')
xlabel('x (meters)'); % X-axis label
h=get(gca,'xlabel');
set(h,'rotation',15)
ylabel('y (meters)'); % y-axis label
h=get(gca,'ylabel');
set(h,'rotation',-25)
zlabel('Amplitude (Volts)') %z-axis label  
set(gcf,'Color','white'); %Set background white 
shading interp
matFileName = sprintf('ECEN6637_hw3_Surf_Difference'); % Create filename
saveas(gcf,[matFileName,'.eps'],'epsc') %Save figure with filename.eps
e = max(Diff(:)); % Maximum Local Error Value
X = sprintf('Maximum difference between two solutions is: %d',e);
disp(X) % Print on screen error value
%
%
