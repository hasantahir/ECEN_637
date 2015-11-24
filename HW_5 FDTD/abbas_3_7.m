% Hasan Tahir Abbas
% ECEN 637
% Homework 5: 2D TMz cylindrical wave propagation
% 10/06/2015
%
%
% This program computes and plots the TM mode fields in 2D rectagnular
% space with PEC boundaries using Yee's FDTD Algorithm. The code is based
% on the method as presented in Taflove,"Compututational Electrodynamics - The Finite
% Difference Time Domain Method", 2nd ed.,
%
%
clear;close all
%
%% Define Parameters
%
%
c = 2.99792458e8; % Speed of light
eps0 = 8.854187817e-12; % Free space permittivity
xmu0 = 4*pi*1e-7; % Free space permeability
%
%
a = 10; b = 10; % Space Dimensions
nx = 100; ny = 100; % Number of grid points in x and y direction
nx2 = ceil(nx/2); ny2 = ceil(ny/2); % Convert to integer in-case of odd number
nt = 200; % Number of time-steps
nskip = 5; % Time-steps to skip while recording
nsnap = ceil(nt/nskip);
%
% Gaussian Source
xndec = 10.0; % Mean
xn0 = 4*xndec; % Variance
%
%
%% Initialize
%
Ez = zeros( nx, ny ); %  Z E-field initialize to zero
Hx = zeros( nx, ny ); %  X H-field initialize to zero
Hy = zeros( nx, ny ); %  Y H-field initialize to zero

dx = a/(nx - 1); % Define spatial step
dy = dx; %
ds = dx;
% **************************
% **************************
dt = ds/(c*sqrt(2)); % Stability Condition
% **************************
% **************************
%
%% Medium Definition
%
% Coefficients as defined in Taflove,"Compututational Electrodynamics - The Finite
% Difference Time Domain Method", 2nd ed., pg 85.
%
% coefficients subject to the medium parameters % % Field Coefficients

dte = ones(nx,ny)*dt/(ds*eps0);
dtm = ones(nx,ny)*dt/(ds*xmu0);
Da = ones(nx,ny);
Db = dtm;
Ca = ones(nx,ny);
Cb = dte;
%
%% Yee's FDTD Algorithm
%
for n = 1:nt
    for i = 1 : nx
        for j = 1 : ny
            if (i == nx2 && j == ny2) % Place source in the center
                
                Ez(i,j) = exp(-((n-xn0)/(xndec))^2); % Gaussian Source
                
            elseif (i == 1 || j == 1 || i == nx || j == ny )
                
                Ez(i,j) = 0; % PEC boundaries at the edges
                
            else
                
                Ez(i,j) = Ez(i,j)*Ca(i,j) + Cb(i,j)*(Hy(i,j) - Hy(i-1,j)...
                    - (Hx(i,j) - Hx(i,j-1)));
                
            end
        end
    end
    
    for i = 1 : nx
        for j = 1 : ny - 1
            
            Hx(i,j) = Hx(i,j)*Da(i,j) - Db(i,j)*(Ez(i,j+1) - Ez(i,j));
            
        end
    end
    
    for i = 1 : nx - 1
        for j = 1 : ny
            
            Hy(i,j) = Hy(i,j)*Da(i,j) + Db(i,j)*(Ez(i+1,j) - Ez(i,j));
            
        end
    end
    %
    %% Plot Routine
    %
    xd = 0:ds:a;
    [xdg, ydg] = meshgrid(xd, xd);
    figure(1)
    EzSurf = surf(xdg,ydg,Ez,'EdgeColor','interp','FaceLighting','gouraud');
    shading interp % Avoid jittered shading
    set(gcf,'Color','white'); % Set background color to white
    set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
    box on %
    colorbar % Show colorbar on the right of the figure
    grid on
    colormap('jet')
    material dull % Set reflectivity of the surface to dull
    % Add three lights below tom improve visuals
    h = light;
    light('Position',[10 0 1]);
    light('Position',[-190 -120 1]);
    light('Position',[0 0 1]);
    title(['E-field at Time Step ',int2str(n)])
    matFileName = sprintf('Abbas_3_7_E_Field_%d', n); % Create Title with mode numbers
    
    if rem(n,nskip) == 0
        
        M_PEC(:,n/nskip) = getframe; % Capture frames to create a movie
        
    end
    
        if n == 20 || n == 60 || n == 80 || n == 100 ||...
                n == 120 || n == 140 || n == 160
            saveas(gcf,[matFileName,'.eps'],'epsc') % Save visualizations
        end
    
end
movie(M_PEC,1); % Play movie