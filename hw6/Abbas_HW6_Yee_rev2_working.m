function Abbas_HW6_Yee()
% Hasan Tahir Abbas
% ECEN 637
% Homework 6: Simulating Yee's Algoirthm
% 10/21/2015
%
%
% This program generates the results presented in the paper:
% "Numerical Solution of Initial Boundary Value Problems involving
% Maxwell's Equations in Isotropic Media", Kane S. Yee (1966)
%
%

clear;close all

% *************************
%% % % % % % % % Parameters
% *************************
global c xmu eps0 asize
global nx  ny  nt mxst mxnd myst mynd
global dt ds
global Ez  Hx  Hy; % Create E and H field components.
% global data type enables the global scope of the variables
% within the code. Unlike global variables they can not be accessed
% outside the code
global mediaEz mediaHx mediaHy; %
global Ca Cb Da Db ; % Define material based coefficients
global Ezs % The excitation signal




c = 2.99792458e8;
xmu = 4*pi*1e-7;
eps0 = 8.854187817e-12;
asize = 5; % Space Dimension in meters
nx = 80;      %% Number of cells in x-direction
ny = 100;     %% Number of cells in y-direction
nt = 100;     %% Number of time steps
mxst = 17;    %% Start of PEC section in x-direction
mxnd = 49;    %% End of PEC section in x-direction
myst = 33;    %% Start of PEC section in y-direction
mynd = 65;    %% End of PEC section in y-direction
strip = 30;   %% y-Position where E-field is recorded
% % % % % % % Initialize




Ez = zeros(nx,ny); %% z-component of E-field
Hx = zeros(nx,ny); %% x-component of H-field
Hy = zeros(nx,ny); %% y-component of H-field

mediaEz = ones(nx,ny); %% z-component of E-field
mediaHx = ones(nx,ny); %% x-component of H-field
mediaHy = ones(nx,ny); %% x-component of H-field

Ca = zeros(2,1); %% x-component of H-field
Cb = zeros(2,1); %% x-component of H-field
Da = zeros(2,1); %% x-component of H-field
Db = zeros(2,1); %% x-component of H-field

% Initialize arrays to create to capture E-field
% line plots at various time instants.
Ez_5 = zeros(1,nx);
Ez_35 = zeros(1,nx);
Ez_65 = zeros(1,nx);
Ez_95 = zeros(1,nx);



ds = asize/(mxnd - mxst - 1); %% Length Increment
% *************************
%% % % % % % % % Stability Condition
% *************************
dt = ds/(c*sqrt(2)); % Stability Condition
% *************************
% *************************


%
iflaga = 2;   %% 1 if Free space; 2 if PEC
%
%
% *************************
% *************************
%% Main Program
% *************************
% *************************
    
    define_media(iflaga); % Structre Definition
    define_coefficients(); % Create Coefficients for equations
    for n = 1:nt
        adv_Ez(n);
        adv_H(n);
            if rem(n,5) == 0
                figure(3)
                title('3D plot of $E_z$','Interpreter','latex')
                my_surface_plot(Ez, iflaga);
                
%                 figure(4)
%                 my_surface_plot(Hx, iflaga);
%                 figure(3)
%                 my_surface_plot(Hy, iflaga);
            end   % end_if
        if n == 5
            Ez_5 = Ez(:,strip); % Record the
        elseif n == 35
            Ez_35 = Ez(:,strip);
        elseif n == 65
            Ez_65 = Ez(:,strip);
        elseif n == 95
            Ez_95 = Ez(:,strip);
        end
        
    end
    my_line_plot(Ez_5, Ez_35, Ez_65, Ez_95,iflaga)
%     iflaga = iflaga+1;








end
% *************************
% *************************
%% Create Structure in the computational space
% *************************
% *************************
function define_media(iflaga)

global nx  ny mxst mxnd myst mynd;
global mediaEz mediaHx mediaHy; %
if (iflaga == 2)
    
    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd)
                if ( j >= myst && j <= mynd)
                    mediaEz(i,j) = 2;
                end
            end
        end
    end
    
    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd)
                if ( j >= myst && j <= mynd-1)
                    mediaHx(i,j) = 2;
                end
            end
        end
    end
    
    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd-1)
                if ( j >= myst && j <= mynd)
                    mediaHy(i,j) = 2;
                end
            end
        end
    end
end

end


% *************************
% *************************
%% Create Coefficients for the equations
% *************************
% *************************
function define_coefficients()

global Ca Cb Da Db ; % Define material based coefficients
global xmu eps0 dt ds
% % % % % % % % Field Coefficients
dte = dt/(ds*eps0);
dtm = dt/(ds*xmu);
Da(1) = 1;
Db(1) = dtm;
Ca(1) = 1;
Cb(1) = dte;
Da(2) = 0;
Db(2) = 0;
Ca(2) = 0;
Cb(2) = 0;
end

% *************************
% *************************
%% Create the excitation signal
% *************************
% *************************
function Ezs = Source(n)

global Ezs
% Creates a half-sinusoidal source between the time increments
% 1 and 10.%
%
%% For Gaussian Source
%     xndec = 10.0;
%     xn0 = 4*xndec;
%     Ezs = exp(-((n-xn0)/(xndec))^2);
    %% For Sinusoidal Source
    
    if ( n >=1 && n <= 10)
        Ezs = sin(n*pi/10);
    end
end


% *************************
% *************************
%% Algorithm for Computing E-field
% *************************
% *************************
function adv_Ez(n)
% Compute z-component of E-field
global Ez  Hx  Hy
global mediaEz
global Ca Cb
global nx ny

for i = 1 : nx
    for j = 1 : ny
        m  = mediaEz(i,j);
        if (i == 1)   %% Incident Field Source Excitation
            Ez(i,j) = Source(n);
        elseif (i >= 2 && j >=2)
            Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1)));
        elseif (j == 1) %% Field at the bottom edge of the boundary
            Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
                - Hx(i,j));
        end
    end
end
end


% *************************
% *************************
%% Algorithm for Computing E-field
% *************************
% *************************
function adv_H(n)
% Compute z-component of E-field
global Ez  Hx  Hy
global mediaHx mediaHy
global Da Db
global nx ny

% % %     Compute x-component of H-field

for i = 1 : nx
    for j = 1 : ny - 1
        m = mediaHx(i,j);
        Hx(i,j) = Hx(i,j)*Da(m) - Db(m)*(Ez(i,j+1) - Ez(i,j));
    end
end

% % %     Compute y-component of H-field

for i = 1 : nx - 1
    for j = 1 : ny
        m = mediaHy(i,j);
        Hy(i,j) = Hy(i,j)*Da(m) + Db(m)*(Ez(i+1,j) - Ez(i,j));
    end
end
end


% *************************
% *************************
%% Plot routine for the 3D surface plots
% *************************
% *************************
function my_surface_plot(field, iflaga)
% This function generates the 3D surface plots for the field in
% the argument of the function
%

colormap('parula')
global mxst mxnd myst mynd n
EzSurf = surf((field),'EdgeColor','interp','FaceLighting','gouraud');     shading interp % Avoid jittered shading
% rectangle('Position',[myst,mxst,mynd-myst,mxnd-mxst],'FaceColor','r')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
box on %
grid on
material dull; % Set reflectivity of the surface to dull
% Add three lights below tom improve visuals
h = light;
light('Position',[10 0 1]);
light('Position',[-190 -120 1]);
light('Position',[0 0 1]);
axis tight
 view([-165 45]);
% view([ 90 90 ])
M(:,n) = getframe(gcf) ;
end

% *************************
% *************************
%% Plot routine for the 2D Line Plot of Ez
% *************************
% *************************
function my_line_plot(Ez_5, Ez_35, Ez_65, Ez_95,iflaga)
% Plots the line plots for E-field

figure('Name', sprintf('Line Plots for E-field_%d',iflaga))
set(gcf,'Color','white');
title('E-field at different Time Steps');
%%%
subplot(4,1,1)
plot(Ez_5,'Color','black','LineWidth',1.2)
a = gca;
a.YMinorGrid = 'on'; % Display Minor Grid lines along the x-axis
a.XTickLabel = {}; % Hide x-labels
a.FontName = 'times new roman';
legend('n=5')
legend('boxoff')
%%%
subplot(4,1,2)
plot(Ez_35,'Color','black','LineWidth',1.2)
a = gca;
a.YMinorGrid = 'on';
a.XTickLabel = {};
a.FontName = 'times new roman';
legend('n=35')
legend('boxoff')
%%%
subplot(4,1,3)
plot(Ez_65,'Color','black','LineWidth',1.2)
a = gca;
a.YMinorGrid = 'on';
a.XTickLabel = {}; % Hide x-labels
a.FontName = 'times new roman';
legend('n=65')
legend('boxoff')
%%%
subplot(4,1,4)
plot(Ez_95,'Color','black','LineWidth',1.2)
a = gca; % Set current axis handle
a.YMinorGrid = 'on';
a.FontName = 'times new roman';
legend('n=95')
legend('boxoff')
xlabel('Cells in x-direction','FontSize',11,...
    'Interpreter','latex','FontName','times new roman')
matlab2tikz('filename',sprintf('ECEN637_HW6_a_%d.tex',iflaga))
saveas(gcf,[sprintf('ECEN637_HW6_a_%d',iflaga'),'.eps'],'epsc') % Save visualizations

end


