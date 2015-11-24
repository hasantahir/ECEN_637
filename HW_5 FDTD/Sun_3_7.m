function [] = Sun_3_7()
clc;
close all;
%% essential constant
global a nx ny nt Hx Hy Ez T
global dx dy ds dt
global Ca Cb Da Db e u c
a = 10; % lenth and width of the numerical space
nx = 100; % number of grid points in the x and y directions
ny = 100;
nt = 100; %number of time steps
e = 8.85e-12; % free space permittivity
u = (4*pi)*10^(-7); %free space permeability
c = 3e8; %speend o light
Hx = zeros(100,100);% 3 main matrix
Hy = zeros(100,100);
Ez = zeros(100,100);
dx = a/(nx-1);
dy = dx ;
ds = dx ;
dt = 1.000*ds/(c*sqrt(2));
%dt = 1.0005*ds/(c*sqrt(2));
Ca = 1;% 4 coefficient
Da = 1;
Cb = dt/e;
Db = dt/u;
%% main function
for T = 1:nt %main time step 't'
    adv_E_field; % compute Ez field
    adv_H_field; % compute Hx & Hy field
    [x,y] = meshgrid(1:100,1:100);
    surf(x,y,Ez);
    F(T) = getframe(gcf); %plot continuous picture
    if T==1 || T==20 || T==30 || T==40 || T==50|| T==60|| T==70|| T==80|| T==90 %capture
        % 10 plots
        figure
        surf(x,y,Ez);
    end
end
end %the main function end
%% subroutines

function [] = adv_E_field ()
% compute the Ez fields by using Hx and Hy(1 time steps before) and Ez(2 time steps
% before)
global nx ny Ez Hx Hy ds
global Ca Cb
for j = 2:ny-1
    for i = 2:nx-1
        if i == 50 && j == 50 % excitation in middle
            Ez(50,50) = Ez_incident;
        else
            Ez(j,i) = Ca*Ez(j,i) + Cb*( Hy(j,i)-Hy(j,i-1)-(Hx(j,i)-Hx(j-1,i)) )/ds ;
        end
    end
end % big for end
end % function end
function [l1] = Ez_incident ()
% input time, returen incident Ez
% l1 means local var 1
% the excitation of Ez field
% gaussian pulse
global T
beta = 10;
l1 = exp( -( (T-4*beta)/beta )^2 ); %guassian function
end
function [] = adv_H_field ()
%compute Hx and Hy by useing Ez(1 time step before) and Hx Hy themselves(2 time steps
% before)
global nx ny Da Db Hx Hy Ez ds
for j = 1:ny-1
    for i = 1:nx
        Hx(j,i) = Da*Hx(j,i) - Db*( Ez(j+1,i)-Ez(j,i) )/ds;
    end
end
for j = 1:ny
    for i = 1:nx-1
        Hy(j,i) = Da*Hy(j,i) + Db*( Ez(j,i+1)-Ez(j,i) )/ds;
    end
end
end