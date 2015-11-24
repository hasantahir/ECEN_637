% just for fun 
% Summation
clear all;close all
a = 2;
b = 3;
V = 100;
nx = 10;
ny = 15;
x = linspace(-a,a,nx);
y = linspace(0,b,ny);
temp = zeros(length(y),length(x));
for i = 1 : nx
    for j = 1 : ny
        sum = 0;
        q = exp(-pi*(b-y(j))/(2*a));
        if q == 1
            q =1-1e-25;
        end
        p = pi/(2*a)*(x(i)+a);
        psi = atan(2*q*sin(p)/(1-q^2));
        for m = 1 : 15
            n = 2*m-1;
            alpha_n = n*pi/(2*a);
            sum = sum + 1/n*((sinh(alpha_n*y(j))/sinh(alpha_n*b) - q.^n).*sin(n*p));
        end
        temp(j,i) = sum +psi/2;
    end
end

Phi = 4*V/pi*temp;

[x,y] = meshgrid(x, y);
surf(x,y,Phi)
figure(1)
shading interp
figure
contour(x,y,Phi, 'ShowText','on')
title('Trough Equipotential Lines')
xlabel('x (meters)') % x-axis label
ylabel('y (meters)') % y-axis label

