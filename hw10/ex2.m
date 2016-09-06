clear;close all

M = 51;
f = 300e6;
c = 2.99792458e8; % Speed of light
xmu = 4*pi*1e-7;  % Permeability of free space
eps0 = 8.854187817e-12; % Permittivity of free space
omega = 2.0*pi*f;
K = omega*sqrt(xmu*eps0);
lambda = c/f;
a = 7.022e-3*lambda;
L = linspace(.5*lambda/pi,3.5*lambda/pi,101);
for i = 1:length(L)
        Y = Pocklington(M, L(i), a, f);
         Y_mine(i) = Y;
end  


figure; plot(pi*L/lambda, real(Y_mine));
hold on; plot(pi*L/lambda, imag(Y_mine));
xlim([1 3.5])

NumTheta = 361;
[E_rad, I] = Pocklington_pattern(M, 4*lambda/5, a, f, NumTheta);
Theta = linspace(0.0, 2.0*pi, NumTheta);
figure
polar(Theta, abs(E_rad));
% figure
% plot(abs(I))

NumTheta = 361;
[E_rad, I] = Pocklington_pattern(M, 5*lambda/5, a, f, NumTheta);
Theta = linspace(0.0, 2.0*pi, NumTheta);
figure
polar(Theta, abs(E_rad));
% figure
% plot(abs(I))

NumTheta = 361;
[E_rad, I] = Pocklington_pattern(M, .7*lambda, .001588*lambda, f, NumTheta);
Theta = linspace(0.0, 2.0*pi, NumTheta);
figure
polar(Theta, abs(E_rad));
figure
x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
plot(x, real(I),x, imag(I));
ax = gca;
ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
axis([ -.3498 .3498 -2.5e-3 2.5e-3])
hold on
xlabel('z')
ylabel('A')