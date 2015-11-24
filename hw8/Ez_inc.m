function Ezs = Ez_inc(n)
% % % Gaussian Source Excitation at the Center of the grid
xndec = 10.0;
xn0 = 4.0*xndec;
Ezs = exp(-(((n-xn0)/xndec))^2);
%  Ezs =1;
% if ( n>= 1 && n <= 10)
%     Ezs = sin(n*pi/10);
% else
%     Ezs = 0;
end