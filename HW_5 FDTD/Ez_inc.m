function Ezs = Ez_inc(n)
% % % Gaussian Source Excitation at the Center of the grid
% % % Format Ezs = Ez_inc(n,f,dt)
% % % where n is the time instant, 
% % % f is the center frequency and 
% % % dt is the time step 
xndec = 10.0;
% w = 40; 
xn0 = 4*xndec;
Ezs = exp(-((n-xn0)/(xndec))^2);
% Ezs = exp( -log(.001) * ( (n - xn0/dx) / w )^2 );
%    Ezs = (exp(-(((n-xn0)/(2*xndec)))^2)).*sin(2*pi*f*(n - xn0)*dt);
%  Ezs =1;
% Ezs = sin(2*pi*f*(n - xn0)*dt);
 if ( n <= xn0)
%         Ezs = (exp(-((n-xn0)/(1.1*xndec))^2)).*sin(2*pi*f*(n - xn0)*dt);
% else
%     Ezs = sin(2*pi*f*(n - xn0)*dt);
% end

end