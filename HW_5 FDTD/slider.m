% Create a common meshgrid for two plots having different z values
[x,y]=meshgrid(-1:0.05:1, -1:0.05:1);
% Two different z-data plots for the same x and y data points
z1=sin(x.^2+y.^2);
z2=cos(x+y);
% Minimum and maximum values of the two plots
% This is useful in setting the limits of the colorbar
bottom = min(min(min(z1)),min(min(z2)));
top  = max(max(max(z1)),max(max(z2)));
% Plotting the first plot
subplot(1,3,1)
h1=surf(x,y,z1);
shading interp;
% This sets the limits of the colorbar to manual for the first plot
caxis manual
caxis([bottom top]);
% Plotting the second plot
subplot(1,3,2)
h2=surf(x,y,z2);
shading interp;
% This sets the limits of the colorbar to manual for the second plot
caxis manual
caxis([bottom top]);
colorbar;