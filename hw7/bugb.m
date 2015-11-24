% Plot a scheme's RGB values:
rgbplot(brewermap(9,'Blues')) % standard
rgbplot(brewermap(9,'*Blues')) % reversed

% View information about a colorscheme:
[~,num,typ] = brewermap(0,'Paired')
num = 12
typ = 'Qualitative'

% Multiline plot using matrices:
N = 6;
axes('ColorOrder',brewermap(N,'Pastel2'),'NextPlot','replacechildren')
X = linspace(0,pi*3,1000);
Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X.', 1:N);
plot(X,Y, 'linewidth',4)

% Multiline plot in a loop:
N = 6;
set(0,'DefaultAxesColorOrder',brewermap(N,'Accent'))
X = linspace(0,pi*3,1000);
Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X.', 1:N);
for n = 1:N
plot(X(:),Y(:,n), 'linewidth',4);
hold all
end

% New colors for the "colormap" example:
load spine
image(X)
colormap(brewermap([],'*YlGnBu'))

% New colors for the "surf" example:
[X,Y,Z] = peaks(30);
surfc(X,Y,Z)
colormap(brewermap([],'RdYlGn'))
axis([-3,3,-3,3,-10,5])

% New colors for the "contourcmap" example:
brewermap('*PuOr'); % preselect the colorscheme.
load topo
load coast
figure
worldmap(topo, topolegend)
contourfm(topo, topolegend);
contourcmap('brewermap', 'Colorbar','on', 'Location','horizontal',...
'TitleString','Contour Intervals in Meters');
plotm(lat, long, 'k')