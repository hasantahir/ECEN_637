clear;
close all
% figure('Name','2','Position',[100 100 850 600])
% % % % % % % % % % PARAMETERS

c = 2.99792458e8; %% Speed of light
eps0 = 8.854187817e-12; %% Free space permittivity
epsr = 4; %% Relative Permitivitty of medium
xmu0 = 4*pi*1e-7; %% Free space permeability
a = 10; b = 10; %% Space Dimensions
nx = 100; ny = 100; 
nx2 = ceil(nx/2); ny2 = ceil(ny/2); 
nt = 100; nskip = 5;
% % % % % % % % % % 

% % % % % % % % % %  INITIALIZE
Ez = zeros( nx, ny ); %%  Z E-field initialize to zero
Hx = zeros( nx, ny ); %%  X H-field initialize to zero
Hy = zeros( nx, ny ); %%  Y H-field initialize to zero

dx = a/(nx - 1); %% 
dy = dx;
ds = dx;
dt = ds/(c*sqrt(2));
nsnap = ceil(nt/nskip);

% % % % % % % % % % MEDIUM DEFINITION
% % % This function sets up the field coefficients subject to the medium
% parameters
% % % Field Coefficients
dte = ones(nx,ny)*dt/(ds*eps0);
dtm = ones(nx,ny)*dt/(ds*xmu0);
Da = ones(nx,ny);
Db = dtm;
Ca = ones(nx,ny);
Cb = dte;

% % % % % % % % % MAIN PROGRAM
for n = 1:nt
    for i = 2 : nx-1
        for j = 2 : ny-1
             if (i == nx2 && j == ny2)
                Ez(i,j) = Ez_inc(n);
             elseif (i == 2 || j == 2 || i == nx-1 || j == ny-1 )
                 Ez(i,j) = 0;
                 
             else
                 Ez(i,j) = Ez(i,j)*Ca(i,j) + Cb(i,j)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1)));
                 
             end
        end
    end
    
    for i = 2 : nx - 1
        for j = 2 : ny - 1
                      
             Hx(i,j) = Hx(i,j)*Da(i,j) - Db(i,j)*(Ez(i,j+1) - Ez(i,j));

             Hy(i,j) = Hy(i,j)*Da(i,j) + Db(i,j)*(Ez(i+1,j) - Ez(i,j));
        end
    end

% % % % % % % % Plot Routune

    V = interp2(Ez,1);
    hSurf = surf(abs(V),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    set(gcf,'Color','white');
    set (gca,'FontName','times new roman')
%     set(hSurf,'Color','black','LineWidth',1.4)
%     xlabel(' f (THz)','FontSize',11,'Interpreter','latex','FontName','times new roman')
%     ylabel('  $log_{10}(\varepsilon''''/\varepsilon'')$ ','FontSize',11,'Interpreter','latex','FontName','times new roman')
%    axis([400 780 -2 2.2])
   box off
   axis tight
    
%      view([0,0,1])
%     axis image
    colorbar
    grid off
%     colormap('hot')
    
    M_PEC(:,n) = getframe ;
    
end
movie(M_PEC,1,30);



