clear all
close all
% writerObj = VideoWriter('Hy_ABC_n_5_yee_dimension.avi');
% open(writerObj);
% figure('Name', 'Field Plot','Position',[0 0 1920 1080])
% % % % % % % Parameters
c = 2.99792458e8;
xmu = 4*pi*1e-7;
eps0 = 8.854187817e-12;
asize = 5;
% nx = 200;      %% Number of cells in x-direction
% ny = 200;     %% Number of cells in y-direction
% nt = 400;     %% Number of time steps
% mxst = 85;    %% Start of PEC section in x-direction
% mxnd = 117;    %% End of PEC section in x-direction
% myst = 85;    %% Start of PEC section in y-direction
% mynd = 117;    %% End of PEC section in y-direction

% Yee's Dimensions
nx = 80;      %% Number of cells in x-direction
ny = 100;     %% Number of cells in y-direction
nt = 600;     %% Number of time steps
mxst = 17;    %% Start of PEC section in x-direction
mxnd = 49;    %% End of PEC section in x-direction
myst = 33;    %% Start of PEC section in y-direction
mynd = 65;    %% End of PEC section in y-direction
% ky = 49;

% % % % % % % Initialize
Ez = zeros(nx,ny); %% z-component of E-field
Hx = zeros(nx,ny); %% x-component of H-field
Hy = zeros(nx,ny); %% y-component of H-field

mediaEz = ones(nx,ny);
mediaHx = ones(nx,ny);
mediaHy = ones(nx,ny);

Ez1 = zeros(nx,ny);
Ez2 = zeros(nx,ny);
Ez3 = zeros(nx,ny);
Ez4 = zeros(nx,ny);
Ez5 = zeros(nx,ny);

% % % % % % % % For 1D plots of Ez
Ezz1 = zeros(1,nx);
Ezz2 = zeros(1,nx);
Ezz3 = zeros(1,nx);
Ezz4 = zeros(1,nx);
Ezz5 = zeros(1,nx);
Ez6 = zeros(1,nx);
Ez7 = zeros(1,nx);

c1=5;	%%% 5th order LIAO ABC Coefficients
c2=10;
c3=10;
c4=5;
c5=1;


% c1=3;	%%% 3rd order LIAO ABC Coefficients
% c2=3;
% c3=1;


Ca = zeros(2,1);
Cb = zeros(2,1);
Da = zeros(2,1);
Db = zeros(2,1);

ds = asize/(mxnd - mxst - 1); %% Length Increment
dt = ds/(2*c);%% Time increment for 2-D

iflaga = 2;   %% 1 if Free space; 2 if PEC

% % % % % % % % Structure Definition
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

% % % % % % % % Main Program
for n = 1:nt
    
% % %     Compute z-component of E-field
    
    for i = 2 : nx - 1      
        for j = 2 : ny - 1 
            m  = mediaEz(i,j);
             if (i == 6)   %% Incident Field Source Excitation
                Ezs = Ez_inc(n);
             else
                Ezs = 0;
             end
             Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1))) + Ezs;
%              elseif (j == 1) %% Field at the bottom edge of the boundary
%                      Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
%                          - Hx(i,j));                      
             
        end
    end
    
% % % % % % % % % % % % The  5th LIAO Absorbing Boundary Condition
    for  j = 1:ny
	  Ez(1,j) = c1*Ez1(2,j)-c2*Ez2(3,j)+c3*Ez3(4,j)...
     		-c4*Ez4(5,j)+c5*Ez5(6,j);		%%%left side
    end
    for j = 1:ny
	  Ez(nx,j) = c1*Ez1(nx-1,j)-c2*Ez2(nx-2,j)+c3*Ez3(nx-3,j) ...
     		-c4*Ez4(nx-4,j)+c5*Ez5(nx-5,j);		%%%right side
    end 
    for i = 2:nx-1
	  Ez(i,1) = c1*Ez1(i,2)-c2*Ez2(i,3)+c3*Ez3(i,4) ...
    		-c4*Ez4(i,5)+c5*Ez5(i,6);	%%%bottom
    end
    for i = 2:nx-1
	  Ez(i,ny) = c1*Ez1(i,ny-1)-c2*Ez2(i,ny-2)+c3*Ez3(i,ny-3) ...
     		-c4*Ez4(i,ny-4)+c5*Ez5(i,ny-5);	%%%top
    end
  
% % % % % % % % % % % % The  2nd LIAO Absorbing Boundary Condition
% for j = 1: ny
%     Ez(1,j) = c1*Ez1(2,j) - c2*Ez2(3,j) + c3*Ez3(4,j); %%% Left Side
% end
% 
% for j = 1: ny
%     Ez(nx,j) = c1*Ez1(nx-1,j) - c2*Ez2(nx-2,j) + c3*Ez3(nx-3,j);
% end
% 
% for i = 1:nx
%     Ez(i,1) = c1*Ez1(i,2) - c2*Ez2(i,3) + c3*Ez3(i,4); %%% Bottom Side
% end
% 
% for i = 1:nx
%     Ez(i,ny) = c1*Ez1(i,ny-1) - c2*Ez2(i,ny-2) + c3*Ez3(i,ny-3); %%% Top Side
% end


%%%   Save previous 5 time step fields
    Ez5=Ez4;
	Ez4=Ez3;
	Ez3=Ez2;
	Ez2=Ez1;
	Ez1=Ez;    
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
 


% % % % % % %  PLOT Routine
 
%   hold on
colormap('viridis')
 hSurf = surf((Ez),'EdgeColor','none','LineStyle','none','BackFaceLighting', 'reverselit',...
      'FaceLighting','phong','FaceColor','interp','SpecularColorReflectance',.1);
% 
axis ([1 nx 1 ny -1 1])
view([-165 45]); % Perspective view
%  set(gcf,'Color','white');
%  set (gca,'FontName','times new roman')
%  xlabel(' Y Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
%  ylabel(' X Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
%  title(['H_y-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','tex','FontName','times new roman')
%  view([ 0 0 1])
%  if n == 10
%     Ezz1 = Ez(:,40);
% elseif n == 30
%     Ezz2 = Ez(:,40);
% elseif n == 60
%     Ezz3 = Ez(:,40);
% elseif n == 90
%     Ezz4 = Ez(:,40);
% elseif n == 120
%     Ezz5 = Ez(:,40);
% elseif n == 150
%     Ezz6 = Ez(:,40);
% elseif n == 180
%     Ezz7 = Ez(:,40);
% end
% grid off

% % % % % % Field Snapshots
%     if n == 20
%        plot2svg('Field_at_20.svg')
%     elseif n == 60
%         plot2svg('Field_at_60.svg')
%     elseif n == 120
%        plot2svg('Field_at_120.svg')
%     elseif n == 160
%         plot2svg('Field_at_160.svg')
%     elseif n == 180
%        plot2svg('Field_at_180.svg')
%     elseif n == 200
%       plot2svg('Field_at_200.svg')
%     elseif n == 240
%       plot2svg('Field_at_240.svg')
%     elseif n == 300
%         plot2svg('Field_at_300.svg')
%     elseif n == 320
%        plot2svg('Field_at_320.svg')  
%     elseif n == 340
%         plot2svg('Field_at_340.svg')
%         
%     end

box off
%      colormap('hot')
  M(:,n) = getframe(gcf);
%   frame = getframe(gcf);
%   writeVideo(writerObj,frame);
        
end
% plot2svg('demo_contour.svg')
%    close(writerObj);
% figure('Name', 'Time Steps of Electric Fields','Position',[100 100 850 600])
% set(gcf,'Color','white');
% title(['E-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','latex','FontName','times new roman')
% subplot(7,1,1)
% plot(Ezz1,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,2)
% plot(Ezz2,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,3)
% plot(Ezz3,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,4)
% plot(Ezz4,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,5)
% plot(Ezz5,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,6)
% plot(Ezz6,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gca,'XTickLabel',{})
% set (gca,'FontName','times new roman')
% subplot(7,1,7)
% plot(Ezz7,'Color','black','LineWidth',1.4)
% axis([0 nx -1.5 1.5])
% set(gcf,'Color','white');
% set (gca,'FontName','times new roman')
% axis([0 nx -1.5 1.5])
% xlabel(' Cells in x-direction','FontSize',11,'Interpreter','latex','FontName','times new roman')
% export_fig E1_ABC_n_2_gaussian '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%  
%    axis([400 780 -2 2.2])
% axis tight%     colorbar
% grid off
% figure('Name', 'Movie','Position',[10 10 850 600])
% set(gcf,'Color','white');
% set (gca,'FontName','times new roman')
% axis ([ 0 100 0 80 -1 1])
% movie(M,1,30);