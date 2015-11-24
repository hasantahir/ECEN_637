clear all;
close all;
% writerObj = VideoWriter('Hy_ABC_n_5_yee_dimension.avi'); %% Movie Plotting file destination
% open(writerObj); 
figure('Name', 'Time Steps of Electric Fields','Position',[0 0 1920 1080])

% ! 1D FDTD dielectric slab reflection coefficient code
% !    recreating Leubber's results in Fig 1 and 6b
% 
% ! Soft Gaussian pulse excitation
% ! observation of reflected field
% ! Calculation of DFT and reflection coefficient
% ! Total absorption at boundaries
% 
% !   nx = number of spatial steps in the x-direction
% !   nxst = start of the dielectric slab
% !   nxnd = end of the dielectric slab
% !   slabwidth = width of the slab in meters
% !   num_freq = number of frequency response sample frequencies
% !   nsnap = number of snapshots of the time domain field
% !   nt = number of time steps
% !   num_med = number of media in the numereical space
% !   coll_pt = time domain field data collection point
% !   mu/epso = free space permeability/permittivity
% !   cpsR = slab dielectric constant 
% !   sigma = slab electric conductivity
% !   rho_prime = slab magnetic conductivity
%  
% !--------------------------------------------------------------------------

% % % % % Initialization   
    f = 4.977e9; %% 500 GHz from the Luebber's paper
    p = 13;
    r = 200;
    nx = 600; 
    nx2 = (nx)/2; 
    nxst = 250; %% Start of the slab 
    nxnd = 309; %% End of the slab
    num_freq = 500; 
    num_samples = 16384;
    nt = 1000; 
    num_med = 2; 
    coll_pt = 240;
	beta = 10; 
    mu = 4 * pi * 1e-7; 
    eps0 = 8.854187817e-12;
    epsR = 4.0;
    sigma = 0;
    rho_prime = 0;
    slabwidth = 0.09; 
	mediaEz = ones(1,nx);
    mediaHy	= ones(1,nx); %%!Field medium index    
    gamma = zeros(1,num_freq);
    d_gamma = zeros(1,num_freq);
    d_num_freq = zeros(1,num_freq);
    Source = zeros(1,nt);
    
    Ca = zeros(2,1);
    Cb = zeros(2,1);
    Da = zeros(2,1);
    Db = zeros(2,1);
% % !*********************************************************************
% % % Initialize
	dx = slabwidth / ( nxnd - nxst + 1 );  %%!????? should it be (nsnd-nxst)=5-1=4
	c = 1 / ( sqrt( mu * eps0 ) );		
	dt = dx/(1 * c );			

% % % !  free space coefficients

	Ca(1) = (1 - sigma * dt / (2 * eps0) ) / ( 1 + sigma * dt / (2 * eps0) );
	Cb(1) = (dt / (eps0 * dx) )/( 1 + sigma * dt /(2 * eps0) );
	Da(1) = (1 - rho_prime * dt/(2*mu) )/(1 + rho_prime *dt/(2*mu) );
	Db(1) = (dt/(mu * dx) )/(1 + rho_prime * dt/(2 * mu) );

% % % !  slab coefficients

	Ca(2) = (1 - sigma * dt / (2 * epsR * eps0) )/(1 + sigma * dt /(2 * epsR * eps0) );
	Cb(2) = (dt /(epsR * eps0 * dx) ) / (1 + sigma * dt/(2 * epsR * eps0) );
	Da(2) = Da(1);
	Db(2) = Db(1);

	Ez = zeros(1,nx);				%%%!zero the incident and total fields
	Hy = zeros(1,nx);	
	Ezinc = zeros(1,nx);	
	Hyinc = zeros(1,nx);	
	
	ftEinc = zeros(1,num_freq);		%%%!zero the DFTs
	ftEref = zeros(1,num_freq);
    
% % !--------------------------------------------------------------------------
%%%!medium 1 is free space, medium 2 is Er=4.0
    				
% % % % % 		mediaEz = 1;		%%%!free space everywhere
% % % % % 		mediaHy = 1;

        for i = nxst:nxnd		%%%!insert the slab into free space
            
                mediaEz(i) = 2;	%%%!Ez exist at ever index point in the slab
                
        end
        
        for i = nxst:nxnd-1	%%%!Hy does not extend to the right index of the slab
            
            mediaHy(i) = 2;
            
        end

        
% % % !********************************Main Program*************************************		
% 	tic
	for  n = 1:nt				%!initiate time stepping procedure
%         tic
        % % Electric Field Computation-------------------------------------------------------------
% % % % ABC 
        Ez(1) = Ez(2);				%%%!left side total field perfect ABC
        Ez(nx) = Ez(nx-1);			%%%!right side total field perfect ABC

        Ezinc(1) = Ezinc(2);		%%%!left side incident field perfect ABC
        Ezinc(nx) = Ezinc(nx-1);	%%%!right side incident field perfect ABC

        for i = 2:nx-1				%%%!Ez Total (field with slab)
            
            if (i == 50)
                
                Ezs = Ez_inc(n);		%%%!Incident field source excitation
                Source(n) = Ezs;
                
            else
                
                Ezs = 0;
                
            end
            m = mediaEz(i);
            Ez(i) = Ez(i) * Ca(m) + ( Hy(i) - Hy(i-1) ) * Cb(m)...
                 + Ezs;   			%%%!soft source
             
            

        end

        for i = 2:nx-1				%%%!Ez Incident (no slab) 
            
            if (i == 50)
                
                Ezs = Ez_inc(n);		%%%!Incident field source excitation
            
            else
                
                Ezs = 0;
            
            end
            m = 1;
            Ezinc(i) = Ezinc(i) * Ca(m)...
                   + ( Hyinc(i) - Hyinc(i-1) ) * Cb(m)...
                   + Ezs; 		%%%!soft source
        end
% % 	  call adv_hfield
      % Magnetic Field Computation--------------------------------------------------------------------------

        for i = 1: nx-1			%%%!advance H total
            
            m = mediaHy(i);
            Hy(i) = Hy(i) * Da(m)...
                + ( Ez(i+1) - Ez(i) ) * Db(m);
        end

        for i = 1:nx-1			%%%!advance H incident
            
            m = 1;			
            Hyinc(i) = Hyinc(i) * Da(m)...
                    + ( Ezinc(i+1) - Ezinc(i) ) * Db(m);
                
        end


%  plot(1:nx,10*log10( abs(Ezinc)))

      	Einc = Ezinc(coll_pt);
        Eref = Ez(coll_pt)-Einc;

        for k = 1 : num_freq
            
              dft_exp = -2 * 1i * pi * k * n / num_samples;
              ftEinc(k) = ftEinc(k) + dt * Einc * exp(dft_exp);	%%%!incident field dft
              ftEref(k) = ftEref(k) + dt * Eref * exp(dft_exp);	%%%!reflected field dft
              
        end
% 	  call dft				!calculate the dft 
%     toc

     hold on
%     axis auto
%     M(:,n) = getframe(gcf);
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);


set(gcf,'Color','white');
title(['E-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','latex','FontName','times new roman')
plot(Ez,'Color','black','LineWidth',1.4)
axis([0 nx -1 1])
% set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
line([nxst nxst],[-1  1],'Color','k') ;
line([nxnd nxnd],[-1  1],'Color','k') ;
plot(Ez)
% % % % % Field Snapshots
   % % 	iplot(1)=215		!My time snaps are different from luebber's 	
% % 	iplot(2)=250		!because my excitation source is different.	
% % 	iplot(3)=320
% % 	iplot(4)=373		
% % 	iplot(5)=395		
% % 	iplot(6)=488


hold off
 pause(.001)
clf
end
    

% % 	call reflect			!calculate the reflection coeff
    for i = 1:num_freq
        
        if ( abs( ftEinc(i) ) > 0.0 )
            
            gamma(i) = ( abs(ftEref(i) ) ) / ( abs( ftEinc(i) ) );
            
        else
            
            gamma(i) = 0;
            
        end
        
        d_num_freq(i) = i * (1  / (num_samples * dt) );
        d_gamma(i) = gamma(i);
        
    end	 
%toc  

% % Plot Routine--------------------------------------------------------------------------

figure(2)
plot(d_num_freq,d_gamma);
ylabels = get(gca, 'YTickLabel');
ylabels = linspace(0,1,length(ylabels));
set(gca,'YTickLabel',ylabels);


xlabels = get(gca, 'XTickLabel');
xlabels = linspace(0,7,length(xlabels));
set(gca,'XTickLabel',xlabels);
title('Reflection Coefficient ');
xlabel('Frequency [GHz]');
ylabel('\Gamma');	
plot2svg('Reflection_coefficient_new.svg')	

figure(3)
plot(d_num_freq,abs(ftEinc));
ylabels = get(gca, 'YTickLabel');
ylabels = linspace(0,1,length(ylabels));
set(gca,'YTickLabel',ylabels);

xlabels = get(gca, 'XTickLabel');
xlabels = linspace(0,7,length(xlabels));
set(gca,'XTickLabel',xlabels);
title('Signal Spectrum');
xlabel('Frequency [GHz]');
ylabel('Normalized Magnitude');	
plot2svg('Source_spectrum_new.svg')

f=figure(4);
plot(1:nt,Source);
axis([ 0 500 -1 1])
title('Source Signal Waveform');
xlabel('Time [ps]');
ylabel('Magnitude');
plot2svg('Source_waveform_new.svg')
	

% !--------------------------------------------------------------------------
