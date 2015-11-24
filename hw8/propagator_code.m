close all;
clear all;
% % Definition of paramaters

nx = 600;   % Numerical Space
nx2 = nx/2;

nxst = 250; % Slab Start
nxnd = 309; % Slab End

num_freq = 1000; % higher the better
num_samples = 16384; % higher --> more resolution

nt = 1000; % longer time

num_med = 2;

coll_pt = 240; % point where we calculate the Reflection Coefficient
beta = 10; % for proper termination of gaussian source

mu = 4*pi*1e-7;
eps0 = 8.854187817e-12;
epsR = 4;

slab_width = .09; % as described in Luebber's paper

media = ones(1,nx); %%!Field medium index    
gamma = zeros(1,num_freq);
d_gamma = zeros(1,num_freq);
d_num_freq = zeros(1,num_freq);
Source = zeros(1,nt);
Recorder = zeros(1,nt);

Ca = zeros(2,1);
Cb = zeros(2,1);
Da = zeros(2,1);
Db = zeros(2,1);    
eta = zeros(2,1);

Ez = zeros(1,nx);				%%%!zero the incident and total fields
temp_Ez = zeros(1,nx);
Hy = zeros(1,nx);	
temp_Hy = zeros(1,nx);
Ezinc = zeros(1,nx);
temp_inc = zeros(1,nx);
Hyinc = zeros(1,nx);	

ftEinc = zeros(1,num_freq);		%%%!zero the DFTs
ftEref = zeros(1,num_freq);

%%
%%%!medium 1 is free space, medium 2 is Er=4.0
    				
% % % % % 		media = 1;		%%%!free space everywhere

        for i = nxst : nxnd		%%%!insert the slab in to free space
            
                media(i) = 2;	%%%!Ez, Hy exist at every index point in the slab
                
        end
        
%% 
eta(1) = sqrt( mu/eps0 ); % free-space impedance
eta(2) = sqrt( mu/(eps0 * epsR) ); % impedance in slab
dx = slab_width / ( nxnd - nxst + 1 );  %%!????? should it be (nsnd-nxst)=5-1=4
c = 1 / ( sqrt( mu * eps0 ) );		
dt = dx/(1 * c );	% Magic ratio

% % % % Coefficients of free-space
Ca(1) = 1;
Cb(1) = eta(1);
Da(1) = 1/eta(1);
Db(1) = 1;

% Ca(1) = 1;
% Cb(1) = dt/(dx*eps0);
% Da(1) = dt/(dx*mu);
% Db(1) = 1;
% 
% % % % !  Slab Coefficients

Ca(2) = 1;
Cb(2) = eta(2);
Da(2) = 1/eta(2);
Db(2) = 1;
% 
% Ca(2) = 1;
% Cb(2) = dt/(dx*eps0*epsR);
% Da(2) = dt/(dx*mu);
% Db(2) = 1;

% Ca = [ 1;1];
% Cb = [3.767303134749689e+02;94.182578368742210];
% Da = [1;1];
% Db = [0.002654418729345;0.002654418729345];
%%
for n = 1 : nt
    
% %     ABC
    Ez(1) = Ez(2);
    Ez(nx) = Ez(nx-1);
    Ezinc(1) = Ez(2);
    Ezinc(nx) = Ezinc(nx-1);
    
    temp_Ez = Ez;
% %     
    for i = 2 : nx - 1
        
        if (i == 50)
            
            Ezs = Ez_inc(n);   % Source Excitation
            Source(n) = Ezs;    % Source Recorder
        
        else
            
            Ezs = 0;
            
        end
        
        
        m = media(i);   % Set the material at current point in space
        
        % Calculate the Electric field
        Ez(i) = .5*( Ca(m)*( temp_Ez(i+1) + temp_Ez(i-1) )...
            + Cb(m)*( Hy(i+1) - Hy(i-1) ) ) + Ezs;
    end
    
    for i = 2 : nx - 1
        
        m = media(i);
        % Calculate the Magnetic Field
        Hy(i) = .5*( -Da(m)*( temp_Ez(i+1) - temp_Ez(i-1) )...
            + Db(m)*( Hy(i+1) + Hy(i-1) ) );
        
    end

%      temp_Ez = Ez;
%      temp_Hy = Hy;
    Recorder(n) = Ez(coll_pt);
    plot(1:nx,(Hy));
    pause(.001);
end
        
        
          
