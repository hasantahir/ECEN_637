% this program will evaluate the fdtd solution for the scattering problem 
% the scattering object is cylinder, the frequency is 1GHz
% the program can be modified to other scattering objects too



%%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Author: ujwol palanchoke, Jacobs University, Bremen,Germany
% application: 2D-FDTD simulation for cylindrical object
% posted date: 13 october 2009

%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


clear all
%declaration of all the vairiables.

lamda = 0.30;       %wavelengths in meter this value can be changed accourding to the specification
c = 3e8;            %speed of light
N = 200;            %no of simulation 
miu = 4*pi*1e-7;    % permiablility of vaccum
e_o = 8.854e-12;    % permitivity of vaccum
f = 1e9;            %frequency in hertz, can be varied accourding to specification 
lx = 10*lamda;       %total length in x-direction: simulation domain
ly = 10*lamda;       %total length in y-direction: simulation domain
dx = lamda/10;      % steps in x-direction: increments in x-direction
dy = dx;            % steps in y-direction: increments in y-direction
dt=1/(10*f);        % time step: should be choosen so that the stability criteria is satisfied
t = 0:dt:N*dt;      % total simulation steps
x = 0:dx:lx;
y = 0:dy:ly;
len_x = length(x);
len_y = length(y);
len_t = length(t);
Ez = zeros(len_x,len_y,len_t);  % initialize the fields to zero
Hx = zeros(len_x,len_y,len_t);
Hy = zeros(len_x,len_y,len_t);
e = ones(len_x,len_y);          % initialize relative-permitivity in grid to 1 (e_0)
e_r = e;
ko = 2*pi*f/c;      % wave number

 

%computation of the  permitivity.

for i = 1:len_x
    for j = 1:len_y
        if sqrt((lx/2 - x(i))^2 + (ly/2 - y(j))^2) <= 0.5*lamda     % using symmetric property of the problem domain
            e(i,j) = 2;     % relative permitivity of scattering object
        else
            e(i,j) = 1;
        end
    end
end


for i = 2:len_x-1
    for j = 2:len_y-1
        e_r(i,j) = 1/4*(e(i+1,j) + e(i-1,j) + e(i,j+1) + e(i,j-1)); % to consider the permitivity in boundary between space and scattering objet
    end
end
imagesc(e_r)

%computaion of Hx,Hy&Ez using central finite difference method

for n = 1:(len_t -1)
    for i = 2:(len_x-1)
        for j = 2:(len_y-1)
            Hx(i,j,n+1) = -(Ez(i,j+1,n) - Ez(i,j,n))*dt/(dy*miu) + Hx(i,j,n);
            Hy(i,j,n+1) =  (Ez(i+1,j,n) - Ez(i,j,n))*dt/(dx*miu) + Hy(i,j,n);
            
            Ez_a        =  (Hy(i,j,n+1) - Hy(i-1,j,n+1))/dx;
            Ez_b        =  (Hx(i,j,n+1) - Hx(i,j-1,n+1))/dy;
            Ez_c        =  (e_r(i,j) - 1)*e_o*2*pi*f*sin(2*pi*f*t(n)-ko*x(i));  %
            Ez_d        =  Ez(i,j,n);
            Ez(i,j,n+1) = (dt/(e_r(i,j)*e_o))*(Ez_a - Ez_b + Ez_c) + Ez_d;
            
            %applying one way boundary condition.
            
            if (i == 2)|(i == (len_x-2))
                Ez(i,j,n+1) = Ez(i-1,j,n) + (c*dt - dx)/(c*dt + dx)*(Ez(i-1,j,n+1) - Ez(i,j,n));                 % wave propagating in x-direction
            elseif (j == 2)|(j == (len_y-2))
                 Ez(i,j,n+1) =   Ez(i,j-1,n) + (c*dt - dx)/(c*dt + dx)*(Ez(i,j-1,n+1) - Ez(i,j,n));               % wave propagation in y-direction
            end  
            
        end
    end
    imagesc(Ez(:,:,n+1)')
     pause
     Display_t = [i,j,n]
end
   
% far field Ez
% converting the time domain solution to frequency domain.

time_resolution = 1;                             %resoloution of the simulation
Ez_s =Ez(:,:,(50:(50+(15*time_resolution))));    % taking field values from the simulation for time steps: 50-65:: stable state solutions.
Hx_s= Hx(:,:,(50:(50+(15*time_resolution))));
Hy_s =Hy(:,:,(50:(50+(15*time_resolution))));

%initialize the fft

Ez_fft=zeros(len_x-1,len_y-1);
Hx_fft=zeros(len_x-1,len_y-1);
Hy_fft=zeros(len_x-1,len_y-1);

%initialize the array for fixed ij and varying time (ie field at different
%time steps for same position)
Ez_s_tem = zeros(1,15*time_resolution);
Hx_s_tem = zeros(1,15*time_resolution); 
Hy_s_tem = zeros(1,15*time_resolution);

% calculate fft

for i=2:len_x-1
    for j=2:len_y-1
        for t=1:(15*time_resolution)
            Ez_s_tem(t)=Ez_s(i,j,t);
            Hx_s_tem(t)=Hx_s(i,j,t);
            Hy_s_tem(t)=Hy_s(i,j,t);
        end

            % evaluate fft at i,j
            Ez_fft_temp=fft(Ez_s_tem);  % fft of Ez field, Ez_fft_tem stores the coefficients
            Hx_fft_temp=fft(Hx_s_tem);     % fft of Hx field
            Hy_fft_temp=fft(Hy_s_tem);      % fft of Hy field
        
            % extracting the coefficients for fundamental component
      
             Ez_fft(i,j)=Ez_fft_temp(2); % second coefficient corresponds to fundamental frequency
             Hx_fft(i,j)=Hx_fft_temp(2);
            Hy_fft(i,j)=Hy_fft_temp(2);
    end 
end 

% evaluating Ez far field for different phi (azimuth angle). We need to make a close counter(loop)
% surrounding the scattering object. We also have to evaluate the magnetic
% current and electric current in the counter. Here phi varies from 0
% degree to 360 degree

% initilizing the constants 

    clear i ; clear j
    rho = 2;        % conductivity : arbitary for given problem
    etta = 377;         %wave impedance in ohms
  
    % defining the close counter (loop)
    I1=ceil(len_x/2-len_x/4);   % left boundary for varying J1
    J1=ceil(len_y/2-len_y/4);   % bottom boundary for varying I1
    I2=ceil(len_x/2+len_x/4);   % right boundary for varying J2
    J2=ceil(len_y/2+len_y/4);   % top boundary for varying I2
    
    % initialize the far field components from the counter of four sides
             dl=dx;              % incremental distance inside integral
             dx=dy;
             scaling_factor =(ko/4)*(sqrt(2/(pi*rho*ko))*exp(-(imag(complex(0,(ko*rho-pi/4))))));
             k=sqrt(-1);     % imaginary number j in lecture
    % evaluate the far fields due to four sides of boundary
        for phi=1:360
             Ez_left=0;
             Ez_right=0;
             Ez_top=0;
             Ez_bottom=0;
             
        for i=I1:I2
            for j=J1:J2
                angle=(phi*pi)/180;    % in radian
                si_ee= (i-1)*dl*cos(angle)+(j-1)*dl*sin(angle);    % power of exponential inside integral
                in_ex=exp(k*ko*si_ee);     % exponent inside integral
           
                    if (i==I1)      % left boundary
                    Jz(i,j)= -Hy_fft(i,j);
                    My(i,j) = -(Ez_fft(i-1,j)+Ez_fft(i+1,j))/2;      % points toward negative y axis
                    Mx(i,j)=0;
                    
                   Ez_left=Ez_left + scaling_factor*(-Mx(i,j)*sin(angle)+My(i,j)*cos(angle)-etta*Jz(i,j))*in_ex*dl;
                    
                elseif (i==I2)      % right boundary
                    Jz(i,j)= Hy_fft(i,j);
                    My(i,j) = ((Ez_fft(i+1,j)+Ez_fft(i-1,j))/2);      % points toward positive y axis
                    Mx(i,j) =0;
                     Ez_right=Ez_right+scaling_factor*(-Mx(i,j)*sin(angle)+My(i,j)*cos(angle)-etta*Jz(i,j))*in_ex*dl;
                    
                elseif (j==J1)      % bottom boundary
                    Jz(i,j)= Hx_fft(i,j);
                    My(i,j) =0; 
                    Mx(i,j)=((Ez_fft(i,j+1)+Ez_fft(i,j-1))/2);      % points toward positive x axis
                     Ez_bottom=Ez_bottom+scaling_factor*(-Mx(i,j)*sin(angle)+My(i,j)*cos(angle)-etta*Jz(i,j))*in_ex*dl;
                    
                elseif (j==J2)      % top boundary
                    Jz(i,j)= -Hx_fft(i,j);
                    My(i,j) =0;  
                    Mx(i,j)= -(Ez_fft(i,j+1)+Ez_fft(i,j-1))/2;  % points toward negative x axis
                    Ez_top=Ez_top+scaling_factor*(-Mx(i,j)*sin(angle)+My(i,j)*cos(angle)-etta*Jz(i,j))*in_ex*dl;
                end
            end
        end
        E_far_field(phi)= Ez_left+Ez_right+Ez_bottom+Ez_top;
 end 
 figure
phi = 1:1:360;
plot(phi,abs(E_far_field(phi)))