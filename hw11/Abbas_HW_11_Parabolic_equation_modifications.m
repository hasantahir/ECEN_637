function Abbas_HW_11_Parabolic_equation_modifications()
% Hasan Tahir Abbas
% ECEN 637
% Homework 10: Computation of Input Admittance and Radiation Pattern of a
% Wire Antenna
% 11/27/2015
%
%
%  Pocklington Integral Equation is used to compute the Admittance and  
%  Radiation Pattern as explained in:
%  Gibson, Walton C. "The Method of Moments in Electromagnetics," Taylor
%                                                  and Francis/CRC, 2008.
%  
%
% Gauss_Quadrature based Numerical Integration (Built-in Function)
% Input Admittance
% Far-field pattern calculation
% Utilization of Toeplitz nature of the impedance matrix (Only one row
% needs to be computed)
%
% PARAMETERS AND VARIABLES LIST
% -----------------------------
% M = Total Wire Segments
% f = Frequency
% c = Speed of light
% xmu = Permeability of free space
% eps0 = Permittivity of free space
% omega = Angular Frequency
% K = Propagation Constant
% lambda = Wavelength
% a = Wire Radius;
% L = Wire Length
% E_rad = Radiated Far-field
% I = Current on the wire
% Y = Input Admittance of the wire
% num_theta = number of angle values for polar plot
% num_length = number of lengths of wire
%
% FUNCTIONS LIST
% -----------------------------
% Pocklington() = Calculates the Pocklington EFIE solution
% Pocklington_Pattern() = Calculates the Pocklington EFIE Far-field pattern
% Pocklington_Current() = Calculates the Current on the wire
% -----------------------------
close all
% *************************
%% Parameters
% *************************
%
% Global variables are used to span across all the functions in this code
% According to MATLAB's documentation, a better and safer option will be
% persistent type variables
%
global N Nx 
global p z pat_p f_z Picture x
% *************************
% *************************

%##################################
%% Main Program
% *************************
% *************************
% *************************
data();
f_z = gen_source();
save f_z.mat
for x = 1 : Nx
        refr_index_n();
        p_space();
        z_space();
        for j = 1 : N
            Picture(j,x) = f_z(j);
        end
% plots(p, pat_p, z, f_z);        
end  
plots(p, pat_p, z, f_z);  
save var.mat
end
%##################################
%% End Main Program
% *************************
% *************************
% *************************

%% Initialize
% *************************
% *************************
function data()
% Set all the variables to zero

global N Nx 
global p_max z_max dp dz
global pe za dx a a_e 
global p z pat_p f_z
global pattern fz temp m
global freq lambda c kappa  theta_max
global beamwidth Length
global epsi_o epsi_r sigma omega epsi_c Picture
global xx

% Data to be used
% *************************
N = 16384;
Nx = 1024;
% Nx = 1024;
freq = 3e9 ; % Antenna operation frequency

c = 2.99792458e8;
epsi_o = 8.854e-12;
epsi_r = 69;
sigma = 6.5;
omega = 2*pi*freq; 
epsi_c = epsi_o * epsi_r + 1i* sigma/omega;


lambda = c/freq; % One wavelength
kappa = 2*pi/lambda; 
Length = 38*lambda; % Length of aperture in wavelengths


za = 31;
dx = 200; % Radius of the earth
a = 6.37e6; % radius of the earth
a_e = a*(4/3); % effective radius of the earth


theta_max = 13;
beamwidth = 22.2718567e-3; % Half-power beamwidth of Taylor pattern
p_max = kappa*sind(theta_max); % Maximum propagation angle %% sind() for degree based arguments
dp = 2*p_max/(N-1);
dz = 1/(p_max/pi); % Delta increment in z-direction
z_max = (N-1)/2*dz; % Physical area analyzed
pe = kappa*sin(beamwidth/2); % Beam tilt angle


z = dz*(0:N-1) - z_max;
p = dp*(0:N-1) - p_max;
xx = (1:Nx)*dx; % Horizontal x-axis profile positions



pattern = zeros(1,N);
fz = zeros(1,N);
temp = zeros(1,N);
m = zeros(1,N);
pat_p = zeros (1,N);
f_z = zeros (1,N);
Picture = zeros(N,Nx);


end
% *************************
% *************************
%
%% Get the refractive index for further calculation
%
% *************************
% *************************
function refr_index_n()

global N a z m

global epsi_o index epsi_c
global dNdz zmax
global xx x choice

if x == 1
    disp('        Choose the media           ');
    disp('                                   ');
    disp('                                   ');
    disp('1. Free space');
    disp('2. for Standard Atmosphere');
    disp('3. for Surface Duct');
    disp('4. for Graded Duct');
    prompt = 'Enter option from 1-4:        ';
    choice = input(prompt);
end

% Set the environment

switch (choice)
    
    
    
    %% Free-space with flat-earth
    case(1)
%         disp('           Choice = 1                  ');
        for i = 1 : N
            if z(i) < 0 % Sea Water
                m(i) = (epsi_c)/(epsi_o) - 1;
            else
                index = 1; %  Flat earth condition (2*z/a = 0.)
                m(i) = (index^2 - 1);
            end
        end
        
    %% Troposphere condition    
    case(2)
%         disp('           Choice = 2                  ');
        for  i = 1 : N
            if z(i) < 0 % Sea Water
                m(i) = (epsi_c)/(epsi_o) - 1;
            else
                index = 1 + (300 - 0.0394 * z(i)) * 1e-6; % !Standard atmosphere condition
                m(i) = index^2 - 1 + 2*z(i)/a;
            end
        end
        
    %% Surface Duct condition    
    case(3)
%         disp('           Choice = 3                  ');
        for i = 1 : N
            if (z(i) < 0) % Sea water
                m(i) = epsi_c/epsi_o - 1;
            elseif (z(i) > 0 && z(i) <= 37) % Surface Duct
                index = 1 + (300 - .5*z(i))*1e-6; 
                m(i) = index^2 - 1 + 2*z(i)/a;
                k = i;
            else
                index = 1 + (300 - .5*z(k) - .0394*( z(i) - z(k)))*1e-6; % Standard atmosphere
                m(i) = index^2 - 1 + 2*z(i)/a;
            end
        end        
        

    %% Graded Surface Duct Condition
    case(4)
%         disp('           Choice = 4                  ');
        if ( xx(x) <= 40000)
            dNdz = (.5- .167)/40000*xx(x) - .5;
            zmax = (150-37)/40000*xx(x) + 37;
        end
        for i = 1 : N
            if (z(i) < 0) % Sea water
                m(i) = epsi_c/epsi_o - 1;
            elseif (z(i) > 0 && z(i) <= zmax) % Graded surface Duct
                index = 1 + (300 + dNdz*z(i))*1e-6; 
                m(i) = index^2 - 1 + 2*z(i)/a;
                k = i;
            else
                index = 1 + (300 + dNdz*z(k) - .0394*( z(i) - z(k)))*1e-6; % Standard atmosphere
                m(i) = index^2 - 1 + 2*z(i)/a;
            end
        end
        
    otherwise
        disp('%%%%%ERROR%%%%%%');
        disp('Wrong choice entered');
end
        
save index.mat
end
% *************************
% *************************
%
% *************************
%% Generate the source
%
% *************************
% *************************
function signal = gen_source()
% The aperture distribution and its far field pattern

global N z Length
global pe za fz f_z zz p
global pattern pat_p;
%% Sinc Source
% x = 0;
% for i = 1 : N
%     pp = p(i) - pe;
%     if pp == 0
%         pattern(i) = Length;
%     else
%         pattern(i) = 2*sin(pp*Length/2)/pp*exp(-1i*pp*za);
%     end
%     pat_p(i) = abs(pattern(i));
% end
%% Create Taylor Line pattern
% SLL = -20 dB 
% n_bar = 10;
% ! Taylor pattern SLL=-20dB, n_bar=10
for  k = 1 : N
    if  z(k) >= (za-Length/2) && z(k) <= (za+Length/2)
            zz = z(k) - za;
            fz(k) = exp(1i*pe*z(k))...
            *( 1.0 ...
            + 2.0*( 0.0977020907)*cos(2*pi*1.0*zz/Length)...
            + 2.0*( 0.045938932)*cos(2*pi*2.0*zz/Length)...
            + 2.0*(-0.0567750651)*cos(2*pi*3.0*zz/Length)...
            + 2.0*( 0.0543140899)*cos(2*pi*4.0*zz/Length)...
            + 2.0*(-0.0472939433)*cos(2*pi*5.0*zz/Length)...
            + 2.0*( 0.0381149255)*cos(2*pi*6.0*zz/Length)...
            + 2.0*(-0.0279769765)*cos(2*pi*7.0*zz/Length)...
            + 2.0*( 0.0177444933)*cos(2*pi*8.0*zz/Length)...
            + 2.0*(-0.0081655839)*cos(2*pi*9.0*zz/Length))*(1.0/Length);
    else
        fz(k) = 0.0;
    end
   
end
f_z = abs(fz(k));
save taylor.mat
% signal =pat_p;\
signal = f_z;
end
% *************************
% *************************
%
%% Fourier Transform of known aperture distribution
%
% *************************
% *************************
function p_space()

    
global N  dz pat_p 
global pattern  ii

ii = -1;
temp1 = replace(); % Store data in temporary array for calculation
% temp2 = four1(temp1, N, ii); % FFT of Data
temp2 = fft(temp1,N);
temp2 = fftshift(temp2);
pattern = temp2*dz;
pat_p = draw_data(pattern);
% save pattern.mat
end

% *************************
% *************************
%
%% Store the data for calculation
%
% *************************
% *************************
function out = replace()

global pattern fz ii

if ii == 1
    out = pattern;
else
    out = fz;
end % end if

% save replace.mat
end % end function
% *************************
% *************************
%
%% Store the data for Plotting in p-space
%
% *************************
% *************************
function out = draw_data(in)

global N 

out = zeros(1,N);
for  i = 1 : N
    if ( i <= floor(N/2+1))
        k = i + floor(N/2) - 1;
        out(k) = abs(in(i));
    else
        k = i - floor(N/2) - 1;
        out(k) = abs(in(i));
    end
end
end
% *************************
% *************************
%
%% Inverse Fourier Transform
%
% *************************
% *************************
function z_space()
% Take inverse fourier transform of the previous pattern at the next step

global N dz
global f_z fz ii

ii = 1;
temp1 = replace();  % Store data in temporary array
temp2 = pre_process(temp1);  % Take care of marching
% temp3 = four1(temp2,N,ii);   % Inverse FFT
temp3 = N*ifft(temp2,N);
temp3 = fftshift(temp3);
temp4 = post_process(temp3); % Take care of index difference in z-direction
fmax = max_fz(temp4);

fz = temp4/(N*dz);
f_z = abs(fz)/fmax;
save z_space.mat
end

%% Marching
%
% *************************
% *************************
function out = pre_process(in)
% Take care of marching to the next step
global N dp dx pp kappa  

for i = 1 : N
    if i <= (N/2 + 1)
        pp = dp*(i-1);
    else
        pp = dp*(i-N-1);
    end
    in(i) = in(i) * exp(-1i*pp^2 * dx/(2 * kappa));
end
out = in;
end % end function
% *************************
% *************************
%
%% Refractive Index difference corretion and Hanning window
%
% *************************
% *************************
function out = post_process(in)
% Take care of refractive index difference in z-direction
global N z_max
global za dx z m kappa han 

out = zeros(1,N);
for i = 1 : N
    han = .5 + .5 * cos( 2*pi*(z(i) - za)/(2*z_max));
    out(i) = in(i)* exp(1i* kappa * m(i) * dx/2) * han;
end
end % end function
% *************************
% *************************
%
%% Maximum of field distribution
%
% *************************
% *************************
function out = max_fz(in)
% 
% Find maximum value of field distribution
global N dz 

fmax = 0;
for i = 1 : N
    ftemp = abs(in(i))/(N*dz);
    if ftemp >= fmax
        fmax = ftemp;
    else
        fmax = fmax;
    end
    out = fmax;
end
save fmax.mat
end % end function
% *************************
% *************************
%
%% FFT from numerical recipes
%
% *************************
% *************************
function dataout = four1(datain, nn, isign)

data = zeros(1,2*nn);
dataout = zeros(1,nn);
for i = 1 : nn
    data(2*i-1) = real(datain(i));
    data(2*i) = imag(datain(i));
end

n = 2*nn;
j = 1;
for i = 1 : 2 : n
    if (j > i)
        
        tempr = data(j);
        tempi = data(j+1);
        data(j) = data(i);
        data(j+1) = data(i+1);
        data(i) = tempr;
        data(i+1) = tempi;
    end
    m = n/2;
    if ( m >= 2 && j > m)
        j = j - m;
        m = m/2;
    end
    j = j + m;
end
mmax = 2;
if (n > mmax)
    istep = 2*mmax;
    theta = 2*pi/(isign * mmax);
    wpr = -2*sin(0.5*theta)^2;
    wpi = sin(theta);
    wr = 1;
    wi = 0;
    for  m = 1 : 2 : mmax
        for i = m : istep : n
            j = i + mmax;
            tempr = (wr) * data(j) - (wi) * data(j+1);
            tempi = (wr) * data(j+1) + (wi) * data(j);
            data(j) = data(i) - tempr;
            data(j+1) = data(i+1) - tempi;
            data(i) = data(i)+tempr;
            data(i+1) = data(i+1) + tempi;
        end
        wtemp = wr;
        wr = wr * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    end
    mmax = istep;
    
end
for i = 1:nn
    dataout(i) = data(2*i-1)+ 1i*data(2*i);
end
end % subroutine four1
% *************************
% *************************
%
%% FFT from numerical recipes
%
% *************************
% *************************
function plots(p, pat_p, z, f_z)

global Nx dx  Picture


dxx = (0:Nx-1)*dx;
% Pic = Picture';
figure(1)
plot(pat_p, p);
pause(.001)

figure(2)
plot(20*log10(f_z),z);
axis([ -50 0 0 500])
pause(.001)
figure(3)
% % % surf(dxx',z',abs(Picture));
% % % shading interp;
% % % axis([0 16000 0 15.625 -1 1]);
% % % view([0 90]);
surf(dxx,z,abs(Picture));
shading interp;
axis([0 16000 -50 300 -1 1]);
view([0 90]);
% 
save plts.mat







% % Plot Input Admittance
% figure(1);
% H = plot(p, pat_p, pi*Lengths/lambda, imag(Yin)*1e3);
% xlim([1 3.5])
% H(1).Color = 'black';
% H(1).LineWidth = 1.4;
% H(2).Color = 'black';
% H(2).LineWidth = 1.4;
% H(2).LineStyle = '--';
% title(['Input Admittance versus length at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
% set(gcf,'Color','white'); % Set background color to white
% set(gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% ax = gca;
% ax.XTick = [0.5 1 2 3];
% xlabel('${\beta L}/2$ ','Interpreter','latex'); % X-axis label
% ylabel('G and B in mmhos ','Interpreter','latex'); % y-axis label
% grid on
% legend('Real Part', 'Imaginary Part');
% % cleanfigure();
% % matlab2tikz('filename',sprintf('ECEN637_HW10_Admittance_plot.tex'));
% 
% % First Polar Plot
% figure(2)
% h1 = polar(theta, abs(E_rad_4by5)./abs(E_rad_4by5(91)));
% h1.Color = 'black';
% h1.LineWidth = 1.4;
% title(['Radiation Pattern of a Wire of length $4\times\lambda/5$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
% set(gcf,'Color','white'); % Set background color to white
% set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% % cleanfigure();
% % matlab2tikz('filename',sprintf('ECEN637_HW10_Polar_plot_h_4by5lambda.tex'));
% 
% % First Polar Plot
% figure(3)
% P = polar(theta, 1 * ones(size(theta))); % this is done to set the polar plot limit to 1.
% set(P, 'Visible', 'off')
% hold on
% h3 = polar(theta, abs(E_rad)./abs(E_rad(91)));
% h3.Color = 'black';
% h3.LineWidth = 1.4;
% title(['Radiation Pattern of a Wire of length $\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
% set(gcf,'Color','white'); % Set background color to white
% set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% % cleanfigure();
% % matlab2tikz('filename',sprintf('ECEN637_HW10_Polar_plot_h_lambda.tex'));
% 
% % Current Plot
% figure(4)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
% H = plot(x, real(I_wire),x, imag(I_wire));
% ax = gca;
% H(1).Color = 'black';
% H(1).LineWidth = 1.4;
% H(2).Color = 'black';
% H(2).LineWidth = 1.4;
% H(2).LineStyle = '--';
% title(['Current on the wire of half-length $ h = .35\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
% set(gcf,'Color','white'); % Set background color to white
% set(gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
% ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
% ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
% ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
% axis([ -.3498 .3498 -2.5e-3 2.5e-3]);
% hold on
% xlabel('z')
% ylabel('A')
% legend('Real Part', 'Imaginary Part');
% grid on
% % cleanfigure();
% % matlab2tikz('filename',sprintf('ECEN637_HW10_Current_on_wire_plot.tex'));
% 
end
