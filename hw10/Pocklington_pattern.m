function [ERad, I] = Pocklington_pattern(TotalSegments, L, a, f, NumTheta)

c = 2.99792458e8; % Speed of light
xmu = 4*pi*1e-7;  % Permeability of free space
eps0 = 8.854187817e-12; % Permittivity of free space
omega = 2.0*pi*f;
K = omega*sqrt(xmu*eps0);
TotalElements = TotalSegments;
deltaz = L / TotalSegments;

z = linspace(-0.5*L , 0.5*L - deltaz, TotalSegments);
Z = zeros(TotalElements);


for j = 1 : 1
    for k = 1 : TotalSegments 
        
        fun = @(zp) exp(-1i*K*sqrt(( z(j)+.5*deltaz - zp).^2 + a^2))...
            ./sqrt(( z(j)+.5*deltaz - zp).^2 + a^2);
        zp_upper = z(k) + deltaz;
        dz_upper = z(j)+.5*deltaz - zp_upper;
        zp_lower = z(k);
        dz_lower = z(j)+.5*deltaz - zp_lower;
%         int_part = integral(fun,zp_lower,zp_upper,'RelTol',1e-8,'AbsTol',1e-12);
        int_part =  quadgk(fun,zp_lower,zp_upper,'RelTol',1e-8,'AbsTol',1e-12);
%                 int_part = 0;
        R_upper = sqrt(dz_upper^2 + a^2);
        R_lower = sqrt(dz_lower^2 + a^2);
        sum_part_upper = (dz_upper)*(1 + 1i*K*R_upper)...
            ./R_upper^3*exp(-1i*K*R_upper);
        sum_part_lower = (dz_lower)*(1 + 1i*K*R_lower)...
            ./R_lower^3*exp(-1i*K*R_lower);
        Z(j,k) = K^2*int_part + sum_part_upper - sum_part_lower;

        if j == k

            R_upper = sqrt(  (-deltaz/2)^2 + a^2);
            R_lower = sqrt((deltaz/2)^2 + a^2);
            sum_part_upper = (-deltaz/2)*(1 + 1i*K*R_upper)...
                ./R_upper^3*exp(-1i*K*R_upper);
            sum_part_lower = (deltaz/2)*(1 + 1i*K*R_lower)...
                ./R_lower^3*exp(-1i*K*R_lower);
            num = sqrt(1 + 4*a^2/deltaz^2) + 1;
            denom = sqrt(1 + 4*a^2/deltaz^2) - 1;
            self = (log(num/denom) - 1i*K*deltaz);
            Z(j,k) = K^2*self + sum_part_upper - sum_part_lower;
         end
    end
end

Z = toeplitz(real(Z(1,:))) + 1i*toeplitz(imag(Z(1,:)));
V = zeros(TotalElements-2,1);
V(floor((TotalSegments-2)/2)+1) = -1i*4*pi*omega*eps0*(1.0/deltaz);   % delta-gap excitation
% solve for currents
% Current =A\rhs;
I = zeros(TotalElements,1);
I(2:TotalElements-1) =Z(2:TotalElements-1,2:TotalElements-1)\V;
Zin = 1.0 / I(floor((TotalSegments-2)/2)+1);
% Yin = 1./Zin;

% array of theta for theta-pol radiated field
Theta = linspace(0.0, 2.0*pi, NumTheta);
for iTheta = 1:NumTheta
    cosTheta = cos(Theta(iTheta));
    sinTheta = sin(Theta(iTheta));
    ERad(iTheta) = 0.0;
    for m = 1:TotalElements
        z_m = z(m) + 0.5*deltaz;
        ERad(iTheta) = ERad(iTheta) + deltaz*I(m)*sinTheta*exp(i*K*z_m*cosTheta);
    end
    ERad(iTheta) = -(1i*omega*xmu/(4.0*pi))*ERad(iTheta);
end


end

