clear;close all
L = .7/pi;
a = 7.022e-3;
TotalSegments = 51;
f = 300e6;
c = 299792458;                       
Mu = 4.0e-7*pi;
Epsilon = 1.0/(Mu*c*c);
omega = 2.0*pi*f;
K = omega*sqrt(Mu*Epsilon);
Eta = sqrt(Mu/Epsilon);
TotalElements = TotalSegments;
A = zeros(TotalElements);
NGauss = 32;
DeltaZ = L/ TotalSegments;
z = linspace(-0.5*L , 0.5*L - DeltaZ, TotalSegments);
Z = zeros(TotalSegments);
for j = 1 : length(z)
    for k = 1 : length(z)
        fun = @(zp) exp(-1i*K*sqrt(( z(j) - zp).^2 + a^2))...
            ./sqrt(( z(j) - zp).^2 + a^2);
        int_part = integral(fun,z(k)-DeltaZ/2,z(k)-DeltaZ/2);
        zp_upper = z(k) + DeltaZ/2;
        zp_lower = z(k) - DeltaZ/2;
        R_upper = sqrt((z(j) - zp_upper)^2 + a^2);
        R_lower = sqrt((z(j) - zp_lower)^2 + a^2);
        sum_part_upper = (z(j) - zp_upper)*(1 + 1i*K*R_upper)...
                           ./R_upper^3*exp(-1i*K*R_upper);
        sum_part_lower = (z(j) - zp_lower)*(1 + 1i*K*R_lower)...
                           ./R_lower^3*exp(-1i*K*R_lower);
        Z(j,k) = int_part + sum_part_upper - sum_part_lower;
    end
end

A = toeplitz(real(Z(1,:))) + 1i*toeplitz(imag(Z(1,:)));
rhs = zeros(TotalElements,1);
rhs(floor(TotalSegments/2)+1) = -1i*4*pi*omega*Epsilon*(1.0/DeltaZ);   % delta-gap excitation
Ainv = inv(A);
% solve for currents
Current =A\rhs;

