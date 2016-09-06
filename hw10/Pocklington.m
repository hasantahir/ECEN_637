function [Yin] = Pocklington(TotalSegments, L, a, f)

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
        
        fun = @(zp) exp(-1i*K*sqrt(( z(j) - zp).^2 + a^2))...
            ./sqrt(( z(j) - zp).^2 + a^2); % Make a symbolic Function of z_prime
        % This calculation is based on the reference cited in the code
        % introduction. 
        % Mathematical Modeling of Eq. 4.63 
        zp_upper = z(k) + deltaz/2; % Upper limit of integration
        dz_upper = z(j) - zp_upper; % Represents (z(m) - z(n)) for the upper limit
        zp_lower = z(k) - deltaz/2; % Lower limit of integration
        dz_lower = z(j) - zp_lower; % Represents (z(m) - z(n)) for the limit limit
        
        %% Numerical Integration Using Gauss_Quadratures
        int_part =  quadgk(fun,zp_lower,zp_upper,'RelTol',1e-8,'AbsTol',1e-12);
        
        R_upper = sqrt(dz_upper^2 + a^2); % Represents R for the upper limit
        R_lower = sqrt(dz_lower^2 + a^2); % Represents R for the lower limit
        sum_part_upper = (dz_upper)*(1 + 1i*K*R_upper)...
            ./R_upper^3*exp(-1i*K*R_upper); % Exact Evalation of the second term (upper limit) in Eq. 4.63 of the reference
        sum_part_lower = (dz_lower)*(1 + 1i*K*R_lower)...
            ./R_lower^3*exp(-1i*K*R_lower); % Exact Evalation of the second term (lower limit) in Eq. 4.63 of the reference
        Z(j,k) = K^2*int_part + sum_part_upper - sum_part_lower; % Pocklington Integral for j not equal to k

        if j == k
        % This is done to avoid very small numbers in the denominator
            % For all the diagonal elements of Z, terms only need to be calculated once
            R_upper = sqrt(  (-deltaz/2)^2 + a^2); % Here z(m) = z(n) 
            R_lower = sqrt((deltaz/2)^2 + a^2);
            sum_part_upper = (-deltaz/2)*(1 + 1i*K*R_upper)...
                ./R_upper^3*exp(-1i*K*R_upper);
            sum_part_lower = (deltaz/2)*(1 + 1i*K*R_lower)...
                ./R_lower^3*exp(-1i*K*R_lower);
            num = sqrt(1 + 4*a^2/deltaz^2) + 1;
            denom = sqrt(1 + 4*a^2/deltaz^2) - 1;
            self = (log(num/denom) - 1i*K*deltaz); % Approximation of the integral term as described in reference
            Z(j,k) = K^2*self + sum_part_upper - sum_part_lower; % Diagonal Elements of the impedance matrix
         end
    end
end

Z = toeplitz(real(Z(1,:))) + 1i*toeplitz(imag(Z(1,:))); % Make a Toeplitz matrix out of a row vector
V = zeros(TotalElements,1); % Initializr Source (RHS)
V(floor(TotalSegments/2)+1) = -1i*4*pi*omega*eps0*(1.0/deltaz);   % Delta Source
%% Calculate the current
I =Z\V;
%% Input impedance and admittance
Zin = 1.0 / I(floor(TotalSegments/2)+1);
Yin = 1./Zin;
end

