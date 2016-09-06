function [Zin,A] = PocklingtonPulseSolution(TotalSegments, NumGauss, WireLength, ...
                            WireRadius, Frequency, ExcitationType, ThetaIn, NumTheta)

% PocklingtonPulseSolution calculates the solution to Pocklington's Integral Equation 
% for a straight wire along z-axis.
%
% CALL:  [Zin, Current, ERadiated] = PocklingtonPulseSolution(TotalSegments, ...
%           NumGauss, WireLength, WireRadius, Frequency, ExcitationType, Theta, NumTheta)
% INPUTS:
%  TotalSegments = total number of thin wire segments to use
%  NumGauss = number of quadrature points for numerical integration
%  WireLength = total wire length (meters)
%  WireRadius = wire radius (meters)
%  Frequency = operating frequency (Hz)
%  ExcitationType = The wire excitation type to use
%       ExcitationType = 0 : Delta-gap source located in center segment
%       ExcitationType = 1 : Theta-polarized incident plane wave
%  ThetaIn = Incident angle of plane wave if ExcitationType = 1 (radians)
%  NumTheta = Number of scattering angles (0 <= theta <= 2*pi)
%
% OUTPUTS:
%  Zin = Input impedance measured at center segment if ExcitationType = 0 (Ohms)
%  Current = Complex-valued current in each segment (Amperes)
%  ERadiated = Complex-valued radiated far-field over NumTheta angles (0 <= theta <= 2*pi)
 
% Reference:
%
% [1] Gibson, Walton C. "The Method of Moments in Electromagnetics," Taylor
% and Francis/CRC, 2008.
%
% Copyright (C) 2007 Tripoint Industries, Inc.
 
                         
                                                    
c = 299792458;                       
Mu = 4.0e-7*pi;
Epsilon = 1.0/(Mu*c*c);
w = 2.0*pi*Frequency;
k = w*sqrt(Mu*Epsilon);
Eta = sqrt(Mu/Epsilon);

TotalElements = TotalSegments;

DeltaZ = WireLength / TotalSegments;

z = linspace(-0.5*WireLength , 0.5*WireLength - DeltaZ, TotalSegments);

% zero matrix
A = zeros(TotalElements);

[GaussZ, GaussW] = lgwt(NumGauss, 0.0, 1);  % Gauss integration weights and locations
GaussZ = GaussZ*DeltaZ;
GaussW = GaussW*DeltaZ;

% Compute the first part of the self term (Eqn. 4.63)
r2 = WireRadius*WireRadius;
num = sqrt(1 + 4*r2/(DeltaZ*DeltaZ)) + 1;
denom = sqrt(1 + 4*r2/(DeltaZ*DeltaZ)) - 1;
self = k*k*log(num/denom) - i*k*k*k*DeltaZ;

% Compute the second part of the self term (Eqn. 4.63)
limit2 = PocklingtonAnalyticTerm(0, DeltaZ/2.0, k, WireRadius);
limit1 = PocklingtonAnalyticTerm(0, -DeltaZ/2.0, k, WireRadius);

% Total self term
self = self + (limit2 - limit1);


% compute first row of matrix only
m = 1;
z_m = z(m) + 0.5*DeltaZ;
for n = 1:TotalElements

    if m == n
        A(m,n) = self;
%         A(m,n) = 0;
    else
        term = 0.0;
        % compute numerical part of Pocklington matrix
        for i_gauss = 1:NumGauss
            z_n = z(n) + GaussZ(i_gauss);
            dz = z_m - z_n;
            dz2 = dz*dz;
            r = sqrt(dz2 + r2); % distance between surface point and interior wire point
            r_inv = 1.0/r;
            expterm = exp(-i*k*r)*r_inv;
            term = term + expterm*GaussW(i_gauss);
        end
        limit2 = PocklingtonAnalyticTerm(z_m, z(n) + DeltaZ, k, WireRadius);
        limit1 = PocklingtonAnalyticTerm(z_m, z(n), k, WireRadius);
        A(m,n) = k*k*term + (limit2 - limit1);

    end
end
 
% construct  rest of matrix using first row
A = toeplitz(real(A(1,:))) + i*toeplitz(imag(A(1,:)));

% invert matrix
Ainv = inv(A);

rhs = zeros(TotalElements,1);
j = zeros(TotalElements,1);

% compute the right-hand side
if ExcitationType == 0
    rhs(floor(TotalSegments/2)+1) = -i*4*pi*w*Epsilon*(1.0/DeltaZ);   % delta-gap excitation
elseif ExcitationType == 1
    for m = 1:TotalElements
        z_m = z(m) + 0.5*DeltaZ;
        rhs(m) = PlanewaveExcitation(z_m, k, ThetaIn);
    end
    rhs = -i*4*pi*w*Epsilon*v;
end

% solve for currents
Current = Ainv*rhs;

% array of theta for theta-pol radiated field
Theta = linspace(0.0, 2.0*pi, NumTheta);
for iTheta = 1:NumTheta
    cosTheta = cos(Theta(iTheta));
    sinTheta = sin(Theta(iTheta));
    ERadiated(iTheta) = 0.0;
    for m = 1:TotalElements
        z_m = z(m) + 0.5*DeltaZ;
        ERadiated(iTheta) = ERadiated(iTheta) + DeltaZ*Current(m)*sinTheta*exp(i*k*z_m*cosTheta);
    end
    ERadiated(iTheta) = -(i*w*Mu/(4.0*pi))*ERadiated(iTheta);
end


Zin = 1.0 / Current(floor(TotalSegments/2)+1);

return


function b = PlanewaveExcitation(z, k, Theta)

b = sin(Theta)*exp(i*k*z*cos(Theta));

return


% This function computes the Pockylinton analytic term following Eqn. (4.63)
function x = PocklingtonAnalyticTerm(zm, zp, k, WireRadius)

dz = zm - zp;
r = sqrt(dz*dz + WireRadius*WireRadius);
r_inv = 1.0/r;

x = dz*(1 + i*k*r)*r_inv*r_inv*r_inv*exp(-i*k*r);

return;


