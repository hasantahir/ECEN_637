N = 16384;
Nx = 800;
freq=3.0e9;
omega = (2.0*pi)*freq;
c=3.0e8;
lambda=c/freq;
kappa=2.0*pi/lambda;
Length = 38*lambda;
epsi_o=8.854e-12;
epsi_r=69.0;
sigma=6.5;
epsi_c=(epsi_o)*(epsi_r)+1i*(sigma/omega);
za=31;
dx=200;
a=6.37e6;
a_e=a*(4.0/3.0);
theta_max=13.0;
beamwidth=22.2718567e-3;
p_max=kappa*sin(theta_max*pi/180.0);
dp=2.*p_max/(N-1);
dz=1.0/(2.0*(p_max/(2.0*pi)));
z_max=((N-1)/2.0)*dz;
pe=kappa*sin((beamwidth/2.0));
for  i =1: N
    z(i)=dz*(i-1)-z_max;
    p(i)=dp*(i-1)-p_max;
end
for i = 1 : Nx
    xx(i)=(i)*dx;
end
save test.mat
