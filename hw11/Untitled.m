parabolic equation method for atmospheric propagation
Module data_modu19 implicit none
integer,parameter :: db1=selected‘rea1_kind(p=13,r=200) integer, parameter :: N=16384,Nx=800 I integer :: i, k, x ,choice, ii, jj ! integer, parameter :: N=131072,Nx=800 ! integer, parameter :: N=512,Nx=20
real(db1) :: p_max,z_max,dp,dz,pe,za,dx,pp,a,a_e,fmax,ftemp,zz complex(dbl) :: j=(0.0,1.0)
real(dbl), dimension(N) :: p, z, pat_p, f_2
complex(db1), dimension(N) :: pattern, fz, temp, m
rea1(dbl) :: freq,pi,lambda,c,kappa,theta_max,beamwidth,Length real(dbl) :: epsi_o,epsi_r,sigma,omega,index
complex(dbl) :: epsi_c
real(dbl) :: ban, de2, zmax
real(dbl), dimension(n,nx) :: Picture lreverse dimensions for Matlab plotting real(dbl), dimension(nx):: xx
end module data_module lti********~k**£i~k*******main programinkn"*i‘irti-i'iiii‘k‘ki‘kiivk‘k‘k'kiriiri‘k‘k‘k
use data_modu1e implicit none
call data
call gen_source Esource pattern = sin(p)/p or Taylor pattern where p=kappasin(theta) >>
do x=l, Nx ESpatial steps in the x direction call refr_index_N lRefraction index call p_space call z_space
do jj=l,n Picture(jj,x)=f_z(jj) !reverse dimension for Matlab plotting end do print *, "X = ",x ! Screen print of iteration end do
call plotter(p,pat_p,z,f_z)
end I**~k**************************Data********‘k***~k***************~k******
! Data for calculations and some values for constants
subroutine data use data_module implicit none
pi=4 . 0*ATAN (1 . 0)
freq=3.0e9 ! Antenna operation frequency omega=(2.0*pi)*freq c=3.0e8
lambda=c/freq ! One wavelength
kappa=2.0*pi/lambda
Length = 38*lambda ! Length of aperture in wavelengths _—'___—_——___—————-———"—‘-——b—“_
Page 1




epsi-o=8.854e-12 3 Dielectric constant for free space ep31_r=69.0 ! Dielectric constant for sea water Sigma=6.5 ! Conductivity of sea water ‘
ePS§IC=(eP3i_O)‘(epsi_r)+j'(sigma/omega) ! Complex permittivity 2a=
za=1520.0 ! Antenna position in z direction
dx=200. ! Distance between horizontal (x-direction) steps a=6.37e6 3 Radius of the earth
a_e=a*(4.0/3.0) 1 Effective radius of the earth
theta_max=l3.0 beamwidth=22.2718567e-3 lHalf—power beamwidth of Taylor pattern (20dB,401ambda)
P_max=kappa*sin(theta_max*pi/180.0) 2Maximum propagation angle
dp=2.’p_max/float(N-l) ! Delta increment in p-space dz=l.0/(2.0‘(p_max/(2.0*pi))) ! Delta increment in z—direction z_max=(float(N—1)/2.0)*dz 5 Physical area analyzed pe=kappa'sin((beamwidth/2.0)) ! Beam tilt angle
do i=1, N
2(i)=dz*(i-l)—z_max p(i)=dp'(i-1)-p_max end do
do i=1, Nx !horizontal (x—axis) profile positions xx(i)=FLOAT(i)*dx end do
end subroutine data Iifl’i‘ki'i******i*i*!i*tlitittiittittiti-A‘itttirt.ttittttilQi***tiit*t**i***
E Get the refractive index for further calculation subroutine refr_index_N
use data_module implicit none
if(x == 1) then write(*,*) "Choose the Media 1" write(*,*) "Type '1' for Free Space” write(*,*) "Type '2‘ for Standard Atmosphere"
write(*,*) "Type '3' {or Surface Duct" —— , write(*,*) "Type '4' For Graded Duct" hynb [)FLEU 5 read(*,*) choice
end if r
! for free space condition with flat earth
if(choice ==1) then ! write (*,*)"choice = 1"
do i=1, N if(z(i)<0.0) then ! Sea Water m(i)= (epsi_c)/(epsi_o)-1.0 else
index=1.0 I Flat earth condition (2*z/a = OJ m(i)=(index**2-1.0) end if end do

epsi-o=8.854e-12 3 Dielectric constant for free space ep31_r=69.0 ! Dielectric constant for sea water Sigma=6.5 ! Conductivity of sea water ‘
ePS§IC=(eP3i_O)‘(epsi_r)+j'(sigma/omega) ! Complex permittivity 2a=
za=1520.0 ! Antenna position in z direction
dx=200. ! Distance between horizontal (x-direction) steps a=6.37e6 3 Radius of the earth
a_e=a*(4.0/3.0) 1 Effective radius of the earth
theta_max=l3.0 beamwidth=22.2718567e-3 lHalf—power beamwidth of Taylor pattern (20dB,401ambda)
P_max=kappa*sin(theta_max*pi/180.0) 2Maximum propagation angle
dp=2.’p_max/float(N-l) ! Delta increment in p-space dz=l.0/(2.0‘(p_max/(2.0*pi))) ! Delta increment in z—direction z_max=(float(N—1)/2.0)*dz 5 Physical area analyzed pe=kappa'sin((beamwidth/2.0)) ! Beam tilt angle
do i=1, N
2(i)=dz*(i-l)—z_max p(i)=dp'(i-1)-p_max end do
do i=1, Nx !horizontal (x—axis) profile positions xx(i)=FLOAT(i)*dx end do
end subroutine data Iifl’i‘ki'i******i*i*!i*tlitittiittittiti-A‘itttirt.ttittttilQi***tiit*t**i***
E Get the refractive index for further calculation subroutine refr_index_N
use data_module implicit none
if(x == 1) then write(*,*) "Choose the Media 1" write(*,*) "Type '1' for Free Space” write(*,*) "Type '2‘ for Standard Atmosphere"
write(*,*) "Type '3' {or Surface Duct" —— , write(*,*) "Type '4' For Graded Duct" hynb [)FLEU 5 read(*,*) choice
end if r
! for free space condition with flat earth
if(choice ==1) then ! write (*,*)"choice = 1"
do i=1, N if(z(i)<0.0) then ! Sea Water m(i)= (epsi_c)/(epsi_o)-1.0 else
index=1.0 I Flat earth condition (2*z/a = OJ m(i)=(index**2-1.0) end if end do
______————————_-——-———————'—‘——’—_. . Page 2 E h
C:\Users\nevels\Desktop\PEM2015.f95
______.__________~__
for trosphere condition else if(choice-12) then ! write (*,*)“choice - 2"
do 1—1, N if(z(i)<0.0) then ! Sea Water m(i)= (epsi_c)/(epsi_o)-1.0 else
index-1v(300—0.0394'z(i))‘1e-6 !Standard atmosphere m(i)=(index'*2-1.0+2.0*z(i)/a) end if end do
! for surface duct condition else if(choice==3) then ! write (*,*)"choice - 3" do i=1, N if(z(i)<=0.0) then ! Sea Water m(i)= (epsi_c)/(epsi_o)-1.0 else if(z(i)>0. .and. z(i)<=37.0) then ! Surface Duct index=1+(300—0.5'z(i))‘1e-6 m(i)=((index'*2)-l.0+2.0'z(i)/a)
k=i
else index=1+(300—0-5*z(k)-0.0394'(z(i)-z(k)))'1e-6! Standard atmosphere m(i)=((index'*2)-1.0+2.0*z(i)la)
end if ‘3) '
end do ! for graded surface duct condition (see your Propagation notes S1) else if(choice==4) then ,vctQV ! write (*,*)"choice = 4" \}h“) ’ " if(xx(x) <= 40000.)then dez=(0.500-0.167)/40000.*xx(x)-0.5 zmax=(150.—37.)/40000.*xx(x)+37.0
endif do i=1, N _ if(z(i)<=0.0) then ! Sea Water *‘
m(i)= (epsi_c)/(epsi_o)—1.0
else if(z(i)>0. .and. z(i)<=zmax) then index=1.+(300.+dez*z(i))*1e-6 m(i)=((index**2)—1.0+2.0*z(i)/a) k=i
else if(z(i)>zmax) then index=1.+(300.+dez*z(k)-0.0394*(z(i)—z(k)))*1e-6! Standard atmosphere
m(i)=((index**2)—1.0+2.0*z(i)/a) end if
end do
end if end subroutine refr_index_N
! graded surface duct
*itiii**k*******i***********i****t**t**
liki*****‘ki’t*iirik‘k‘ki-kiz‘kiiriiiri'i‘k'ki‘i'tiiiii: distribution and its far field pattern
!Generate the source, 1.9. the aperture subroutine gen_source
use data_module implicit none
x=0 1 sin(p)/p pattern 1 do i=l,N ' PP=P(1)'P€
. a: h j” I 1“ PP 0 ) t en r dorm“. ,‘ W 0M“;


