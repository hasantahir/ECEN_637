clear;close
[a1, dph1] = taylor1p(0.5, 90, 50, 20);
[g, ph] = array(0.5, a1, 200); % compute array gain
dbz(ph, g, 45, 40); % plot gain in dB with 40-dB scale
addcirc(3, 40, '--'); % add 3-dB grid circle
addray(90 + dph1/2, '-'); % add rays at 3-dB angles
addray(90 - dph1/2, '-');
nbar = 10;
R = 40;
d = 22.27e-3;
ph0 = 13;
N = 10;
[a,dph] = taylornb(d,ph0,N,R,nbar)