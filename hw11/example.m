n = [0:29];
x1 = cos(2*pi*n/10); % 3 periods
x2 = [x1 x1]; % 6 periods
x3 = [x1 x1 x1]; % 9 periods
N = 2048;
X1 = abs(fft(x1,N));
X2 = abs(fft(x2,N));
X3 = abs(fft(x3,N));
F = [0:N-1]/N;
figure(1)
subplot(3,1,1)
plot(F,X1),title('3 periods'),axis([0 1 0 50])
subplot(3,1,2)
plot(F,X2),title('6 periods'),axis([0 1 0 50])
subplot(3,1,3)
plot(F,X3),title('9 periods'),axis([0 1 0 50])


n = [0:29];
x1 = cos(2*pi*n/10); % 3 periods
x2 = [x1 x1]; % 6 periods
x3 = [x1 x1 x1]; % 9 periods
N = 30;
N2 = 2*N;
N3 = 3*N;
X1 = abs(my_fft(x1,N,1));
X2 = abs(my_fft(x2,N2,1));
X3 = abs(my_fft(x3,N3,1));
F1 = [0:N-1]/N;
F2 = [0:N2-1]/N2;
F3 = [0:N3-1]/N3;
figure(2)
subplot(3,1,1)
plot(F1,X1),title('3 periods')
subplot(3,1,2)
plot(F2,X2),title('6 periods')
subplot(3,1,3)
plot(F3,X3),title('9 periods')