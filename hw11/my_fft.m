function dataout = my_fft(datain, nn, isign)

data = zeros(1,2*nn);
dataout = zeros(1,nn);
% data = real(data);
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

