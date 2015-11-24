% Hasan Ta h
% clear;
a = 4;
b = 3;
M = 3;
N = 3;
k = zeros(M,N);
row = zeros(M*N,1);
col = zeros(M*N,1);
for i = 1 : M
    for j = 1 : N
        k(i,j) = sqrt( (i*pi/a)^2 + (j*pi/b)^2 );
    end
end
% Sort the modes in ascending order and make a vector out of it
sorted = reshape(sort(k(:)),[],1); 
for i = 1 : length(sorted)
    [row(i),col(i)] = find(k==sorted(i));
end
Mode_order = [row,col];
  
