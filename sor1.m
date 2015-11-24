function x= sor1(A, b, N)

%SOR(A, b, N) solve iteratively a system of linear equations whereby
%A is the coefficient matrix, and b is the right-hand side column vector.
%N is the maximum number of iterations.
%The method implemented is that of Successive Over Relaxation. 
%The starting vector is the null vector, but can be adjusted to one's needs.
%The iterative form is based on the S.O.R. transition/iteration matrix 
%Tw = inv(D-w*L)*((1-w)*D+w*U) and the constant vector cw = w*inv(D-w*L)*b.
%The optimal parameter w is calculated using the spectral raduis of the Jacobi
%transition matrix Tj = inv(D)*(L+U).
%The output is the solution vector x.

%This file follows the algorithmic guidelines given in the book 
%Numerical Analysis, 7th Ed, by Burden & Faires

%Author: Alain G. Kapitho
%Date  : Aug. 2007


n = size(A,1);
%splitting matrix A into the three matrices L, U and D
D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);

Tj = inv(D)*(L+U);				%Jacobi iteration matrix
rho_Tj = max(abs(eig(Tj)));     %spectral radius of Tj	
w = 2./(1+sqrt(1-rho_Tj^2));	%optimal overrelaxation parameter
disp('w =');disp(w);
Tw = inv(D-w*L)*((1-w)*D+w*U);	%SOR iteration matrix	
cw = w*inv(D-w*L)*b;			%constant vector needed for iterations

tol = 1e-05;
k = 1;
x = zeros(n,1);                 %starting vector

while k <= N
   x(:,k+1) = Tw*x(:,k) + cw;
   if norm(x(:,k+1)-x(:,k)) < tol
      disp('The procedure was successful')
      disp('Condition ||x^(k+1) - x^(k)|| < tol was met after k iterations')
      disp(k); disp('x = ');disp(x(:,k+1));
      break
   end
   k = k+1;
end

if norm(x(:,k+1)- x(:,k)) > tol || k > N
   disp('Maximum number of iterations reached without satisfying condition:')
   disp('||x^(k+1) - x^(k)|| < tol'); disp(tol);
   disp('Please, examine the sequence of iterates')
   disp('In case you observe convergence, then increase the maximum number of iterations')
   disp('In case of divergence, the matrix may not be diagonally dominant')
   disp(x');
end