function [ x, error_norm, iter, flag ]  = sorr ( A, x, b, w, max_it, tol )

%*****************************************************************************80
%
%% SOR solves the linear system Ax=b using the Successive Over-Relaxation Method.  
%
%  Discussion:
%
%    When the parameter W
%    is set to 1, this is equivalent to the Gauss-Seidel iteration.
%
%  Reference:
%
%    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
%      June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
%      Charles Romine, Henk van der Vorst
%    Templates for the Solution of Linear Systems: Building Blocks for 
%      Iterative Methods, 
%    SIAM Publications, 1993.
%
%  Parameters:
%
%    Input, real A(N,N), the symmetric positive definite matrix.
%
%    Input, real X(N), the initial guess vector.
%
%    Input, real B(N), the right hand side vector.
%
%    Input, real W, the relaxation scalar, between 0 and 2.
%
%    Input, integer MAX_IT, the maximum number of iterations.
%
%    Input, real TOL, an error tolerance.
%
%    Output, real X(N), the solution.
%
%    Output, real ERROR_NORM, the norm of the error.
%
%    Output, integer ITER, the number of iterations performed.
%
%    Output, integer FLAG, the return flag.
%    0 = the solution was found to within the specified tolerance.
%    1 = a satisfactory solution was not found.  The iteration limit
%        was exceeded.
%

%
%  Initialization.
%

  bnrm2 = norm ( b );

  r = b - A * x;
  error_norm = norm ( r ) / bnrm2;
  errorhist = [ ];
  errorhist(1) = error_norm;

%
%  Split the matrix.
%
    b = w * b;
    D = diag(diag(A));
    U = -triu(A,1);
    L = -tril(A,-1);
    M =(D - w*L);
    N = w*U + (1-w)*D;
%     M =  w * tril ( A, -1 ) + diag ( diag ( A ) );
%     N = -w * triu ( A,  1 ) + ( 1.0 - w ) * diag ( diag ( A ) );
  for iter = 1 : max_it

    x_1 = x;
%
%  Update the approximation.
%
%     x = M \ ( N * x + b );
      x = M \ (N*x + b);
%
%  Compute the error.
%
    error_norm = norm ( x - x_1 ) / norm ( x );
    errorhist(iter+1) = error_norm;
%
%  Check for convergence.
%
    if ( error_norm <= tol )
      flag = 0;
      break 
    end

  end
%
%  Restore the right hand side.
%
  b = b / w;

  error_norm = errorhist;

  return
end