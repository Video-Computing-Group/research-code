function c = nchoosek(v,k)
%NCHOOSEK Binomial coefficient or all combinations.
%   NCHOOSEK(N,K) where N and K are non-negative integers returns N!/K!(N-K)!.
%   This is the number of combinations of N things taken K at a time.
%   When a coefficient is large, a warning will be produced indicating 
%   possible inexact results. In such cases, the result is only accurate 
%   to 15 digits for double-precision inputs, or 8 digits for single-precision
%   inputs.
%
%   NCHOOSEK(V,K) where V is a vector of length N, produces a matrix 
%   with N!/K!(N-K)! rows and K columns. Each row of the result has K of 
%   the elements in the vector V. This syntax is only practical for 
%   situations where N is less than about 15.
%
%   Class support for inputs N,K,V:
%      float: double, single
%
%   See also PERMS.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.21.4.11 $  $Date: 2010/11/17 11:29:57 $

if ~isscalar(k) || k < 0 || ~isreal(k) || k ~= round(k)
  error(message('MATLAB:nchoosek:InvalidArg2'));
end

[m, n] = size(v);

if min(m,n) ~= 1
   error(message('MATLAB:nchoosek:InvalidArg1'));
end

% the first argument is a scalar integer
if isscalar(v) && isfloat(v) && isreal(v) && v==round(v) && v >= 0  
   % if the first argument is a scalar, then, we only return the number of 
   % combinations. Not the actual combinations.
   % We use the Pascal triangle method. No overflow involved. c will be 
   % the biggest number computed in the entire routine.
   %

   % Do the computation in doubles. Determine the output type first.
   classOut = superiorfloat(v,k);

   n = double(v); % rename v to be n. the algorithm is more readable this way.
   k = double(k);

   if k > n 
     error(message('MATLAB:nchoosek:KOutOfRange')); 
   end 
   if k > n/2, k = n-k; end
   if k <= 1
      c = n^k;
   else
      if strcmp(classOut,'single')
         tolerance = 1e8;
      else % double case
         tolerance = 1e15;
      end
      nums = (n-k+1):n;
      dens = 1:k;
      nums = nums./dens;
      c = round(prod(nums));
      if c > tolerance
         warning(message('MATLAB:nchoosek:LargeCoefficient', sprintf( '%e', tolerance ), log10( tolerance )));
      end
   end
   % Convert answer back to the correct type
   c = cast(c,classOut);
   
else
   % the first argument is a vector, generate actual combinations.
   
   if n == 1, n = m; v = v.'; end
   
   if n == k
      c = v(:).';
   elseif n == k + 1
      tmp = v(:).';
      c   = tmp(ones(n,1),:);
      c(1:n+1:n*n) = [];
      c = reshape(c,n,n-1);
   elseif k == 1
      c = v.';
   elseif k == 0
      c = zeros(1,0);
   elseif n < 17 && (k > 3 || n-k < 4)
      rows = 2.^(n);
      ncycles = rows;
      
      for count = 1:n
         settings = [1 0];
         ncycles = ncycles/2;
         nreps = rows./(2*ncycles);
         settings = settings(ones(1,nreps),:);
         settings = settings(:);
         settings = settings(:,ones(1,ncycles));
         x(:,n-count+1) = settings(:);
      end
      
      idx = x(sum(x,2) == k,:);
      nrows = size(idx,1);
      [rows,~] = find(idx');
      c = reshape(v(rows),k,nrows).';
   else 
      c = [];
      if k < n && k > 1
         for idx = 1:n-k+1
            Q = combs(v(idx+1:n),k-1);
            c = [c; [v(ones(size(Q,1),1),idx) Q]];
         end
      end
      
   end
end

%----------------------------------------
function P = combs(v,m)
%COMBS  All possible combinations.
%   COMBS(1:N,M) or COMBS(V,M) where V is a row vector of length N,
%   creates a matrix with N!/((N-M)! M!) rows and M columns containing
%   all possible combinations of N elements taken M at a time.
%
%   This function is only practical for situations where M is less
%   than about 15.

if nargin~=2, error(message('MATLAB:nchoosek:WrongInputNum')); end

v = v(:).'; % Make sure v is a row vector.
n = length(v);
if n == m
   P = v;
elseif m == 1
   P = v.';
else
   P = [];
   if m < n && m > 1
      for k = 1:n-m+1
         Q = combs(v(k+1:n),m-1);
         P = [P; [v(ones(size(Q,1),1),k) Q]];
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
