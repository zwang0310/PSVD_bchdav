function  [w, mvcput] = user_Ax(v)   
%
% Usage: [w] = user_Hx( v,  Mmat, varargin )
%
% compute matrix-vector products.
% the "Mmat" is optional, when it is omitted,
% a global matrix named "A" is necessary.
%
% 
% all matrix-vector products are performed through calls to
% this subroutine so that the mat-vect count can be accurate.
% a global veriable "MVprod" is needed to count the 
% matrix-vector products
%  
% if Mmat is a string function, in many applications one needs to 
% input more variables to @Mmat than just v. the additional 
% variables, if exist, are passed through the varargin cell array.
%
  
  global A_operator
  global MVprod       %count the number of matrix-vector products
  global MVcpu        %count the cputime for matrix-vector products
  
  mvcput = cputime;
  
  [nr,nc] = size(A_operator);
  % v is min(nr,nc)-by-k
  if nr >= nc   % calculate (A'A)v, where v is nc-by-k     
      w_temp = A_operator * v;      % w_temp is nr-by-k  right multiplication
      w = (w_temp' * A_operator)';  % w is nc-by-k       left multiplication
  else          % calculate (AA')v, where v is nr-by-k
      w_temp = (v' * A_operator)';  % w_temp is nc-by-k  left multiplication
      w = A_operator * w_temp;      % w is nr-by-k       right multiplication
  end
  
  mvcput = cputime - mvcput;
%  
% increase the global mat-vect-product count accordingly
%   
MVcpu  = MVcpu + mvcput;
MVprod = MVprod + 2*size(v,2); % count twice: a left and a right multiplication to A

  
%end function user_Ax


