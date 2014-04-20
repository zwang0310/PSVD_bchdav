function [eval, eigV, kconv, history] = bchdav2(varargin)
%
% bchdav2.m implements the block chebyshev-davidson method for 
% computing the largest eigenpairs of A'A or AA'.
%
% Usage:
%      [eval, eigV, kconv, history] = bchdav(varargin)
%
% where   bchdav(varargin)  is one of the following:
%------------------------------------------------------------------------    
% 1.  bchdav(A);
% 2.  bchdav(A, nwant);
% 3.  bchdav(A, nwant, opts);
% 4.  bchdav(Astring, dim, nwant, opts);  %when A is a script function,
%                                         %need to specify the dimension
%------------------------------------------------------------------------
% 
% A is the matrix for the eigenproblem, (dim=size(A,1)), 
%   (the A can be input as a function for matrix-vector products)
%
% nwant is the number of wanted eigen pairs,
%   by default: nwant=6.
%
% opts is a structure containing the following fields: (order of field names 
% does not matter. there is no need to specify all fields since most of the time
% the default values should suffice.
% the more critical fields that can affect performance are blk, polm, vimax)
%
%        blk -- block size, 
%               by default: blk = 3.
%      polym -- the degree of the Chebyshev polynomial; 
%               by default:  polym=20.
%      vimax -- maximum inner-restart subspace (the active subspace) dimension, 
%               by default:  vimax = max(max(5*blk, 30), ceil(nwant/4)).
%
%   do_outer -- do outer restart (if no outer restart, inner restart may gradually 
%               increase subspace dim until convergence if necessary)
%               by default:  do_outer = logical(1).
%      vomax -- maximum outer-restart subspace dimension;
%               by default:  vomax= nwant+30.  
%
%     filter -- filter method to be used;
%               currently only two Chebshev filters are implemented,
%               by default: filter=1. (the one without scaling)
%        tol -- convergence tolerance for the residual norm of each eigenpair;
%               by default:  tol =1e-8.
%      itmax -- maximum iteration number;
%               by default:  itmax= max(floor(dim/2), 300).
%      ikeep -- number of vectors to keep during inner-restart,
%               by default:  ikeep= max(floor(vimax/2), vimax-3*blk).
%
%         v0 -- the initial vectors (can be of any size);
%               by default: v0 = rand(dim,blk).
%      displ -- information display level; 
%               (<=0 --no output; 1--some output; >=2 --more output, 
%                when displ>5, some expensive debugging calculations will be performed)
%               by default: displ=1.
%     chksym -- check the symmetry of the input matrix A.
%               if chksym==1 and A is numeric, then isequal(A,A') is called.
%               the default is not to check symmetry of A. 
%      kmore -- additional number of eigen-pairs to be computed.
%               by default, kmore=3.
%        upb -- upper bound of all the eigenvalues of the input matrix A.
%               (provide this bound only when you know a good bound; otherwise, 
%                the code has an estimator to compute this upper bound.)
%    augment -- choose how many filtered vectors to keep in the basis,  
%               by default augment=1,  only the last blk filtered vectors are kept;
%               if augment=2,  the last 2*blk  filtered vectors are kept;
%               if augment=3,  the last 3*blk  filtered vectors are kept.
%
%========== Output variables:
%
%       eval:  converged eigenvalues (optional).
%
%       eigV:  converged eigenvectors (optional, but since eigenvectors are
%              always computed, not specifying this output does not save cputime).
%
%      kconv:  number of converged eigenvalues (optional).
%
%    history:  log information (optional)
%              log the following info at each iteration step:
%
%              history(:,1) -- iteration number (the current iteration step)
%              history(:,2) -- cumulative number of matrix-vector products 
%                              at each iteration, history(end,2) is the total MVprod.
%              history(:,3) -- residual norm at each iteration
%              history(:,4) -- current approximation of the wanted eigenvalues
%
%---------------------------------------------------------------------------
%
% As an example:
%
%    A = delsq(numgrid('D',90));   A = A - 1.6e-2*speye(size(A));
%    k = 10;  blk=3;  v=ones(size(A,1), blk);
%    opts = struct('vmax', k+5, 'blk', blk, 'v0', v, 'displ', 0);
%    [eval, eigV] = bchdav(A, k, opts);  
%
% will compute the k smallest eigenpairs of A, using the specified
% values in opts and the defaults for the other unspecified parameters.
%

%  
%---coded by  y.k. zhou,   yzhou@smu.edu
%   Jan 2010 
%
%(revision june 2012:  1. remove a redundant check of symmetry of A
%                      2. change some displ setting so that when disp=1, 
%                         converging histry will be displayed as iteration continues)

  %
  % use a global variable 'MVprod' to count the number of matrix-vector products.
  % this number is adjusted in the user provided mat-vect-product script 'user_Ax', 
  % therefore it is automatically incremented  whenever 'user_Ax' is called.
  % (by this the mat-vect-product count will be accurate, there is no need 
  % to manually increase the count within this code in case one may miss
  % increasing the count somewhere by accident)
  % use a global variable 'MVcpu' to store the cputime spent on mat-vect products.
  %
  global  MVprod   MVcpu
  % initialize mat-vect-product count and mat-vect-product cpu to zero
  MVprod = 0;    MVcpu = 0;

  global  filt_non_mv_cput   filt_mv_cput 
  filt_non_mv_cput = 0;   filt_mv_cput = 0;
  
  cputotal = cputime;   %init the total cputime count

  %
  % Process inputs and do error-checking
  %
  
  % if no input arguments, return help.
  if nargin == 0, help bchdav2, return, end
  
  if isnumeric(varargin{1})
    global A_operator
    A_operator = varargin{1}; 
    dim = min(size(A_operator));
%     [dim]=size(A_operator,1);
%     if (dim ~= size(A_operator,2)),
%       error('The input numeric matrix A must be a square matrix')
%     end
%     if (dim <= 300), 
%       warning('small dimension problem, use eig instead')
%       [eigV, eval]=eig(full(A_operator));
%       eval = diag(eval);
%       if (nargout >2), kconv = dim; end
%       if (nargout >3), history=[]; end
%       return
%     end
    Anumeric = 1;
  else
    A_operator = fcnchk(varargin{1});
    Anumeric = 0;
    dim = varargin{2};  % when A is a string, need to explicitly
                        % input the dimension of A
    if (~isnumeric(dim) || ~isreal(dim) || round(dim)~=dim || dim <1)
      error('A valid matrix dimension is required for the input string function')
    end
  end 
  
  if (nargin < 3-Anumeric),
    nwant = min(dim,6);
  else
    nwant = varargin{3-Anumeric};
    if (~isnumeric(nwant) || nwant < 1 || nwant > dim),
      warning('invalid # of wanted eigenpairs. use default value')
      nwant = min(dim,6);
    end
  end

  %
  % Set default values and apply existing input options:
  % default opt values will be overwritten by user input opts (if exist).
  % 
  % there are a few unexplained parameters which are mainly used to
  % output more info for comparision purpose, these papameters can be 
  % safely neglected by simply using the defaults.
  %
  blk = 3;
  opts=struct('blk', blk, 'filter', 1, 'polym', 20, 'tol', 1e-8, ...
              'vomax', nwant+30,  'do_outer', true, ...
              'vimax', max(max(5*blk, 30), ceil(nwant/4)), ...
              'adjustvimax', false, ... 
              'itmax', max(floor(dim/2), 300), 'augment', 1, ...
              'chksym', false,  'displ', 1,  'kmore', 3);

  if (nargin >= (4-Anumeric))
    user_opts = varargin{4-Anumeric};
    if ~isa(user_opts,'struct')
        error('Options must be a structure. (note bchdav does not need ''mode'')')
    else
        % overwrite default options by user input options
        opts = struct_merge(user_opts, opts);        
    end
  end  

  if (opts.chksym)
      if (Anumeric),
          if(~isequal(A_operator, A_operator')),
              error('input matrix is not symmetric/Hermitian');
          end
      end
  end
  
  % save opt values in local variables 
  blk=opts.blk;  filter=opts.filter;  polym=opts.polym;  tol=opts.tol;
  vomax=opts.vomax;  vimax=opts.vimax;  itmax=opts.itmax;  
  augment=opts.augment;  kmore=opts.kmore;  displ=opts.displ; 

  if (isfield(opts, 'v0')),
      sizev0 = size(opts.v0,2);      %fprintf('use predefined initial vectors\n');
      if (sizev0 < blk), 
          warning('*** input size(v0,2)=%i, blk=%i, augment %i random vectors\n', ...
                  sizev0, blk, blk-sizev0),
          opts.v0(:,sizev0+1:blk)=rand(dim,blk-sizev0);
      end
  else
      opts.v0 = rand(dim, blk);  sizev0 = blk;
  end

  if (opts.do_outer), 
      vomaxtmp = max( min(nwant + 6*blk, nwant+30), ceil(nwant*1.2) );
      if ( vomax < vomaxtmp ),
          fprintf('--> Warning: vomax=%i, nwant=%i, blk=%i\n',...
              vomax, nwant, blk),
          vomax = vomaxtmp;
          fprintf('--> Warnning: increase vomax to %i\n',vomax);
      end  
      if (vimax > vomax),
          fprintf('--> Warning:  (vimax > vomax)  vimax=%i, vomax=%i\n', vimax, vomax),
          vimax = max(min(6*blk, nwant), ceil(nwant/4));  %note vomax > nwant
          fprintf('--> reduce vimax to %i\n', vimax);
      end
  end
  if (vimax < 5*blk), 
      fprintf('--> Warning:  (vimax < 5*blk)  vimax=%i, blk=%i\n', vimax, blk),
      if (opts.adjustvimax), 
          vimax = 5*blk;
          fprintf('--> increase vimax to %i\n', vimax);        
%      else
%          if (3*blk>vimax), vimax=3*blk; fprintf('--> adjust vimax to %i\n', vimax); end
      end
  end
  if (vimax > vomax), vomax=vimax+2*blk;  end
  ikeep = max(floor(vimax/2), vimax-3*blk);
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Now start the main algorithm:
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Comment: when the matrix A is large, passing A explicitly as a 
  % variable in the function interface is not as efficient as using 
  % A as a global variable. 
  % In this code A is passed as a global variable when A is numeric. 
  %
  
  if (nargout > 3),  
      longlog = 1;  
  else            
      longlog = 0;
  end 
  
  %--Preallocate memory, useful if dim and vomax are large--
  V = zeros(dim, vomax); 
  W = zeros(dim, vimax); %note vimax<vomax, inner restart also saves memory
  Hn= zeros(vimax, vimax);
  eval = zeros(nwant,1);    
  resnrm=zeros(nwant,1);

  %
  % compute the upper bound of unwanted eigenvalues (up_nwb) and an
  % estimated upper bound of the spectrum (upb). upb is used in the scaled
  % filter.
  %
  [up_nwb, upb, maxritz] = lancz_high(dim);
 
  %
  % add some variables to measure computational cost
  %
  iter_total = 0;           % init the total iteration number count
  orth_cputotal = 0;        % init the total cputime for orthogonalization
  orth_flopstotal = 0;      % init the total flops for orthogonalization
  refinement_cputotal = 0;  % init the total cputime for rayleigh-ritz refinement
  nlog = 1;                 % only used if longlog = 1 (output history)
  
  kinner = 0;  kouter = 0;  % count the number of inner and outer restarts
  Hn_cpu=0;    conv_cpu=0;  % temp variables, used only for measuring cpu
  
  %-----------------------------------------------------------------------------
  % the two-indeces trick: use n and kact,  where n = kact+kconv, 
  % n is the same as ksub in the paper, note that n is NOT the dimension 'dim'
  %-----------------------------------------------------------------------------
  n = 0;        % n stores current dimension of subspace V 
  kact = 0;     % kact counts the dimension of the active subspace
  kconv = 0;    % init number of converged eigenvalues
  kconv1 = 1;   % auxiliary, kconv1 always stores kconv+1 
  
  % initialize V: V(:,1:blk) is orthonormal
  [V(:,1:blk), rt] = qr(opts.v0(:,1:blk), 0); 
  ec = 0;      % # of converged eigenpair at current step
  kinit = blk;    
  
  %==================================
  % start the main iteration
  %==================================
  
  while (iter_total <= itmax) 
      
      iter_total = iter_total +1;

      if (displ >= 4)
          fprintf('in iteration %i ,',iter_total)
          fprintf('upper bound of the unwanted eigvalues =%e\n', up_nwb);
      end 

      % get blk new vectors for filtering and augmenting the subspace
      if (ec > 0  &&  kinit + ec <= sizev0)
          if (displ > 4)
            fprintf('---> using column [%i : %i] of initial vectors\n',kinit+1, kinit+ec); 
          end
          Vtmp = [opts.v0(:,kinit+1:kinit+ec), V(:,kconv1:kconv+blk-ec)];
          kinit = kinit + ec;
      else
          Vtmp = V(:,kconv1:kconv+blk);
      end
      
      switch (filter)
          case 1 % the default chebyshev filter (non-scaled)
              Vtmp = Cheb_filter_lg(Vtmp, polym, up_nwb, augment);
          case 2 % the scaled chebyshev filter
              Vtmp = Cheb_filter_scal_lg(Vtmp, polym, up_nwb, upb, augment);
          otherwise
              error('selected filter is not available')
      end
 
      %
      % make the filtered vectors orthogonal to the current subspace V
      %
      n = kconv + kact;
      orth_cpu = cputime;
      [Vtmp, orth_flops] = DGKS_blk(V(:,1:n), Vtmp);
      orth_cputotal = orth_cputotal + (cputime - orth_cpu);  
      orth_flopstotal = orth_flopstotal + orth_flops;
      
      %
      % augment V by adding the orthogonalized vectors
      %
      kb = size(Vtmp,2);    % kb could be smaller than blk
      n1 = n+1;  kact = kact + kb;  n = kconv + kact;   
      V(:, n1:n) = Vtmp;
      
      %
      % compute new matrix-vector product.
      %
      if (Anumeric),
          W(:, kact-kb+1:kact) =  user_Ax(V(:, n1:n));
      else
          W(:, kact-kb+1:kact) =  user_Ax(V(:, n1:n), A_operator);  
      end

      
      Hn_cpu0=cputime;
      %
      % update Hn, compute only the active part (Hn does not include any deflated part)
      %
      Hn(1:kact, kact-kb+1:kact) = V(:, kconv1:n)'* W(:, kact-kb+1:kact); 
      %
      % symmetrize the diagonal block of Hn, this is very important since
      % otherwise Hn may not be numerically symmetric from last step 
      % (the orthonormality of eigenvectors of Hn will be poor without this step)
      %
      Hn(kact-kb+1:kact,kact-kb+1:kact)= ...
          (Hn(kact-kb+1:kact,kact-kb+1:kact)+(Hn(kact-kb+1:kact,kact-kb+1:kact))')/2;
      if (kact > kb),  % symmetrize Hn
          Hn(kact-kb+1:kact, 1:kact-kb) = Hn(1:kact-kb, kact-kb+1:kact)'; 
      end
      Hn_cpu = Hn_cpu + (cputime - Hn_cpu0);

      refinement_cpu = cputime;
      %
      % Rayleigh-Ritz refinement (at each ietration)
      % First compute the eigen-decomposition of the rayleigh-quotient matrix
      % (sorting is unnecessary since eig(Hn) already sorts the Ritz values). 
      % Then rotate the subspace basis.
      %
      [Eig_vec, Eig_val] = eig(Hn(1:kact,1:kact));  
      d_eig = diag(Eig_val); 
      % d_eig is in non-decreasing order      
      % sort d_eig in non-increasing order, and rearrange Eig_vec Eig_val
      % accordingly
      [d_eig, indx] = sort(d_eig, 'descend');
      Eig_vec = Eig_vec(:,indx);
      Eig_val = diag(d_eig);
      
      kold = kact;
      if (kact + blk > vimax),
          %
          % inner-restart, it can be easily controlled by the two-indeces (kact,n) trick
          %
          if (displ > 4), 
              fprintf('==> Inner-restart: vimax=%i, n=%i, kact from %i downto %i \n', ...
                      vimax, n, kact, ikeep),
          end
          kact = ikeep;     % truncate the active subspace to ikeep-dim.
          n = kconv+kact;   % should always keep n=kconv+kact
          kinner = kinner+1;  %used only for counting
      end 
      
      % subspace rotation (final step of RR refinement)
      V(:,kconv1:kconv+kact)=V(:,kconv1:kconv+kold)*Eig_vec(1:kold,1:kact);
      if (displ > 5),  %can be expensive 
          fprintf('Refinement: n=%i, kact=%i, kconv=%i,\n', n, kact, kconv),  
          orth_err = norm(V(:,kconv1:kconv+kact)'*V(:,kconv1:kconv+kact) ...
                          - eye(kact));
          if (orth_err > 1e-10),
              error(' After RR refinement: orth-error = %e\n', orth_err),
          end
      end
      W(:,1:kact)=W(:,1:kold)*Eig_vec(1:kold,1:kact);
      refinement_cputotal = refinement_cputotal + (cputime-refinement_cpu);  
      
      beta1 = max(abs(d_eig));
      %--deflation and restart may make beta1 too small, (the active subspace
      %--dim is small at the end), so use beta1 only when it is larger.     
      if (beta1 > maxritz),
          maxritz = beta1; 
      end
      tolr = tol*maxritz;     

      % test convergence     
      ec = 0;    %ec conunts # of converged eigenpair at current step

      conv_cpu0=cputime;
      % check convergence only for the smallest blk # of Ritz pairs. 
      % i.e., check the first blk Ritz values (since d_eig is in non-increasing order)
      kconv0 = kconv;
      for jc = 1 : blk  
          
          rho = d_eig(jc);
          r = W(:, jc)  - rho*V(:,kconv0+jc);   
          normr = norm(r);
          if (displ >= 4),
              fprintf (' n = %i,  rho= %e,  rnorm= %e\n', n, rho, normr);
          end
          swap = false;

          if (longlog == 1),
              history(nlog, 1) = iter_total;
              history(nlog, 2) = MVprod;
              history(nlog, 3) = normr;
              history(nlog, 4) = rho;
              nlog = nlog+1;
          end
          
          if  ( normr < tolr )
              kconv = kconv +1;   kconv1 = kconv +1;   ec=ec+1;
              if (displ >= 1),
                  fprintf ('#%i converged in %i steps, ritz(%3i)=%e, rnorm= %6.4e\n', ...
                           kconv,  iter_total, kconv, rho,  normr)
              end
              eval(kconv) = rho; 
              resnrm(kconv) = normr;
              
              %
              %%--compare with converged eigenvalues (sort in non-increasing order)
              %
              % determine to which position we should move up the converged eigenvalue
              imove = kconv - 1;
              while (imove >0),
                  if rho > eval(imove)
                      imove = imove -1; 
                  else
                      break
                  end
              end
              imove = imove+1;  %the position to move up the current converged pair
              if ( imove < kconv ),  
                  swap = true;  
                  if (displ > 3),
                      fprintf(' ==> swap %3i  upto %3i\n', kconv, imove);
                  end
                  vtmp =  V(:,kconv);
                  for i = kconv : -1 : imove+1
                      V(:,i)=V(:,i-1);  eval(i)=eval(i-1);
                  end
                  V(:,imove)=vtmp; eval(imove)=rho;
              end

              if (kconv >= nwant && ~swap && blk>1 ) || (kconv >= nwant+kmore),
                  if (displ > 1),
                      fprintf('  The converged eigenvalues and residual_norms are:\n')
                      for i = 1 : kconv
                          fprintf( '  eigval(%3i) = %11.8e,   resnrm(%3i) = %8.5e \n', ...
                                   i,  eval(i), i,  resnrm(i))
                      end
                  end
                  % change to output V= V(:, 1:kconv); later 
                  eigV = V(:, 1:kconv);
                  
                  if (displ > 0), %these output info may be useful for profiling
                      fprintf('#converged eig=%i,  #wanted eig=%i,  #iter=%i, kinner=%i, kouter=%i\n', ...
                              kconv,  nwant,  iter_total, kinner, kouter),

                      fprintf(' info of the eigenproblem and solver parameters :  dim=%i\n', dim),
                      fprintf(' polym=%i,  blk=%i,  vomax=%i,  vimax=%i,  n=%i, augment=%i, tol=%4.2e\n',...
                              polym, blk, vomax, vimax, n, augment, tol),
                      fprintf(' ORTH-cpu=%6.4e, ORTH-flops=%i,  ORTH-flops/dim=%6.4e\n',...
                              orth_cputotal, orth_flopstotal,  orth_flopstotal/dim),
                      fprintf(' mat-vect-cpu=%6.4e  #mat-vect-prod=%i,  mat-vec-cpu/#mvprod=%6.4e\n',...
                              MVcpu, MVprod, MVcpu/MVprod),
                      cputotal = cputime - cputotal; 
                      fprintf(' filt_MVcpu=%6.4e, filt_non_mv_cpu=%6.4e, refinement_cpu=%6.4e\n', ...
                              filt_mv_cput, filt_non_mv_cput, refinement_cputotal);
                      fprintf(' CPU%%: MV=%4.2f%%(filt_MV=%4.2f%%), ORTH=%4.2f%%, refinement=%4.2f%%\n', ...
                              MVcpu/cputotal*100, filt_mv_cput/cputotal*100, orth_cputotal/cputotal*100, ...
                              refinement_cputotal/cputotal*100);      
                      fprintf('       other=%4.2f%% (filt_nonMV=%4.2f%%, Hn_cpu=%4.2f%%, conv_cpu=%4.2f%%)\n', ...
                              (cputotal-MVcpu-orth_cputotal-refinement_cputotal)/cputotal*100,...
                              filt_non_mv_cput/cputotal*100, Hn_cpu/cputotal*100, conv_cpu/cputotal*100);
                      fprintf(' TOTAL CPU seconds = %e\n', cputotal),
                  end
                  return
              end
          else
              break; % exit when the first non-convergent Ritz pair is detected
          end
      end

      if (ec > 0), 
          W(:,1:kact-ec) = W(:,ec+1:kact); 
          % update the current active subspace dimension 
          kact = kact - ec;   
      end
      
      % save only the non-converged Ritz values in Hn
      Hn(1:kact,1:kact) = Eig_val(ec+1:kact+ec,ec+1:kact+ec);     
      
      %
      % update up_nwb to be the mean value of d_eig. 
      % from many practices this choice turn out to be efficient,
      % but there are other choices for updating the up_nwb.
      % (the convenience in adapting this bound without extra computation
      % shows the reamrkable advantage in integrating the Chebbyshev 
      % filtering in a Davidson-type method)
      %
      %low_nwb = median(d_eig(1:max(1, length(d_eig)-1)));
      up_nwb = median(d_eig);
     

      %
      % determine if need to do outer restart (only n need to be updated)
      %
      if (n + blk > vomax && opts.do_outer),
          nold = n;
          n = max(kconv+blk, vomax - 2*blk);
          if (displ > 4 ),
              fprintf('--> Outer-restart: n from %i downto %i, vomax=%i, vimax=%i, kact=%i\n', ...
                      nold, n, vomax, vimax, kact),
          end
          kact = n-kconv;
          kouter = kouter+1;  %used only for counting
      end
      conv_cpu = conv_cpu + (cputime - conv_cpu0);
     
   end
   
   if ( iter_total > itmax),
       %
       % the following should rarely happen unless the problem is
       % extremely difficult (highly clustered eigenvalues)
       % or the vomax or vimax is set too small
       %
       fprintf('***** itmax=%i, it_total=%i\n', itmax, iter_total)
       warning('***** bchdav.m:  Maximum iteration exceeded\n')
       fprintf('***** nwant=%i, kconv=%i, vomax=%i\n', nwant, kconv, vomax)
       warning('***** it could be that your vimax/blk, or vomax, or polym is too small')
   end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function [ynew] = Cheb_filter_lg(x, polm, high, augment)
%   
%  [ynew] = Cheb_filter_lg(x, polm, high, augment) 
%
%  Chebyshev iteration, normalized version.
%
%...input variables:
%
%        x -- the input vector to be filtered
%     polm -- polynomial degree  
%     high -- upper bound of the unwanted eigenvalues
%  augment -- by default augment=1, it can also be 2 or 3.
%
%...output:
%     ynew -- by default (augment=1), ynew = P_m(A)*x (m = polm)
%             augment=2, ynew = [P_(m-1)(A)*x, P_m(A)*x]
%             augment=3, ynew = [P_(m-2)(A)*x, P_(m-1)(A)*x, P_m(A)*x] 
%----------------------------------------------------------------------
%  This is the simplest Chebyshev filter for computing the largest 
%  eigenvalues. It does not apply internal scaling.
%  Lower bound of the spectrum is 0 since the matrix is SPSD.
%  Eigenvalues inside [0, high] are to be dampened and eigenvalues
%  greater than high are to be magnified: The interval [0, high] 
%  is mapped to [-1, 1] using f(t) = (t-c)/e, where 
%  c is the center and e is half the length of [low, high].
%

global filt_non_mv_cput
global filt_mv_cput

filt_cput = cputime;
mvcput = 0;

c = high/2;    % e = c;

% initialize the first two vectors:
% x_0 = x; 
% x_1 = (A*x-c*x)/e;
[y, mvcpu0] = user_Ax(x);
mvcput = mvcput + mvcpu0;
y = y/c - x;          

% start Chebyshev iteration: 
% x_(j+1) = (A*x_j-c*x_j)*2/e - x_(j-1),j=1,...,polm-1
for ii = 2 : polm-1
    
    [ynew, mvcpu0] = user_Ax(y);
    mvcput = mvcput + mvcpu0;
    ynew = 2 * (ynew/c - y) - x;
    
    x = y;
    y = ynew;
end

[ynew, mvcpu0] = user_Ax(y);
mvcput = mvcput + mvcpu0;
% default return (augment=1)
ynew = 2 * (ynew/c - y) - x;

filt_cput = cputime - filt_cput;
filt_mv_cput = filt_mv_cput + mvcput;
filt_non_mv_cput = filt_non_mv_cput + (filt_cput - mvcput);

if (nargin > 4)
    if (augment == 2)
        ynew = [y, ynew];
    elseif (augment == 3)
        ynew = [x, y, ynew];
    end
end

%end % end of function Cheb_filter_lg

%-------------------------------------------------------------------------
function [ynew] = Cheb_filter_scal_lg(x, polm, high, rightb, augment)
%   
%  [ynew] = Cheb_filter_scal_lg(x, polm, high, rightb, augment) 
%
%  Chebshev iteration, normalized version.
%
%...input variables:
%
%      x  -- the input vector to be filtered
%   polm  -- polynomial degree 
%   high  -- upper bound of the (unwanted) eigenvalues 
% rightb  -- a number larger than high for the purpose of scaling
%  augment -- by default augment=1, it can also be 2 or 3.
%
%...output:
%     ynew -- by default (augment=1), ynew = P_m(A)*x (m = polm)
%             augment=2, ynew = [P_(m-1)(A)*x, P_m(A)*x]
%             augment=3, ynew = [P_(m-2)(A)*x, P_(m-1)(A)*x, P_m(A)*x]
%
%------------------------------------------------------------------
%  This is the scaled Chebyshev filter for computing the largest 
%  eigenvalues. Lower bound of the spectrum is 0 since the matrix 
%  is SPSD.
%  Eigenvalues inside [0, high] are to be dampened and eigenvalues
%  greater than high are to be magnified: The interval [0, high] 
%  is mapped to [-1, 1] using f(t) = (t-c)/e, where 
%  c is the center and e is half the length of [low, high].
%  In order to enhance stability, Chebyshev polynomial C_k(x) is 
%  scaled as (1/rho_k)*C_k(x), where rho_k = C_k((rightb-c)/e).
%  In the computation, use scaling factor sigma_k = rho_(k-1)/rho_k.

global filt_non_mv_cput
global filt_mv_cput

filt_cput = cputime;
mvcput = 0;

e = high/2;
c = e;

% initialize the first scaling factor: sigma_1 = e/(rightb - c).
sigma1 = e/(rightb - c);

sigma = sigma1;
tau = 2/sigma1;

% initialize the first two vectors x,y:
% x_0 = x;
% x_1 = sigma1*(A*x-c*x)/e;
[y, mvcpu0] = user_Ax(x);
mvcput = mvcput + mvcpu0;
y = (y - c*x)*(sigma1/e);

% start (scaled) Chebyshev iteration:
% sigma_(j+1) = 1/(2/sigma1 - sigma_j)
% x_(j+1) = (A*x_j-c*x_j)*2/e*sigma_(j+1) - x_(j-1)*sigma_j*sigma_(j+1),j=1,...,polm-1
for jj = 2 : polm-1
    
    sigma_new = 1/(tau - sigma);
    [ynew, mvcpu0] = user_Ax(y);
    mvcput = mvcput + mvcpu0;
    ynew = (ynew - c*y) * (2/e * sigma_new) - x * (sigma * sigma_new);
    
    sigma = sigma_new;
    x = y;
    y = ynew;
    
end

[ynew, mvcpu0] = user_Ax(y);
mvcput = mvcput + mvcpu0;
sigma_new = 1/(tau - sigma);
ynew = (ynew - c*y) * (2/e * sigma_new) - x * (sigma * sigma_new);

filt_cput = cputime - filt_cput;
filt_mv_cput = filt_mv_cput + mvcput;
filt_non_mv_cput = filt_non_mv_cput + (filt_cput - mvcput);

if (nargin > 4)
    if (augment == 2),
        ynew = [y, ynew];
    elseif (augment == 3),
        ynew = [x, y, ynew];
    end
end

%end % end of function Cheb_filter_scal_lg

%-------------------------------------------------------------------------
function  [high, upb, maxritz] = lancz_high(n, k)
%
% Usage: [high, upb, maxritz] = lancz_high(n, k)
%
% apply k steps safeguarded Lanczos to get the upper bound of the unwanted 
% eigenvalues and the upper bound of spectrum of A.
%  
% Input: 
%        n  --- dimension
%        k  --- (optional) perform k steps of Lanczos
%               if not provided, k =6 (a relatively small k is enough)
%  
% Output:
%     high  --- estimated upper bound of the unwanted eigenvalues
%      upb  --- rough estimate of the upper bound of the spectrum of A.
%               (upb is used for the scaling purpose, it could be smaller
%               than the largest eigenvalue of A)
%


   if (nargin < 2), 
     k = 6; 
   else
     k = min(max(k,6), 12);    %do not go over 12 steps
   end 

   T = zeros(k);
   
   % start with a normal random vector to get the desired bound
   v = rand(n,1);     
   v = v/norm(v);
   
   tol = 2.5e-16;
   
   f = user_Ax(v);   
   alpha = v'*f;
   f = f - alpha*v;
   T(1,1) = alpha;
   
   isbreak = 0;
   for j = 2:k      % run k steps
       beta = norm(f);
       if (beta > tol)
           v0 = v; 
           v = f/beta;
           f = user_Ax(v);
           f = f - beta*v0;
           alpha = v'*f;
           f = f - alpha*v;
           T(j,j-1) = beta; T(j-1,j) = beta; T(j,j) = alpha;
       else
           isbreak = 1;
           break
       end
   end
   
   if (isbreak ~=1)
       [Q, D] = eig(T(1:j,1:j));
   else
       [Q, D] = eig(T(1:j-1,1:j-1));
   end
   
   if (beta > 1e-1)
       beta = beta * max(Q(end,:));
   end
   
   d = diag(D);
   maxritz = max(d);
   minritz = min(d);
   high = (maxritz + minritz)/2; % upper bound for unwanted eigenvalues
   upb = maxritz + beta;         % upper bound of spectrum(A)
   
%end % end of function lancz_high

%-------------------------------------------------------------------------
function [V, orth_flops] = DGKS_blk(X, V, polym, up_nwb, augment, orthtest)
% 
% Usage:  V = DGKS_blk( X, V ); 
%
% Apply DGKS to ortho-normalize V against X.
% V can have more than 1 column. the returned V is orthonormal to X.
%
% It is assumed that X is already ortho-normal, this is very important 
% for the projection  P = I - X*X^T  to be orthogonal.
%
% For debugging purpose, a 3rd variable can be provided to test if
% X is ortho-normal or not. when orthtest is 1, test will be performed.
%
%--y.k. zhou,  
%  Jan. 2010
%

%  if (nargin==3),
%    if (orthtest==1),
%      xorth=norm(X'*X - eye(size(X,2)));
%      if (xorth > 1e-8),
%          fprintf('--> Input orthgonality test:  ||X^t*X - I|| = %e\n', xorth),          
%          error('input X is not ortho-normal'); 
%      end
%    end
%  end
  
  epsbig = 2.22045D-16;  
  reorth=0.717;
  one=1.0D0; 
  [ndim, colv] = size(V);
  colx = size(X,2);
  orth_flops = 0; 
  vrank = 0;
  
  for k = 1 : colv
    
    nrmv = norm(V(:,k),2);   orth_flops=orth_flops+ndim;
    if (nrmv < epsbig*sqrt(ndim)), continue;  end
    
    %
    % normalize for the sake of numerical stability (important)
    %
    if (nrmv <= 2*epsbig || nrmv >= 300),
      V(1:ndim,k) = V(1:ndim,k)/nrmv;  orth_flops=orth_flops+ndim;
      nrmv = one;
    end
    
    h = [X, V(1:ndim,1:vrank)]'*V(1:ndim,k);
    V(1:ndim,k) = V(1:ndim,k) -  [X, V(1:ndim,1:vrank)]*h;
    orth_flops = orth_flops + ndim*((colx+vrank)*2 + 1);
    
    nrmproj = norm(V(1:ndim,k),2);  orth_flops=orth_flops+ndim;
    if (nrmproj > reorth*nrmv),
      %
      % pass the reorthogonalization test, no need to refine
      %
      vrank = vrank +1;
      if (abs(nrmproj - one) > epsbig),
	V(1:ndim,k)=V(1:ndim,k)/nrmproj;   orth_flops=orth_flops+ndim;
      end
      if (vrank ~= k),
	V(1:ndim,vrank)=V(1:ndim,k);
      end
%fprintf('V(:,%i) orthogonalized in step 1.\n', k)
    else
      %
      % fail the reorthogonalization test, refinement necessary
      %
      nrmv = nrmproj;      

      h = [X, V(1:ndim,1:vrank)]'*V(1:ndim,k);
      V(1:ndim,k) = V(1:ndim,k) -  [X, V(:,1:vrank)]*h;
      orth_flops = orth_flops + ndim*((colx+vrank)*2 + 1);
      
      nrmproj = norm(V(1:ndim,k),2);   orth_flops=orth_flops+ndim;     
      if (nrmproj > reorth*nrmv  && nrmproj >= sqrt(ndim)*epsbig),
      %if (nrmproj > reorth*nrmv),
	vrank = vrank +1;
	if (abs(nrmproj - one) > epsbig),
	  V(1:ndim,k)=V(1:ndim,k)/nrmproj;   orth_flops=orth_flops+ndim;
	end
	if (vrank ~= k),
	  V(1:ndim,vrank)=V(1:ndim,k);
	end
%fprintf('V(:,%i) orthogonalized in step 2.\n', k)	
      else
	%
	% fail the 2nd reorthogonalization test,
	%    V(:,k) is numerically in V(:, 1:vrank),
	% do not increase vrank, but go for the next k 
	%
%fprintf('V(:,%i) lies in the subspace. cannot orthogonalize it.\n',k)
      end
    
    end    
  end
  
  if (vrank > 0),
       V = V(:,1:vrank);
  else %recursively call DGKS_blk to make output V not a zero vector
      fprintf('DGKS: # of columns replaced by random vectors =%i\n', colv-vrank); 
      [V(:,vrank+1:colv)] = DGKS_blk([X,V(:, 1:vrank)], rand(ndim,colv-vrank));
  end
  
  if (1==0), %set it to be true only for debugging purpose
      if (nargin==3),
          if (orthtest==1),
              %orthVX = norm(V'*X);
              xorth=norm([X,V]'*[X,V] - eye(size([X,V],2)));
              if (xorth > 1e-8), 
                  fprintf('--> Output orthgo-normal error = %e\n', xorth),
                  error('output [X, V] is not ortho-normal'); 
              end
          end
      end
  end
%end % end of function DGKS_blk

%-------------------------------------------------------------------------
function [structm] = struct_merge(varargin)
% 
% Usage: [structm] = struct_merge(struct1, struct2, struct3, ...) 
%
% Merges any number of structures such that the output structm 
% will have all the fields in the listed input structures.
%
% Redundant field values will be removed. 
% Field values from the earlier input structures will overwrite
% field values from latter input structures (in case there are conflicts).
%

%
% (credit: struct_merge is based on a nice idea from David Gleich's pagerank.m)
% --y.k. zhou,    11/2009
%

structm =  varargin{1};
if ~isstruct(structm) , 
    error('the 1st input to struct_merge should be a structure');
end

for i = 2 : nargin
    struct2 = varargin{i};
    if ~isstruct(struct2) , 
       error(['the #',num2str(i),' input to struct_merge should be a structure']);
    end
    fn = fieldnames(struct2);

    for ii = 1 : length(fn)
        if (~isfield(structm, fn{ii}))
            structm.(fn{ii}) = struct2.(fn{ii});
        end
    end
end
%end % end of function struct_merge
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
   
%end %end of function bchdav2
