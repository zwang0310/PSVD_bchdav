function [U, S, V] = bchdav_psvd(A, k, opts)

% A ~ U*diag(S)*V'

    [nr, nc] = size(A);
    if nr >= nc % bchdav2 computes eigenvectors of A'A, which is V
        [eval, V] = bchdav2(A, k, opts);
        S = sqrt(eval);
        U = zeros(nr,k);
        for j = 1:k
            U(:,j) = (A * V(:,j))/S(j); % u_j = A*v_j/s_j
        end
    else       % bchdav2 computes eigenvectors of AA', which is U
        [eval, U] = bchdav2(A, k, opts);
        S = sqrt(eval);
        V = zeros(nc,k);
        for j = 1:k
            V(:,j) = (U(:,j)' * A)'/S(j); % v_j = A'*u_j/s_j
        end
    end
        
end