function [lam, vec] = full_spec(H0, H1, m)
    k = size(H0)(1);
    [lam0, vecr0, vecl0] = sda_3_left(H0, H1);
    lam = zeros(1, 2*k); 
    vec = zeros(m*k, 2*k);
    tmp = (0:m-1)';
    for jj = 1:k
        lam(jj) = lam0(jj)^m;
        lam(k+jj) = 1/lam(jj);
        tmp2 = lam0(jj) .^ tmp;
        vec(:, jj) = kron(tmp2, vecr0(:, jj));
        tmp2 = 1./ tmp2;
        vec(:, k+jj) = kron(tmp2, vecl0(:, jj));
    end
end
