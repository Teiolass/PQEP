function [lam, vecr, vecl] = sda_fast(H0, H1)
    tic;
    PHI = doubling(H1, H0); 
    [G, T, U, V, vecr, vecl] = qz(H1, -PHI);
    lam = diag(G) ./ diag(T);
    vecl = conj(vecl);
    vecl = T * V' * vecl;
    kk = size(vecl)(2);
    for ii = 1:kk
        vecl(:, ii) = (-lam(ii) * G + T) \ vecl(:, ii);
    end
    vecl = V * vecl;
    toc;
end

