function [lam, vecr, vecl] = sda_3_left(H0, H1)
% [lam, vecr, vecl] = sda_3_left(H0, H1)
    tic;
    PHI = doubling(H1, H0); 
    disp('ok doubling');
    [vecr, lam, vecl] = eig(H1, -PHI, 'vector');
    vecl = conj(vecl);
    kk = size(vecl)(2);
    for ii = 1:kk
        vecl(:, ii) = (PHI + lam(ii) * H1) \ (PHI * vecl(:, ii));
    end
    toc;
end
