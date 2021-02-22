k = 303;
m = 19;
omegas = [100 1000 3000 5000];
rtol = 1e-16;

load new_mtx_K_M.mat;

% K0 = rand(k, k);
% K0 = 0.5 * (K0 + K0.');
% K1 = rand(k, k);
% M0 = rand(k, k);
% M0 = 0.5 * (M0 + M0.');
% M1 = rand(k, k);

c1 = 0.8;
c2 = 0.2;

for ii = 1:4
    omega = omegas(ii);
    H1 = K1 + i * omega * (c1 * M1 + c2 * K1) - omega * omega * M1;
    H0 = K0 + i * omega * (c1 * M0 + c2 * K0) - omega * omega * M0;

    out = [];

    AA1 = H1;
    XX1 = H0;
    YY1 = zeros(k,k);
    while 1
        AA = AA1;
        XX = XX1;
        YY = YY1;
        [L,U,P] = lu(XX - YY);
        SOL_AA = U \ (L \ (P*AA));
        AA1 = AA * SOL_AA;
        XX1 = XX - AA.' * SOL_AA;
        YY1 = YY + AA * (U \ (L \ P*AA.'));
        err = norm(XX1-XX) / norm(XX);
        if err > 0
            out = [out err];
        end
        if  err < rtol
            break;
        end
    end
    gamma = abs(eigs(inv(XX1)*H1, 1));
    disp(strcat(num2str(omega), ': ', num2str(gamma)));
    semilogy(out, strcat(";", num2str(omega), ";"));
    hold on;
end
hold off;
