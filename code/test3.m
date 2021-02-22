k = 303;
m = 19;
omegas = [100 1000 3000 5000];
rtol = 1e-16;

% load new_mtx_K_M.mat;
K0 = rand(k, k);
K0 = 0.5 * (K0 + K0.');
K1 = rand(k, k);
M0 = rand(k, k);
M0 = 0.5 * (M0 + M0.');
M1 = rand(k, k);

c1 = 0.8;
c2 = 0.2;

for ii = 1:4
    omega = omegas(ii);
    H1 = K1 + i * omega * (c1 * M1 + c2 * K1) - omega * omega * M1;
    H0 = K0 + i * omega * (c1 * M0 + c2 * K0) - omega * omega * M0;
    
    PHI = doubling(H1, H0);
    R = PHI + H1.' * inv(PHI) * H1 - H0;
    raw = norm(PHI) + norm(H1)^2*norm(inv(PHI)) + norm(H0);
    err = norm(R) / raw;
    disp(strcat(num2str(omega), ': ', num2str(err)));
end

