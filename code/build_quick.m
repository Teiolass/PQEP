load new_mtx_K_M.mat;
k = 303;
m = 19;
omega = 1000;

rtol = 1e-16;

c1 = 0.8;
c2 = 0.2;

% K0 = rand(k, k);
% K0 = 0.5 * (K0 + K0.');
% K1 = rand(k, k);
% M0 = rand(k, k);
% M0 = 0.5 * (M0 + M0.');
% M1 = rand(k, k);

% filt1 = rand(k,k);
% filt2 = rand(k,k);
% filt3 = rand(k,k);
% filt4 = rand(k,k);
% H0 = H0 .* filt1;
% H1 = H1 .* filt2;
% M0 = M0 .* filt3;
% M1 = M1 .* filt4;

H1 = K1 + i * omega * (c1 * M1 + c2 * K1) - omega * omega * M1;
H0 = K0 + i * omega * (c1 * M0 + c2 * K0) - omega * omega * M0;

