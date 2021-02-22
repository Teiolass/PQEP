k = 30;
m = 10;
omega = 1000;

rtol = 1e-16;

n = m*k;

c1 = 0.8;
c2 = 0.2;

D = eye(m);
U = diag(ones(m-1,1),  1);
L = diag(ones(m-1,1), -1);
T = zeros(m);
T(1,m)=1;

% load new_mtx_K_M;
K0 = rand(k, k);
K0 = 0.5 * (K0 + K0.');
K1 = rand(k, k);
M0 = rand(k, k);
M0 = 0.5 * (M0 + M0.');
M1 = rand(k, k);

Kt = kron(D,K0) + kron(U,K1.') + kron(L,K1);
Mt = kron(D,M0) + kron(U,M1.') + kron(L,M1);
Kc = kron(T, K1);
Mc = kron(T, M1);
Dt = c1 * Mt + c2 * Kt;
Dc = c1 * Mc + c2 * Kc;

Q = Kt + i * omega * Dt - omega * omega * Mt;
A = Kc + i * omega * Dc - omega * omega * Mc;

H1 = Q(k+1:2*k, 1:k);
H0 = Q(1:k,1:k);

% gamma = abs(eigs(inv(PHI)*H1, 1, "lm"))


