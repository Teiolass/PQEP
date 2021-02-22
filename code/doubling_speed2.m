load new_mtx_K_M.mat;
k = 303;
m = 4;
omegas = [100 1000 3000 5000];
rtol = 1e-16;
n = k*m;

c1 = 0.8;
c2 = 0.2;
D = eye(m);
U = diag(ones(m-1,1),  1);
L = diag(ones(m-1,1), -1);
T = zeros(m);
T(1,m)=1;

for ii = 4:4
    omega = omegas(ii);

    Kt = kron(D,K0) + kron(U,K1.') + kron(L,K1);
    Mt = kron(D,M0) + kron(U,M1.') + kron(L,M1);
    Kc = kron(T, K1);
    Mc = kron(T, M1);
    Dt = c1 * Mt + c2 * Kt;
    Dc = c1 * Mc + c2 * Kc;
    Q = Kt + i * omega * Dt - omega * omega * Mt;
    A = Kc + i * omega * Dc - omega * omega * Mc;
    Q = sparse(Q);
    A = sparse(A);

    disp('begin');

    out = [];

    AA1 = A;
    XX1 = Q;
    YY1 = zeros(n,n);
    steps = 0;
    while 1
        steps = steps + 1;
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
        disp('step');
    end
    gamma = abs(eigs(inv(XX1)*A, 1));
    disp(strcat(num2str(omega), ': ', num2str(gamma)));
    disp(strcat('steps: ', num2str(steps)));
    semilogy(out, strcat(";", num2str(omega), ";"));
    hold on;
end
hold off;
