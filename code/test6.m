build;
disp('built');

outx = [];
outy = [];
outy2 = [];

[lam, vec] = full_spec(H0, H1, m);
disp('spectrum found');

bar = waitbar(0, 'comuputing errors');

na = norm(A, 'fro');
nq = norm(Q, 'fro');
n1 = norm(H1, 'fro');
A = sparse(A);
Q = sparse(Q);

disp('begin calc');

for ii = 1:2*k
    mu = lam(ii);
    num = mu*(mu*A.' + Q) + A;
    num = norm(num*vec(:, ii));
    mod = abs(lam(ii));
    den = mod*mod*na + mod*nq + na;
    den = den * norm(vec(:,ii));
    den2 = mod*mod* n1 * norm(vec(1:k,ii)) + mod* nq * norm(vec(:,ii)) + n1 * norm(vec(n-k:k, ii));
    rres = num/den;
    rres2 = num/den2;
    outx = [outx mod];
    outy = [outy rres];
    outy2 = [outy2 rres2];
    waitbar(ii/2/k, bar);
end
hold on;
scatter(outx, outy, 'b', 'x');
scatter(outx, outy2, 'r', 'x');
set(gca,'Xscale','log', 'Yscale', 'log');
