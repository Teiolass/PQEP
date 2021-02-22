outx2 = [];
outy2 = [];

for iii = 1:1
    build_quick;
    disp('built');

    [vec2, lam2] = polyeig(H1, H0, H1.');
    disp('polyeig completed');

    na = norm(H1, 'fro');
    nq = norm(H0, 'fro');
    for ii = 1:size(lam2)(1)
        mu = lam2(ii);
        v = vec2(:,ii);
        num = norm(mu*mu*H1.'*v + mu*H0*v + H1*v);
        mod = abs(mu);
        den = mod*mod*na + mod*nq + na;
        den = den * norm(v);

        rres = num/den;
        outx2 = [outx2 mod];
        outy2 = [outy2 rres];
    end
end
hold on;
scatter(outx2, outy2, 'r', 'x');
set(gca,'Xscale','log', 'Yscale', 'log');
