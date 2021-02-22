% build_quick;

[lam, vecr, vecl] = sda_fast(H0, H1);
na = norm(H1, 'fro');
nq = norm(H0, 'fro');
outx = [];
outy = [];

for ii = 1:k
        mu = 1/lam(ii);
        v = vecl(:,ii);
        num = norm(mu*mu*H1.'*v + mu*H0*v + H1*v);
        mod = abs(mu);
        den = mod*mod*na + mod*nq + na;
        den = den * norm(v);
        rres = num/den;
        outx = [outx mod];
        outy = [outy rres];

        mu = lam(ii);
        v = vecr(:,ii);
        num = norm(mu*mu*H1.'*v + mu*H0*v + H1*v);
        mod = abs(mu);
        den = mod*mod*na + mod*nq + na;
        den = den * norm(v);
        rres = num/den;
        outx = [outx mod];
        outy = [outy rres];
end

scatter(outx, outy, 'b', '^');
set(gca,'Xscale','log', 'Yscale', 'log');
