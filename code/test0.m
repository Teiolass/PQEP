hold on
set(gca,'Xscale','log', 'Yscale', 'log')

outx = [];
outy = [];

for iii = 1:100
    build_quick;

    disp('built');

    [lam, vec] = sda_3(H0, H1);

    disp('sda completed');

    k = size(lam)(1);
    for ii = 1:k
        mu = lam(ii,ii);
        v = vec(:,ii);
        num = norm(mu*mu*H1.'*v + mu*H0*v + H1*v);
        nn = norm(H1, 'fro');
        mod = abs(mu);
        den = mod*mod*nn + mod*norm(H0, 'fro') + nn;
        den = den * norm(v);

        rres = num/den;
        outx = [outx mod];
        outy = [outy rres];
    end
end
scatter(outx, outy, 'r', 'x');
