build_quick;
disp('built');

[lam, vecr, vecl] = sda_3_left(H0, H1);
disp('sda completed');
[vec2, lam2] = polyeig(H1, H0, H1.');
disp('polyeig completed');

% lam = [lam; 1./lam];
% vec = [vecr vecl];
vec = vecr;

vec2 = vec2.';
nro = size(vec2)(2); % dovrebbe essere 2*k
nro = nro+1;
vec2(:, nro) = abs(lam2);
vec2(:, nro+1) = lam2;
vec2 = sortrows(vec2, nro).';
lam2 = vec2(end,:);
vec2(end,:) = [];
vec2(end,:) = [];
vec2(:,end/2+1:end) = [];
lam2(end/2+1:end) = [];

vec = vec.';
nro = size(vec)(2); % dovrebbe essere 2*k
nro = nro+1;
vec(:, nro) = abs(lam);
vec(:, nro+1) = lam;
vec = sortrows(vec, nro).';
lam = vec(end,:);
vec(end,:) = [];
vec(end,:) = [];

err_vec = norm(vecnorm(vec-vec2)./vecnorm(vec))
err_lam = norm(abs(lam-lam2)./abs(lam2))

return;

outx = [];
outy = [];
outx2 = [];
outy2 = [];
na = norm(H1, 'fro');
nq = norm(H0, 'fro');
for ii = 1:2*k
    mu = lam(ii);
    v = vec(:,ii);
    num = norm(mu*mu*H1.'*v + mu*H0*v + H1*v);
    mod = abs(mu);
    den = mod*mod*na + mod*nq + na;
    den = den * norm(v);
    rres = num/den;
    outx = [outx mod];
    outy = [outy rres];
end
for ii = 1:2*k
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
hold on;
plt1 = scatter(outx, outy, 'b', 'x');
plt2 = scatter(outx2, outy2, 'r', 'x');
legend([plt1; plt2], ['sda';'polyieg']);
set(gca,'Xscale','log', 'Yscale', 'log');
