outx = [];
outy = [];
outx2 = [];
outy2 = [];

for iii = 1:1
    build;
    disp('built');

    [lam, vec] = full_spec(H0, H1, m);
    disp('sda completed');
    [vec2, lam2] = polyeig(A, Q, A.');
    disp('polyeig completed');

    fin = and(isfinite(lam2),(lam2 != 0));
    lam2 = lam2(fin);
    vec2 = vec2(:, fin).';
    nro = size(vec2)(2); % dovrebbe essere 2*k
    nro = nro+1;
    vec2(:, nro) = abs(lam2);
    vec2 = sortrows(vec2, nro).';
    %vec2(end,:) = [];

    vec = vec.';
    nro = size(vec)(2); % dovrebbe essere 2*k
    nro = nro+1;
    vec(:, nro) = abs(lam);
    vec = sortrows(vec, nro).';
    %vec(end,:) = [];
    vecnorm(abs(vec-vec2)/abs(vec)).'
    %abs(vec-vec2)/abs(vec)
end
