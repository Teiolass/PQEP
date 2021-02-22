omegas = [100; 1000; 3000; 5000];
set(gca, 'Yscale', 'log')
hold on
plots = [];
for iii = 1:4
    omega = omegas(iii);
    mx = [];
    my = []
    for times = 1:10
        build_quick;
        [PHI, errs] = doubling(H1, H0);
        mx = [mx 1:size(errs)(2)];
        my = [my errs];
    end
    plots = [plots; plot(mx, my)];
end

legend(plots, num2str(omegas))
