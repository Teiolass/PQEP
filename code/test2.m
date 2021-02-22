% Riscrivo questo test
% Mi serve trovare il valore di gamma cup che e` definito come 
% gamma_cup = rho ( Phi^-1 H1)
% con rho raggio spettrale. Forse per calcoli precisi piu` che usare
% abs(eigs(inv(PHI)*H1, 1)) che risolve inv(PHI)*H1*v = lam*v, meglio
% H1*v=PHI*lam*v, ie avere abs(eigs(-H1, PHI), 1)
% In realta` non cambia molto, il problema non era quello.

for ii = 1:5
    clear;
    build_quick;
    PHI = doubling(H1, H0);
    disp('ok doubling');
    gamma_cup = abs(eigs(inv(PHI)*H1, 1))
end

