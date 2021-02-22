function [lam, vec] = sda_3(H0, H1)
    % Ritorna lam, vec nel formato di eig (occhio che sono invertiti perche`
    % si`. Argomenti: H0 e H1
    PHI = doubling(H1, H0); 
    disp('ok doubling');
    [vec, lam] = eig(H1, -PHI);
end
