function XX1 =  doubling(A, Q, rtol=1e-16)
% returns a solution of the equation
% X + A.' X^-1 A = Q
    AA1 = A;
    XX1 = Q;
    YY1 = zeros(size(A));
    while 1
        AA = AA1;
        XX = XX1;
        YY = YY1;
        [L,U,P] = lu(XX - YY);
        SOL_AA = U \ (L \ (P*AA));
        AA1 = AA * SOL_AA;
        XX1 = XX - AA.' * SOL_AA;
        YY1 = YY + AA * (U \ (L \ P*AA.'));
        err = norm(XX1-XX) / norm(XX);
        if  err < rtol
            break;
        end
    end
end
