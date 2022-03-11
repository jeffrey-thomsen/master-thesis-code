function[P,F] = p0MultipleMatrix(p0,N,M)

    for n=1:N
        for m=1:M
            P(n,m) = n*p0/m;
        end
    end
    F = 1./P;

end