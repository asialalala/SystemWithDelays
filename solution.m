function x = solution(A, phi0, f, t)
    integrandFunction = @(s) expm(A * (t-s)) * f(s);
    x = expm(A*t)*phi0 + integral(integrandFunction,0,t,'ArrayValued',true);
end