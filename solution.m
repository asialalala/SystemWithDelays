function x = solution(A, phi0, f, t)
    term1 = expm(A*t)*phi0;
    integrand = @(s) exmp(A*(t-s)*f(s));
    term2 = integral(@(s) integrand(s), 0, t, 'ArrayValued',true);
    x = term1 + term2;
end