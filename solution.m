function X = solution(A, phi0, f, t)
% A      : n×n matrix
% phi0   : n×1 initial state
% f      : input function
% t      : vector of times 1×m or mx1
% X      : solution - n×m matrix, each column is x(t_k)

    phi0 = phi0(:); % Ensure phi0 is a column vector
    t = t(:).';
    
    if size(A, 1) ~= size(phi0,1)
        error("Incorrect dimentions for phi0 and A ")
    end

    m = numel(t);
    n = size(A, 1);
    X = zeros(n, m);

    for k = 1:m
        tk = t(k);
        integrandFunction = @(s) expm( A * (tk-s)) * f(s);
        I = integral(integrandFunction, 0, tk, 'ArrayValued', true);
        X(:,k) = expm(A * tk) * phi0 + I; 
    end
end