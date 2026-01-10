function x = solution(A,  A1, h, phi, f, t)
% A      : n×n matrix
% A1     : nxn matrix
% f      : function handle f(s) -> n×1
% t      : vector of times 1×m or mx1
% x      : solution -> n×m matrix, each column is x(t_k)
% h      : delay
% phi    : nx1 function describing [-h,0] behaviour
    
    X = @(s) (s < 0).*0 + (s == 0).*1 + (s > 0).*expm( A * (s)); % CZY TAK POWINNO BYC?

    % TODO ensure phi is corect dimention
    t = t(:).'; % Ensure t is a row vector

    m = numel(t);
    n = size(A, 1);
    x = zeros(n, m);

    for k = 1:m
        tk = t(k);
        integrandFunction1 = @(teta) X(tk - teta - h) * A1 * phi(teta);
        integrandFunction2 = @(s) X(tk-s)*f(s);
        I1 = integral(integrandFunction1, -h, 0, 'ArrayValued', true);
        I2 = integral(integrandFunction2, 0, tk, 'ArrayValued', true);

        x(:,k) = X(tk) * phi(0) + I1 + I2; 
    end
end