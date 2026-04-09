function [numHandle, formulaSym] = getImplicitAdamsFormula(k)
% Derive formula for the implicit adams method
% with order k.
%   k - rzad metody
%
%   Returns numHandle czyli ochwyt, funkcję numeryczna, która może zostać
%   wykorzystana do dalszych obliczeń i formulaSym - symboliczny wzór.

    % 1. Definicja symboli
    h = sym('h');
    y_n = sym('y_n');
    syms f(n);

    gamma =  getImpliciteCoefficients(k); % Zbiór współczynników metody

    % 2. Wyznaczenie symbolicznego wzoru
    sum = 0;
    for j = 0:k
        sum = sum + gamma(j + 1) * backwardDifference(f, j, n+1);
    end
    formulaSym = simplify(y_n + h * sum);

    % 3. Zdefiniowanie funkcji numerycznej
    f_calls = arrayfun(@(idx) f(n+1-idx), 0:k) % [f(n), f(n-1), ..., f(n-k+1)]
    f_names = sym('f', [1 k+1]) % [f1, ..., f_k]
    cleanFormula = subs(formulaSym, f_calls, f_names)
    
    numHandle = matlabFunction(cleanFormula, 'Vars', [h, y_n, f_names]);
end