function [numHandle, formulaSym] = getFormula(k)
% W sposób symboliczy wyznacza rownanie dla kroku metod¹ explicite adams o
% podanym rzedzie k.
%   k - rzad metody
%
%   Zwraca numHandle czyli ochwyt, funkcjê numeryczna, która mo¿e zostaæ
%   wykorzystana do dalszych obliczeñ i formulaSym - symboliczny wzór.

    % 1. Definicja symboli
    h = sym('h');
    y_n = sym('y_n');
    syms f(n);

    gamma =  getCoefficients(k); % Zbiór wspó³czynników metody

    % 2. Wyznaczenie symbolicznego wzoru
    sum = 0;
    for j = 0:k-1
        sum = sum + gamma(j + 1) * backwardDifference(f, j, n);
    end
    formulaSym = simplify(y_n + h * sum); 

    % 3. Zdefiniowanie funkcji numerycznej
    f_calls = arrayfun(@(idx) f(n-idx), 0:k-1); % [f(n), f(n-1), ..., f(n-k+1)]
    f_names = sym('f', [1 k]); % [f1, ..., f_k]
    cleanFormula = subs(formulaSym, f_calls, f_names);
    
    numHandle = matlabFunction(cleanFormula, 'Vars', [h, y_n, f_names]);
end
