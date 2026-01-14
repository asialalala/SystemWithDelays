function [numHandle, formulaSym] = getFormulaNewApproach(k)
% W sposób symboliczy wyznacza rownanie dla kroku metod¹ explicite adams o
% podanym rzedzie k.
%   k - rzad metody
%   tau - delay from the delay differential equation
%
%   Zwraca numHandle czyli ochwyt, funkcjê numeryczna, która mo¿e zostaæ
%   wykorzystana do dalszych obliczeñ i formulaSym - symboliczny wzór.

% 1. Definicja symboli
h = sym('h');
y_n = sym('y_n');
syms f(n);
Ntau = sym('Ntau');

gamma =  getCoefficients(k); % Zbiór wspó³czynników metody

% 2. Determine symbolic formula
sum = 0;
for j = 0:k-1
    sum = sum + gamma(j + 1) * backwardDifference(f, j, n);
end

tauSum = 0;
for j = 0:k-1
    tauSum = tauSum + gamma(j + 1) * backwardDifference(f, j, n - Ntau);
end
formulaSym = simplify(y_n + h * sum + h*tauSum);

% 3. Define numerical function
f_calls = arrayfun(@(idx) f(n-idx), 0:k-1);
f_tau_calls = arrayfun(@(idx) f(n-Ntau-idx), 0:k-1); % Poprawiony indeks dla tau

f_names = sym('f', [1 k]);
f_tau_names = sym('ft', [1 k]);

old_vars = [f_calls, f_tau_calls];
new_vars = [f_names, f_tau_names];

cleanFormula = subs(formulaSym, old_vars, new_vars);

% 4. Definition of numerical handle
all_vars = [h, y_n, f_names, f_tau_names];

numHandle = matlabFunction(cleanFormula, 'Vars', all_vars);
end
