function [numHandle, formulaSym] = getExplicitAdamsFormula(k)
% Symbolically determines the equation for a step using the explicit Adams method
% of a given order k.
%   k - order of the method
%
%   Returns numHandle, which is a function handle that can be
%   used for further numerical calculations, and formulaSym - the symbolic formula.

    % 1. Define symbols
    h = sym('h');
    y_n = sym('y_n');
    f_names = arrayfun(@(x) sprintf('f%d', x), 0:k-1, 'UniformOutput', false);
    f_vars = sym(f_names);
    
    % 2. Formulate symbolic formula
    betas = getBetas(k); % Zbiór współczynników metody
    step_sum = sum(betas(:)' .* f_vars);
    formulaSym = simplify(y_n + h * step_sum); 

    % 3. Define numeric function
    numHandle = matlabFunction(formulaSym, 'Vars', [h, y_n, f_vars]);
end
