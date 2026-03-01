function [numHandle, formulaSym] = getFormulaNA(k, h, tau)
% Symbolically determines the equation for a step using the explicit Adams method
% of a given order k.
%   k - order of the method
%
%   Returns numHandle, which is a function handle that can be
%   used for further numerical calculations, and formulaSym - the symbolic formula.
    
    Ntau = round(tau/h);

    % Define fn
    f_names1 = arrayfun(@(x) sprintf('f%d', x), 0:ceil(k/2) - 1, 'UniformOutput', false);
    f_names2 = arrayfun(@(x) sprintf('f%d', x), Ntau:(Ntau + floor(k/2) - 1), 'UniformOutput', false);
    f_names = [f_names1, f_names2];

    % 1. Define symbols
    y_n = sym('y_n');
    f_vars = sym(f_names);
    
    % 2. Formulate symbolic formula
    betas = getBetas(k,h,tau); % Zbiór współczynników metody
    step_sum = sum(betas(:)' .* f_vars);
    formulaSym = simplify(y_n + h * step_sum); 

    % 3. Define numeric function
    numHandle = matlabFunction(formulaSym, 'Vars', [y_n, f_vars]);
end
