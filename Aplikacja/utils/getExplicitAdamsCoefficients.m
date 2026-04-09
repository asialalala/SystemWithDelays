function gamma = getCoefficients(j)
% Calculates the values of the gamma coefficients for the Explicit Adams method.
%   j - number of coefficients, corresponding to the order of the method.
%
%  Returns gamma, which is the vector of coefficients.

gamma = sym(zeros(1, j));
gamma(1) = 1;

for i = 2:j
    sum = 0;
    for  k = 2:i
        sum = sum + gamma(i - k + 1)/k;
    end
    gamma(i) = 1 - sum;
end
end