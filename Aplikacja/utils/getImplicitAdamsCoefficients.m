function gamma = getImpliciteCoefficients(j)
% Calculate values of alpha in Implicit Adams Method
%   j - number of coefficients - 1, corresponded to the method's order k
%
%  Returns gammas, as coefficients vector.

gamma = sym(zeros(1, j+1));
gamma(1) = 1;

for m = 2:(j+1)
    sum = 0;
    for  i = 1:m
        sum = sum + gamma(m - i + 1)/(i);
    end
    gamma(m) =  - sum;
end
end