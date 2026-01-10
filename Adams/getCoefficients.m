function gamma = getCoefficients(j)
% Oblicza wartoœci wspó³czynników gamma dla metody Explicit Adams 
%   j - libczan wspó³czynników, odpowiadaj¹ca rzêdowi metody
%
%  Zwraca gamma, czyli wektor wspó³czynników.

gamma = zeros(1, j);
gamma(1) = 1;

for i = 2:j
    sum = 0;
    for  k = 2:i
        sum = sum + gamma(i - k + 1)/k;
    end
    gamma(i) = 1 - sum;
end
end