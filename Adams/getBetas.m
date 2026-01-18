function betas = getBetas(k)
% Calusates beta coefficiences based on the Vandermod matrix and 
%  get seps
s = 0:-1:-(k-1);

getC = @(j) 1./(j+1);
c = getC(0:1:(k-1));


% get Vandermod matrix transposition
VT = zeros(k, k); 
for i = 1:k
    VT(i, :) = s.^(i-1); 
end

betas = VT^(-1)*c';
end


