function betas = getBetas(k, h, tau)
% GETBETAS Calculates beta coefficients based on the Vandermonde matrix.
% This function supports two modes of operation based on the number of input arguments:
% 1 argument (k): Standard calculation with uniform steps.
% 3 arguments (k, h, tau): Calculation with a custom time delay (tau) and step (h).

    if nargin == 1
        % --- Logic for standard Adams Method ---
        % Define uniform steps: [0, -1, -2, ..., -(k-1)]
        s = 0:-1:-(k-1);
        
    elseif nargin == 3
        % --- Logic for dynamic step (delay) ---
        L = round(k/2);
        hl = 1;             % Local step size
        hd = round(tau/h);  % Discrete delay factor
        
        % Initialize step vector
        s = zeros(1,k);
        % Fill uniform part
        s(1: L) = 0:-hl:-(L-1)*hl;
        % Apply jump (delay)
        s(L+1) = s(L) - hd;
        % Fill the remaining part after the jump
        s(L+2:end) = (s(L+1)-hl):-hl:(-(k-2)*hl - hd);
    else
        error('Invalid number of arguments. Use getBetas(k) or getBetas(k, h, tau).');
    end

    % --- Common Calculation Core ---
    
    % Define the target vector c based on the integration of monomials
    % c_j = integral from 0 to 1 of x^j dx = 1/(j+1)
    getC = @(j) 1./(j+1);
    c = getC(0:1:(k-1));

    % Construct the transposed Vandermonde matrix
    % VT(i, j) = s(j)^(i-1)
    VT = zeros(k, k); 
    for i = 1:k
        VT(i, :) = s.^(i-1); 
    end

    % Solve the linear system VT * betas = c'
    betas = VT \ c';
end