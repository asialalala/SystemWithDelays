% Examples of functions from "A Set of Functional Differential Equation"
% In numerical form
classdef exampleNddeFunctions
    methods (Static)
        % Equation 1.1 (Neutral DDE)
        % Model: $x'(t) = x(t) + x(t - \tau) + c \cdot x'(t - \tau)$
        % Here $\tau = 1$, and $x'(t - \tau)$ is approximated by finite difference:
        % $x'(t - \tau) \approx \frac{x(t - \tau) - x(t - \tau - \delta)}{\delta}$
        % Solution ( Only in range $t \in [0,2]$):
        % For $t \in [0, 1]: x(t) = -\frac{1}{4} + t + \frac{1}{4}e^t$
        % For $t \in (1, 2]: x(t) = \frac{1}{2} - t + \frac{1}{4}e^t + \frac{17}{16}e^{t-1} + \frac{3}{16}t e^{t-1}$
        
        function f_handle = get_f_ode1(h)
            tau = 1;
            Ntau = round(tau/h);
            Ndelta = 1;
            delta = Ndelta*h;
            f_handle = @(idxT, x) x(idxT) + x(idxT - Ntau) - 1/4 * (( x(idxT - Ntau) - x(idxT - Ntau - Ndelta) )/delta);
        end
        
        function phi_handle = get_phi1()
            phi_handle = @(t) -t;
        end
        
        function sol = sol1(t)
            sol = zeros(size(t));
            
            for i = 1:length(t)
                if (t(i) >= 0 && t(i) <= 1)
                    sol(i) = -1/4 + t(i) + 1/4*exp(t(i));
                elseif  (t(i) > 1 && t(i) <=2)
                    sol(i) = 1/2 - t(i) + 1/4*exp(t(i)) + 17/16*exp(t(i) - 1) + 3/16*t(i)*exp(t(i) - 1);
                else
                    % For t > 2 not supported
                    sol(i) = NaN; 
                end
            end

        end 

        function handler = get_sol1()
            handler = @(t) exampleNddeFunctions.sol1(t);
        end
      
        % TODO - implement more functions
    end
end
