% Examples of functions from "A Set of Functional Differential Equation"
% In numerical form

classdef exampleRddeFunctions
    methods (Static)
        % Equation 1.1
        % Model: $x'(t) = A \cdot x(t - \tau)$, where $\tau = B$
        % Initial condition: $\phi(t) = C$
        % Solution: $x(t) = C \cdot \sum_{n=0}^{\lfloor t/B \rfloor + 1} \frac{A^n (t - (n-1)B)^n}{n!}$
        function f_handle = get_f_ode1(A, B, C)
            f_handle = @(t, x_curr, x_del, dx_del) A*x_del;
        end
        
        function phi_handle = get_phi1(C)
            phi_handle = @(t) C.* ones(size(t));
        end
        
        function sol = sol1(t, A, B, C)
            sol = zeros(size(t));
            
            for i = 1:length(t)
                current_t = t(i);
                current_sum = 0;
            
                maxN = floor(current_t/B) + 1;
            
                for n = 0:maxN
                    term = (A^n * (current_t - (n-1)*B)^n) / factorial(n);
                    current_sum = current_sum + term;
                end
            
                sol(i) = current_sum*C;
            end
        end 

        function handler = get_sol1(A, B, C)
            handler = @(t) exampleRddeFunctions.sol1(t, A, B, C);
        end
        
        % Equation 1.2 
        % Model: $x'(t) = -x(t - \frac{\pi}{2})$
        % Initial condition: $\phi(t) = \sin(t)$
        % Solution: $x(t) = \sin(t)$

        function f_handle = get_f_ode2()
            f_handle = @(t, x_curr, x_del, dx_del) -x_del;
        end

        function phi_handle = get_phi2()
            phi_handle = @(t) sin(t);
        end

        function sol_handle = get_sol2()
            sol_handle = @(t)sin(t);
        end 
     
        % TODO - implement more functions
    end
end
