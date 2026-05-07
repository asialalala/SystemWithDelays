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
        
        function f_handle = get_f_ode1()
            f_handle = @(t, x_curr, x_del, dx_del) x_curr + x_del - 1/4*dx_del; 
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
      

        % Equation 1.2 (Neutral DDE)
        % Model: $x'(t) = x(t) + x(t - \tau) + c \cdot x'(t - \tau)$
        % Parameters: $\tau = 1, c = -2$
        % Initial condition: $\phi(t) = -t$
        % Solution (Only in range $t \in [0,2]$):
        % For $t \in [0, 1]: x(t) = -2 + t + 2 \cdot e^t$
        % For $t \in (1, 2]: x(t) = 4 - t + 2 \cdot e^t - 2(t + 1) \cdot e^{t-1}$
        function f_handle = get_f_ode2()
            f_handle = @(t, x_curr, x_del, dx_del) x_curr + x_del - 2 * dx_del;
        end

        function phi_handle = get_phi2()
            phi_handle = @(t) -t;
        end
        
        function sol = sol2(t)
            sol = zeros(size(t));
            
            for i = 1:length(t)
                ti = t(i);
                if (ti >= 0 && ti <= 1)
                    sol(i) = -2 + ti + 2*exp(ti);
                elseif (ti > 1 && ti <= 2)
                    sol(i) = 4 - ti + 2*exp(ti) - 2*(ti + 1)*exp(ti - 1);
                else
                    sol(i) = NaN; 
                end
            end
        end 

        function handler = get_sol2()
            handler = @(t) exampleNddeFunctions.sol2(t);
        end

        % Equation 1.3 (Neutral DDE)
        % Model: $x'(t) = x(t) + x(t - \tau) + c \cdot x'(t - \tau) + \sin(t)$
        % Parameters: $\tau = 1, c = -1/4$
        % Initial condition: $\phi(t) = -t$
        % Solution (Only in range $t \in [0,2]$):
        % For $t \in [0, 1]: x(t) = -1/4 + t + 3/4 \cdot e^t - 1/2 \cdot \cos(t) - 1/2 \cdot \sin(t)$
        function f_handle = get_f_ode3(h)
            tau = 1;
            Ntau = round(tau/h);
            Ndelta = 1;
            delta = Ndelta*h;
            f_handle = @(idxT, x) x(idxT) + x(idxT - Ntau) ...
                - 1/4 * (x(idxT - Ntau) - x(idxT - Ntau - Ndelta))/delta + sin(idxT*h);
        end

        function phi_handle = get_phi3()
            phi_handle = @(t) -t;
        end
        
        function sol = sol3(t)
            sol = zeros(size(t));
     
            for i = 1:length(t)
                ti = t(i);
                if (ti >= 0 && ti <= 1)
                    sol(i) = -1/4 + ti + 3/4*exp(ti) - 1/2*cos(ti) - 1/2*sin(ti);
                elseif (ti > 1 && ti <= 2)
                    sol(i) = 1/2 - ti + 3/4*exp(ti) + 3/16*(3*ti + 1)*exp(ti - 1) ...
                             + 1/2*(cos(ti - 1) - cos(ti) - sin(ti)) + 1/8*sin(ti - 1);
                else
                    sol(i) = NaN; 
                end
            end
        end 

        function handler = get_sol3()
            handler = @(t) exampleNddeFunctions.sol3(t);
        end

        % TODO - implement more functions
    end
end
