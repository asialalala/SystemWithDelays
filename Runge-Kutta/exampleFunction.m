classdef exampleFunction
    methods (Static)
        function f_handle = get_f_ode1()
            f_handle = @(x, y) x+y;
        end

        function phi_handle = get_phi1()
            phi_handle = @(t) 0;
        end

        function sol = sol1(t)
            sol = exp(t) - t - 1;
        end 

        function handler = get_sol1()
            handler = @(t) exampleFunction.sol1(t);
        end
    end
end