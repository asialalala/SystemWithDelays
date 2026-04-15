classdef ImplicitAdamsSolverTest < matlab.unittest.TestCase
    
    properties
        % Tolerance for comparative testing
        TolSmall = 1e-1;
        TolHigh =  1e-4;
    end
    
    methods (Test)
        
        function testLinearOde(testCase)
            % Test 1: Wimple equation without delay (dx/dt = x)
            % Solution: x(t) = exp(t)
            k = 2; h = 0.01; tk = 1;
            f_ode = @(t, x, xtau, dxtau) x;
            phi = @(t) exp(t);
            Xtau = 0; dXtau = 0;
            
            [time, sol] = implicitAdamsSolver(k, h, tk, f_ode, Xtau, dXtau, phi);
            
            analytical = exp(time);
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolSmall, ...
                'Solver can not handle a simple exponential equation - small tolerance.');
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolHigh, ...
                'Solver can not handle a simple exponential equation - high tolerance.');
        end
        
        function testDdeDelay(testCase)
            % Test 2: Ewuation with a delay (DDE): dx/dt = -x(t-1)
            % for t in [0, 1] and phi(t)=1, solution to x(t) = 1 - t
            k = 1; h = 0.001; tk = 1;
            Xtau = 1; dXtau = 0;
            f_ode = @(t, x, xtau, dxtau) -xtau;
            phi = @(t) ones(size(t));
            
            [time, sol] = implicitAdamsSolver(k, h, tk, f_ode, Xtau, dXtau, phi);
            
            analytical = 1 - time;
            testCase.verifyEqual(sol, analytical, 'AbsTol', testCase.TolSmall, ...
                'Error in handle with delay (Xtau) - small tolerance.');
            testCase.verifyEqual(sol, analytical, 'AbsTol', testCase.TolHigh, ...
                'Error in handle with delay (Xtau) - high tolerance.');
        end
        
        function testNeutralDde(testCase)
            % Test 3: Neutral equation (NDDE): dx/dt = -0.5*dx/dt(t-1)
            % Very simple test for delay element
            k = 2; h = 0.01; tk = 0.5;
            Xtau = 1; dXtau = 1;
            % f(t, x, x_tau, dx_tau) = -0.5 * dx_tau
            f_ode = @(t, x, xtau, dxtau) -0.5 * dxtau;
            phi = @(t) t; % dx/dt = 1 dla t <= 0
            
            [time, sol] = implicitAdamsSolver(k, h, tk, f_ode, Xtau, dXtau, phi);
            
            % Sprawdzenie czy rozwiązanie nie eksploduje
            testCase.verifyTrue(all(isfinite(sol)), 'Solver lost stability na NDDE.');
        end

        function testStepSizeConsistency(testCase)
            % Test 4: Check if decreasing the h decreases the error
            f_ode = @(t, x, xtau, dxtau) -x;
            phi = @(t) exp(-t);
            tk = 1; k = 2;
            
            [~, sol1] = implicitAdamsSolver(k, 0.1, tk, f_ode, 0, 0, phi);
            [~, sol2] = implicitAdamsSolver(k, 0.01, tk, f_ode, 0, 0, phi);
            
            err1 = abs(sol1(end) - exp(-tk));
            err2 = abs(sol2(end) - exp(-tk));
            
            testCase.verifyTrue(err2 < err1, 'Erro, decreasing h shoud decrease the error.');
        end
    end
end