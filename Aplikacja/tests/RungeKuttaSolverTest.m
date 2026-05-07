classdef RungeKuttaSolverTest < matlab.unittest.TestCase
    
    properties
        h = 0.01;
        tk = 2;

        TolSmall = 1e-1;
        TolHigh =  1e-4;
    end
    
    methods (Test)
        
        function testExponentialODE(testCase)
            % Test 1: dx/dt = x ( ODE withaout delay)
            % Solution: x(t) = exp(t)
            f_ode = @(t, x, xtau, dxtau) x;
            phi = @(t) exp(t);
            [t, sol] = rungeKuttaSolver(4, testCase.h, testCase.tk, f_ode, 0.1, 0, phi);
            
            analytical = exp(t);
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolSmall, ...
                'RK4 should solve it correctly - tolerance small.');
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolHigh, ...
                'RK4 should solve it correctly - tolerance high.');
        end
        
        function testDelayDDE(testCase)
            % Test 2: dx/dt = -x(t-1) for t w [0, 1]
            % If phi(t) = 1,  x(t) = 1 - t
            h_small = 0.001;
            f_ode = @(t, x, xtau, dxtau) -xtau;
            phi = @(t) ones(size(t));
            [t, sol] = rungeKuttaSolver(2, h_small, 1.0, f_ode, 1.0, 0, phi);
            
            analytical = 1 - t;
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolSmall, ...
                'DDE solved incorrectly - tolerance small.');
            testCase.verifyEqual(sol, analytical, 'RelTol', testCase.TolHigh, ...
                'DDE solved incorrectly - tolerance high.');
        end
        
        function testOrderConsistency(testCase)
            % Test 3: Check if higher k decreases error
            f_ode = @(t, x, xtau, dxtau) -x;
            phi = @(t) ones(size(t));
            
            [~, sol1] = rungeKuttaSolver(1, 0.1, 1.0, f_ode, 0, 0, phi);
            [~, sol4] = rungeKuttaSolver(4, 0.1, 1.0, f_ode, 0, 0, phi);
            
            err1 = abs(sol1(end) - exp(-1));
            err4 = abs(sol4(end) - exp(-1));
            
            testCase.verifyTrue(err4 < err1, 'RK4 should be more precise than RK1 (Euler) for h=0.1.');
        end
        
        function testInvalidOrder(testCase)
            % Test 4: Check error's throws for incorrect k
            f_ode = @(t, x, xtau, dxtau) x;
            phi = @(t) 1;
            testCase.verifyError(@() rungeKuttaSolver(3, 0.1, 1, f_ode, 0, 0, phi), ...
                'RKSolver:OrderNotSupported'); 
        end
    end
end