% Lambert Function W

% Model parameters dy(t)+alfa*y(t-T)+beta(t)=0, T>0
alfa = 1;
beta = 2;
T = 1;
phi = @(teta) -exp(teta); % initial 

N = 1000; % Lambert modes used in the method

s = @(k) (1/T * lambertw(k,-alfa*T*exp(beta*T)) - beta);
epsilon = @(k,t) exp(s(k)*t);

k_vec = -N:N;
% [0, N] divided on 2N elements
t_vec = linspace(T, 0, 2*N + 1)';

% Omega initialization
Omega2 = zeros(2*N + 1, 2*N + 1);
for i = 1:length(t_vec)
    for j = 1:length(k_vec)
        Omega2(i,j) = epsilon(k_vec(j), t_vec(i));
    end
end

[K, T_grid] = meshgrid(k_vec, t_vec);
S_vals = arrayfun(s, K);
Omega = exp(S_vals .* T_grid);


Phi = phi(t_vec-T);

C = Omega\Phi;

% Calculate the response for each branch k
syms t_sym real
x_sym = 0;
for idx = 1:length(k_vec)
    Sk = s(k_vec(idx)); 
    x_sym = x_sym + C(idx) * exp(Sk * t_sym); 
end

% disp('Analityczny wzór na x(t):');
% pretty(vpa(x_sym, 4))

t_plot = linspace(0, 5, 100);
x_numeric = double(subs(x_sym, t_sym, t_plot));

figure;
plot(t_plot, real(x_numeric), 'LineWidth', 2);
grid on;
title('Rozwiązanie x(t) wyznaczone metodą Lamberta W');
xlabel('t'); ylabel('x(t)');