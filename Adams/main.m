% 1. Definicja parametrów modelu
k = 3;
h = 0.01;
tk = 3;

% 2. Definicja problemu
tau = 2;
f_ode = @(t, x, xTau) xTau + x + 1; % f(t, x(t), x(t-tau)) = dx = x(t-tau) + x(t) + 1 UWAGA! Zmieniasz tu, zmien tez 16 linie
phi = @(teta) -exp(teta); % Warunek poczatkowy w postaci funkcji

% 3. Rozwiazanie Metoda Explicite Adams
[t_adams, x_adams] = expliciteAdamsSolver(k, h, tk, f_ode, tau, phi);

lags = [tau];

% 4. Rozwizanie Solveram Matlaba
tspan = [0 tk];
ddefun = @(t, x, Z) Z(:,1) + x + 1;  % f(t, x(t), x(t-tau))
sol_matlab = dde23(ddefun, lags, phi, tspan);

% 5. Ewaluacja dde23 w punktach t_adams
x_matlab = deval(sol_matlab, t_adams);

% 7. Obliczenie b³êdu bezwzglêdnego
error = abs(x_adams - x_matlab);

% 8. Statystyki b³êdu
max_err = max(error);
rmse = sqrt(mean(error.^2)); % B³±d ¶redniokwadratowy

fprintf('Maksymalny b³ad bezwzglêdny: %e\n', max_err);
fprintf('B³ad RMSE: %e\n', rmse);

% 9. Wizualizacja porównania i b³êdu
figure('Position', [100, 100, 1000, 400]);

% Wykres rozwiazañ
subplot(1, 2, 1);
plot(t_adams, x_matlab, 'k-', 'LineWidth', 2); hold on;
plot(t_adams, x_adams, 'r--', 'LineWidth', 1.5);
grid on;
legend('dde23 (Reference)', ['Adams k=', num2str(k)], 'Location', 'best');
title('Porownanie rozwiazan');
xlabel('t'); ylabel('x(t)');

% Wykres b³êdu
subplot(1, 2, 2);
semilogy(t_adams, error + eps, 'b', 'LineWidth', 1.5); % eps zapobiega log(0)
grid on;
title('B³ad bezwzglêdny (skala logarytmiczna)');
xlabel('t'); ylabel('|x_{adams} - x_{dde}|');

