% 1. Definicja parametrów modelu
k = 3;
h = 0.01;
tk = 10;

% 2. Definicja problemu
dXtau = 2;
Xtau = 2;
f_ode = exampleNddeFunctions.get_f_ode1(); % f(t, x(t), x(t-tau)) = dx = x(t-tau) + x(t) + 1 UWAGA! Zmieniasz tu, zmien tez 16 linie
phi = exampleNddeFunctions.get_phi1(); % Warunek poczatkowy w postaci funkcji


% 3. Rozwiazanie Metoda Explicite Adams
[t_adams, x_adams] = StepMethodSolver(tk, f_ode, Xtau, dXtau, phi);

% 4. Rozwizanie Solveram Matlaba
sol_matlab = exampleNddeFunctions.get_sol1()

% 5. Ewaluacja dde23 w punktach t_adams
x_matlab = sol_matlab(t_adams);

% 7. Obliczenie błędu bezwzględnego
error = abs(x_adams - x_matlab);

% 8. Statystyki błędu
max_err = max(error);
rmse = sqrt(mean(error.^2)); % Błąd średniokwadratowy

fprintf('Maksymalny bład bezwzględny: %e\n', max_err);
fprintf('Bład RMSE: %e\n', rmse);

% 9. Wizualizacja porównania i błędu
figure('Position', [100, 100, 1000, 400]);

% Wykres rozwiazań
subplot(1, 2, 1);
plot(t_adams, x_matlab, 'y-', 'LineWidth', 2); hold on;
plot(t_adams, x_adams, 'r--', 'LineWidth', 1.5);
grid on;
legend('dde23 (Reference)', ['Adams k=', num2str(k)], 'Location', 'best');
title('Porownanie rozwiazan');
xlabel('t'); ylabel('x(t)');

% Wykres błędu
subplot(1, 2, 2);
semilogy(t_adams, error + eps, 'b', 'LineWidth', 1.5); % eps zapobiega log(0)
grid on;
title('Bład bezwzględny (skala logarytmiczna)');
xlabel('t'); ylabel('|x_{adams} - x_{dde}|');

