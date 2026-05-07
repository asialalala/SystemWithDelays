% 1. Definicja parametrów
tau = 2;

% 2. Definicja zmiennych symbolicznych
syms y(t)

% 3. Definicja warunku początkowego (phi) i równania
% Skoro phi(teta) = -exp(teta), to dla t-tau mamy:
phi_t_tau = -exp(t - tau); 

% Równanie: dy/dt = x(t-tau) + x(t) + 1
% Podstawiamy za x(t-tau) naszą funkcję phi (dla przedziału t < tau)
eq = diff(y, t) == phi_t_tau + y + 1;

% 4. Warunek początkowy w punkcie t=0
% phi(0) = -exp(0) = -1
cond = y(0) == -1;

% 5. Rozwiązanie analityczne
sol = dsolve(eq, cond);

% Wyświetlenie wyniku
disp('Rozwiązanie analityczne y(t):')
disp(sol)

% 6. Wizualizacja
figure;
fplot(sol, [0 6]) % Rysujemy do czasu tk = 6
grid on
xlabel('Czas t')
ylabel('Rozwiązanie y(t)')
title('Analityczne rozwiązanie równania (dla t < \tau)')

% 7. Konwersja na wektor (jeśli potrzebujesz danych do innych obliczeń)
t_vec = linspace(0, 6, 100);
y_vec = double(subs(sol, t, t_vec));