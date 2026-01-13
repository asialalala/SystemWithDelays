clear all;
close all;

% 1. Model parameters definition
h = 0.01;
tk = 2;

A = 1;
B = 2;
C = 3;

% 2. Problem definition
phi = exampleRddeFunctions.get_phi1(C);
f_ode = exampleRddeFunctions.get_f_ode1(h, A, B, C);
f_sol = exampleRddeFunctions.get_sol1(A, B, C);
tau = B;

for k = 1:18
    
    % 3. solution with Explicite Adams Method
    [t_adams, x_adams] = expliciteNddeAdamsWithoutStartSolver(k, h, tk, f_ode, tau, phi);
    x_exact = f_sol(t_adams); % Dokładne rozwiązanie dla porównania
    
    % 7. Obliczenie błędu bezwzględnego
    error_val = abs(x_adams - x_exact);
    
    % 8. Statystyki błędu
    max_err = max(error_val);
    rmse = sqrt(mean(error_val.^2));
    
    fprintf('Maksymalny błąd bezwzględny dla k=%d: %e\n', k, max_err);
    fprintf('Błąd RMSE dla k=%d: %e\n', k, rmse);
    
    % --- ZAPIS DO PLIKU CSV ---
    % Nagłówki: Rząd, Max_Error, RMSE
    csv_filename = 'Results/bledy.csv';
    data_row = [k, max_err, rmse];
    
    % Jeśli plik nie istnieje, zapisujemy z nagłówkami (opcjonalnie)
    if ~exist(csv_filename, 'file')
        header = {'Rzad_k', 'Max_Error', 'RMSE'};
        writecell(header, csv_filename);
    end
    
    % Dopisywanie danych do pliku CSV
    writematrix(data_row, csv_filename, 'WriteMode', 'append');
    
    % --- WIZUALIZACJA I ZAPIS DO PNG ---
    fig = figure('Position', [100, 100, 1000, 400], 'Visible', 'on'); % 'Visible', 'off' jeśli nie chcesz wyświetlać okien
    
    % Wykres rozwiązań
    subplot(1, 2, 1);
    plot(t_adams, x_exact, 'b-', 'LineWidth', 2); hold on;
    plot(t_adams, x_adams, 'r--', 'LineWidth', 1.5);
    grid on;
    legend('exact sol (Reference)', ['Adams k=', num2str(k)], 'Location', 'best');
    title(['Porównanie rozwiązań k=', num2str(k)]);
    xlabel('t'); ylabel('x(t)');
    
    % Wykres błędu
    subplot(1, 2, 2);
    semilogy(t_adams, error_val + eps, 'b', 'LineWidth', 1.5);
    grid on;
    title('Błąd bezwzględny (log)');
    xlabel('t'); ylabel('|x_{adams} - x_{dde}|');
    
    % Zapis wykresu do pliku PNG
    img_name = sprintf('Results/wykres_k%d.png', k);
    saveas(fig, img_name);
    
    % Zamknięcie figury, aby nie zaśmiecać pamięci przy wielu iteracjach
    close(fig);
end