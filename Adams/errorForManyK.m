for k = 10:19
    % 1. Definicja parametrów modelu
    h = 0.01;
    tk = 5;
    
    % 2. Definicja pd
    % roblemu
    tau = 1;
    f_ode = @(t, x, xTau) xTau + x + 1; 
    phi = @(teta) -exp(teta);
    
    % 3. Rozwiazanie Metoda Explicite Adams
    [t_adams, x_adams] = expliciteAdamsSolver(k, h, tk, f_ode, tau, phi);
    lags = [tau];
    
    % 4. Rozwizanie Solverem Matlaba
    tspan = [0 tk];
    ddefun = @(t, x, Z) Z(:,1) + x + 1;
    sol_matlab = dde23(ddefun, lags, phi, tspan);
    
    % 5. Ewaluacja dde23 w punktach t_adams
    x_matlab = deval(sol_matlab, t_adams);
    
    % 7. Obliczenie błędu bezwzględnego
    error_val = abs(x_adams - x_matlab);
    
    % 8. Statystyki błędu
    max_err = max(error_val);
    rmse = sqrt(mean(error_val.^2));
    
    fprintf('Maksymalny błąd bezwzględny dla k=%d: %e\n', k, max_err);
    fprintf('Błąd RMSE dla k=%d: %e\n', k, rmse);
    
    % --- ZAPIS DO PLIKU CSV ---
    % Nagłówki: Rząd, Max_Error, RMSE
    csv_filename = 'Results2/bledy.csv';
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
    plot(t_adams, x_matlab, 'k-', 'LineWidth', 2); hold on;
    plot(t_adams, x_adams, 'r--', 'LineWidth', 1.5);
    grid on;
    legend('dde23 (Reference)', ['Adams k=', num2str(k)], 'Location', 'best');
    title(['Porównanie rozwiązań k=', num2str(k)]);
    xlabel('t'); ylabel('x(t)');
    
    % Wykres błędu
    subplot(1, 2, 2);
    semilogy(t_adams, error_val + eps, 'b', 'LineWidth', 1.5);
    grid on;
    title('Błąd bezwzględny (log)');
    xlabel('t'); ylabel('|x_{adams} - x_{dde}|');
    
    % Zapis wykresu do pliku PNG
    img_name = sprintf('Results2/wykres_k%d.png', k);
    saveas(fig, img_name);
    
    % Zamknięcie figury, aby nie zaśmiecać pamięci przy wielu iteracjach
    close(fig);
end