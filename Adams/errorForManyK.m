clear all;
close all;

% 1. Model parameters definition
h = 0.01;
tk = 6;

A = 3;
B = 2;
C = 1;
% 2. Problem definition
phi = exampleRddeFunctions.get_phi1(C);
f_ode = exampleRddeFunctions.get_f_ode1(h, A, B, C);
f_sol = exampleRddeFunctions.get_sol1(A, B, C);
tau = B;
t_span = 0:h:tk;

if ~exist('Results', 'dir'), mkdir('Results'); end

% Plik zbiorczy na błędy (nadpisywany przy nowym uruchomieniu)
summary_filename = 'Results/summary_errors.csv';
header = {'k', 'Method', 'Max_Error', 'RMSE'};
writecell(header, summary_filename);

for k = 2:3

    % exact solution
    x_exact = f_sol(t_span); % Reference

    % 3. solution with Explicite Adams Method and Start
    [t_adamsS, x_adamsS] = expliciteAdamsSolver(k, h, tk, f_ode, tau, phi);
    error_valS = abs(x_adamsS - x_exact);

    % 4. solution with Explicite Adams Method without Start
    [t_adamsWS, x_adamsWS] = expliciteAdamsWithoutStartSolver(k, h, tk, f_ode, tau, phi);
    error_valWS = abs(x_adamsWS - x_exact);

    
    % 4. solution with Explicite Adams Method New Approachwithout Start
    [t_adamsWSNA, x_adamsWSNA] = expliciteAdamsWithoutStartSolverNA(k, h, tk, f_ode, tau, phi);
    error_valWSNA = abs(x_adamsWSNA - x_exact);

    % 6. ddends (Matlab) solution
    % sol_matlab = ddensd(@ddefun,@dely,@delyp,@history,[0,tk]);
    ddefun = @(t, y, lags) A * lags;
    sol_matlab = dde23(ddefun, [B], C, t_span);
    t_matlab = t_span;
    x_matlab = deval(sol_matlab,t_matlab);
    error_valM = abs(x_matlab - x_exact);


   % --- ZAPIS SZCZEGÓŁOWY DLA DANEGO K ---
    csv_filename = sprintf('Results/wynik_k%d.csv', k);
    T = table(t_span(:), x_exact(:), x_adamsS(:), x_adamsWS(:), x_adamsWSNA(:), x_matlab(:), ...
        error_valS(:), error_valWS(:), error_valWSNA(:), error_valM(:), ...
        'VariableNames', {'Time', 'Exact', 'Adams_S', 'Adams_WS', 'Adams_WSNA', 'ddends', ...
                          'Err_S', 'Err_WS', 'Err_WSNA', 'Err_ddends'});
    writetable(T, csv_filename);

    fprintf('Wyniki zostały zapisane do pliku: %s\n', csv_filename);

    % --- ZAPIS ZBIORCZY BŁĘDÓW ---
    methods = {'Adams_Start', 'Adams_WithoutStart', 'ddends'};
    errors = {error_valS, error_valWS, error_valM};
    
    for m = 1:3
        m_err = max(errors{m});
        m_rmse = sqrt(mean(errors{m}.^2));
        row = {k, methods{m}, m_err, m_rmse};
        writecell(row, summary_filename, 'WriteMode', 'append');
    end


    % --- WIZUALIZACJA I ZAPIS DO PNG ---
    fig = figure('Position', [100, 100, 1000, 400], 'Visible', 'off'); ...
        % 'Visible', 'off' jeśli nie chcesz wyświetlać okien
    set(fig, 'Color', 'w');

    % Subplot 1: Porównanie rozwiązań
    mSize = 4;
    ax1 = subplot(1, 2, 1);
    plot(t_span, x_exact, 'k-', 'LineWidth', 1.5); hold on;
    plot(t_span, x_adamsS, 'r-o', 'MarkerSize', mSize, 'LineWidth', 0.8);
    plot(t_span, x_adamsWS, 'g-s', 'MarkerSize', mSize, 'LineWidth', 0.8);
    plot(t_span, x_adamsWSNA, 'b-d', 'MarkerSize', mSize, 'LineWidth', 0.8);
    plot(t_span, x_matlab, 'm-x', 'MarkerSize', mSize, 'LineWidth', 0.8);

    grid on;
    set(ax1, 'XColor', 'k', 'YColor', 'k', 'Color', 'w'); 
    legend('Exact', 'Adams Start', 'Adams Without Start', 'Adams WSNA', ...
        'ddends', 'Location', 'best', 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
    title(['Comparison of Methods (k=', num2str(k), ')'], 'Color', 'k');
    xlabel('t', 'Color', 'k'); ylabel('x(t)', 'Color', 'k');

    ax2 = subplot(1, 2, 2); 
    semilogy(t_span, error_valS, 'r-o', 'MarkerSize', mSize); hold on;
    semilogy(t_span, error_valWS, 'g-s', 'MarkerSize', mSize);
    semilogy(t_span, error_valWSNA, 'b-d', 'MarkerSize', mSize);
    semilogy(t_span, error_valM, 'm-x', 'MarkerSize', mSize);

    grid on;
    set(ax2, 'XColor', 'k', 'YColor', 'k', 'Color', 'w'); 
    legend('Err Adams S', 'Err Adams WS', 'Err Adams WSNA', ... 
        'Err ddends', 'Location', 'best', 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
    title('Absolute Error (log scale)', 'Color', 'k');
    xlabel('t', 'Color', 'k'); ylabel('Error', 'Color', 'k');

    % Dodatkowe upewnienie się, że napisy na osiach (tyknięcia) są czarne
    set([ax1, ax2], 'GridColor', [0.8 0.8 0.8], 'MinorGridColor', ... 
        [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k');

    % Zapis wykresu do pliku PNG
    img_name = sprintf('Results/wykres_k%d.png', k);
    saveas(fig, img_name);

    % Zamknięcie figury, aby nie zaśmiecać pamięci przy wielu iteracjach
    close(fig);
end