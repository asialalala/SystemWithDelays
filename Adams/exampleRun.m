fileID = fopen('Results/wzory.txt', 'w');
if fileID == -1
    error('Nie mo¿na otworzyæ pliku do zapisu.');
end

for k = 1:20
    % 1.Pobranie fukcji do obliczania korku metoda explicite adams ze wskazanym
    % rzedem k
    [adamsFunc, formula] = getFormula(k);
    
    % 2. Wyswietlenie wzoru
    fprintf('Wzór dla k=%d:\n', k);
    disp(formula);

    % 3. Zapisz do pliku
    fprintf(fileID, '--- Rz±d k = %d ---\n', k);
    fprintf(fileID, '%s\n\n', formula);
end

fclose(fileID);

fprintf('Zapisywanie zakoñczone pomy¶lnie.\n');

% 3. Definicja prametrow
% h_val = 0.1;
% y_val = 2.0;
% f_history = [0.8, 0.75]; % Przyk³adowe próbki dla n, n-1, n-2

% result = adamsFunc(h_val, y_val, f_history(1), f_history(2));
% disp(['Wynik kroku: ', num2str(result)]);