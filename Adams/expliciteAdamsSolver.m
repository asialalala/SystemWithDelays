function [time, solution] = expliciteAdamsSolver(k, h, tk, f_ode, tau, phi)
% Rozwiazuje przekazane rownanie rozniczkowe metoda explicite adams
%   k - rzad metody
%   h - krok (jako czas), co jaki czas obliczany jest kolejny punkt
%   tk - czas koñcowy, do którego obliczane jest rozwiazanie
%   f_ode - funkcja opisujaca rownanie rozniczkowe
%   tau - opoznienie wystêpujace w równaniu ró¿niczkowym
%   phi - funkcja opisujaca warunek poczatkowy na przedziale [-tau, 0]
%
%   Zwraca solution - rozwiazanie od czasu 0 do tk

    % 1. Parametry metody i generowanie wzoru
    [adamsFunc, formulaSym] = getFormula(k);

    % 2. Inicjalizacja siatki i pamiêci
    Ntau = round(tau/h); % Liczba próbek z opó¼nienia
    N = round(tk/h);     % Docelowa liczba próbek rozwiazania   
    Nx = N + Ntau + 1;   % Liczba wszystkich próbek (rozwi±zanie od 0 + próbki opó¼nione + 0)

    x = zeros(1, Nx);
    t = (0:Nx-1)*h - tau; % Os czasu przesuniêta o tau
    
    % 3. Warunek poczatkowy (historia)
    tetaSpan = -tau : h : 0;
    x(1:Ntau+1) = phi(tetaSpan);

    % 4. G£ÓWNA PÊTLA OBLICZENIOWA
    for n = Ntau+1 : Nx-1 % Start od t=0 (probki Ntau+1), obliczamy wartosc dla probki n+1 probki

        % Sprawdzamy, czy mamy wystarczajaco du¿o probek wstecz dla rzêdu k
        % Jesli nie (rozruch), u¿yj ni¿szego rzêdu  (Euler)
        k_eff = min(k, n - Ntau); 

        if k_eff < k 
            % ROZRUCH Metoda eulera (k = 1)
            fn = f_ode(n, x(1:n));
            x(n+1) = x(n) + h * fn;
        else
            % Przygotowanie wektora wartosci pochodnych: [f(n), f(n-1), ..., f(n-k+1)]
            f_vals = zeros(1, k); % Dla k-tego rzedu zawsze k pochodnych
            for i = 0:k-1
                idx = n - i;
                f_vals(i+1) = f_ode(idx, x(1:idx));
            end

            % Wywo³anie: adamsFunc(h, y_n, f1, f2, ...)
            args = num2cell([h, x(n), f_vals]);
            x(n+1) = adamsFunc(args{:});
        end
    end
    % Zwróæ solution od 0 do tk
    solution = x(Ntau+1:end);
    time = t(Ntau+1:end);
end
