clear all;
close all;

% 1. Definicja zmiennej s jako obiektu Transfer Function
s = tf('s'); 

% Przykład bardziej złożony (dwa szeregowo połączone człony)
Td = 2 * pi;
% G2 = (10 / (s + 1)) * (3 / (s + 5)) * exp(-s * Td);
 G2 = (1 - exp(-s * Td))/(2 * pi * s);
% G2 = 1/(1+exp(-s * Td));
% G2 = 1/(1-exp(-s * Td));
% G2 = 1/(s+exp(-s * Td));
% G2 = 1/(s-exp(-s * Td));

% 3. Sprawdzenie wyniku i wykres Bode
disp('Transmitancja G2:');
disp(G2);
figure; 
bode(G2);
grid on;
title('Wykres Bode dla G2(s)');


figure; 
nyquist(G2);
grid on;
title('Wykres Nyquista dla G2(s)');