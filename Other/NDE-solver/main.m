% NDE
% czy on to może robić macierzowo???
% czy jest coś co potrafiło by to opakować w sys?

tspan = [0 2];

% define all necesarry functions for the solver as functions

sol = ddensd(@ddefun,@dely,@delyp,@history,tspan);

tn = linspace(0,2);
yn = deval(sol,tn);

plot(tn,yn);
xlim([0 pi]);
ylim([-1.2 1.2]);
xlabel('time t');
ylabel('solution y');