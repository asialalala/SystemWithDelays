tspan = [0 2];
delays = 1;
sol = dde23(@ddefun,delays,@history,tspan);
plot(sol.x,sol.y)