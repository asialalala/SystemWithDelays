% Initial condition
h = 5;
tau = [-h:0];
x = @(t) (t < 0).*0 + (t == 0).*1;



