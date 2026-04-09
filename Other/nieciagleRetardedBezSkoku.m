hold off;
t0 = linspace(-1,0, 10);
t1 = linspace(0,1, 10);
t2 = linspace(1,2, 10);
t3 = linspace(2,3, 10);
t4 = linspace(3,4, 10);

f0 = @(t) ones(size(t));
f1 = @(t) 1 - t;
f2 = @(t) 0.5*t.^2 - 2*t + 1.5;
f3 = @(t) -1/6*t.^3 + 1.5*t.^2 - 4*t + 17/6;
f4 = @(t) (1/24)*t.^4 - (2/3)*t.^3 + (15/4)*t.^2 - (17/2)*t + 149/24;

y0 = f0(t0);
y1 = f1(t1);
y2 = f2(t2);
y3 = f3(t3);
y4 = f4(t4);

pts_t = [0, 1, 2, 3];
pts_y = [1, 0, -1/2, -0.166666666666];

% concatenate y1 and y2
y = [y0, y1, y2, y3, y4];
t = [t0, t1, t2, t3, t4];
plot(t, y, "k", 'LineWidth', 0.8)
hold on;
grid on;
scatter(pts_t, pts_y, "k")
xlabel("t")
ylabel("x")

