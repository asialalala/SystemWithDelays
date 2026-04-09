hold off;
t0 = linspace(-1,0, 10);
t1 = linspace(0,1, 10);
t2 = linspace(1,2, 10);
t3 = linspace(2,3, 10);
t4 = linspace(3,4, 10);

f0 = @(t) t+1;
f1 = @(t) -t+1;
f2 = @(t) t-1;
f3 = @(t) -t+3;
f4 = @(t) t-3;

y0 = f0(t0);
y1 = f1(t1);
y2 = f2(t2);
y3 = f3(t3);
y4 = f4(t4);

pts_t = [0, 1, 2, 3];
pts_y = [1, 0, 1, 0];

% concatenate y1 and y2
y = [y0, y1, y2, y3, y4];
t = [t0, t1, t2, t3, t4];
plot(t, y, "k", 'LineWidth', 1.2)
hold on;
grid on;
scatter(pts_t, pts_y, "k")
xlabel("t")
ylabel("x")

