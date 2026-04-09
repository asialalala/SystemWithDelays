t0 = linspace(-1,0, 10);
t1 = linspace(0,1, 10);
t2 = linspace(1,2, 10);
t3 = linspace(2,3, 10);
t4 = linspace(3,4, 10);

f0 = @(t) zeros(size(t))
f1 = @(t) ones(size(t));
f2 =  @(t) -t+2;
f3 = @(t) 1/2*t.^2 - 3.*t + 4;
f4 = @(t) -1/6.*t.^3 + 2*t.^2 -15/2*t + 17/2;

y0 = f0(t0);
y1 = f1(t1);
y2 = f2(t2);
y3 = f3(t3);
y4 = f4(t4);

pts_t = [0, 1, 2, 3];
pts_y = [1, 1, 0, -0.5];

% concatenate y1 and y2
y = [y0, y1, y2, y3, y4];
t = [t0, t1, t2, t3, t4];
plot(t, y, "k", 'LineWidth', 0.8)
hold on;
grid on;
scatter(pts_t, pts_y, "k")
xlabel("t")
ylabel("x")

