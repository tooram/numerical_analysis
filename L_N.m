% 插值范围和插值点数量的设置
a = -5; % 插值区间的左端点
b = 5;  % 插值区间的右端点
n = 10; % 插值点数量-1 （插值多项式的次数）

% 定义原函数
f = @(x) exp(x);

% 生成插值节点和函数值
x_nodes = linspace(a, b, n + 1);
y_nodes = f(x_nodes);

% 定义绘图范围
x_plot = linspace(a, b, 1000);
y_true = f(x_plot);

%% 牛顿插值法
% 计算差商表
divided_diff = y_nodes;
for j = 2:n+1
    for i = n+1:-1:j
        divided_diff(i) = (divided_diff(i) - divided_diff(i-1)) / (x_nodes(i) - x_nodes(i-j+1));
    end
end
fprintf('牛顿插值法系数（差商表的第一行）:\n');
disp(divided_diff');

% 使用差商表进行牛顿插值
y_interp_newton = divided_diff(n+1) * ones(size(x_plot));
for k = n:-1:1
    y_interp_newton = y_interp_newton .* (x_plot - x_nodes(k)) + divided_diff(k);
end

%% 拉格朗日插值法
y_interp_lagrange = zeros(size(x_plot));
for i = 1:n+1
    L = ones(size(x_plot));
    for j = 1:n+1
        if j ~= i
            L = L .* (x_plot - x_nodes(j)) / (x_nodes(i) - x_nodes(j));
        end
    end
    y_interp_lagrange = y_interp_lagrange + y_nodes(i) * L;
end

fprintf('拉格朗日插值法系数（插值点权重系数）:\n');
disp(y_nodes');

%% 绘图
figure;
plot(x_plot, y_true, 'k', 'LineWidth', 1.5); hold on;
plot(x_plot, y_interp_newton, 'r--', 'LineWidth', 1.2);
plot(x_plot, y_interp_lagrange, 'b-.', 'LineWidth', 1.2);
plot(x_nodes, y_nodes, 'ko', 'MarkerFaceColor', 'k'); % 插值点
legend('原函数', '牛顿插值', '拉格朗日插值', '插值点', 'Location', 'Best');
xlabel('x');
ylabel('f(x)');
title(sprintf('原函数与插值多项式比较 (n = %d)', n));
grid on;
hold off;
