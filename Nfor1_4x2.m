% 定义函数
f = @(x) 1 ./ (1 + 4 * x.^2);

% 定义区间
a = -5;
b = 5;

% 定义 x 的范围，用于绘制真实函数和插值结果
x_plot = linspace(a, b, 1000);
y_true = f(x_plot);

% 定义不同的 n 值
n_values = [2, 4, 6, 8, 10];

% 绘图
figure;
plot(x_plot, y_true, 'k', 'LineWidth', 1.5); % 绘制真实函数
hold on;

% 遍历不同的 n 值
for n = n_values
    % 生成 n+1 个等距节点
    x_nodes = linspace(a, b, n + 1);
    y_nodes = f(x_nodes);
    
    % 输出插值点
    fprintf('当 n = %d 时，插值点如下：\n', n);
    disp(table(x_nodes', y_nodes', 'VariableNames', {'x', 'f(x)'}));
    
    % 计算牛顿插值多项式的差商表
    divided_diff = y_nodes;
    for j = 2:n+1
        for i = n+1:-1:j
            divided_diff(i) = (divided_diff(i) - divided_diff(i-1)) / (x_nodes(i) - x_nodes(i-j+1));
        end
    end
    
    % 使用差商表进行插值
    y_interp = divided_diff(n+1) * ones(size(x_plot));
    for k = n:-1:1
        y_interp = y_interp .* (x_plot - x_nodes(k)) + divided_diff(k);
    end
    
    % 绘制插值结果
    plot(x_plot, y_interp, 'LineWidth', 1, 'DisplayName', sprintf('n = %d', n));
end

% 设置图例和图形属性
legend('f(x)', 'n = 2', 'n = 4', 'n = 6', 'n = 8', 'n = 10', 'Location', 'Best');
xlabel('x');
ylabel('f(x)');
title('f(x) = 1 / (1 + 4x^2)的牛顿插值函数图像');
grid on;
hold off;
