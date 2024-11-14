% Newton_Cotes_Numerical_Integration.m 
% 实现不同阶数的牛顿-科特斯公式对函数 y = pi*cos(pi*x) 在区间 [0,1] 上进行数值积分
% 并与精确积分结果 y = sin(pi*x) 进行比较

% 清除环境
clear; close all; clc;

%% 定义被积函数和精确积分函数
integrand = @(x) pi * cos(pi * x);      % 被积函数 y = pi*cos(pi*x)
exact_integral = @(x) sin(pi * x);      % 精确积分结果 y = sin(pi*x)

%% 定义积分区间和绘图点
a_initial = -1;                          % 积分下限
b_initial = 1;                          % 积分上限
num_points = 200;                       % 绘图点的数量
x_plot = linspace(a_initial, b_initial, num_points);  % 绘图的 x 轴点
exact_values = exact_integral(x_plot);  % 精确积分值

%% 定义节点数量列表
n_list = [2, 4, 6, 8, 10, 12, 16, 18, 20, 22];          % 牛顿-科特斯公式的阶数

%% 初始化图形
figure('Name', 'Numerical Integration using Newton-Cotes Formulas', 'NumberTitle', 'off');
hold on;
plot(x_plot, exact_values, 'k-', 'LineWidth', 2, 'DisplayName', 'F(x) = sin(\pi x) ');
colormap(jet);
colors = lines(length(n_list));         % 为不同的n分配不同的颜色

%% 初始化误差数组
errors = zeros(length(n_list), 1);

%% 循环遍历每个n值
for idx = 1:length(n_list)
    n = n_list(idx);                    % 当前的n值
    numerical_result = zeros(size(x_plot));  % 存储数值积分结果
    error_vector = zeros(size(x_plot));      % 存储误差结果

    %% 循环遍历每个x_i点
    for k = 1:length(x_plot)
        x_i = x_plot(k);                % 当前的x_i值
        
        if x_i == 0
            numerical_integral = 0;      % 当x_i等于下限0时，积分结果为0
            exact_val = 0;
        else
            % 设置积分区间为 [0, x_i]
            a = 0;
            b = x_i;
            
            % 计算牛顿-科特斯权重
            weights = compute_newton_cotes_weights(a, b, n);
            
            % 定义节点
            x_nodes = linspace(a, b, n);
            
            % 评估被积函数在节点处的值
            f_values = integrand(x_nodes);
            
            % 计算数值积分
            numerical_integral = sum(weights .* f_values);
            
            % 计算精确积分值
            exact_val = exact_integral(x_i);
        end
        
        % 存储结果
        numerical_result(k) = numerical_integral;
        error_vector(k) = numerical_integral - exact_val;
    end
    
    %% 绘制数值积分结果曲线
    plot(x_plot, numerical_result, 'LineWidth', 1.5, 'DisplayName', sprintf('n=%d', n), 'Color', colors(idx,:));
    
    %% 计算二范数误差
    errors(idx) = norm(error_vector, 2);
end

%% 完善图形
title('牛顿-柯特斯公式实现数值积分');
xlabel('x');
ylabel('积分值');
legend('Location', 'best');
grid on;
hold off;

%% 输出二范数误差
fprintf('二范数误差:\n');
for idx = 1:length(n_list)
    fprintf('n=%d, 二范数误差: %.6e\n', n_list(idx), errors(idx));
end

%% 定义计算牛顿-科特斯权重的辅助函数
function weights = compute_newton_cotes_weights(a, b, n)
    % 计算牛顿-科特斯公式的权重
    % 输入:
    %   a - 积分下限
    %   b - 积分上限
    %   n - 节点数量
    % 输出:
    %   weights - 牛顿-科特斯权重向量

    x_nodes = linspace(a, b, n);         % 等距节点
    weights = zeros(1, n);               % 初始化权重向量

    for i = 1:n
        % 计算分母: product(x_i - x_j) for j ~= i
        denominator = 1;
        for j = 1:n
            if j ~= i
                denominator = denominator * (x_nodes(i) - x_nodes(j));
            end
        end

        % 计算分子多项式的系数: 多项式的根为x_j, j ~= i
        if i == 1
            other_nodes = x_nodes(2:n);
        elseif i == n
            other_nodes = x_nodes(1:n-1);
        else
            other_nodes = [x_nodes(1:i-1), x_nodes(i+1:n)];
        end
        numerator_coeffs = poly(other_nodes);  % 获取多项式系数 (降幂排列)

        % 计算L_i(x)的多项式系数: 分子系数 / 分母
        L_i_coeffs = numerator_coeffs / denominator;

        % 将系数从降幂转换为升幂
        L_i_coeffs = fliplr(L_i_coeffs);

        % 计算多项式在 [a, b] 上的定积分
        integral = 0;
        for k = 0:length(L_i_coeffs)-1
            integral = integral + L_i_coeffs(k+1)/(k+1) * (b^(k+1) - a^(k+1));
        end

        weights(i) = integral;               % 存储权重
    end
end
