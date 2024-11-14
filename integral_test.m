close all; % 关闭所有图形窗口

% 参数设置：积分区间分段次数（必须是偶数，以满足 Simpson 公式的要求）
n_segments = 4; % 可以调整此参数控制分段次数

% 定义被积函数
f = @(x) pi * cos(pi * x);

% 定义积分区间
a = -1;
b = 1;

% 定义离散点
x_vals = linspace(a, b, 100); % 主要用于绘制的采样点

% 计算原函数 F(x) = sin(pi * x)
F_exact = sin(pi * x_vals);

% 初始化用于存储每种方法的积分结果
F_integral = zeros(size(x_vals));
F_matlab_trapz = zeros(size(x_vals));
F_custom_trapz = zeros(size(x_vals));
F_cumtrapz = zeros(size(x_vals));
F_simpson = zeros(size(x_vals));
F_segmented_simpson = zeros(size(x_vals));
F_segmented_trapz = zeros(size(x_vals));
F_gauss_legendre = zeros(size(x_vals));

% 循环计算每个点的积分
for i = 1:length(x_vals)
    xi = x_vals(i);
    sub_x = linspace(0, xi, n_segments + 1); % 将 [0, xi] 分为 n_segments 等分
    sub_f = f(sub_x);
    h = (xi - 0) / n_segments; % 子区间的宽度

    % MATLAB 的 integral 方法
    F_integral(i) = integral(f, 0, xi); 

    % 梯形方法
    F_custom_trapz(i) = h * (sum(sub_f) - (sub_f(1) + sub_f(end)) / 2); % 自定义梯形公式
    F_matlab_trapz(i) = trapz(sub_x, sub_f);    % MATLAB 自带 trapz 方法
    F_cumtrapz_vals = cumtrapz(sub_x, sub_f);   % MATLAB 的 cumtrapz 方法
    F_cumtrapz(i) = F_cumtrapz_vals(end);       % 取 cumtrapz 的最终积分值
    F_segmented_trapz(i) = 0;  % 分段复合梯形
    for j = 1:n_segments
        F_segmented_trapz(i) = F_segmented_trapz(i) + (h / 2) * (sub_f(j) + sub_f(j+1));
    end

    % Simpson 方法
    F_simpson(i) = (h / 3) * (sub_f(1) + 4 * sum(sub_f(2:2:end-1)) + 2 * sum(sub_f(3:2:end-2)) + sub_f(end)); % 自定义 Simpson
    F_segmented_simpson(i) = 0;  % 分段复合 Simpson
    for j = 1:2:n_segments-1
        F_segmented_simpson(i) = F_segmented_simpson(i) + (h/3) * (sub_f(j) + 4*sub_f(j+1) + sub_f(j+2));
    end

    % 高斯-勒让德积分（quadgk）
    F_gauss_legendre(i) = quadgk(f, 0, xi);
end

% 计算每种方法相对于原函数的误差（使用二范数）
error_integral = norm(F_integral - F_exact, 2);
error_matlab_trapz = norm(F_matlab_trapz - F_exact, 2);
error_custom_trapz = norm(F_custom_trapz - F_exact, 2);
error_cumtrapz = norm(F_cumtrapz - F_exact, 2);
error_simpson = norm(F_simpson - F_exact, 2);
error_segmented_simpson = norm(F_segmented_simpson - F_exact, 2);
error_segmented_trapz = norm(F_segmented_trapz - F_exact, 2);
error_gauss_legendre = norm(F_gauss_legendre - F_exact, 2);

% 分类输出误差
fprintf('--- MATLAB 的 integral 方法 ---\n');
fprintf('Integral 方法误差（二范数）：%.10f\n', error_integral);

fprintf('\n--- 梯形方法 ---\n');
fprintf('自定义梯形公式误差（二范数）：%.10f\n', error_custom_trapz);
fprintf('MATLAB trapz 方法误差（二范数）：%.10f\n', error_matlab_trapz);
fprintf('MATLAB cumtrapz 方法误差（二范数）：%.10f\n', error_cumtrapz);
fprintf('分段复合梯形方法误差（二范数）：%.10f\n', error_segmented_trapz);

fprintf('\n--- Simpson 方法 ---\n');
fprintf('自定义 Simpson 方法误差（二范数）：%.10f\n', error_simpson);
fprintf('分段复合 Simpson 方法误差（二范数）：%.10f\n', error_segmented_simpson);

fprintf('\n--- MATLAB 的 高斯-勒让德方法 ---\n');
fprintf('Gauss-Legendre 方法误差（二范数）：%.10f\n', error_gauss_legendre);

% 绘制数值积分结果与原函数的比较
figure;
plot(x_vals, F_exact, 'k-', 'LineWidth', 2); % 原函数
hold on;
plot(x_vals, F_integral, 'b-', 'LineWidth', 1.5); % integral 方法
plot(x_vals, F_custom_trapz, 'g--', 'LineWidth', 1.5); % 自定义梯形公式
plot(x_vals, F_matlab_trapz, 'y-.', 'LineWidth', 1.5); % MATLAB 自带 trapz 方法
plot(x_vals, F_cumtrapz, 'c-', 'LineWidth', 1.5); % MATLAB cumtrapz 方法
plot(x_vals, F_segmented_trapz, 'b--', 'LineWidth', 1.5); % 分段复合梯形
plot(x_vals, F_simpson, 'm-.', 'LineWidth', 1.5); % Simpson 方法
plot(x_vals, F_segmented_simpson, 'r-', 'LineWidth', 1.5); % 分段 Simpson
plot(x_vals, F_gauss_legendre, 'c:', 'LineWidth', 1.5); % Gauss-Legendre 方法 (quadgk)

title('数值积分方法与原函数 y = \sin(\pi \cdot x) 的比较');
xlabel('x');
ylabel('F(x)');
legend('原函数 y = \sin(\pi \cdot x)', 'integral 方法', '自定义梯形公式', ...
       'MATLAB trapz 方法', 'MATLAB cumtrapz 方法', '分段复合梯形', ...
       'Simpson 方法', '分段 Simpson', 'Gauss-Legendre 方法', 'Location', 'Best');
grid on;
hold off;

% 绘制误差比较（用于可视化每种方法的相对误差）
figure;
plot(x_vals, abs(F_integral - F_exact), 'b-', 'LineWidth', 1.5);
hold on;
plot(x_vals, abs(F_custom_trapz - F_exact), 'g--', 'LineWidth', 1.5);
plot(x_vals, abs(F_matlab_trapz - F_exact), 'y-.', 'LineWidth', 1.5);
plot(x_vals, abs(F_cumtrapz - F_exact), 'c-', 'LineWidth', 1.5);
plot(x_vals, abs(F_segmented_trapz - F_exact), 'b--', 'LineWidth', 1.5);
plot(x_vals, abs(F_simpson - F_exact), 'm-.', 'LineWidth', 1.5);
plot(x_vals, abs(F_segmented_simpson - F_exact), 'r-', 'LineWidth', 1.5);
plot(x_vals, abs(F_gauss_legendre - F_exact), 'c:', 'LineWidth', 1.5);

title('数值积分方法的误差比较');
xlabel('x');
ylabel('误差 |F_{numerical} - F_{exact}|');
legend('integral 方法误差', '自定义梯形公式误差', 'MATLAB trapz 方法误差', ...
       'MATLAB cumtrapz 方法误差', '分段复合梯形方法误差', 'Simpson 方法误差', ...
       '分段 Simpson 方法误差', 'Gauss-Legendre 方法误差', 'Location', 'Best');
grid on;
hold off;
