clear; clc; close all;
%% 1. 自动生成数据点

rng(0);

% 生成自变量 x
n = 50;
x = linspace(0, 10, n)';

% 定义真实模型，例如 y = 3 + 2*x + 0.5*x^
true_coefs = [3; 2; 0.5];
y_true = true_coefs(1) + true_coefs(2)*x.^(1/2) + true_coefs(3)*x.^(1/3);

% 添加高斯噪声
noise_sigma = x.*0.05;
noise = noise_sigma .* randn(n,1);
y = y_true + noise;

% 绘制生成的数据点和真实模型
figure;
scatter(x, y, 15,'k', 'filled');
hold on;
plot(x, y_true, 'r-', 'LineWidth', 2);
legend('数据点', '真实模型');
xlabel('x');
ylabel('y');
title('生成的数据点与真实模型');

%% 2. 手动设置拟合函数组

m = 3; % 多项式的最高阶数

% % 构建矩阵 A
% A = zeros(n, m+1);
% for i = 0:m
%     A(:,i+1) = x.^i;
% end

% 使用其他基函数
A = [ones(n,1), x.^(1/2), x.^(1/3)];
m = 3; % 相应调整阶数

%% 3. 给出法方程的矩阵形式

% 法方程：A' * A * c = A' * y
AtA = A' * A;
Aty = A' * y;

% 显示法方程的矩阵形式
disp('法方程的矩阵形式：');
disp('A'' * A =');
disp(AtA);
disp('A'' * y =');
disp(Aty);

% 求解法方程，得到拟合系数 c
my_coefs = AtA \ Aty;

% 显示拟合系数
disp('手动计算的拟合系数 c =');
disp(my_coefs);

% 计算拟合曲线
y_fit_manual = A * my_coefs;

% 绘制手动拟合曲线
plot(x, y_fit_manual, 'b--', 'LineWidth', 2);

%% 4. 使用 MATLAB 内置的最小二乘法函数进行拟合

% 方法1：使用 polyfit
degree = m; % 多项式的阶数
coefficients_polyfit = polyfit(x, y, degree);

% 调整polyfit系数结果
coefficients_polyfit_corrected = flip(coefficients_polyfit)';

% 计算 polyfit 的拟合曲线
y_fit_polyfit = polyval(coefficients_polyfit, x);

% 显示 polyfit 拟合系数
disp('使用 polyfit 的拟合系数 =');
disp(coefficients_polyfit_corrected);

% 绘制 polyfit 拟合曲线
plot(x, y_fit_polyfit, 'g-.', 'LineWidth', 2);

% 方法2：使用 lscov
[coefficients_lscov, ~] = lscov(A, y);

% 计算 lscov 的拟合曲线
y_fit_lscov = A * coefficients_lscov;

% 显示 lscov 拟合系数
disp('使用 lscov 的拟合系数 =');
disp(coefficients_lscov);

% 绘制 lscov 拟合曲线
plot(x, y_fit_lscov, 'm:', 'LineWidth', 2);


%% 5. 比较不同方法的拟合结果

legend('数据点', '真实模型', '手动拟合曲线', 'polyfit 拟合曲线', 'lscov 拟合曲线'); 

title('最小二乘法曲线拟合对比');

%% 6. 计算并显示拟合误差

% 手动拟合误差
residual_manual = y - y_fit_manual;
MSE_manual = sqrt(sum(residual_manual.^2));
disp(['手动拟合误差 MSE = ', num2str(MSE_manual)]);

% polyfit 拟合误差
residual_polyfit = y - y_fit_polyfit;
MSE_polyfit = sqrt(sum(residual_polyfit.^2));
disp(['polyfit 拟合误差 MSE = ', num2str(MSE_polyfit)]);

% lscov 拟合误差
residual_lscov = y - y_fit_lscov;
MSE_lscov = sqrt(sum(residual_lscov.^2));
disp(['lscov 拟合误差 MSE = ', num2str(MSE_lscov)]);

