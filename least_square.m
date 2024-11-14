clear; clc; close all;
%% 1. 自动生成数据点

rng(0);

% 生成自变量 x
n = 50;
x = linspace(0, 10, n)';
% x = [0.2, 0.5, 0.7, 0.85, 1]';

% y = [1.221, 1.649, 2.014, 2.340, 2.718]';

% 定义函数模型
true_coefs = [3; 2; 0.5];
y_true = true_coefs(1) + true_coefs(2)*x.^(1/2) + true_coefs(3)*x.^(1/3);

% 添加高斯噪声
noise_sigma = x.*0.05;
noise = noise_sigma .* randn(n,1);
y = y_true + noise;


figure;
scatter(x, y, 15,'k', 'filled');
hold on;
plot(x, y_true, 'r-', 'LineWidth', 2);
legend('数据点', '函数原型');
xlabel('x');
ylabel('y');
title('生成的数据点与函数模型');

%% 2. 手动设置拟合函数组

% m = 2; % 多项式的最高阶数

% % 构建矩阵 A
% A = zeros(n, m+1);
% for i = 0:m
%     A(:,i+1) = x.^i;
% end

% 使用其他基函数
A = [ones(n,1), x.^(1/2), x.^(1/3)];
m = 3; 

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
degree = m; 
coefs_polyfit = polyfit(x, y, degree);
coefs_polyfit_corrected = flip(coefs_polyfit)';
y_fit_polyfit = polyval(coefs_polyfit, x);

% 显示 polyfit 拟合系数
disp('使用 polyfit 的拟合系数 =');
disp(coefs_polyfit_corrected);

% 绘制 polyfit 拟合曲线
plot(x, y_fit_polyfit, 'g-.', 'LineWidth', 2);

% 方法2：使用 lscov
[coefs_lscov, ~] = lscov(A, y);
y_fit_lscov = A * coefs_lscov;

% 显示 lscov 拟合系数
disp('使用 lscov 的拟合系数 =');
disp(coefs_lscov);

% 绘制 lscov 拟合曲线
plot(x, y_fit_lscov, 'm:', 'LineWidth', 2);


%% 5. 比较不同方法的拟合结果

legend('数据点', '函数原型', '手动拟合曲线', 'polyfit 拟合曲线', 'lscov 拟合曲线'); 

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

