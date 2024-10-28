clc;
clear all;
close all;
% 定义区间和函数
a = -1;
b = 1;
f = @(x) abs(x);

% 定义 M1 和 M2 中的基函数
M1 = {@(x) 1, @(x) x, @(x) x.^3};
M2 = {@(x) 1, @(x) x.^2, @(x) x.^4};

% 计算内积函数，加入 'ArrayValued', true
inner_product = @(g, h) integral(@(x) g(x) .* h(x), a, b, 'ArrayValued', true);

% 计算在空间 M1 中的最佳平方逼近
A1 = zeros(3, 3);
b1 = zeros(3, 1);
for i = 1:3
    for j = 1:3
        A1(i, j) = inner_product(M1{i}, M1{j});
    end
    b1(i) = inner_product(f, M1{i});
end
c1 = A1 \ b1;

% 计算在空间 M2 中的最佳平方逼近
A2 = zeros(3, 3);
b2 = zeros(3, 1);
for i = 1:3
    for j = 1:3
        A2(i, j) = inner_product(M2{i}, M2{j});
    end
    b2(i) = inner_product(f, M2{i});
end
c2 = A2 \ b2;

% 显示结果
fprintf('M1 中的最佳平方逼近多项式系数为：\n');
disp(c1);
fprintf('M2 中的最佳平方逼近多项式系数为：\n');
disp(c2);

% 绘图
x = linspace(a, b, 100);
f_approx_M1 = c1(1) + c1(2)*x + c1(3)*x.^3;
f_approx_M2 = c2(1) + c2(2)*x.^2 + c2(3)*x.^4;

figure;
plot(x, f(x), 'k', 'LineWidth', 1.5); hold on;
plot(x, f_approx_M1, 'r--', 'LineWidth', 1.5);
plot(x, f_approx_M2, 'b-.', 'LineWidth', 1.5);
legend('|x|', 'M1的最佳平方逼近', 'M2的最佳平方逼近');
xlabel('x');
ylabel('f(x)');
title(' M1 和 M2 的最佳平方逼近');
grid on;
