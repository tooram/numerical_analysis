t = [0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8]';
C = [19.21, 18.15, 15.36, 14.10, 12.89, 9.32, 7.45, 5.24, 3.01]'; 

lnC = log(C);
p_lib = polyfit(t, lnC, 1);
b_lib = p_lib(1);
a_lib = exp(p_lib(2));
n = length(t);
sum_t = sum(t);
sum_lnC = sum(lnC);
sum_t2 = sum(t.^2);
sum_t_lnC = sum(t .* lnC);

b_manual = (n * sum_t_lnC - sum_t * sum_lnC) / (n * sum_t2 - sum_t^2);
ln_a_manual = (sum_lnC - b_manual * sum_t) / n;
a_manual = exp(ln_a_manual);

fprintf('使用库函数的拟合结果:\n');
fprintf('a = %.4f, b = %.4f\n', a_lib, b_lib);
fprintf('\n不使用库函数的拟合结果:\n');
fprintf('a = %.4f, b = %.4f\n', a_manual, b_manual);

figure;
scatter(t, C, 'o', 'filled');
hold on;
t_fit = linspace(min(t), max(t), 100);
C_fit_lib = a_lib * exp(b_lib * t_fit);
C_fit_manual = a_manual * exp(b_manual * t_fit);

plot(t_fit, C_fit_lib, '-r', 'LineWidth', 2, 'DisplayName', '库函数拟合');
plot(t_fit, C_fit_manual, '--b', 'LineWidth', 2, 'DisplayName', '手动最小二乘法');
xlabel('时间 t (h)');
ylabel('血药浓度 C (微克/mL)');
title('血药浓度与时间的拟合结果');
legend('原始数据', '库函数拟合', '手动最小二乘法');
grid on;
hold off;
