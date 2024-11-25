% MATLAB代码：Romberg积分法计算 I = ∫(sin(x)/x)dx 的近似值

% 定义被积函数
f = @(x) (x == 0) .* 1 + (x ~= 0) .* (sin(x)./x);

% 积分上下限
a = 0;
b = 1;

% Romberg积分法参数
tol = 1e-7; % 误差容限
max_iter = 20; % 最大迭代次数

% 初始化日志文件
log_file = fopen('romberg_integration_log.txt', 'w');
fprintf(log_file, 'Romberg Integration Table Log:\n\n');

% 动态生成表格标题
header = sprintf('%-10s', 'h_k'); % 初始化标题字符串
for m = 0:max_iter-1
    header = strcat(header, sprintf('%-15s', sprintf('R(k,%d)', m)));
end
fprintf(log_file, '%s\n', header); % 打印标题
fprintf(log_file, '%s\n', repmat('-', 1, length(header))); % 分隔线

% Romberg表格初始化
R_table = zeros(max_iter, max_iter); % 初始化表
h = b - a; % 初始步长

% 初始梯形公式计算
R_table(1, 1) = (h / 2) * (f(a) + f(b));
fprintf(log_file, '%-10.5f %-15.10f\n', h, R_table(1, 1));

% Romberg积分递推计算
for k = 2:max_iter
    % 更新步长
    h = (b - a) / (2^(k-1));

    % 计算当前层的梯形积分值
    sum_mid = 0;
    for j = 1:2^(k-2)
        sum_mid = sum_mid + f(a + (2*j-1)*h);
    end
    R_table(k, 1) = 0.5 * R_table(k-1, 1) + h * sum_mid;

    % Romberg递推公式
    for m = 2:k
        R_table(k, m) = R_table(k, m-1) + ...
            (R_table(k, m-1) - R_table(k-1, m-1)) / (4^(m-1) - 1);
    end

    % 打印当前行数据到日志文件
    row = sprintf('%-10.5f', h); % 步长
    for m = 1:k
        row = strcat(row, sprintf('%-15.10f', R_table(k, m))); % 当前行值
    end
    fprintf(log_file, '%s\n', row);

    % 判断收敛条件
    if abs(R_table(k, k) - R_table(k-1, k-1)) < tol
        fprintf(log_file, 'Converged at step %d with R(%d,%d) = %.10f\n', ...
            k, k, k, R_table(k, k));
        break;
    end
end

% 关闭日志文件
fclose(log_file);

% 显示最终结果
final_result = R_table(k, k);
fprintf('积分结果为: %.10f\n', final_result);
fprintf('详细计算过程已保存到文件: romberg_integration_log.txt\n');
