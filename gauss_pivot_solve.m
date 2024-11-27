function x = gauss_pivot_solve(A, b, log_filename)
    % 高斯列主元消去法求解线性方程组 Ax = b，并输出日志
    % 输入:
    %   A: 系数矩阵 (n x n)
    %   b: 常数向量 (n x 1)
    %   log_filename: 日志文件名
    % 输出:
    %   x: 解向量 (n x 1)

    [n, m] = size(A);
    if n ~= m
        error('矩阵 A 必须是方阵');
    end

    if length(b) ~= n
        error('向量 b 的维度必须与矩阵 A 的行数一致');
    end

    % 打开日志文件
    log_file = fopen(log_filename, 'w');
    fprintf(log_file, '高斯列主元消去法\n');
    fprintf(log_file, '初始增广矩阵:\n');
    Aug = [A, b];
    log_matrix(log_file, Aug);

    % 消元过程
    for k = 1:n-1
        fprintf(log_file, '\n第 %d 步:\n', k);

        % 选择列主元
        [~, max_row] = max(abs(Aug(k:n, k)));
        max_row = max_row + k - 1;

        % 交换行
        if max_row ~= k
            fprintf(log_file, '交换第 %d 行和第 %d 行\n', k, max_row);
            Aug([k, max_row], :) = Aug([max_row, k], :);
            fprintf(log_file, '交换后的增广矩阵:\n');
            log_matrix(log_file, Aug);
        end

        % 消去
        for i = k+1:n
            factor = Aug(i, k) / Aug(k, k);
            fprintf(log_file, '消去第 %d 行的第 %d 列: 乘以 %.2f\n', i, k, factor);
            Aug(i, k:end) = Aug(i, k:end) - factor * Aug(k, k:end);
        end
        fprintf(log_file, '消元后的增广矩阵:\n');
        log_matrix(log_file, Aug);
    end

    % 回代过程
    fprintf(log_file, '\n回代过程:\n');
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (Aug(i, end) - Aug(i, i+1:n) * x(i+1:n)) / Aug(i, i);
        fprintf(log_file, 'x(%d) = %.2f\n', i, x(i));
    end

    % 关闭日志文件
    fclose(log_file);
end

function log_matrix(file, M)
    % 辅助函数：将矩阵 M 写入日志文件，保留两位小数
    [rows, cols] = size(M);
    for i = 1:rows
        for j = 1:cols
            fprintf(file, '%8.2f', M(i, j));
        end
        fprintf(file, '\n');
    end
end
