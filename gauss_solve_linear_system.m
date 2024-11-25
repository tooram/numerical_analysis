function x = gauss_solve_linear_system(A, b, log_filename)
    % 检查输入是否有效
    [n, m] = size(A);
    if n ~= m
        error('系数矩阵 A 必须是方阵');
    end
    if length(b) ~= n
        error('右侧向量 b 的长度必须与矩阵 A 的维度一致');
    end

    % 打开日志文件
    log_file = fopen(log_filename, 'w');
    fprintf(log_file, '高斯消元法求解线性方程组\n\n');
    fprintf(log_file, '初始矩阵 A 和向量 b：\n');
    Aug = [A, b]; % 构造增广矩阵
    log_matrix(log_file, Aug, 10, 2);

    % 高斯消元法
    for i = 1:n
        % 主元为零时交换行
        if Aug(i, i) == 0
            for j = i+1:n
                if Aug(j, i) ~= 0
                    Aug([i, j], :) = Aug([j, i], :);
                    fprintf(log_file, '\n交换第 %d 行和第 %d 行：\n', i, j);
                    log_matrix(log_file, Aug, 10, 2);
                    break;
                end
            end
            if Aug(i, i) == 0
                error('矩阵不可逆，线性方程组无解或无唯一解');
            end
        end

        % 将主元归一化
        fprintf(log_file, '\n将第 %d 行的主元归一化：\n', i);
        Aug(i, :) = Aug(i, :) / Aug(i, i);
        log_matrix(log_file, Aug, 10, 2);

        % 消去其他行的对应列
        for j = 1:n
            if j ~= i
                fprintf(log_file, '\n将第 %d 行中消去列 %d 的元素：\n', j, i);
                Aug(j, :) = Aug(j, :) - Aug(j, i) * Aug(i, :);
                log_matrix(log_file, Aug, 10, 2);
            end
        end
    end

    % 提取解向量
    x = Aug(:, end);

    % 输出最终结果
    fprintf(log_file, '\n最终增广矩阵 [I | x]：\n');
    log_matrix(log_file, Aug, 10, 2);

    % 关闭日志文件
    fclose(log_file);
end

function log_matrix(log_file, matrix, width, precision)
    % 将矩阵写入日志文件，控制列宽和小数位数
    [rows, cols] = size(matrix);
    format_str = sprintf('%%%d.%df ', width, precision); % 构造格式字符串
    for i = 1:rows
        for j = 1:cols
            fprintf(log_file, format_str, matrix(i, j));
        end
        fprintf(log_file, '\n');
    end
    fprintf(log_file, '\n');
end