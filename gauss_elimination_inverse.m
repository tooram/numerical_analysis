function A_inv = gauss_elimination_inverse(A, log_filename)
    % 检查输入是否为方阵
    [n, m] = size(A);
    if n ~= m
        error('输入的矩阵必须是方阵');
    end

    % 打开日志文件
    log_file = fopen(log_filename, 'w');
    fprintf(log_file, '高斯消元法计算矩阵逆\n\n');
    fprintf(log_file, '初始矩阵 A：\n');
    log_matrix(log_file, A, 8, 2); % 调整宽度为8，精度为2

    % 初始化单位矩阵
    I = eye(n);

    % 构造增广矩阵 [A | I]
    Aug = [A, I];
    fprintf(log_file, '\n初始增广矩阵 [A | I]：\n');
    log_matrix(log_file, Aug, 8, 2);

    % 高斯消元法
    for i = 1:n
        % 检查主元是否为零，若为零则交换行
        if Aug(i, i) == 0
            for j = i+1:n
                if Aug(j, i) ~= 0
                    % 交换第 i 行和第 j 行
                    Aug([i, j], :) = Aug([j, i], :);
                    fprintf(log_file, '\n交换第 %d 行和第 %d 行：\n', i, j);
                    log_matrix(log_file, Aug, 8, 2);
                    break;
                end
            end
            if Aug(i, i) == 0
                error('矩阵不可逆');
            end
        end

        % 将主元归一化
        fprintf(log_file, '\n将第 %d 行的主元归一化：\n', i);
        Aug(i, :) = Aug(i, :) / Aug(i, i);
        log_matrix(log_file, Aug, 8, 2);

        % 对其他行进行消元
        for j = 1:n
            if j ~= i
                fprintf(log_file, '\n将第 %d 行中消去列 %d 的元素：\n', j, i);
                Aug(j, :) = Aug(j, :) - Aug(j, i) * Aug(i, :);
                log_matrix(log_file, Aug, 8, 2);
            end
        end
    end

    % 提取逆矩阵
    A_inv = Aug(:, n+1:end);
    fprintf(log_file, '\n最终增广矩阵 [I | A_inv]：\n');
    log_matrix(log_file, Aug, 8, 2);

    % 检查结果
    if any(isnan(A_inv(:))) || any(isinf(A_inv(:)))
        error('矩阵不可逆');
    end

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
