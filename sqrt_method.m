function x = sqrt_method(A, b)
    if ~isequal(A, A')
        error('矩阵 A 不是对称矩阵，无法使用平方根法');
    end
    n = size(A, 1);
    G = zeros(n, n);
    for i = 1:n
        sum_diag = 0;
        for k = 1:i-1
            sum_diag = sum_diag + G(i, k)^2;
        end
        G(i, i) = sqrt(A(i, i) - sum_diag);
        if G(i, i) <= 0
            error('矩阵 A 不是正定矩阵，无法进行 Cholesky 分解');
        end
        for j = i+1:n
            sum_off_diag = 0;
            for k = 1:i-1
                sum_off_diag = sum_off_diag + G(i, k) * G(j, k);
            end
            G(j, i) = (A(j, i) - sum_off_diag) / G(i, i);
        end
    end
    fprintf('Cholesky 分解得到的下三角矩阵 G 为:\n');
    disp(G);
    y = zeros(n, 1);
    for i = 1:n
        sum_y = 0;
        for k = 1:i-1
            sum_y = sum_y + G(i, k) * y(k);
        end
        y(i) = (b(i) - sum_y) / G(i, i);
    end
    x = zeros(n, 1);
    for i = n:-1:1
        sum_x = 0;
        for k = i+1:n
            sum_x = sum_x + G(k, i) * x(k);
        end
        x(i) = (y(i) - sum_x) / G(i, i);
    end
    fprintf('线性方程组的解 x 为:\n');
    disp(x);
end
