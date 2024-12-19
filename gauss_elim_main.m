clc; clear;

% 定义系数矩阵 A 和向量 b
A = [ 2 -1  0 ;
     -1  2  1 ;
      0 -1  2 ];

b = [0 ; 1 ; 0];


% 求解线性方程组 Ax = b
disp('1. 高斯若当法求解线性方程组 Ax = b');
log_filename1 = 'gauss_solution_fixed_log.txt';
x = gauss_solve_linear_system(A, b, log_filename1);
disp('线性方程组 Ax = b 的解为：');
disp(x);

% 求解线性方程组 Ax = b
disp('2. 高斯列主元法求解线性方程组 Ax1 = b');
log_filename2 = 'gauss_pivot_log.txt';
x1 = gauss_pivot_solve(A, b, log_filename2);
disp('线性方程组 Ax = b 的解为：');
disp(x1);

% 求矩阵 A 的逆
disp('2. 求矩阵 A 的逆');
log_filename3 = 'gauss_inverse_fixed_log.txt';
A_inv = gauss_elimination_inverse(A, log_filename3);
disp('矩阵 A 的逆为：');
disp(A_inv);

% 比较结果验证
disp('验证结果：');
fprintf('A * x 应等于 b：\n');
disp(A * x);
fprintf('A * x1 应等于 b：\n');
disp(A * x1);
fprintf('A * A_inv 应为单位矩阵：\n');
disp(A * A_inv);
