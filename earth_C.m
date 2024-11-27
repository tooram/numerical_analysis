% 地球参数
R = 6371; % 地球半径 (km)
H = 2384; % 远地点高度 (km)
h = 439;  % 近地点高度 (km)

% 椭圆参数计算
a = (2 * R + H + h) / 2; % 长半轴
c = (H - h) / 2;         % 焦距

% 被积函数
f = @(theta) 4.*a.*sqrt(1 - (c / a)^2 * sin(theta).^2);

% Romberg 积分参数设置
max_iter = 10;    % 最大迭代次数
theta_start = 0;  % 积分下限
theta_end = pi/2; % 积分上限

% 初始化梯形积分和 Romberg 表
T = zeros(max_iter, max_iter); % 存储梯形积分和加速值

% 第一步：计算 T_{2^0}
h = theta_end - theta_start; % 初始步长
T(1, 1) = (h / 2) * (f(theta_start) + f(theta_end)); % 梯形公式初值

% 输出表格头
fprintf('%-5s %-15s %-15s %-15s %-15s\n', 'k', 'T_2^k', 'S_2^(k-1)', 'C_2^(k-2)', 'R_2^(k-3)');

% 输出第 1 行，仅包含 T_0
fprintf('%-5d %-15.7f\n', 0, T(1, 1));

% Romberg 递推
for k = 2:max_iter
    % k 对应公式中的 k-1
    formula_k = k - 1;

    % 计算 T_{2^k}
    h = h / 2; % 步长减半
    sum_mid = 0;

    % 计算中点累加和
    for j = 1:2^(formula_k-1)
        x_mid = theta_start + (2 * j - 1) * h;
        sum_mid = sum_mid + f(x_mid);
    end

    % 更新 T_k
    T(k, 1) = 0.5 * T(k-1, 1) + h * sum_mid;

    % 使用加速公式逐步计算 S, C, R
    for m = 2:k
        T(k, m) = T(k, m-1) + (T(k, m-1) - T(k-1, m-1)) / (4^(m-1) - 1);
    end

    % 输出当前行（显示公式中的 k-1）
    fprintf('%-5d %-15.7f', formula_k, T(k, 1)); % 输出 T_k
    for m = 2:k
        fprintf(' %-15.7f', T(k, m)); % 输出 S, C, R
    end
    fprintf('\n');

    % 停止计算：在第一个 Romberg 值 \( R_{2^{k-3}} \) 后终止
    if formula_k >= 3
        break;
    end
end

% 使用最后计算的 Romberg 值计算轨道周长
C_romberg =T(k, k);

% 使用 MATLAB 内置数值积分函数
C_builtin =integral(f, theta_start, theta_end, 'AbsTol', 1e-12);

% 输出结果，保留小数点后10位有效数字
fprintf('\n比较结果：\n');
fprintf('Romberg 方法计算的轨道周长为：%.7f km\n', C_romberg);
fprintf('MATLAB 内置数值积分计算的轨道周长为：%.7f km\n', C_builtin);
fprintf('误差为：%.7e km\n', abs(C_romberg - C_builtin));
