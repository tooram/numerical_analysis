import numpy as np
from scipy.integrate import quad, quadrature
from numpy import trapz, cumsum  # 从 numpy 导入 trapz，用 cumsum 模拟 cumtrapz
import matplotlib.pyplot as plt

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题

# 参数设置
n_segments = 4  # 分段数，可以调整来控制精度
a, b = -1, 1      # 积分区间
x_vals = np.linspace(a, b, 100)  # 用于绘制的采样点

# 定义被积函数
f = lambda x: np.pi * np.cos(np.pi * x)

# 定义原函数 (理论积分结果)
F_exact = np.sin(np.pi * x_vals)

# 初始化积分结果
F_integral = np.zeros(len(x_vals))
F_custom_trapz = np.zeros(len(x_vals))
F_scipy_trapz = np.zeros(len(x_vals))
F_cumtrapz = np.zeros(len(x_vals))
F_simpson = np.zeros(len(x_vals))
F_segmented_simpson = np.zeros(len(x_vals))
F_segmented_trapz = np.zeros(len(x_vals))
F_gauss_legendre = np.zeros(len(x_vals))

# 自定义 cumtrapz 实现（使用 numpy 的 cumsum）
def custom_cumtrapz(y, x):
    dx = np.diff(x)
    return np.concatenate([[0], cumsum(dx * (y[:-1] + y[1:]) / 2)])

# 简单梯形方法：在整个区间 [0, xi] 上计算梯形积分
def simple_trapz(y, x):
    h = (x[-1] - x[0])
    return h * (y[0] + y[-1]) / 2

# 简单 Simpson 方法：在整个区间 [0, xi] 上计算 Simpson 积分
def simple_simpson(y, x):
    h = (x[-1] - x[0]) / 2
    return (h / 3) * (y[0] + 4 * y[1] + y[2])

# 计算积分
for i, xi in enumerate(x_vals):
    sub_x = np.array([0, xi / 2, xi])  # Simpson 方法需要三个点
    sub_f = f(sub_x)

    # SciPy 的 quad 方法 (相当于 MATLAB 的 integral)
    F_integral[i], _ = quad(f, 0, xi)

    # 简单梯形公式
    F_custom_trapz[i] = simple_trapz(sub_f, sub_x)

    # NumPy 的 trapz 方法
    F_scipy_trapz[i] = trapz(sub_f, sub_x)

    # 自定义 cumtrapz 实现
    sub_x_full = np.linspace(0, xi, n_segments + 1)
    sub_f_full = f(sub_x_full)
    F_cumtrapz_vals = custom_cumtrapz(sub_f_full, sub_x_full)
    F_cumtrapz[i] = F_cumtrapz_vals[-1]

    # 简单 Simpson 方法
    F_simpson[i] = simple_simpson(sub_f, sub_x)

    # 自定义分段 Simpson 方法（复合 Simpson 仅作为对比）
    F_segmented_simpson[i] = 0
    for j in range(0, n_segments, 2):
        F_segmented_simpson[i] += (xi / n_segments / 3) * (sub_f_full[j] + 4 * sub_f_full[j + 1] + sub_f_full[j + 2])

    # 自定义分段梯形方法（复合梯形仅作为对比）
    F_segmented_trapz[i] = 0
    for j in range(n_segments):
        F_segmented_trapz[i] += (xi / n_segments / 2) * (sub_f_full[j] + sub_f_full[j + 1])

    # SciPy 的 Gauss-Legendre (quadrature 方法)
    F_gauss_legendre[i], _ = quadrature(f, 0, xi)

# 计算每种方法的误差 (二范数)
error_integral = np.linalg.norm(F_integral - F_exact)
error_custom_trapz = np.linalg.norm(F_custom_trapz - F_exact)
error_scipy_trapz = np.linalg.norm(F_scipy_trapz - F_exact)
error_cumtrapz = np.linalg.norm(F_cumtrapz - F_exact)
error_simpson = np.linalg.norm(F_simpson - F_exact)
error_segmented_simpson = np.linalg.norm(F_segmented_simpson - F_exact)
error_segmented_trapz = np.linalg.norm(F_segmented_trapz - F_exact)
error_gauss_legendre = np.linalg.norm(F_gauss_legendre - F_exact)

# 分类输出误差
print('--- SciPy 的 integral 方法 (quad) ---')
print(f'Integral 方法误差（二范数）：{error_integral:.10f}')

print('\n--- 简单梯形方法 ---')
print(f'简单梯形公式误差（二范数）：{error_custom_trapz:.10f}')
print(f'NumPy trapz 方法误差（二范数）：{error_scipy_trapz:.10f}')
print(f'自定义 cumtrapz 方法误差（二范数）：{error_cumtrapz:.10f}')
print(f'分段复合梯形方法误差（二范数）：{error_segmented_trapz:.10f}')

print('\n--- 简单 Simpson 方法 ---')
print(f'简单 Simpson 方法误差（二范数）：{error_simpson:.10f}')
print(f'分段复合 Simpson 方法误差（二范数）：{error_segmented_simpson:.10f}')

print('\n--- 高斯-勒让德方法 ---')
print(f'Gauss-Legendre 方法误差（二范数）：{error_gauss_legendre:.10f}')

# 绘制积分结果与原函数的比较
plt.figure()
plt.plot(x_vals, F_exact, 'k-', linewidth=2, label='原函数 y = sin(πx)')
plt.plot(x_vals, F_integral, 'b-', linewidth=1.5, label='Integral 方法 (quad)')
plt.plot(x_vals, F_custom_trapz, 'g--', linewidth=1.5, label='简单梯形公式')
plt.plot(x_vals, F_scipy_trapz, 'y-.', linewidth=1.5, label='NumPy trapz 方法')
plt.plot(x_vals, F_cumtrapz, 'c-', linewidth=1.5, label='自定义 cumtrapz 方法')
plt.plot(x_vals, F_segmented_trapz, 'b--', linewidth=1.5, label='分段复合梯形')
plt.plot(x_vals, F_simpson, 'm-.', linewidth=1.5, label='简单 Simpson 方法')
plt.plot(x_vals, F_segmented_simpson, 'r-', linewidth=1.5, label='分段 Simpson')
plt.plot(x_vals, F_gauss_legendre, 'c:', linewidth=1.5, label='Gauss-Legendre 方法 (quadrature)')
plt.title('数值积分方法与原函数 y = sin(πx) 的比较')
plt.xlabel('x')
plt.ylabel('F(x)')
plt.legend(loc='best')
plt.grid()
plt.show()

# 绘制误差比较
plt.figure()
plt.plot(x_vals, np.abs(F_integral - F_exact), 'b-', linewidth=1.5, label='Integral 方法 (quad)')
plt.plot(x_vals, np.abs(F_custom_trapz - F_exact), 'g--', linewidth=1.5, label='简单梯形公式')
plt.plot(x_vals, np.abs(F_scipy_trapz - F_exact), 'y-.', linewidth=1.5, label='NumPy trapz 方法')
plt.plot(x_vals, np.abs(F_cumtrapz - F_exact), 'c-', linewidth=1.5, label='自定义 cumtrapz 方法')
plt.plot(x_vals, np.abs(F_segmented_trapz - F_exact), 'b--', linewidth=1.5, label='分段复合梯形')
plt.plot(x_vals, np.abs(F_simpson - F_exact), 'm-.', linewidth=1.5, label='简单 Simpson 方法')
plt.plot(x_vals, np.abs(F_segmented_simpson - F_exact), 'r-', linewidth=1.5, label='分段 Simpson')
plt.plot(x_vals, np.abs(F_gauss_legendre - F_exact), 'c:', linewidth=1.5, label='Gauss-Legendre 方法 (quadrature)')
plt.title('数值积分方法的误差比较')
plt.xlabel('x')
plt.ylabel('误差 |F_numerical - F_exact|')
plt.legend(loc='best')
plt.grid()
plt.show()
