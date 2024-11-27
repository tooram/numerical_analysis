import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题
# 定义被积函数和精确积分函数
def integrand(x):
    return np.pi * np.cos(np.pi * x)

def exact_integral(x):
    return np.sin(np.pi * x)

# 计算牛顿-科特斯公式的权重
def compute_newton_cotes_weights(a, b, n):
    x_nodes = np.linspace(a, b, n)      # 等距节点
    weights = np.zeros(n)               # 初始化权重数组

    for i in range(n):
        # 计算Lagrange基函数的分母: product(x_i - x_j) for j != i
        xi = x_nodes[i]
        denominator = 1
        for j in range(n):
            if j != i:
                denominator *= (xi - x_nodes[j])

        # 计算分子多项式的系数
        numerator = Polynomial([1])
        for j in range(n):
            if j != i:
                numerator *= Polynomial([-x_nodes[j], 1])

        # 计算分母上的基函数L_i(x)
        L_i = numerator / denominator

        # 对L_i(x)在区间[a, b]上积分
        L_i_integ = L_i.integ()
        weights[i] = L_i_integ(b) - L_i_integ(a)

    return weights

# 设置积分区间和绘图点
a_initial = -1     # 积分下限
b_initial = 1     # 积分上限
num_points = 200  # 绘图点数量
x_plot = np.linspace(a_initial, b_initial, num_points)  # x 轴点
exact_values = exact_integral(x_plot)                   # 精确积分值

# 定义不同的节点数量 n
n_list = [2, 4, 6, 8, 10, 12]

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(x_plot, exact_values, 'k-', linewidth=2, label=r'$F(x) = \sin(\pi x)$')

# 存储二范数误差
errors = []

# 循环遍历每个 n 值
for n in n_list:
    numerical_result = []    # 存储数值积分结果
    error_vector = []        # 存储误差

    # 对于每个 x_i，计算从 0 到 x_i 的积分
    for x_i in x_plot:
        if x_i == 0:
            numerical_integral = 0  # 当 x_i 为 0 时积分结果为 0
            exact_val = 0
        else:
            # 设置积分区间为 [0, x_i]
            a = 0
            b = x_i

            # 计算牛顿-科特斯权重
            weights = compute_newton_cotes_weights(a, b, n)

            # 定义节点
            x_nodes = np.linspace(a, b, n)

            # 评估被积函数在节点处的值
            f_values = integrand(x_nodes)

            # 计算数值积分
            numerical_integral = np.dot(weights, f_values)

            # 计算精确积分值
            exact_val = exact_integral(x_i)

        # 保存当前积分结果和误差
        numerical_result.append(numerical_integral)
        error_vector.append(numerical_integral - exact_val)

    # 绘制数值积分结果曲线
    plt.plot(x_plot, numerical_result, linewidth=1.5, label=f'n={n}')

    # 计算二范数误差
    error_norm = np.linalg.norm(error_vector, ord=2)
    errors.append(error_norm)

# 完善图形
plt.title('牛顿-柯特斯公式实现数值积分')
plt.xlabel('x')
plt.ylabel('积分值')
plt.legend()
plt.grid(True)
plt.show()

# 输出二范数误差
print("二范数误差:")
for n, error in zip(n_list, errors):
    print(f'n={n}, 二范数误差: {error:.6e}')
