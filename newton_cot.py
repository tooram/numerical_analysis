import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

def compute_newton_cotes_weights(a, b, n):
    """
    计算牛顿-科特斯公式的系数。

    参数：
    a, b: 积分区间
    n: 节点数量

    返回：
    x_nodes: 等距节点
    weights: 牛顿-科特斯公式系数
    """
    x_nodes = np.linspace(a, b, n)
    weights = np.zeros(n)

    for i in range(n):
        # 构建Lagrange基函数的多项式
        xi = x_nodes[i]
        numerator = Polynomial([1])
        denominator = 1

        for j in range(n):
            if j != i:
                xj = x_nodes[j]
                # 更新分子
                numerator *= Polynomial([-xj, 1])
                # 更新分母
                denominator *= (xi - xj)

        # 计算权重为Lagrange基函数在区间上的积分
        Li = numerator / denominator
        Li_integ = Li.integ()
        weights[i] = Li_integ(b) - Li_integ(a)

    return x_nodes, weights

def numerical_integration(f, x_nodes, weights):
    """
    使用牛顿-科特斯公式进行数值积分。

    参数：
    f: 被积函数
    x_nodes: 积分节点
    weights: 牛顿-科特斯系数

    返回：
    积分结果
    """
    return np.dot(weights, f(x_nodes))

# 定义被积函数和原函数
def integrand(x):
    return np.pi * np.cos(np.pi * x)

def original_function(x):
    return np.sin(np.pi * x)

# 积分区间
a = -1
b = 1

# 节点数量列表
n_list = [2, 4, 6, 8, 10, 12]

# 用于绘图的数据
x_plot = np.linspace(a, b, 200)
original_values = original_function(x_plot)

plt.figure(figsize=(12, 8))

# 保存误差
errors = []

# 绘制原函数曲线
plt.plot(x_plot, original_values, 'k-', label='Exact Function $\\sin(\\pi x)$')

for n in n_list:
    numerical_result = []
    exact_result = []
    error = []

    for x_end in x_plot:
        if x_end == a:
            numerical_integral = 0
            exact_integral = 0
        else:
            x_nodes, weights = compute_newton_cotes_weights(a, x_end, n)
            numerical_integral = numerical_integration(integrand, x_nodes, weights)
            exact_integral = original_function(x_end) - original_function(a)

        numerical_result.append(numerical_integral)
        exact_result.append(exact_integral)
        error.append(numerical_integral - exact_integral)

    # 绘制数值积分结果曲线
    plt.plot(x_plot, numerical_result, label=f'Numerical Integration n={n}')

    # 计算二范数误差
    error_norm = np.linalg.norm(error, ord=2)
    errors.append(error_norm)

# 设置图形参数
plt.title('Comparison of Exact Function and Numerical Integration Results')
plt.xlabel('x')
plt.ylabel('Integral Value')
plt.legend()
plt.grid(True)
plt.show()

# 输出二范数误差
for n, err in zip(n_list, errors):
    print(f'n={n}, 二范数误差: {err:.6e}')
