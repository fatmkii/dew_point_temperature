from sympy import diff, symbols, log


def f(x, y_0=0.0):
    """
    待求解的原函数
    :param x: float类型
    :param y_0: float类型。待求的x_0值，使得f(x_0)=y_0。
    :return: float类型
    """

    # Goff-Gratch方程，x=温度/K，y=饱和蒸汽压/hPa
    y = pow(10,
            (10.79574 * (1 - 273.16 / x)
             - 5.028 * log(x / 273.16, 10.0)
             + 1.50475 * 1e-4 * (1 - pow(10, -8.2969 * (x / 273.16 - 1)))
             + 0.4287 * 1e-3 * (pow(10, 4.76955 * (1 - 273.16 / x)) - 1)
             + 0.78614
             )
            )

    return y - y_0


def df(x):
    """
    原函数的导数
    :param x: float类型
    :return: float类型，返回df(x)。
    """

    t = symbols('t', real=True)  # 用于f(x)求导的sympy.symbol类型
    y = diff(f(t), t, 1)  # 对x求导

    return y.evalf(subs={t: x})


def newtonMethod(x_n, y_0, time=0):
    """
    牛顿迭代法
    :param x_n: float类型，本次迭代的x
    :param y_0: float类型，使f(x_n)值尽可能逼近y_0
    :param time: int类型，迭代次数显示
    :return: 迭代退出后返回x_n，使得f(x_n)=y_0
    """

    fx_n = f(x_n, y_0=y_0)  # 避免后续过程多次计算fx
    dfx_n = df(x_n)  # 避免后续过程多次计算dfx

    print('第%s次迭代：\tx=%s,\tf(x) = %s,\tdf(x) =%s' %
          (str(time), str(x_n), str(fx_n), str(dfx_n),))

    if f(x_n, y_0=y_0) == 0.0:  # 如果f(x)恰好为0，则完成迭代
        print('解为：\tx = %s,\tf(x) = %s' % (str(x_n), str(f(x_n))))
        return x_n
    else:
        x_n1 = x_n - fx_n / dfx_n

    dx_n1 = f(x_n1, y_0=y_0)
    if abs(fx_n - dx_n1) < 1e-6:  # 迭代完成条件
        print('解为：\tx = %s,\tf(x) = %s' %
              (str(x_n1), str(dx_n1)))
        return x_n1
    else:
        return newtonMethod(x_n1, y_0, time=time + 1)


if __name__ == '__main__':
    dry_bulb = float(input('请输入干球温度/℃：'))
    relative_humidity = float(input('请输入相对湿度/%：'))
    print('________________')

    saturated_vapor_pressure = f(dry_bulb + 273.15)  # 饱和蒸汽压/hPa
    recent_vapor_pressure = saturated_vapor_pressure * relative_humidity / 100  # 当前蒸汽压/hPa
    absolute_humidity = 0.622 * relative_humidity / (1013.25 - relative_humidity)  # 绝对湿度  kg/kg绝干气。1013.25hPa是标准大气压

    # 露点温度/℃，以牛顿迭代式求解（以干球温度为初始逼近值）
    dew_point_temperature = newtonMethod(dry_bulb + 273.15, recent_vapor_pressure) - 273.15

    print('绝对湿度：\t%s\t%s' % (round(absolute_humidity, 4), 'kg/kg绝干气'))
    print('露点温度：\t%s\t%s' % (round(dew_point_temperature, 1), '℃'))
