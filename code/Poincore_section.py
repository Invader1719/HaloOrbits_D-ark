import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Constants import Const


# Определяем функцию для правых частей системы
def Poincore(t, variables):
    '''
     q_1 = q
     q_2 = q'
    dq2_dt = q''
    '''
    x1, x2, y1, y2, z1, z2 = variables

    dx1_dt = x2
    dx2_dt = 2 * y2 + x1 - ((1 - Const.Mu) * (x1 + Const.Mu)) / ((x1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2) ** 1.5 - (Const.Mu * (x1 - 1 + Const.Mu)) / ((x1 - 1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2)**1.5

    dy1_dt = y2
    dy2_dt = - 2 * x2 + y1 - ((1 - Const.Mu) * y1) / ((x1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2) ** 1.5 - (Const.Mu * y1) / ((x1 - 1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2) ** 1.5

    dz1_dt = z2
    dz2_dt = -((1 - Const.Mu) * z1) / ((x1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2) ** 1.5 - (Const.Mu * z1) / ((z1 - 1 + Const.Mu) ** 2 + y1 ** 2 + z1 ** 2) ** 1.5

    return [dx1_dt, dx2_dt, dy1_dt, dy2_dt, dz1_dt, dz2_dt]


# Начальные условия: x(0), x'(0), y(0), y'(0), z(0), z'(0)
initial_conditions = [1-Const.Mu + 0.16784756, 0, 0, 0, 0, 0]

t_span = (0, 1)
t_eval = np.linspace(t_span[0], t_span[1], 10000)

solution = solve_ivp(Poincore, t_span, initial_conditions, t_eval=t_eval)



# Добавляем ЛУны и ЗЕмли
#Земля
xe = -Const.Mu
ye = 0
ze = 0

xm = 1 -Const.Mu
ym = 0
zm = 0


# График решения
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(solution.y[0], solution.y[2], solution.y[4], label="Траектория (x, y, z)")
ax.set_xlabel("x(t)")
ax.set_ylabel("y(t)")
ax.set_zlabel("z(t)")
ax.legend()
plt.title("Ограниченная задача трех тел")


# Рисуем точки
ax.scatter([xe, xm], [ye, ym], [ze, zm], color="red", label="Точки")

# Подписываем точки
ax.text(xe, ye, ze, "Земля", color="red")
ax.text(xm, ym, zm, "Луна", color="red")



plt.show()
