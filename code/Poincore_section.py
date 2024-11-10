from token import COMMA

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import solve_undetermined_coeffs

from Constants import Const

Sampling_frequency = 1000

# Определяем функцию для правых частей системы
def Poincore(t, variables):
    x, y, z, vx, vy, vz = variables

    r1 = np.array([-Const.Mu, 0, 0])
    r2 = np.array([1 - Const.Mu, 0, 0])

    # Ось вращения системы
    w = np.array([0, 0, 1])

    # Позиция и скорость маломассивного тела
    r = np.array([x, y, z])
    dr_dt = np.array([vx, vy, vz])

    dr2_dt = (-(1 - Const.Mu) * (r - r1) / np.linalg.norm(r - r1) ** 3
              - Const.Mu * (r - r2) / np.linalg.norm(r - r2) ** 3
              - np.cross(w, np.cross(w, r))
              - 2 * np.cross(w, dr_dt))

    return [dr_dt[0], dr_dt[1], dr_dt[2], dr2_dt[0], dr2_dt[1], dr2_dt[2]]

def Jacoby_Integral(solution):
    integral = np.zeros(len(solution[0]))
    for i in range(len(solution[0])):
        r = np.array([solution[0][i], solution[1][i], solution[2][i]])
        v = np.array([solution[3][i], solution[4][i], solution[5][i]])

        r1 = np.array([-Const.Mu, 0, 0])
        r2 = np.array([1 - Const.Mu, 0, 0])
        z = np.array([0, 0, 1])

        integral[i] = (-(r[0] ** 2 + r[1] ** 2) -
                       2 * (1 - Const.Mu) / np.linalg.norm(r - r1) -
                       2 * Const.Mu / np.linalg.norm(r - r2) +
                       np.linalg.norm(v) ** 2)

    return integral


# Начальные условия: x(0), y(0), z(0), x'(0), y'(0), z'(0)
initial_conditions = [1 - Const.Mu + 0.01, 0, 0, 0, np.sqrt(Const.Mu * 100), 0] #0.16784756630621525

t_span = (0, 5)
t_eval = np.linspace(t_span[0], t_span[1], Sampling_frequency)

solution = solve_ivp(Poincore, t_span, initial_conditions, t_eval=t_eval,  rtol=1e-9, atol=1e-9)



# Добавляем ЛУны и ЗЕмли
#Земля
xe = -Const.Mu
ye = 0
ze = 0

xm = 1 -Const.Mu
ym = 0
zm = 0


# График решения
fig = plt.figure("Траектория")
ax = fig.add_subplot(111, projection='3d')

fig2 = plt.figure("Интеграл Якоби")
ax2 = fig2.add_subplot()

ax.plot(solution.y[0], solution.y[1], solution.y[2], label="Траектория (x, y, z)")
ax2.plot(solution.t, Jacoby_Integral(solution.y))

ax.set_xlabel("x(t)")
ax.set_ylabel("y(t)")
ax.set_zlabel("z(t)")
ax.legend()
plt.title("Ограниченная круговая задача трех тел")


# Рисуем точки
ax.scatter([xe, xm], [ye, ym], [ze, zm], color="red", label="Точки")

# Подписываем точки
ax.text(xe, ye, ze, "Земля", color="red")
ax.text(xm, ym, zm, "Луна", color="red")

ax.set_xlim([-0.2, 1])
ax.set_ylim([-0.6, 0.6])
ax.set_zlim([-0.6, 0.6])


plt.show()
