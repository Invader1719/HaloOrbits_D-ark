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
    x1, y1, z1, x2, y2, z2 = variables


    r1 = np.array([-Const.Mu, 0, 0])
    r2 = np.array([1 - Const.Mu, 0, 0])
    z = np.array([0, 0, 1])

    r = np.array([x1, y1, z1])
    dr_dt = np.array([x2, y2, z2])
    dr2_dt = np.array([0, 0, 0])

    dr2_dt = (-(1 - Const.Mu) * pow(np.linalg.norm(r - r1), -3) * (r - r1) -
              Const.Mu * pow(np.linalg.norm(r - r2), -3) * (r - r2) -
              np.cross(z, np.cross(z, r)) -
              2 * np.cross(z, dr_dt))


    return [dr_dt[0], dr_dt[1], dr_dt[2], dr2_dt[0], dr2_dt[1], dr2_dt[2]]


# Начальные условия: x(0), y(0), z(0), x'(0), y'(0), z'(0)
initial_conditions = [1 - Const.Mu + 0.16784756630621525, 0, 0, 0, 0, 0] #0.16784756

t_span = (0, 22)
t_eval = np.linspace(t_span[0], t_span[1], 10000)

solution = solve_ivp(Poincore, t_span, initial_conditions, t_eval=t_eval, method='RK45')



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
ax.plot(solution.y[0], solution.y[1], solution.y[2], label="Траектория (x, y, z)")

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

ax.set_xlim([-0.2, 1])
ax.set_ylim([-0.6, 0.6])
ax.set_zlim([-0.6, 0.6])


plt.show()
