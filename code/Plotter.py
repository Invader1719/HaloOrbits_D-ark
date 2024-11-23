import matplotlib.pyplot as plt

from Constants import Const

class Plotter:
    # График решения
    def __init__(self, three_body_problem):
        self.__t_fig = plt.figure("Траектория")
        self.__t_ax = self.__t_fig.add_subplot(111, projection='3d')

        self.__t_ax.set_xlabel("x(t)")
        self.__t_ax.set_ylabel("y(t)")
        self.__t_ax.set_zlabel("z(t)")
        #self.__t_ax.legend()
        plt.title("Ограниченная круговая задача трех тел")

        xe = -Const.Mu
        ye = 0
        ze = 0

        xm = 1 - Const.Mu
        ym = 0
        zm = 0

        # Рисуем точки
        self.__t_ax.scatter([xe, xm], [ye, ym], [ze, zm], color="red", label="Точки")

        # Подписываем точки
        self.__t_ax.text(xe, ye, ze, three_body_problem.bodies[0], color="red")
        self.__t_ax.text(xm, ym, zm, three_body_problem.bodies[1], color="red")

        self.__t_ax.set_xlim([1 - Const.Mu - 0.6, 1 - Const.Mu + 0.6])
        self.__t_ax.set_ylim([-0.6, 0.6])
        self.__t_ax.set_zlim([-0.6, 0.6])

        self.__i_fig = plt.figure("Интеграл Якоби")
        self.__i_ax = self.__i_fig.add_subplot()

    def Plot(self, three_body_problem):
        self.__t_ax.plot(three_body_problem.solution.y[0],
                         three_body_problem.solution.y[1],
                         three_body_problem.solution.y[2],
                         label="Траектория (x, y, z)")

        self.__i_ax.plot(three_body_problem.solution.t,
                         three_body_problem.jacoby_integral)

        plt.show()


