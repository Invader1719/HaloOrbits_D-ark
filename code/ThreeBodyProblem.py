from token import COMMA

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from Constants import Const

class ThreeBodyProblem:

    def __init__(self, type):
        if type == "Sun Earth":
            self.__mu = Const.Mass_Earth / (Const.Mass_Earth + Const.Mass_Sun)
            self.bodies = ["Солнце", "Земля"]
        elif type == "Earth Moon":
            self.__mu = Const.Mass_Moon / (Const.Mass_Earth + Const.Mass_Moon)
            self.bodies = ["Земля", "Луна"]
        else:
            raise "wrong system type"

        self.__sampling_frequency = 1000
        self.solution = np.zeros(self.__sampling_frequency)
        self.jacoby_integral = np.zeros(self.__sampling_frequency)
        self.__rtol = 1e-9
        self.__atol = 1e-9

    def SetPrecision(self, rtol, atol):
        self.__rtol = rtol
        self.__atol = atol

    # Определяем функцию для правых частей системы
    def CalculateForce(self, t, variables):
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

    def SolveDI(self, initial_conditions, time_limit):
        # Начальные условия: x(0), y(0), z(0), x'(0), y'(0), z'(0)
        #initial_conditions = [1 - Const.Mu + 0.01, 0, 0, 0, np.sqrt(Const.Mu * 120), 0]  # 0.16784756630621525

        t_span = (0, time_limit)
        t_eval = np.linspace(t_span[0], t_span[1], self.__sampling_frequency)

        self.solution = solve_ivp(self.CalculateForce,
                                    t_span,
                                    initial_conditions,
                                    t_eval=t_eval,
                                    rtol=self.__rtol,
                                    atol=self.__atol)
        return self.solution

    def CaclulateJacobyIntegral(self):
        self.jacoby_integral = np.zeros(self.__sampling_frequency)
        for i in range(self.__sampling_frequency):
            r = np.array([self.solution.y[0][i], self.solution.y[1][i], self.solution.y[2][i]])
            v = np.array([self.solution.y[3][i], self.solution.y[4][i], self.solution.y[5][i]])

            r1 = np.array([-Const.Mu, 0, 0])
            r2 = np.array([1 - Const.Mu, 0, 0])
            z = np.array([0, 0, 1])

            self.jacoby_integral[i] = (-(r[0] ** 2 + r[1] ** 2) -
                                        2 * (1 - Const.Mu) / np.linalg.norm(r - r1) -
                                        2 * Const.Mu / np.linalg.norm(r - r2) +
                                        np.linalg.norm(v) ** 2)

        return self.jacoby_integral






