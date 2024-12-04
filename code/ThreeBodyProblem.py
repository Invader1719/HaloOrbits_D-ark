from token import COMMA

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Constants import Const, countlibra


class ThreeBodyProblem:

    def __init__(self, type):
        if type == "Sun Earth":
            self.mu = Const.Mass_Earth / (Const.Mass_Earth + Const.Mass_Sun)
            self.a = Const.Dist_Sun
            self.bodies = ["Солнце", "Земля"]

        elif type == "Earth Moon":
            self.mu = Const.Mass_Moon / (Const.Mass_Earth + Const.Mass_Moon)
            self.a = Const.Dist_Moon
            self.bodies = ["Земля", "Луна"]
        else:
            raise "wrong system type"

        self.Libration = countlibra(self.mu)

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

    def halo_r3oa(self, Az_km, lp_name, n_phase=1):
        '''
        :param Az_km: Стартовая точка по оси z в ккилометрах, мысчитаем что начинаем максимально высоко
        :param lp_name: Вокруг какой точки лагража начинаем
        :param n_phase:
        :return: параметра старта для гало орбиты и период
        '''
        def lambda_eqw(lmb):
            return lmb**4 + (c2 - 2.) * lmb**2 - (c2 - 1.) * (1. + 2. * c2)

        # Get L points
        if lp_name == "L1":
            gamma_l = (1. - self.mu) - self.Libration[0][0] # Расстояние от Луны до точки лагранжа
            gamma_expr = gamma_l / (1. - gamma_l)
            c2 = (self.mu + (1. - self.mu) * gamma_expr**3) / gamma_l**3
            c3 = (self.mu - (1. - self.mu) * gamma_expr**4) / gamma_l**3
            c4 = (self.mu + (1. - self.mu) * gamma_expr**5) / gamma_l**3
        elif lp_name == "L2":
            gamma_l = self.Libration[0][1] - (1. - self.mu)
            gamma_expr = gamma_l / (1. + gamma_l)
            c2 = (self.mu + (1. - self.mu) * gamma_expr**3) / gamma_l**3
            c3 = (-self.mu - (1. - self.mu) * gamma_expr**4) / gamma_l**3
            c4 = (self.mu + (1. - self.mu) * gamma_expr**5) / gamma_l**3
        else:
            raise Exception('Unsupported libration point specification')

        unit_length = gamma_l * self.a

        root = fsolve(lambda_eqw, np.array(1000))
        lmbda = root[(np.imag(root) == 0 and root > 0)][0]

        k = 2. * lmbda / (lmbda**2 + 1. - c2)


        # Все константы вида буква_номер смотреть ричардсона
        d1 = 3. * lmbda**2 * (k * (6. * lmbda**2 - 1.) - 2. * lmbda) / k
        d2 = 8. * lmbda**2 * (k * (11. * lmbda**2 - 1.) - 2. * lmbda) / k

        a21 = 3. * c3 * (k**2 - 2) / 4. / (1. + 2. * c2)
        a22 = 3. * c3 / 4. / (1. + 2 * c2)
        a23 = - 3. * c3 * lmbda * (3. * k ** 3 * lmbda - 6. * k * (k - lmbda) + 4.) / (4. * k * d1)
        a24 = - 3. * c3 * lmbda * (2. + 3. * k * lmbda) / (4. * k * d1)

        d21 = -c3 / 2. / lmbda**2
        d31 = 3. * (4. * c3 * a24 + c4) / 64. / lmbda**2
        d32 = 3. * (4. * c3 * (a23 - d21) + c4 * (4. + k ** 2)) / 64. / lmbda**2

        b21 = -3. * c3 * lmbda * (3. * k * lmbda - 4.) / 2. / d1
        b22 = 3. * c3 * lmbda / d1
        b31 = 3. * (8. * lmbda * (3. * c3 * (k * b21 - 2. * a23) - c4 * (2. + 3. * k**2)) +
                    (9. * lmbda**2 + 1. + 2. * c2) * (4. * c3 * (k * a23 - b21) + k * c4 * (4. + k**2))) / 8. / d2
        b32 = (9. * lmbda * (c3 * (k * b22 + d21 - 2. * a24) - c4) +
               3. * (9. * lmbda**2 + 1. + 2. * c2) * (4. * c3 * (k * a24 - b22) + k * c4) / 8.) / d2

        a31 = -9. * lmbda * (4. * c3 * (k * a23 - b21) + k * c4 * (4. + k**2)) / 4. / d2 + \
            (9. * lmbda**2 + 1. - c2) * (3. * c3 * (2. * a23 - k * b21) + c4 * (2. + 3. * k**2)) / 2. / d2
        a32 = -(9. * lmbda * (4. * c3 * (k * a24 - b22) + k * c4) / 4. +
                1.5 * (9. * lmbda**2 + 1. - c2) * (c3 * (k * b22 + d21 - 2. * a24) - c4)) / d2

        denom = 2. * lmbda * (lmbda * (1. + k**2) - 2. * k)
        s1 = (1.5 * c3 * (2. * a21 * (k ** 2 - 2.) - a23 * (k**2 + 2.) - 2. * k * b21) -
              3. * c4 * (3. * k**4 - 8 * k**2 + 8.) / 8.) / denom
        s2 = (1.5 * c3 * (2. * a22 * (k**2 - 2.) + a24 * (k**2 + 2.) + 2. * k * b22 + 5. * d21) +
              3. * c4 * (12. - k**2) / 8.) / denom

        a1 = -1.5 * c3 * (2. * a21 + a23 + 5. * d21) - 3. * c4 * (12. - k**2) / 8.
        a2 = 1.5 * c3 * (a24 - 2. * a22) + 9. * c4 / 8.

        l1 = a1 + 2. * lmbda**2 * s1
        l2 = a2 + 2. * lmbda**2 * s2

        Delta = lmbda**2 - c2

        Az = Az_km / (unit_length)

        Ax = np.sqrt((-l2 * Az**2 - Delta) / l1)
        Ay = k * Ax

        omega = 1 + s1 * Ax**2 + s2 * Az**2
        period = 2 * np.pi / omega / lmbda

        theta = 0.
        delta_n = 2. - n_phase  # n_phase can be 1 or 3 (input argument)
        tau_1 = 0. + theta  # tau_z = 0. + (theta + np.pi * n / 2) this is taken care of by delta_n

        x0 = a21 * Ax**2 + a22 * Az**2 - Ax * np.cos(tau_1) + (a23 * Ax**2 - a24 * Az**2) * np.cos(2. * tau_1) + (a31 * Ax**3 - a32 * Ax * Az**2) * np.cos(3. * tau_1)
        y0 = k * Ax * np.sin(tau_1) + (b21 * Ax**2 - a22 * Az**2)*np.sin(2 * tau_1) + (b31 * Ax**3 - b32 * Ax * Az**2) * np.sin(3. * tau_1)
        z0 = delta_n * (Az * np.cos(tau_1) + d21 * Ax * Az * (np.cos(2. * tau_1) - 3.) + (d32 * Az * Ax**2 - d31 * Az**3) * np.cos(3. * tau_1))
        vy0 = (k * Ax * np.cos(tau_1) + (b21 * Ax**2 - b22 * Az**2) * (2. * np.cos(2. * tau_1)) + (b31 * Ax**3 - b32 * Ax * Az**2) * 3. * np.cos(3. * tau_1))

        return np.array([x0, 0, z0, 0, vy0, 0]), period

