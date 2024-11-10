from Functions import Solver5


class Const:
    Mass_Earth = 5.9722e24
    Mass_Moon = 7.3477e22
    Grav = 6.67430e-11
    Dist_Moon = 384399000

    # Обезразмереннные величины

    Mu = Mass_Moon / (Mass_Earth + Mass_Moon)

    # self.L1 =
    L2 = Solver5([1, (3 - Mu), (3 - 2 * Mu), -Mu, -2 * Mu, -Mu]).real_roots()[0]
    # self.L3 =
    # self.L4 =
    # self.L5 =
