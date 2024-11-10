import numpy as np


class Solver5:
    """
    Коэффициенты многочлена в формате [a, b, c, d, e, f]
    для многочлена a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f = 0
    Мы работаем np.array и храним там compex128
    """

    def __init__(self, array):
        self.coefficients = array
        self.roots = np.roots(self.coefficients)

    def roots(self):
        return self.roots

    def real_roots(self):
        ans = []
        for root in self.roots:
            if root.imag == 0:
                ans.append(root.real)
        return ans
