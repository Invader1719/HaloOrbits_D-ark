import numpy as np
class Solver5:
    '''
    Коэффициенты многочлена в формате [a, b, c, d, e, f]
    для многочлена a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f = 0
    Мы работаем np.array и храним там compex128
    '''
    def __init__(self, Array):
        self.coefficients = Array
        self.roots = np.roots(self.coefficients)

    def Roots(self):
        return self.roots

    def RealRoots(self):
        ans = []
        for root in self.roots:
            if (root.imag == 0):
                ans.append(root.real)
        return ans


"""
coefficients = [1, (3-Const.Mu), (3-2*Const.Mu), -Const.Mu, -2*Const.Mu, -Const.Mu]

roots = Solver5(coefficients).RealRoots()
print(roots)
"""