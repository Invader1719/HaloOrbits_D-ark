from Constants import Const
from  Functions import Solver5

coefficients = [1, (3-Const.Mu), (3-2*Const.Mu), -Const.Mu, -2*Const.Mu, -Const.Mu]
coef = [1, -4, 3, -2, 1, 5]

roots = Solver5(coefficients).RealRoots()
print(roots)