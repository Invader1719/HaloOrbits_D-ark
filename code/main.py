import numpy as np

from Constants import Const
from ThreeBodyProblem import ThreeBodyProblem
from Plotter import Plotter

system = ThreeBodyProblem("Sun Earth") ## type: "Sun Earth" or "Earth Moon"
plotter = Plotter(system)

A, period = system.halo_r3oa(1.25 * 1e5, "L2") # У ричардсона эта цифра

print(A)

# Вот тут нужно перевести координаты ричардсона в наши. Сейчас ничего не работает
system.SolveDI(A, 10)
system.CaclulateJacobyIntegral()

plotter.Plot(system)
