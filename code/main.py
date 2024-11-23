import numpy as np

from Constants import Const
from  Functions import Solver5
from ThreeBodyProblem import ThreeBodyProblem
from Plotter import Plotter

system = ThreeBodyProblem("Sun Earth") ## type: "Sun Earth" or "Earth Moon"
plotter = Plotter(system)

system.SolveDI([1 - Const.Mu + 0.01, 0, 0, 0, np.sqrt(Const.Mu * 120), 0], 5)
system.CaclulateJacobyIntegral()

plotter.Plot(system)