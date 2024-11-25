import numpy as np

from Constants import Const
from ThreeBodyProblem import ThreeBodyProblem
from Plotter import Plotter

system = ThreeBodyProblem("Sun Earth") ## type: "Sun Earth" or "Earth Moon"
plotter = Plotter(system)

system.SolveDI([1 - system.mu + 0.01, 0, 0, 0, np.sqrt(system.mu * 120), 0], 5)
system.CaclulateJacobyIntegral()

plotter.Plot(system)