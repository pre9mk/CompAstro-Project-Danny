import numpy as np

G = 4.498e-12

class Potential:
    """Class for potential components"""

    def evaluate(x, y, z):
        """Returns potential energy at (x, y, z)"""
        raise NotImplementedError("not implemented")

    def acceleration():
        """Returns acceleration vector [ax, ay, az]"""
        raise NotImplementedError("not implemented")


class MiyamotoNagai(Potential):
    """Miyamoto-Nagai disk potential"""

    def init():


    def evaluate(x, y, z):


    def acceleration():

        