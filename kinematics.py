import numpy as np
from functools import reduce

class Kinematics:
    def __init__(self):
        return

    @classmethod
    def TCPFrame(cls,a,alfa,d,angles):
        z1 = zip(a,alfa,d,angles)
        matrices = map(lambda tup: cls.DHmatrix(tup[0], tup[1], tup[2], tup[3]), z1)
        finalPos = reduce(lambda m1,m2: m1 * m2, matrices)
        return finalPos

    @staticmethod
    def DHmatrix(a,alfa,d,theta):
        lam = np.cos(np.array(alfa))
        mu = np.sin(np.array(alfa))
        DH = np.matrix([
            [np.cos(theta), -lam*np.sin(theta), mu*np.sin(theta), a*np.cos(theta)],
            [np.sin(theta), lam*np.cos(theta), -mu*np.cos(theta), a*np.sin(theta)],
            [0, mu, lam, d ],
            [0, 0, 0, 1]])
        return DH

        