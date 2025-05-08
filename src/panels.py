import numpy as np
from utilities import read_aerofoil

class PanelObj:
    def __init__(self, xa, ya, xb, yb):
        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb
        self.xc = 0.5*(xa + xb)
        self.yc = 0.5*(ya + yb)
        self.length = np.hypot(xb - xa, yb - ya)
        self.beta = np.arctan2(yb - ya, xb - xa)
        self.gamma_0 = 0.0
        self.gamma_1 = 0.0
        self.v_tangent = 0.0
        self.Cp = 0.0


def discretise_panels(x, y):
        """
        Defines panels array from (closed loop) data
        """
        panels = []
        N = len(x) - 1      # As last point == first
        for i in range(N):
             panels.append(PanelObj(x[i], y[i], x[i+1], y[i+1]))
        return panels


def influence_coefficient(panel_j, panel_i):
     """
     Calculates influence coefficient for given panel_j
     """

     dx = panel_i.xc - panel_j.xa
     dy = panel_i.yc - panel_j.ya
     phi = -panel_j.beta
     x_local = dx*np.cos(phi) - dy*np.sin(phi)
     y_local = dx*np.sin(phi) + dy*np.cos(phi)
     L = panel_j.length

     a1_ij = (
          (1 / (4 * np.pi)) *
          np.log((x_local**2 + y_local**2) / ((x_local - L_**2 + y_local**2)))
     )
     a2_ij = (
          (1 / (2*np.pi)) * (
               (x_local / 2) *
               np.log((x_local**2 + y_local**2) / ((x_local - L_**2 + y_local**2))) -
               L +
               y_local * (np.atan2(x_local / y_local) - np.atan2((x_local - L) / y_local))
          )
     )

     return a1_ij, a2_ij


def influence_matrix(panels):
     """
     Creates the full matrix of influence coefficients
     """

     N = len(panels)
     A = np.zeros((2*N, 2*N))
     b = np.zeros(2*N)

     for i, panel_i in enumerate(panels):
          for j, panel_j in enumerate((panels)):
               I0, I1 = 