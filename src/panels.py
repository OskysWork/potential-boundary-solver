import numpy as np

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

     a0_ij = (
          -(1 / (4 * np.pi)) *
          np.log((x_local**2 + y_local**2) / ((x_local - L)**2 + y_local**2))
     )
     a1_ij = (
          -(1 / (2*np.pi)) * (
               (x_local / 2) *
               np.log((x_local**2 + y_local**2) / ((x_local - L)**2 + y_local**2)) -
               L +
               y_local * (np.arctan2(x_local, y_local) - np.arctan2((x_local - L), y_local))
          )
     )

     return a0_ij, a1_ij


def influence_matrix(panels):
    """
    Creates the full matrix of influence coefficients
    """

    N = len(panels)
    A = np.zeros((2*N, 2*N))
    b = np.zeros(2*N)
    for i, panel_i in enumerate(panels):
        for j, panel_j in enumerate(panels):
            if i == j:      # Added if else statement for self-influence singularity
                a0_ij = 0.0
                a1_ij = panels[j].length/(2*np.pi)
                A[i, 2*j] = a0_ij
                A[i, 2*j+1] = a1_ij
            else:
                a0_ij, a1_ij = influence_coefficient(panel_j, panel_i)
                A[i, 2*j] = a0_ij
                A[i, 2*j+1] = a1_ij
    
    # Panel strength continuity
    for j in range(N-1):
        row = N + j
        A[row, 2*j] = 1.0
        A[row, 2*j+1] = panels[j].length
        A[row, 2*(j+1)] = -1.0
        b[row] = 0.0
    
    row = N + (N-1)
    A[row, 2*(N-1)] = 1.0
    A[row, 2*(N-1)+1] = panels[-1].length
    A[row, 0] = 0.0#-1.0
    b[row] = 0.0

    cond_A = np.linalg.cond(A)
    print(cond_A)

    return A, b


def solve_panels(panels, U_free=1.0, alpha=0.0):
    """
    Solves system of equations to give panel vortext strengths
    """

    N = len(panels)
    A, b = influence_matrix(panels)

    for i, panel_i in enumerate(panels):
        nx = np.sin(panel_i.beta)
        ny = -np.cos(panel_i.beta)
        b[i] = -U_free * (nx * np.cos(alpha) + ny * np.sin(alpha))

    gamma = np.linalg.solve(A, b)

    for j, panel in enumerate(panels):
        panel.gamma_0 = gamma[2*j]
        panel.gamma_1 = gamma[2*j+1]

    return panels

def compute_cp(panels, U_free=1.0, alpha=0.0):
    """
    Computes xp and Cp from panel strengths
    """
    
    xp = []
    Cp = []

    for panel in panels:
        xp.append(panel.xc)

        gamma_mean = panel.gamma_0 + panel.gamma_1 * (panel.length / 2)
        Vt = U_free * np.sin(alpha - panel.beta) + 0.5 * gamma_mean
        Cp_panel = 1.0 - (Vt / U_free)**2
        Cp.append(Cp_panel)

    return xp, Cp