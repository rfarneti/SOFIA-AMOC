import numpy as np
from scipy import integrate
#from matplotlib import pyplot as plt
from momsofia.make_func import make_func


class Psi_Thermwind(object):
  r"""
  Thermal Wind Closure

  Instances of this class represent the overturning circulation between two columns,
  given density profiles in those columns.

  The model assumes a thermal-wind based equation for the overturning circulation
  as in Nikurashin and Vallis (2012):

  .. math::
    d_{zz}\left(\Psi\right) = g / ( f \rho_0 ) (\rho_b - \rho_n)
    d_{zz}\left(\Psi\right) = f^{-1} (b_n - b_b)

  This equation is solved subject to the boundary conditions:

  .. math::
    \Psi(0) = \Psi(-H) = 0 

  (these BCs are different than NV2012)

  [From PyMOC - M. Jansen]

  USAGE:
  twb = Psi_Thermwind(z=z, rhob=rho_basin, rhon=rho_north)
  twb.solve()

  """
  def __init__(
      self,
      g=9.8,     # gravitational acceleration (input)
      rho0=1025, # reference density (input) 
      f=1.2e-4,  # Coriolis parameter (input)
      z=None,    # grid (input)
      sol_init=None,    # Initial conditions for ODE solver (input)
      bb=None,    # Buoyancy in the basin (input, output)
      bn=None,    # Buoyancy in the deep water formation region (input, output)
  ):

    self.f = f
    self.g = g
    self.rho0 = rho0
    # initialize grid:
    if isinstance(z, np.ndarray):
      self.z = z
      nz = np.size(z)
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    self.bb = make_func(bb, self.z, 'bb')
    self.bn = make_func(bn, self.z, 'bn')

    # Set initial conditions for BVP solver
    if sol_init is None:
      self.sol_init = np.zeros((2, nz))
    else:
      self.sol_init = sol_init

  def bc(self, ya, yb):
    return np.array([ya[0], yb[0]])

  def ode(self, z, y):
    #Buoyancy 
    #return np.vstack((y[1], 1. / self.f * (self.bn(z) - self.bb(z))))
    #Density 
    return np.vstack((y[1], self.g / (self.f * self.rho0) * (self.bb(z) - self.bn(z))))

  def solve(self):
    r"""
    Solve for the thermal wind overturning streamfunction as a boundary value problem
    based on the system of equations defined in :meth:`pymoc.modules.Psi_Thermwind.ode`.
    """

    # Note: The solution to this BVP is a relatively straightforward integral
    # it would probably be faster to just code it up that way.
    res = integrate.solve_bvp(self.ode, self.bc, self.z, self.sol_init)
    # interpolate solution for overturning circulation onto original grid (and change units to SV)
    self.Psi = res.sol(self.z)[0] / 1e6

