# -*- coding: utf-8 -*-

import unittest
import numpy as np
import mod1
from numpy.testing import assert_allclose

class PolarangleTest(unittest.TestCase):
  
    def test1(self):
      """
      Polar angle of r=[0,1,0] equals to pi/2
      """
      r = np.array([0,1,0])
      self.assertEqual(mod1.polarangle(r), np.pi/2, msg="""polar angle is 
                       incorrect, error1""")
      
    def test2(self):
      """
      Polar angle of r=[1,0,0] equals to 0
      """
      r = np.array([1,0,0])
      self.assertEqual(mod1.polarangle(r), 0, msg="""polar angle is incorrect, 
                       error2""")
      
    def test3(self):
      """
      Polar angle of r=[-1,1,0] equals to 3pi/4
      """
      r = np.array([-1,1,0])
      self.assertEqual(mod1.polarangle(r), 3*np.pi/4, msg="""polar angle is 
                       incorrect, error3""")
      
      
class DistributionTest(unittest.TestCase):
    """
    Check of equality of the distribution function to Maxwellian
    distribution function on the heliosphere boundary
    """
  
    def test1(self):
      Vin, dV, Vfin = -2.45, 0.06, 2.47
      N = int((Vfin-Vin)/dV)
      Vfin = Vin + N*dV
      arr_maxw = np.zeros((N, N))
      arr_func = np.zeros((N, N))
      Rs = 100
      beta = 1e-7
      Vs = np.array([0, 0, -1])
      for i in range(N):
        u = Vin + i*dV
        for j in range(N):
            w = Vin + j*dV
            v = np.array([u, 0, w])
            arr_func[i][j] = mod1.distrfunc2(np.array([0, 0, Rs]), v, Vs, 1.2, beta, Rs)
            par = (Vs[0] - v[0])**2 + (Vs[1] - v[1])**2 + (Vs[2] - v[2])**2 
            arr_maxw[i][j] = mod1.norm*np.exp(-mod1.m*par*mod1.V**2/2/mod1.kB/mod1.Ts)

      assert_allclose(arr_func, arr_maxw, atol=1e-10)      
