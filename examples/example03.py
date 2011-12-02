"""
pypgapack/examples/example03.py  --  maxreal
"""
from pypgapack import PGA
import numpy as np
import sys

class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def maxreal(self, p, pop) :
        """
        The maximum real sum problem.
    
        The alleles are doubles, and we solve
        .. math:: 
          \max f(x) &= \sum^N_{n=1} x_n \\
              s.t.  &= |x_i| \leq 100
        That maximum is  :math:`f_{\textrm{max}}(x) = 100n` obtained for
        :math:`x_i = 100, i = 1\ldots N`.
        """
        c = self.GetRealChromosome(p, pop) # Get pth string as Numpy array 
        val = np.sum(c)                    #   and sum it up.     
        del c                              # Delete "view" to internals.
        return val                         # Already a float.   
        
n   = 10              # String length.
# (Command line arguments, doubles, string length, and maximize it)
opt = MyPGA(sys.argv, PGA.DATATYPE_REAL, n, PGA.MAXIMIZE)
opt.SetRandomSeed(1)  # Set random seed for verification.
u_b = 100*np.ones(n)  # Define lower bound.
l_b =-100*np.ones(n)  # Define upper bound.
# Set the bounds.  Default floats are handled without issue.
opt.SetRealInitRange(l_b, u_b)
# Force mutations to keep values in the initial range, a useful
#   feature for bound constraints.
opt.SetMutationType(PGA.MUTATION_RANGE)
opt.SetMaxGAIterValue(50) # 50 generations for short output.
opt.SetUp()               # Internal allocations, etc.
opt.Run(opt.maxreal)      # Set the objective.
opt.Destroy()             # Clean up PGAPack internals.
