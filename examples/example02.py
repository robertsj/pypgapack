"""
pypgapack/examples/example02.py  --  maxint
"""
from pypgapack import PGA
import numpy as np
import sys
class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def maxint(self, p, pop) :
        """
        The maximum integer sum problem.
    
        The alleles are integers, and we solve
            max f(x) = x_1 + x_2 + ... + x_N
        subject to
            |x_i| <= 100 .
        That maximum is f(x) = 100n obtained for x_i = 100 for all i.
        """
        c = self.GetIntegerChromosome(p, pop) # Get pth string as Numpy array 
        val = np.sum(c)                       #   and sum it up.     
        del c                                 # Delete "view" to internals.
        return float(val)                     # Always return a float.   
        
n   = 10                # String length.
# (Command line arguments, integers, string length, and maximize it)
opt = MyPGA(sys.argv, PGA.DATATYPE_INTEGER, n, PGA.MAXIMIZE)
opt.SetRandomSeed(1)    # Set random seed for verification.
u_b =  100*np.ones(n)   # Define lower bound.
l_b = -100*np.ones(n)   # Define upper bound.
# Set the bounds.  Note, need to cast as C-combatible integers.
opt.SetIntegerInitRange(l_b.astype('intc'), u_b.astype('intc')) 
opt.SetMaxGAIterValue(50)    # 50 generations for short output.
opt.SetUp()                  # Internal allocations, etc.
opt.Run(opt.maxint)          # Set the objective.
opt.Destroy()                # Clean up PGAPack internals.
