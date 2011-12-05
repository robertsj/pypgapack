"""
pypgapack/examples/example06.py  --  maxint with user mutation.
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
        """
        c = self.GetIntegerChromosome(p, pop) # Get pth string as Numpy array 
        val = np.sum(c)                       #   and sum it up.     
        del c                                 # Delete "view" to internals.
        return float(val)                     # Always return a float.   

    def mutate(self, p, pop, pm) :
        """
        Mutate randomly within -n to n
        """
        n     = self.GetStringLength()
        c     = self.GetIntegerChromosome(p, pop)
        count = 0
        for i in range(0, n) :
            if self.RandomFlip(pm) :
                k = self.RandomInterval(1, 2*n)-n
                c[i] = k
                count += 1
        del c 
        return count 
        
n   = 10                # String length.
opt = MyPGA(sys.argv, PGA.DATATYPE_INTEGER, n, PGA.MAXIMIZE)
opt.SetRandomSeed(1)    # Set random seed for verification.
np.random.seed(1)       # Do the same with Numpy.
u_b =  100*np.ones(n)   # Define lower bound.
l_b = -100*np.ones(n)   # Define upper bound.
# Set the bounds.  Note, need to cast as C-combatible integers.
opt.SetIntegerInitRange(l_b.astype('intc'), u_b.astype('intc')) 
opt.SetMaxGAIterValue(50)    # 50 generations for short output.
opt.SetMutation(opt.mutate)  # Set a custom mutation.
opt.SetUp()                  # Internal allocations, etc.
opt.Run(opt.maxint)          # Set the objective.
opt.Destroy()                # Clean up PGAPack internals.
