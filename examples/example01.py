"""
pypgapack/examples/example01.py  --  maxbit
"""
from pypgapack import PGA
import sys
class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def maxbit(self, p, pop) :
        """
        Maximum when all alleles are 1's, and that maximum is n. 
        """
        val = 0
        # Size of the problem
        n = self.GetStringLength()               
        for i in range(0, n) :               
            # Check whether ith allele in string p is 1
            if self.GetBinaryAllele(p, pop, i) : 
                val = val + 1
        # Remember that fitness evaluations must return a float
        return float(val)      
                
# (Command line arguments, 1's and 0's, string length, and maximize it)
opt = MyPGA(sys.argv, PGA.DATATYPE_BINARY, 100, PGA.MAXIMIZE)
opt.SetRandomSeed(1)       # Set random seed for verification.
opt.SetMaxGAIterValue(50)  # 50 generations (default 1000) for short output.
opt.SetUp()                # Internal allocations, etc.
opt.Run(opt.maxbit)        # Set the objective.
opt.Destroy()              # Clean up PGAPack internals

