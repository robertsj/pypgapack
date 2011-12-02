"""
pypgapack/examples/example08.py  --  parallel maxbit
"""
from pypgapack import PGA
from mpi4py import MPI
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

comm = MPI.COMM_WORLD           # Get the communicator.
rank = comm.Get_rank()          # Get my rank.  
if rank == 0 :                  # Just to show it works, have node 0
    seed = 1                    #   set seed=1 and n=500 and have all
    n    = 500
else :                          #   other nodes set seed=n=0.  Then,
    seed = 0                    #   broadcast them to all nodes with
    n    = 0
seed = comm.bcast(seed, root=0) #   node 0 as the root process,
n    = comm.bcast(n,    root=0) #   and verify by printing.
print " node=", rank, " seed=", seed, " n=", n

if rank == 0 :
    t_start = MPI.Wtime()  # Start the clock.
# (Command line arguments, 1's and 0's, string length, and maximize it)
opt  = MyPGA(sys.argv, PGA.DATATYPE_BINARY, n, PGA.MAXIMIZE)
opt.SetPopSize(1000)
opt.SetRandomSeed(seed)        # Set random seed for verification.
opt.SetMaxGAIterValue(50)      # 50 generations for short output.
opt.SetUp()                    # Internal allocations, etc.
opt.Run(opt.maxbit)            # Set the objective.
opt.Destroy()                  # Clean up PGAPack internals
if rank == 0 :
    t_end = MPI.Wtime()
    print "Elapsed time = ", t_end-t_start, " seconds."
MPI.Finalize() # Should be called automatically, but good practice.
