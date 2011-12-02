"""
pypgapack/examples/example09.py  --  traveling salesman in parallel
"""
from pypgapack import PGA
from mpi4py import MPI
import numpy as np
import sys

class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def tsm(self, p, pop) :
        """
        Hamiltonian traveling salesman problem.
        """
        c = self.GetIntegerChromosome(p, pop) # Get pth string as Numpy array 
        val = self.distance(c)
        del c
        return val
    
    def distance(self, c) :
        """
        Compute the total distance for a set of cities.
        """
        # x and y coordinates by city
        x = np.array([54.0,54.0,37.0,41.0, 2.0, 7.0,25.0,22.0,18.0, 4.0,\
                      13.0,18.0,24.0,25.0,44.0,41.0,45.0,58.0,62.0,82.0,\
                      91.0,83.0,71.0,64.0,68.0,83.0,87.0,74.0,71.0,58.0])
        y = np.array([67.0,62.0,84.0,94.0,99.0,64.0,62.0,60.0,54.0,50.0,\
                      40.0,40.0,42.0,38.0,35.0,26.0,21.0,35.0,32.0, 7.0,\
                      38.0,46.0,44.0,60.0,58.0,69.0,76.0,78.0,71.0,69.0])
        n   = self.GetStringLength()   
        val = 0.0
        for i in range(0, n-1) :
            val += np.sqrt( (x[c[i]]-x[c[i+1]])**2 + (y[c[i]]-y[c[i+1]])**2 )
        val += np.sqrt( (x[c[0]]-x[c[n-1]])**2 + (y[c[0]]-y[c[n-1]])**2 )
        assert(val > 423.70) # Debug.
        return val
    
    def tbx(self, p1, p2, pop1, c1, c2, pop2) :
        """
        Tie-breaking cross-over.  See Poon and Carter for details.
        """
        # Grab the city id's by direct access to the memory.  This way,
        #   we can modify the contents directly, avoiding the "set"
        #   functions.
        paren1 = self.GetIntegerChromosome(p1, pop1)
        paren2 = self.GetIntegerChromosome(p2, pop1)
        child1  = self.GetIntegerChromosome(c1, pop2)
        child2  = self.GetIntegerChromosome(c2, pop2) 
        assert(np.sum(paren1)==435) # Ensure we haven't
        assert(np.sum(paren2)==435) # in
         
        # String length.
        n = self.GetStringLength()
        parent1 = np.zeros(n)
        parent2 = np.zeros(n)
                     
        for i in range(0, n) :
            parent1[i] = paren1[i]
            parent2[i] = paren2[i]               
        
        # Code the parents using "position listing".  For example, a string
        #   (2,1,3,0) 
        # is represented as 
        #   [4,2,1,3] .
        # Be careful, the algorithm is most straightforward using 1-indexing 
        # rather than 0-indexing, but the id's are stored starting at 0.
        code1   = np.zeros(n)
        code2   = np.zeros(n)
        for i in range(0, n) :
            code1[parent1[i]] = i + 1
            code2[parent2[i]] = i + 1
        
        # Randomly choose two cross-over points.  (Is there a Numpy 
        #   function to select two that aren't the same?)
        perm   = np.random.permutation(n)
        point1 = np.min(perm[0:2])
        point2 = np.max(perm[0:2])+1 # a:b does a up *to* b; want inclusive
        
        # Exchange all alleles between the two points.
        temp = np.zeros(point2-point1)
        for i in range(point1, point2) :
            temp[i-point1] = parent1[i]
            parent1[i]     = parent2[i]
            parent2[i]     = temp[i-point1] 
        
        # Generate a cross-over map, a random ordering of the 0,1,...,n-1
        crossovermap = np.random.permutation(n) # could use old one
        
        # Multiply each allele of the strung by n and add the map.
        parent1 = parent1*n + crossovermap
        parent2 = parent2*n + crossovermap
        
        # Replace the lowest allele by 0, the next by 1, up to n-1.  Here,
        #   we sort the parents first, and then for each element, find
        #   where the increasing values are found in the original.  There
        #   is probably a simpler set of functions built in somewhere.
        sort1 = np.sort(parent1)
        sort2 = np.sort(parent2)
        for i in range(0, n) :
            index = np.where(parent1 == sort1[i])
            parent1[index[0][0]] = i
            index = np.where(parent2 == sort2[i])
            parent2[index[0][0]] = i 
        
        tmpc1 = np.zeros(n)
        tmpc2 = np.zeros(n)
        # Map the string back to elements.  These are the offspring.
        for i in range(0, n) :
            tmpc1[parent1[i]] = i
            tmpc2[parent2[i]] = i            
        for i in range(0, n) :
            child1[i] = tmpc1[i]
            child2[i] = tmpc2[i]       
  
    def swap(self, p, pop, pm) :
        """
        Random swap of allele pairs.  Note, nswap must be set!
        """
        n     = self.GetStringLength()
        c     = self.GetIntegerChromosome(p, pop)
        index = np.random.permutation(n)
        for i in range(0, self.nswap) :
            i1 = index[2*i  ]
            i2 = index[2*i+1]
            tmp1      = c[i1]
            tmp2      = c[i2]
            c[i1]     = tmp2  
            c[i2]     = tmp1
        del c 
        return 0 
            
    def init(self, p, pop) :
        """
        Random initial states.  We do this so that we can enforce the same
        initial guesses for all runs to compare against the Poon and Carter.
        """
        n = self.GetStringLength()
        c = self.GetIntegerChromosome(p, pop)
        np.random.seed(p)
        perm = np.random.permutation(n)
        for i in range(0, n) :
            c[i] = perm[i]
        del c
         
comm = MPI.COMM_WORLD           # Get the communicator.
rank = comm.Get_rank()          # Get my rank.  
t_start = MPI.Wtime()           # Start the clock. 
n        = 30                   # Number of cities.              
numrun   = 25                   # Number of runs to average.
besteval = np.zeros(numrun)
for i in range(0, numrun) :
    opt = MyPGA(sys.argv, PGA.DATATYPE_INTEGER, n, PGA.MINIMIZE)
    opt.SetInitString(opt.init) # Set an initialization operator.
    opt.SetCrossoverProb(0.8)   #
    opt.SetPopSize(22)          # 22 rather than 21
    opt.SetNumReplaceValue(21)  # Keep the best 22-21 = 1 string = elitist.
    opt.SetMaxGAIterValue(300)  # 300 generations, like the reference.
    opt.SetCrossover(opt.tbx)   # Set a cross-over operation.
    opt.SetMutation(opt.swap)   # Set a mutate operation. 
    opt.nswap = 3               # Number of pairs to swap in mutation.
    opt.SetUp()                 # Internal allocations, etc.
    opt.Run(opt.tsm)            # Set the objective and run.
    if rank == 0 :
        best        = opt.GetBestIndex(PGA.OLDPOP)
        besteval[i] = opt.GetEvaluation(best, PGA.OLDPOP)
    opt.Destroy()              # Clean up PGAPack internals.
    
if rank == 0 :
    print "  MEAN: ", np.mean(besteval) # Print out the mean
    print " SIGMA: ", np.std(besteval)  # and standard deviation.
    print "   MIN: ", np.min(besteval)  # Print out the mean
    print "   MAX: ", np.max(besteval)  # Print out the mean
    t_end = MPI.Wtime()
    print "Elapsed time = ", t_end-t_start, " seconds."
    
MPI.Finalize() # Should be called automatically, but good practice.