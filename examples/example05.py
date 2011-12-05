"""
pypgapack/examples/example05.py  --  traveling salesman
"""
from pypgapack import PGA
import numpy as np
import sys

class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def tsm(self, p, pop) :
        """
        Oliver's 30 city traveling salesman problem.
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
        assert(val > 423.70) # DEBUG.
        return val
    
    def tbx(self, p1, p2, pop1, c1, c2, pop2) :
        """
        Tie-breaking cross-over.  See Poon and Carter for details. 
        """
        # Grab the city id's.
        paren1 = self.GetIntegerChromosome(p1, pop1)
        paren2 = self.GetIntegerChromosome(p2, pop1)
        child1 = self.GetIntegerChromosome(c1, pop2)
        child2 = self.GetIntegerChromosome(c2, pop2) 
        assert(np.sum(paren1)==435) # DEBUG
        assert(np.sum(paren2)==435) # DEBUG
         
        # Copy the parents to temporary vector for manipulation.
        n = self.GetStringLength()
        parent1 = np.zeros(n)
        parent2 = np.zeros(n)
        for i in range(0, n) :
            parent1[i] = paren1[i]
            parent2[i] = paren2[i]               
        
        # Code the parents using "position listing".
        code1   = np.zeros(n)
        code2   = np.zeros(n)
        for i in range(0, n) :
            code1[parent1[i]] = i + 1
            code2[parent2[i]] = i + 1
        
        # Randomly choose two cross-over points. 
        perm   = np.random.permutation(n)
        point1 = np.min(perm[0:2])
        point2 = np.max(perm[0:2])+1 
        
        # Exchange all alleles between the two points.  (It's unclear to me
        #   whether these points should be inclusive or not; here, they are.)
        temp = np.zeros(point2-point1)
        for i in range(point1, point2) :
            temp[i-point1] = parent1[i]
            parent1[i]     = parent2[i]
            parent2[i]     = temp[i-point1] 
        
        # Generate a cross-over map, a random ordering of the 0,1,...,n-1
        crossovermap = np.random.permutation(n) 
        
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
            
        # Map the string back to elements.  These are the offspring.
        tempchild1 = np.zeros(n)
        tempchild2 = np.zeros(n)
        for i in range(0, n) :
            tempchild1[parent1[i]] = i
            tempchild2[parent2[i]] = i            
        for i in range(0, n) :
            child1[i] = tempchild1[i]
            child2[i] = tempchild2[i]       
  
# Number of cities.                    
n   = 30    

opt = MyPGA(sys.argv, PGA.DATATYPE_INTEGER, n, PGA.MINIMIZE) 
    
# One possible benchmark solution
reference = np.array([ 0, 2, 3, 4, 5, 6, 7, 8, 9,10, \
                      11,12,13,14,15,16,17,18,19,20, \
                      21,22,24,23,25,26,27,28,29, 1] )
print "Reference distance is: ", opt.distance(reference)

opt.SetRandomSeed(1)                  # Set seed for verification.
np.random.seed(1)                     # Do the same with Numpy.
opt.SetIntegerInitPermute(0, n-1)     # Start with random permutations. 
opt.SetPopSize(400)                   # Large enough to see some success.
opt.SetMaxGAIterValue(100)            # Small number for output.
opt.SetCrossover(opt.tbx)             # Set a cross-over operation.
opt.SetMutation(PGA.MUTATION_PERMUTE) # Mutate by permutation.
opt.SetUp()                           # Internal allocations, etc.
opt.Run(opt.tsm)                      # Set the objective and run.
opt.Destroy()                         # Clean up PGAPack internals.
