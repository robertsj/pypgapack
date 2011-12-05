"""
pypgapack/examples/example10.py  --  optimize a reflected 8 slab reactor
"""

from pypgapack import PGA
import numpy as np
import sys

class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def objective(self, p, pop) :
        """
        Minimize peaking and maximize keff using weighted objective.
        """
        pattern = self.GetIntegerChromosome(p, pop) 
        val = self.distance(pattern)
        del pattern
        return val
    
    def flux(self, c) :
        """
        Solve for the flux.
        """
        # Coarse mesh boundaries.
        coarse = np.array([  0.0, 20.0, 40.0, 60.0, 80.0,100.0, \
                           120.0,140.0,160.0,180.0,200.0,220.0] )
        # Fine meshes per coarse mesh.
        fine = 20
        
        # Coefficient Matrix (just diagonals)
        AD = np.zeros(n,numg);
        AL = np.zeros(n,numg);
        AU = np.zeros(n,numg);
        
        for g in range(0, 1) :
            
            b = 1 # reflective left boundary
            AU[1,g]   = -2*dd(mt(1),g)*dd(1,g)/(dd(mt(1),g)*dx+dd(mt(1),g)*dx);
            AD[1,g]   = -AU(1,g) + 2*dd(mt(1),g)*(1-b)/  \
                        (4*dd(mt(1),g)*(1+b)+dx(1)*(1-b)) + dx(1)*aa(mt(1),g); 
            b = 0 # vacuum right boundary
            AL[n-1,g] = -2*dd(mt(end),g)*dd(mt(end),g)/(dd(mt(end),g)*dx+dd(mt(end),g)*dx);
            AD[n,g]   = -AL(n-1,g) + 2*dd(mt(end),g)*(1-b)/  \
                         (4*dd(mt(end),g)*(1+b)+dx*(1-b)) + dx*aa(mt(end),g);                                    
            for i in range(1, N-2) :
                AL(i-1,g) = -2*dd(mt(i-1),g)*dd(i,g)/(dd(mt(i-1),g)*dx(i)+dd(i,g)*dx(i-1));
                AU(i,g)   = -2*dd(i+1,g)*dd(i,g)/(dd(i+1,g)*dx(i)+dd(i,g)*dx(i+1));   
                AD(i,g)   = -( AL(i-1,g) + AU(i,g) ) + dx(i,1)*aa(i,g);

        
    def tridiag(self, U, L, C, f):
        """ 
        tridiagonal solver. c=upper, b=lower, a=central.  Must be same length.
        """
        N = len(C)
        w = np.zeros(N)
        v = w
        z = w 
        y = w 
        w[0] = C[0]
        v[0] = U[0]/w[0]
        z[0] = f[0]/w[0]
        for i in range(1, N) :
            w[i] = C[i] - L[i] * v[i-1]
            v[i] = U[i] / w[i]
            z[i] = (f[i] - L[i] * z[i])/w[i]

        y[N-1] = z[N-1]
        for i in range(N-1, 0) :
            y[i] = z[i] - v[i] * y[i+1]
            
    
    def materials(self):
        """
        8 fuels and 1 reflector by row.  Represents burnup of 0, 5, 10, 15, 
        20,25,30, and 35 MWd/kg for one assembly type.
        """
        data = np.array([\
            [1.4402,0.37939,0.025800,1.1817e-01,7.9653e-03,1.6359E-01,1.5204e-02],\
            [1.4429,0.37516,0.025751,1.2301e-01,7.6255e-03,1.7301E-01,1.5152e-02],\
            [1.4453,0.37233,0.025755,1.2306e-01,7.2724e-03,1.7681e-01,1.4958e-02],\
            [1.4467,0.37045,0.025840,1.2277e-01,6.9344e-03,1.7634e-01,1.4828e-02],\
            [1.4476,0.36913,0.025958,1.2223e-01,6.6169e-03,1.7355e-01,1.4744e-02],\
            [1.4483,0.36818,0.026090,1.2136e-01,6.3189e-03,1.6931e-01,1.4691e-02],\
            [1.4489,0.36749,0.026226,1.2017e-01,6.0453e-03,1.6424e-01,1.4652e-02],\
            [1.4496,0.36699,0.026362,1.1871e-01,5.7908e-03,1.5869e-01,1.4626e-02],\
            [1.3200,0.2672, 0.025700,0.05150000,0.00000000,0.00000000,0.02310000]])
        return data
    
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
    
    def htbx(self, p1, p2, pop1, c1, c2, pop2) :
        """
        Heuristic tie-breaking cross-over.  See Carter for details.
        """
        # Grab the city id's.
        paren1 = self.GetIntegerChromosome(p1, pop1)
        paren2 = self.GetIntegerChromosome(p2, pop1)
        child1  = self.GetIntegerChromosome(c1, pop2)
        child2  = self.GetIntegerChromosome(c2, pop2) 
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
        
        # Exchange all alleles between the two points.
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

opt.SetRandomSeed(1)        # Set seed for verification.
np.random.seed(1)           # Do the same with Numpy.
opt.SetIntegerInitPermute(0, n-1) # Start with random permutations. 
opt.SetPopSize(400)         # Large enough to see some success.
opt.SetMaxGAIterValue(100)  # Small number for output.
opt.SetCrossover(opt.tbx)   # Set a cross-over operation.
opt.SetMutation(PGA.MUTATION_PERMUTE) # Mutate by permutation.
opt.SetUp()                 # Internal allocations, etc.
opt.Run(opt.tsm)            # Set the objective and run.
opt.Destroy()               # Clean up PGAPack internals.

