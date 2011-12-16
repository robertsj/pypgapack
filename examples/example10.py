"""
pypgapack/examples/example10.py  --  optimize a reflected 8 slab reactor
"""

from pypgapack import PGA
from mpi4py import MPI
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import factorial
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

class MyPGA(PGA) :
    """  
    Derive our own class from PGA.
    """
    def f1(self, p, pop) :
        """
        Minimize peaking and maximize keff using weighted objective.
        """
        pattern = self.GetIntegerChromosome(p, pop) 
        keff, peak = self.flux(pattern)
        del pattern
        delta = 0
        if peak > 1.8 :
            delta = peak - 1.8
        val = 1.0*keff - 10.0*delta# / 1.1 + peak / 1.5 
        self.evals+=1
        return val
       
    def f2(self, p, pop) :
        """
        Minimize peaking and maximize keff using weighted objective.
        """
        pattern = self.GetIntegerChromosome(p, pop) 
        keff, peak = self.flux(pattern)
        del pattern
        delta = 0
        if peak > 1.5 :
            delta = peak - 1.5
        val = 2.0*keff - 1.0*delta# / 1.1 + peak / 1.5 
        self.evals+=1
        return val
    
    def flux(self, pattern, plot=0) :
        """
        Solve for the flux via simple mesh-centered finite differences.
        
        Returns keff and the maximum-to-average assembly-averaged fission density ratio. 
        """
        # Coarse mesh boundaries.
        #print np.arange(0, 4)
        #coarse = np.arange(0, 3 + self.number_slabs) * 20.0
        coarse = np.array([  0.0, 20.0, 40.0, 60.0, 80.0,100.0,120.0,140.0,160.0,180.0,200.0] )
        # Fine meshes per coarse mesh.
        fine = 20
        dx = 20.0 / fine
        # Material map (simplifies coefficients)
        mat = np.zeros(fine * (len(coarse) - 1))
        mat[0:fine] = 10 # reflector
        j = fine
        for i in range(0, self.number_slabs) :
            mat[j:(j + fine)] = pattern[i] # one of the slabs
            j += fine 
        mat[j:(j + fine)] = 10 # reflector
        
        # System size 
        n = fine * (len(coarse) - 1)
        
        # Materials
        D, R, F, S = self.materials()
        
        # Coefficient Matrix (just diagonals) and Vectors
        AD = np.zeros((n, 2))
        AL = np.zeros((n, 2))
        AU = np.zeros((n, 2))
        nufission = np.zeros((n, 2))
        scatter12 = np.zeros(n)
        
        for g in range(0, 2) :
            # reflective left boundary
            b = 1 
            AU[0, g] = -2.0 * D[g, mat[0]] * D[g, mat[0]] / \
                       (D[g, mat[0]] * dx + D[g, mat[1]] * dx)
            AD[0, g] = -AU[0, g] + 2.0 * D[g, mat[0]] * (1.0 - b) / \
                       (4.0 * D[g, mat[0]] * (1.0 + b) + dx * (1.0 - b)) + \
                       dx * R[g, mat[0]] 
            nufission[0, g] = dx * F[g, mat[0]]
            scatter12[0] = dx * S[mat[0]]
            # vacuum right boundary
            b = 0 
            AL[n - 2, g] = -2.0 * D[g, mat[n - 1]] * D[g, mat[n - 2]] / \
                           (D[g, mat[n - 1]] * dx + D[g, mat[n - 2]] * dx)
            AD[n - 1, g] = -AL[n - 2, g] + 2.0 * D[g, mat[n - 1]] * (1.0 - b) / \
                       (4.0 * D[g, mat[n - 1]] * (1.0 + b) + dx * (1.0 - b)) + \
                       dx * R[g, mat[n - 1]]
            nufission[n - 1, g] = dx * F[g, mat[n - 1]]
            scatter12[n - 1] = dx * S[mat[n - 1]]
            # internal cells                                              
            for i in range(1, n - 1) :
                AL[i - 1, g] = -2.0 * D[g, mat[i]] * D[g, mat[i - 1 ]] / \
                               (D[g, mat[i]] * dx + D[g, mat[i - 1]] * dx)
                AU[i, g] = -2.0 * D[g, mat[i]] * D[g, mat[i + 1 ]] / \
                           (D[g, mat[i]] * dx + D[g, mat[i + 1]] * dx)   
                
                AD[i, g] = -(AL[i - 1, g] + AU[i, g]) + dx * R[g, mat[i]]
                nufission[i, g] = dx * F[g, mat[i]]
                scatter12[i] = dx * S[mat[i]]

        # Initiate the fluxes.
        phi1 = np.zeros(n)
        phi2 = np.zeros(n)
        # Use a guess for fission density based on fission cross section.
        fission_density = nufission[:, 1]
        # and normalize it.
        fission_density = fission_density / \
                           np.sqrt(np.sum(fission_density ** 2))                          
        fission_density0 = np.zeros(n)
        # Initialize the downscatter source.
        scatter_source = np.zeros(n)
        # Initial eigenvalue guess.
        keff = 1
        keff0 = 0
        # Set errors.
        errorfd = 1.0
        errork = 1.0
        it = 0
        while (errorfd > 1e-5 and errork > 1e-5 and it < 200) :
            # Solve fast group.
            self.tridiag(AU[:, 0], AL[:, 0], AD[:, 0], fission_density / keff, phi1)
            # Compute down scatter source.
            scatter_source = phi1 * scatter12
            # Solve thermal group.
            self.tridiag(AU[:, 1], AL[:, 1], AD[:, 1], scatter_source, phi2)
            # Keep old values.
            fission_density0[:] = fission_density
            keff0 = keff
            # Update density and eigenvalue
            fission_density[:] = phi1 * nufission[:, 0] + \
                                 phi2 * nufission[:, 1] 
            keff = keff0 * np.sum(fission_density) / \
                             np.sum(fission_density0)
            # Update errors.  Use Linf norm on density.
            errorfd = np.max(np.abs(fission_density - fission_density0))
            errork = np.abs(keff - keff0)
            it += 1
  
        # Now we average the fission density over each fueled coarse mesh.
        slab_fission_density = np.zeros(self.number_slabs)
        j = fine
        for i in range(0, self.number_slabs) :
            slab_fission_density[i] = np.mean(fission_density[j:(j + fine) - 1])
            j += fine 
        mean_fission_density = np.mean(fission_density)  
        peaking = slab_fission_density / mean_fission_density  
        max_peaking = np.max(peaking)   
        # plot if desired
        if plot == 1 :
            xf = np.linspace(0.5 * dx, coarse[len(coarse) - 1] - 0.5 * dx, n)
            xc = coarse[1:len(coarse) - 1]
            plt.plot(xf, fission_density, 'b', lw=2)
            plt.step(xc, np.hstack((slab_fission_density[0], \
                     slab_fission_density)), 'g', lw=2)
            plt.title('Fission density: peak = ' + str(max_peaking))
            plt.grid(True)
            plt.show()
        return keff, max_peaking
        
    def tridiag(self, U, L, D, f, y):
        """ 
        Tridiagonal solver.  
        
        This assumes vectors U, L, and D are of the same length.  The right
        hand side is f and the solved unknowns are returned in y.
        """
        N = len(D)
        w = np.zeros(N)
        v = np.zeros(N)
        z = np.zeros(N)
        w[0] = D[0]
        v[0] = U[0] / w[0]
        z[0] = f[0] / w[0]   
        for i in range(1, N) :
            w[i] = D[i] - L[i - 1] * v[i - 1]
            v[i] = U[i] / w[i]
            z[i] = (f[i] - L[i - 1] * z[i - 1]) / w[i]
        y[N - 1] = z[N - 1]
        for i in range(N - 2, -1, -1) :
            y[i] = z[i] - v[i] * y[i + 1]
    
    def materials(self):
        """
        10 fuels with one 1 reflectors by row.  Represents burnup of 0, 5, 10, 15, 
        20, 25, 30, 35, 40, and 45 MWd/kg for one assembly type.
        """
        assert(self.GetPopSize() <= 10) # Only have 10 slab materials.
        D = np.array([[1.4402e+00, 1.4429e+00, 1.4453e+00, 1.4467e+00, 1.4476e+00, \
                       1.4483e+00, 1.4489e+00, 1.4496e+00, 1.4507e+00, 1.4525e+00, 1.3200e+00], \
                      [3.7939e-01, 3.7516e-01, 3.7233e-01, 3.7045e-01, 3.6913e-01, \
                       3.6818e-01, 3.6749e-01, 3.6699e-01, 3.6649e-01, 3.6615e-01, 2.6720e-01]])
        R = np.array([[2.5800e-02, 2.5751e-02, 2.5755e-02, 2.5840e-02, 2.5958e-02, \
                       2.6090e-02, 2.6226e-02, 2.6362e-02, 2.6559e-02, 2.6802e-02, 2.5700e-02], \
                      [1.1817e-01, 1.2301e-01, 1.2306e-01, 1.2277e-01, 1.2223e-01, \
                       1.2136e-01, 1.2017e-01, 1.1871e-01, 1.1622e-01, 1.1275e-01, 5.1500e-02]])
        F = np.array([[7.9653e-03, 7.6255e-03, 7.2724e-03, 6.9344e-03, 6.6169e-03, \
                       6.3189e-03, 6.0453e-03, 5.7908e-03, 5.4413e-03, 5.0379e-03, 0.00000000], \
                      [1.6359e-01, 1.7301e-01, 1.7681e-01, 1.7634e-01, 1.7355e-01, \
                       1.6931e-01, 1.6424e-01, 1.5869e-01, 1.5005e-01, 1.3894e-01, 0.00000000]])
        S = np.array([1.5204e-02, 1.5152e-02, 1.4958e-02, 1.4828e-02, 1.4744e-02, \
                      1.4691e-02, 1.4652e-02, 1.4626e-02, 1.4601e-02, 1.4587e-02, 2.3100e-02])
                      
        return D, R, F, S
    
    def kinf(self, pattern) :
        """
        Compute kinf for each slab.
        """
        D, R, F, S = self.materials()
        k = np.zeros(len(pattern))   
        for i in range(0, len(pattern)) :
            k[i] = (F[0, i] + F[1, i] * S[i] / R[1, i]) / R[0, i] 
        return k
            
    def htbx(self, p1, p2, pop1, c1, c2, pop2) :
        """
        Heuristic tie-breaking cross-over.  See Carter for details. Actually, for this problem,
        htbx is equivalent to tbx, since the materials are already ordered by reactivity.
        """
        # Grab the city id's.
        paren1 = self.GetIntegerChromosome(p1, pop1)
        paren2 = self.GetIntegerChromosome(p2, pop1)
        child1 = self.GetIntegerChromosome(c1, pop2)
        child2 = self.GetIntegerChromosome(c2, pop2) 
         
        # Copy the parents to temporary vector for manipulation.
        n = self.GetStringLength()
        parent1 = np.zeros(n)
        parent2 = np.zeros(n)
        for i in range(0, n) :
            parent1[i] = paren1[i]
            parent2[i] = paren2[i]               
        
        # Code the parents using "position listing".
        code1 = np.zeros(n)
        code2 = np.zeros(n)
        for i in range(0, n) :
            code1[parent1[i]] = i + 1
            code2[parent2[i]] = i + 1
        
        # Randomly choose two cross-over points. 
        perm = np.random.permutation(n)
        point1 = np.min(perm[0:2])
        point2 = np.max(perm[0:2]) + 1 
        
        # Exchange all alleles between the two points.
        temp = np.zeros(point2 - point1)
        for i in range(point1, point2) :
            temp[i - point1] = parent1[i]
            parent1[i] = parent2[i]
            parent2[i] = temp[i - point1] 
        
        # Generate a cross-over map, a random ordering of the 0,1,...,n-1
        crossovermap = np.random.permutation(n) 
        
        # Multiply each allele of the strung by n and add the map.
        parent1 = parent1 * n + crossovermap
        parent2 = parent2 * n + crossovermap
        
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
                
    def eog(self):
        """
        Randomly set a bit to 1 in each string 
        """
        best                = self.GetBestIndex(PGA.OLDPOP)
        bestpattern         = self.GetIntegerChromosome(best, PGA.OLDPOP)
        iter                = self.GetGAIterValue()
        keff, peak          = opt.flux(bestpattern)
        self.besteval[iter-1] = self.GetEvaluation(best, PGA.OLDPOP)                                                
        self.bestkeff[iter-1] = keff
        self.bestpeak[iter-1] = peak        
                                   
comm = MPI.COMM_WORLD           # Get the communicator.
rank = comm.Get_rank()          # Get my rank.  
t_start = MPI.Wtime()           # Start the clock.                                             
n = 8                           # Number of fueled slabs (1-10)
opt = MyPGA(sys.argv, PGA.DATATYPE_INTEGER, n, PGA.MAXIMIZE)
#opt.SetRandomSeed(1)                    # Set seed for verification.
#np.random.seed(1)                       # Do the same with Numpy.
opt.SetIntegerInitPermute(0, n - 1)     # Start with random permutations. 
opt.SetPopSize(50)                      # Large enough to see some success.
opt.SetMaxGAIterValue(300)              # Small number for output.
opt.SetCrossover(opt.htbx)              # Set a cross-over operation.
opt.SetEndOfGen(opt.eog)                # End of generation info
opt.SetMutation(PGA.MUTATION_PERMUTE)   # Mutate by permutation.
opt.SetNoDuplicatesFlag(PGA.TRUE)       # Keep no duplicate patterns.
opt.SetUp()                             # Internal allocations, etc.
opt.number_slabs = n
opt.evals = 0
opt.besteval = np.zeros(301)
opt.bestkeff = np.zeros(301)
opt.bestpeak = np.zeros(301)
opt.Run(opt.f2)                         # Set the objective and run.
# For f1, the reference solution is
    
MPI.Finalize()  
opt.Destroy()


