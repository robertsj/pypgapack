.. _sec_examples:

Examples
========

.. highlight:: python
   :linenothreshold: 25

All examples are located in the ``pypgapack/examples`` and the 
reference output for all examples is in ``pypgapack/examples/output``.
Aside from small floating point differences, the values should be
the same give the use of a fixed random number generator seed in
all the examples.  A utility script ``run_examples.py`` is included
to test user output to the included reference cases.

Also, the maximum generation count is limited to
50 for all cases to produce short output.  Experiment with that limit
to see better solutions.


.. _sec_basicexamples:

Basic Examples
--------------

The following are some simple examples that illustrate
the basic pgapack functionality.



.. _sec_maxbitexample:

Example 1: MAXBIT
"""""""""""""""""

:mod:`pypgapack` is pretty easy to use, and to demonstrate, we'll solve
the maxbit problem, the first example in the pgapack documentation.

.. literalinclude:: ../../examples/example01.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example01_ref
  


.. _sec_maxintexample:

Example 2: MAXINT
"""""""""""""""""

This is a similar problem, but the alleles are integers ranging
from -100 to 100.  Note that when the integer ranges are
set, a cast to ``intc`` is used.  Python uses high precision
datatypes, and there doesn't seem to be a safe implicit conversion
between the Python integer type and the C integer type behind
the scenes (in SWIG land).  Casting
explicitly circumvents the issue.

.. literalinclude:: ../../examples/example02.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example02_ref



.. _sec_maxrealexample:

Example 3: MAXREAL
""""""""""""""""""

This is the same problem, but for real alleles.  Here, note use
of :meth:`SetMutationalType` with the option of :const:`PGA.MUTATION_RANGE`.
This forces mutated allele values to remain within the initial range
specified, useful for cases with constrained inputs.  The default adds
some (small) random amount, but over many iterations, this can cause allele
values to go significantly beyond the initial range.

.. literalinclude:: ../../examples/example03.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example03_ref



.. _sec_useroperatorexamples:

Examples of User Defined Operators
----------------------------------

These examples explore one of the strengths of PGAPack,
namely user-defined operators.  



.. _sec_initstringexamples:

Example 4: User-defined String Initialization
"""""""""""""""""""""""""""""""""""""""""""""

We redo :ref:`sec_maxintexample` by initializing the strings
with our own routine.  Here, that's done by generating a 
permutation using Numpy.

.. literalinclude:: ../../examples/example04.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example04_ref



.. _sec_crossoverexamples:

Example 5: User-defined Crossover Operator
""""""""""""""""""""""""""""""""""""""""""

This example solves "Oliver's 30-city Hamiltonian cycle Traveling
Salesman Problem", as described in Poon and Carter, *Computer Ops
Res*., **22**, (1995).  More importantly, it demonstrates use of a 
user-defined crossover operator, namely the "Tie-Breaking 
Crossover" (TBX1) of the same work.  

The problem has 30 cities in a plane, and the goal is to minimize
the distance traveled when visiting each city just once.
The problem has 40 equivalent optima, each with a distance of
423.741 units.  The coordinates of the cities are in the code, and
are from Oliver's original paper by way of Steve Dower's 
`site <http://www.stevedower.id.au/other/oliver30>`_.  See 
:ref:`sec_parallelexamples` for a parallel version that tries
matching the cited results.

.. literalinclude:: ../../examples/example05.py

Running it once yields the following output:

.. literalinclude:: ../../examples/output/example05_ref



.. _sec_mutationexamples:

Example 6: User-defined Mutation Operator
"""""""""""""""""""""""""""""""""""""""""

We redo :ref:`sec_maxintexample` using a custom mutation
operator, largely following the PGAPack example.

.. literalinclude:: ../../examples/example04.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example04_ref



.. _sec_endofgenexamples:

Example 7: User-defined End of Generation Operator
""""""""""""""""""""""""""""""""""""""""""""""""""

This example is almost the same as :ref:`sec_maxbitexample` but it
adds an end-of-generation operator.  Here, we're using it to 
flip a random bit to 1.  Of course, since we're maximizing the 
bit sum, this is "climbing the hill" to a better answer.  This is
a trivial example of such heuristics; in other situations, there are
more complex, physically-motivated approaches.  Another use of an
end-of-generation operator would be for post-generation processing, 
such as plotting fitnesses, writing to file, etc.

.. literalinclude:: ../../examples/example07.py

Running it yields the following output:

.. literalinclude:: ../../examples/output/example07_ref

Notice that for the same settings, the heuristic improved the best
solution a little bit.  Of course, flipping just one bit out of
100 shouldn't be expected to work miracles.



.. _sec_parallelexamples:

Parallel Examples
-----------------

The following examples illustrate use of :mod:`pypgapack` in
a parallel setting using the ``mpi4py`` package.



Example 8: Parallel MAXBIT
""""""""""""""""""""""""""

We adapt our favorite :ref:`sec_maxbitexample` using MPI.  We up
the string length and population to bring out timing differences.

.. literalinclude:: ../../examples/example08.py

Running it using ``mpirun -np 1 python example08.py`` yields 
the following output:

.. literalinclude:: ../../examples/output/example08_ref_np_1

Running it using ``mpirun -np 2 python example08.py`` yields 
the following output:

.. literalinclude:: ../../examples/output/example08_ref_np_2

This was using a dual core laptop with probably more browser 
windows open than needed, but the results look good.  Overall,
pgapack focuses parallelism on the object function evaluation. 
Hence, if your objective function is expensive to evaluate, you
can expect relatively good scaling up to the number of strings
replaced every generation (the default is 1/10 of the total).



Example 9: Parallel Traveling Salesman
""""""""""""""""""""""""""""""""""""""

We adapt the Traveling Salesman Problem to parallel and try
matching the results Poon and Carter found for the TBX1 
cross-over operator.  The GA 
parameters used by Poon and Carter are not entirely clear.  They 
cite results for a population of 21 over 300 iterations with a
cross-over to mutation probability ratio of 0.8/0.2.  The exact
nature of the "swap" mutation is cited from another work I don't have at
hand, and the selection appears to be standard elitist, i.e all but 
the best is replaced.  I set a swap operator that always swaps
a user-set number of pairs.  Using 3 pairs (i.e. 10%) seems
to get close to the cited results.  Also, Because PGAPack doesn't like
odd population sizes, I use 22.  

.. literalinclude:: ../../examples/example09.py

Running it using ``mpirun -np 1 python example09.py`` yields 
the following output (last several lines):

.. literalinclude:: ../../examples/output/example09_ref_np_1
   :lines: 1741-
Running it using ``mpirun -np 6 python example09.py`` yields 
the following output (last several lines):

.. literalinclude:: ../../examples/output/example09_ref_np_6
   :lines: 1990-

This was using a 6-core machine, but the speedup was barely
25%---why?  Well, evaluating the distance of 30 cities is 
cheap compared to the sorting occuring on the master process
for the cross-over operation.  As noted earlier, PGAPack's 
parallelism is definitely meant for expensive evaluations
relative to everything else.  Here, we see *some* speedup, 
just not very good.  Also compare the mean and standard 
deviation to the cited 829.00 and 72.69.  The parallel results
are much closer, which may be a fluke, but it may be worth 
investigating.
