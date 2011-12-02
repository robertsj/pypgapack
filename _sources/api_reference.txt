.. _sec_reference:

API Reference
=============

Introduction
------------


This section provides a reference for all functions
defined in the :mod:`pypgapack` module that are *extensions* of the basic 
pgapack API.  All pgapack functions are contained in the 
:class:`pypgapack.PGA` class.  The pgapack library is typically 
used as follows:

.. code-block:: c

 double evaluate(PGAContext *ctx, int p, int pop);
 PGAContext *ctx; 
 ctx = PGACreate(&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
 PGASetUp(ctx);
 PGARun(ctx, evaluate);
 PGADestroy(ctx);

The ``ctx`` object is created explicitly by the user and then passed as 
the first argument to all subsequent function calls, with function 
names taking the form ``PGAxxx``.  For :mod:`pypgapack`, ``ctx`` is
a *private* member of :class:`PGA` created during construction, and all 
:class:`PGA` members drop the `PGA` prefix and the initial
``ctx`` argument. So, for example,

.. code-block:: c

  ctx = PGACreate(&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
  PGASetUp(ctx);

in C/C++ becomes

.. code-block:: python

  obj = pypgapack.PGA(sys.argv, PGA.DATATYPE_BINARY, 100, PGA.MAXIMIZE)
  obj.PGASetUp()

in Python. For all functions included in pgapack, the user is directed to the 
pgapack documentation.  What follows is a description of the few new
methods added for :mod:`pypgapack` that make life in Python a 
bit easier.


pypgapack API
-------------

The easiest way to see what :mod:`pypgapack` offers is to
do the following:

>>> import pypgapack as pga
>>> dir(pga)
['PGA', 'PGA_swigregister', '__builtins__', '__doc__', '__file__', 
 '__name__', '__package__', '_newclass', '_object', '_pypgapack', 
 '_swig_getattr', '_swig_property', '_swig_repr', '_swig_setattr', 
 '_swig_setattr_nondynamic']


This command works with any Python module.  Our interest is in the
:class:`PGA` class.  We do the same for this:

>>> dir(pga.PGA)
['BinaryBuildDatatype', 'BinaryCopyString', 'BinaryCreateString', 
 'BinaryDuplicate', 'BinaryHammingDistance', 'BinaryInitString', 
 'BinaryMutation', 'BinaryOneptCrossover', 'BinaryPrint', 
 'BinaryPrintString', 'BinaryTwoptCrossover', 'BinaryUniformCrossover',...

and find a really long list of class members, most of which are
directly from pgapack.  In the following, we document only those
not included in pgapack, as use of the pgapack functionality is
covered above (i.e. drop the ``ctx`` argument and ``PGA`` prefix). 

.. class:: PGA

   PGA wrapper class.

   .. method:: __init__(argv, datatype, n, direction)

      Construct the PGA context.  This essentially wraps the PGACreate
      function, so see the pgapack documentation.

      :param argv:      system argument
      :param datatype:  allele dataype; can be :const:`PGA.DATATYPE_XXX`,
                        where ``XXX`` is ``BINARY``, ``INTEGER``, and so on. 
      :param n:         size of the unknown, i.e. number of alleles of type
                        ``datatype``
      :param direction: either :const:`PGA.MAXIMIZE` or :const:`PGA.MINIMIZE`

   .. method::  GetIntegerChromosome(p, pop)
       
      Get direct access to the *p*-th integer chromosome string in 
      population *pop* .

      :param p:         string index
      :param pop:       population index
      :returns:         string as numpy array of integers

   .. method::  GetRealChromosome(p, pop)
       
      Get direct access to the *p*-th double chromosome string in population *pop*.

      :param p:         string index
      :param pop:       population index
      :returns:         string as numpy array of floats

   .. method:: SetInitString(f)
      
      Set a function for initializing strings.  The function ``f`` provided 
      **must** have the signature ``f(p, pop)``, but should almost certainly
      be an inerited class member with the signature ``f(self, p, pop)``.  
      See pgapack documentation for more about user functions.     

      :param f:         Python function 

      .. seealso:: 

        :ref:`sec_initstringexamples` for an example on string initialization.

   .. method:: SetCrossover(f)
      
      Set a function for the crossover operation.  The function ``f`` provided 
      **must** have the signature ``f(a,b,c,d,e,f)``, but should almost 
      certainly be an inerited class member with the signature 
      ``f(self,a,b,c,d,e,f)``. 
      See pgapack documentation for more about user functions.

      :param f:         Python function 

      .. seealso:: 

        :ref:`sec_crossoverexamples` for an example on setting the crossover 
        operator.

   .. method:: SetMutation(f)
      
      Set a function for the mutation operator.  The function ``f`` provided 
      **must** have the signature ``f(p, pop, prob)``, but should almost
      certainly be an inerited class member with the signature 
      ``f(self, p, pop, prob)``.  
      See pgapack documentation for more about user functions.

      :param f:         Python function 

      .. seealso:: 

        :ref:`sec_mutationexamples` for an example on setting the mutation 
        operator.

   .. method:: SetEndOfGen(f)
      
      Set a function for an operator to be performed at the end of each
      generation.  The function ``f`` provided 
      **must** have the signature ``f(pop)``, but should almost certainly
      be an inerited class member with the signature ``f(self, pop)``.  
      Such an operator can be used to implement hill-climbing heuristics.
      See pgapack documentation for more about user functions.

      :param f:         Python function 

      .. seealso:: 

        :ref:`sec_endofgenexamples` for an example on setting the an end of 
        generation operator.


