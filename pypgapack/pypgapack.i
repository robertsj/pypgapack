//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   pypgapack.i
 * \author Jeremy Roberts
 * \date   12/02/2011
 * \brief  SWIG interface file for PGAPack
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

%module pypgapack
%{
#include "pgapack.hh"
%}

//////////////////////////////////////////////////////////////////////////////
// Structures from pgapack.h

%ignore PGAIndividual;

//%{
//typedef struct {                    /* primary population data structure   */
//  double evalfunc;                  /* evaluation function value           */
//  double fitness;                   /* fitness    function value           */
//  int    evaluptodate;              /* flag whether evalfunc is current    */
//  void   *chrom;                    /* pointer to the GA string            */
//} PGAIndividual;
//%}

//////////////////////////////////////////////////////////////////////////////
// Include numpy support for arrays.

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
  import_array();
%}

// Support for floats and integers. 

%numpy_typemaps(double, NPY_DOUBLE, int)
//%numpy_typemaps(int,    NPY_INT,    int)
%apply (int DIM1, double* IN_ARRAY1) 
      {(int len1, double* x1),
       (int len2, double* x2)}
%apply (int DIM1,  int* IN_ARRAY1)
      {(int len1,  int* x1),
       (int len2,  int* x2)}

//////////////////////////////////////////////////////////////////////////////
// Python interfaces for initializing variables.


%ignore SetRealInitRange;
%ignore SetIntegerInitRange;
%rename(SetRealInitRange)    PySetRealInitRange;
%rename(SetIntegerInitRange) PySetIntegerInitRange ;
/*%ignore GetIntegerAllele;*/
/*%rename(GetIntegerAllele)    PyGetIntegerAllele ;*/

%extend PGA {
  // Real range
  %exception PySetRealInitRange 
  {
      $action
      if (PyErr_Occurred()) SWIG_fail;
  }
  void PySetRealInitRange(int len1, double* x1, int len2, double* x2) 
  {
    if (len1 != len2) 
    {
      PyErr_Format(PyExc_ValueError,
                   "Arrays of lengths (%d,%d) given", len1, len2);
      return;
    }
    $self->SetRealInitRange(x1, x2);
    return;
  }
  // Integer range
  %exception PySetIntegerInitRange 
  {
      $action
      if (PyErr_Occurred()) SWIG_fail;
  }
  void PySetIntegerInitRange(int len1, int* x1, int len2, int* x2) 
  {
    if (len1 != len2) 
    {
      PyErr_Format(PyExc_ValueError,
                   "Arrays of lengths (%d,%d) given", len1, len2);
      return;
    }
    $self->SetIntegerInitRange(x1, x2);
    return;
  }
}; 

//////////////////////////////////////////////////////////////////////////////
// Getting a chromosome (aka string), which represents a possible solution.
//
// PGA has functions to access all the other members of
// an individual, but the void* nature of individual->chrom is a little bit
// tricky.  We use numpy.i to get an "argoutview" of the data; that is, we
// get direct, "could mess it up" access to the memory.  However, this is 
// pretty much the only way to do it without using loops around the 
// GetXXXAllele functions.  On the Python side, the data shows up as a 
// Numpy array.  *It is up to the client* not to do anything stupid like
// deleting the object in which the data resides while trying to view that
// data.  This should largely be avoided if the client sticks to inheriting
// the PGA class.

// Map bools to a numpy unsigned int (c bool (4 bit) != numpy bool (1 bit))
//%numpy_typemaps(bool, NPY_UINT, int)

// Apply the applicable numpy.i typemaps to the native allele types
//%apply (int* DIM1, bool**    ARGOUTVIEW_ARRAY1) 
//{
//  (int* n, bool**    chromo)
//}
//%apply (int* DIM1, unsigned char**    ARGOUTVIEW_ARRAY1) 
//{
//  (int* n, unsigned char**    chromo)
//}
%apply (int* DIM1, long signed int**    ARGOUTVIEW_ARRAY1) 
{
  (int* n, long signed int**    chromo)
}
%apply (int* DIM1, double** ARGOUTVIEW_ARRAY1) 
{
  (int* n, double** chromo)
}


// Define the chromosome accessors
//%ignore GetChromosome; // Ignore the empty stub definition
%extend PGA {
    void GetIntegerChromosome(int p, int pop, int* n, long signed int** chromo) 
    {
      *n = $self->GetStringLength();
      *chromo = (long signed int*)$self->GetIndividual(p, pop)->chrom;
      return;
    }
    void GetRealChromosome   (int p, int pop, int* n, double** chromo) 
    {
      *n = $self->GetStringLength();
      *chromo = (double *)$self->GetIndividual(p, pop)->chrom;
      return;
    }
    
//  \todo Try getting this templated version of GetChromosom to work.
//  template <class T> void GetChromosome(int p, int pop, int* n, T** chromo) 
//  {
//    *n = $self->GetStringLength();
//    *chromo = (T *)$self->GetIndividual(p, pop)->chrom;
//    return;
//  }
// Overload for the native types.
//  %template(GetChromosome) GetChromosome<bool>;    
//  %template(GetChromosome) GetChromosome<char>;    
//  %template(GetChromosome) GetChromosome<int>;    
//  %template(GetChromosome) GetChromosome<double>;

}; 


//////////////////////////////////////////////////////////////////////////////
// Python wrapper for system arguments 

%typemap(in) (int argc, char *argv[]) {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
	$2[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError,"list must contain strings");
	free($2);
	return NULL;
      }
    }
    $2[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(freearg) (int argc, char **argv) {
  free((char *) $2);
}

////////////////////////////////////////////////////////////////////////////////
//// Python wrappers for function callbacks
%{

// Pointers to Python functions.
static PyObject *objective_   = NULL;
static PyObject *initstring_  = NULL;
static PyObject *crossover_   = NULL;
static PyObject *mutation_    = NULL;
static PyObject *endofgen_    = NULL;

// "glue" function for objective
static double objective(PGAContext* c, int p, int pop)
{
  // Ensure callback has been defined.
  if (objective_ == NULL) 
    throw std::runtime_error("Error: objective not set!");
  // Build the argument list.
  PyObject *arglist = Py_BuildValue("ii", p, pop);
  // and get the result.
  PyObject *result  = PyEval_CallObject(objective_, arglist);
  // Decrement the reference counts.
  Py_DECREF(arglist);
  // Deal with result.
  double val = HUGE_VAL;
  if (PyErr_Occurred()) 
  {
    Py_XDECREF(result);
    PyErr_Print();
    throw std::runtime_error("Forced stop in objective.");
  }
  else if (result && PyFloat_Check(result)) 
  {
    val = PyFloat_AsDouble(result);
    Py_DECREF(result);
  }
  else 
  {
    Py_XDECREF(result);
    throw std::invalid_argument("Invalid result passed to PGA in objective.");
  }
  return val;
}

// "glue" function for initstring
static void initstring(PGAContext* c, int p, int pop)
{
  // Ensure callback has been defined.
  if (initstring_ == NULL) 
    throw std::runtime_error("Error: initstring not set!");
  // Build the argument list.
  PyObject *arglist = Py_BuildValue("ii", p, pop);
  // and get the result.
  PyObject *result  = PyEval_CallObject(initstring_, arglist);
  // Decrement the reference counts.
  Py_DECREF(arglist);
  // Deal with result.
  if (PyErr_Occurred()) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    throw std::runtime_error("Forced stop in initstring.");
  }
  else if (result == Py_None) 
  {
    Py_DECREF(result);
  }
  else 
  {
    Py_XDECREF(result);
    throw std::invalid_argument("Invalid result passed to PGA in initstring.");
  }
  return;
}

// "glue" function for crossover
static void crossover(PGAContext* c, int p1, int p2, int pop1, int pop2, int c1, int c2)
{
  
  // Ensure callback has been defined.
  if (crossover_ == NULL) 
    throw std::runtime_error("Error: crossover not set!");
  // Build the argument list
  PyObject *arglist = Py_BuildValue("iiiiii", p1, p2, pop1, pop2, c1, c2);
  // and get the result.
  PyObject *result  = PyEval_CallObject(crossover_, arglist);
  // Decrement the reference counts.
  Py_DECREF(arglist);
  // Deal with result.
  if (PyErr_Occurred()) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    throw std::runtime_error("Forced stop in crossover.");
  }
  else if (result == Py_None) // must return void
  {
    Py_DECREF(result);
  }
  else 
  {
    Py_XDECREF(result);
    throw std::invalid_argument("Invalid result passed to PGA in crossover.");
  }
  return;
}

// "glue" function for mutation
static int mutation(PGAContext* c, int p, int pop, double pm)
{
  // Ensure callback has been defined.
  if (mutation_ == NULL) 
    throw std::runtime_error("Error: mutation not set!");
  // Build the argument list
  PyObject *arglist = Py_BuildValue("iif", p, pop, pm);
  // and get the result.
  PyObject *result  = PyEval_CallObject(mutation_, arglist);
  // Decrement the reference counts.
  Py_DECREF(arglist);
  // Deal with result.
  int val = 0;
  if (PyErr_Occurred()) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    throw std::runtime_error("Forced stop in mutation.");
  }
  else if (result && PyInt_Check(result)) // must return int
  {
    val = PyInt_AsLong(result);
    Py_DECREF(result);
  }
  else 
  {
    Py_XDECREF(result);
    throw std::invalid_argument("Invalid result passed to PGA in mutation.");
  }
  return val;
}

// "glue" function for endofgen
static void endofgen(PGAContext* c)
{
  // Ensure callback has been defined.
  if (endofgen_ == NULL) 
    throw std::runtime_error("Error: endofgen not set!");
  // Build the argument list.  None for this one.
  PyObject *arglist = Py_BuildValue("()");
  // and get the result.
  PyObject *result  = PyEval_CallObject(endofgen_, arglist);
  // Decrement the reference counts.
  Py_DECREF(arglist);
  // Deal with result.
  if (PyErr_Occurred()) 
  {
    PyErr_Print();
    Py_XDECREF(result);
    throw std::runtime_error("Forced stop in endofgen.");
  }
  else if (result == Py_None) 
  {
    Py_DECREF(result);
  }
  else 
  {
    Py_XDECREF(result);
    throw std::invalid_argument("Invalid result passed to PGA in endofgen.");
  }
  return;
}

%}

// Interface for Run(...).  
%ignore Run;
%rename (Run) PyRun;
%extend PGA {
  void PyRun(PyObject *PyFunc)
  {
    if (!PyCallable_Check(PyFunc)) 
    {
      PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
      return;
    }
    Py_XDECREF(objective_); // Dispose of previous objective callback
    Py_XINCREF(PyFunc);     // Add a reference to new objective callback
    objective_ = PyFunc;    // Remember new callback
    $self->Run(objective);  // Run via interface the "glue" objective function
  }
}


// Interface for Evaluate(...).  

#ifdef PARALLEL
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%ignore Evaluate;
%rename (Evaluate) PyEvaluate;
%extend PGA {
  void PyEvaluate(MPI_Comm comm)
  {
    if (objective_ == NULL) 
      throw std::runtime_error("Error: objective not set!  Can't use Evaluate.");
    $self->Evaluate($self->NEWPOP, objective, comm);  // Run via interface the "glue" objective function
  }
}


// Interface for SetUserFunction(...,PGA.USERFUNCTION_INITSTRING,...).
%extend PGA {
  void SetInitString(PyObject *PyFunc)
  {
    if (!PyCallable_Check(PyFunc)) 
    {
      PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
      return;
    }
    Py_XDECREF(initstring_);   // Dispose of previous initstring callback
    Py_XINCREF(PyFunc);        // Add a reference to new initstring callback
    initstring_ = PyFunc;      // Remember new callback
    $self->SetUserFunction(PGA::USERFUNCTION_INITSTRING, (void *) initstring);
  }
}

// Interface for SetUserFunction(...,PGA.USERFUNCTION_CROSSOVER,...)
%extend PGA {
  void SetCrossover(PyObject *PyFunc)
  {
    if (!PyCallable_Check(PyFunc)) 
    {
      PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
      return;
    }
    Py_XDECREF(crossover_);   // Dispose of previous crossover callback
    Py_XINCREF(PyFunc);       // Add a reference to new crossover callback
    crossover_ = PyFunc;      // Remember new callback
    $self->SetUserFunction(PGA::USERFUNCTION_CROSSOVER,  (void *) crossover);
  }
}

// Interface for SetUserFunction(...,PGA.USERFUNCTION_MUTATION,...).
%extend PGA {
  void SetMutation(PyObject *PyFunc)
  {
    if (!PyCallable_Check(PyFunc)) 
    {
      PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
      return;
    }
    Py_XDECREF(mutation_);   // Dispose of previous mutation callback
    Py_XINCREF(PyFunc);      // Add a reference to new mutation callback
    mutation_ = PyFunc;      // Remember new callback
    $self->SetUserFunction(PGA_USERFUNCTION_MUTATION,  (void *) mutation);
  }
}

// Interface for SetUserFunction(...,PGA.USERFUNCTION_ENDOFGEN,...).
%extend PGA {
  void SetEndOfGen(PyObject *PyFunc)
  {
    if (!PyCallable_Check(PyFunc)) 
    {
      PyErr_SetString(PyExc_TypeError, "Parameter must be callable!");
      return;
    }
    Py_XDECREF(endofgen_);   // Dispose of previous initstring callback
    Py_XINCREF(PyFunc);      // Add a reference to new endofgen callback
    endofgen_ = PyFunc;      // Remember new callback
    $self->SetUserFunction(PGA::USERFUNCTION_ENDOFGEN, (void *) endofgen);
  }
}

%include "pgapack.hh"

