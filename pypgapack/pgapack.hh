//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   pgapack.hh
 * \author Jeremy Roberts
 * \date   10/28/2011
 * \brief  PGAPack library C++ wrapper class definition for use with SWIG.
 * \note   Copyright (C) 2011 Jeremy Roberts
 */
//---------------------------------------------------------------------------//
// $Rev:: 152                                           $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date:: 2011-09-27 06:34:56 -0400 (Tue, 27 Sep 2011) $:Date of last commit
//---------------------------------------------------------------------------//

#ifndef PGAPACK_HH
#define PGAPACK_HH



#include "pgapack.h"
#include <stdexcept>
#include <iostream>

/*!
 *  \class PGA
 *  \brief A C++ wrapper for the PGAPack genetic algorithm library.
 *
 *  This class wraps most of PGAPack's function calls, particularly
 *  those that seem most useful to the general user.  Adding the rest should
 *  be straightforward if desired.  The purpose of this class is to make
 *  SWIG interfaces easier to generate (as opposed to wrapping pgapack.h 
 *  directly).
 */
class PGA
{

private:

  // The "context".  This is hidden from the client.
  PGAContext *ctx;

public:

  //------------------------------------------------------------------------//
  // Useful Typedefs
  //------------------------------------------------------------------------//

  // function pointer with context argument
  typedef double (*func)(PGAContext* c, int pop, int p);
  typedef PGAContext Context;
  // types
  typedef PGABinary       Binary;
  typedef PGAInteger      Integer;
  typedef PGAReal         Real;
  typedef PGACharacter    Character;

  //------------------------------------------------------------------------//
  // Wrappers for Basic PGAPack Methods
  //------------------------------------------------------------------------//

  // Constructor (replaces direct PGACreate call).  
  PGA(int argc, char *argv[], int datatype, int stringlength, int maxormin)
    : ctx(PGACreate(&argc, argv, datatype, stringlength, maxormin))
  {
    std::cout << "***Constructing PGA***" << std::endl;
  }

  // Destructor (replaces direct PGADestroy call)
  ~PGA()
  {
  }

  // Setup the context (uses defaults unless settings are already altered)
  void SetUp()
  {
    PGASetUp(ctx);
  }

  // Solve the problem.
  void Run(func f)
  {
    PGARun(ctx, f);
  }

  // Delete the context.  Must be called explicitly.
  void Destroy()
  {
    if (ctx != NULL)
    {
      std::cout << "***Destroying PGA context***" << std::endl;
      PGADestroy(ctx);
    }
    else
    {
      std::cout << "***Warning: PGA context was already NULL***" << std::endl;
    }
  }

// stub function to be replace in swig interface
//  template <class T> 
//  void GetChromosome(int p, int pop, int* n, T** chromo){}

  //------------------------------------------------------------------------//
  // PGAPack consts
  //------------------------------------------------------------------------//

  /// \name PGAPack Datatypes
  //\{
  static const int DATATYPE_BINARY      = PGA_DATATYPE_BINARY;
  static const int DATATYPE_INTEGER     = PGA_DATATYPE_INTEGER;
  static const int DATATYPE_REAL        = PGA_DATATYPE_REAL;
  static const int DATATYPE_CHARACTER   = PGA_DATATYPE_CHARACTER;
  static const int DATATYPE_USER        = PGA_DATATYPE_USER;
  //\}
  /// \name Directions
  //\{
  static const int MAXIMIZE = PGA_MAXIMIZE; 
  static const int MINIMIZE = PGA_MINIMIZE;
  //\}
  /// \name Termination
  //\{
  static const int STOP_MAXITER = PGA_STOP_MAXITER; 
  static const int STOP_NOCHANGE = PGA_STOP_NOCHANGE; 
  static const int TOP_TOOSIMILAR = PGA_STOP_TOOSIMILAR; 
  //\}
  /// \name Crossover
  //\{
  static const int CROSSOVER_ONEPT = PGA_CROSSOVER_ONEPT; 
  static const int CROSSOVER_TWOPT = PGA_CROSSOVER_TWOPT; 
  static const int CROSSOVER_UNIFORM = PGA_CROSSOVER_UNIFORM;
  //\}
  /// \name Selection
  //\{
  static const int SELECT_PROPORTIONAL = PGA_SELECT_PROPORTIONAL; 
  static const int SELECT_SUS = PGA_SELECT_SUS; 
  static const int SELECT_TOURNAMENT = PGA_SELECT_TOURNAMENT; 
  static const int SELECT_PTOURNAMENT = PGA_SELECT_PTOURNAMENT;
  //\}
  /// \name Fitness
  //\{
  static const int FITNESS_RAW = PGA_FITNESS_RAW;
  static const int FITNESS_NORMAL = PGA_FITNESS_NORMAL;
  static const int FITNESS_RANKING = PGA_FITNESS_RANKING; 
  //\}
  /// \name Mutation
  //\{
  static const int MUTATION_CONSTANT = PGA_MUTATION_CONSTANT; 
  static const int MUTATION_RANGE    = PGA_MUTATION_RANGE;
  static const int MUTATION_UNIFORM  = PGA_MUTATION_UNIFORM;  
  static const int MUTATION_GAUSSIAN = PGA_MUTATION_GAUSSIAN; 
  static const int MUTATION_PERMUTE  = PGA_MUTATION_PERMUTE;
  //\}
  /// \name Population Replacement
  //\{
  static const int POPREPL_BEST         = PGA_POPREPL_BEST;
  static const int POPREPL_RANDOM_NOREP = PGA_POPREPL_RANDOM_NOREP;
  static const int POPREPL_RANDOM_REP   = PGA_POPREPL_RANDOM_REP;
  //\}
  /// \name User Functions
  //\{
  static const int USERFUNCTION_INITSTRING = PGA_USERFUNCTION_INITSTRING;
  static const int USERFUNCTION_CROSSOVER = PGA_USERFUNCTION_CROSSOVER;
  static const int USERFUNCTION_MUTATION = PGA_USERFUNCTION_MUTATION;
  static const int USERFUNCTION_ENDOFGEN = PGA_USERFUNCTION_ENDOFGEN;
  //\}
  /// \name User Functions
  //\{
  static const int OLDPOP = PGA_OLDPOP;
  static const int NEWPOP = PGA_NEWPOP;
  //\}
  /// \name PGA Bools
  //\{
  static const int FALSE = PGA_FALSE;
  static const int TRUE  = PGA_TRUE;
  //\}

  //------------------------------------------------------------------------//
  // More specialized PGAPack Methods
  //------------------------------------------------------------------------//

  //
  // binary.c
  //
  void SetBinaryAllele(int p, int pop, int i, int val)
  {
    PGASetBinaryAllele(ctx, p, pop, i, val);
  }
  int GetBinaryAllele(int p, int pop, int i)
  {
    return PGAGetBinaryAllele(ctx, p, pop, i);
  }
  void SetBinaryInitProb(double probability)
  {
    PGASetBinaryInitProb(ctx, probability);
  }
  double GetBinaryInitProb()
  {
    return PGAGetBinaryInitProb(ctx);
  }
  void BinaryCreateString(int p, int pop, int initflag)
  {
    PGABinaryCreateString(ctx, p, pop, initflag);
  }
  int BinaryMutation(int p, int pop, double mr)
  {
    return PGABinaryMutation(ctx, p, pop, mr);
  }
  void BinaryOneptCrossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGABinaryOneptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void BinaryTwoptCrossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGABinaryTwoptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void BinaryUniformCrossover(int p1,
                              int p2,
                              int pop1,
                              int c1,
                              int c2,
                              int pop2)
  {
    PGABinaryUniformCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void BinaryPrintString(FILE* fp, int p, int pop)
  {
    PGABinaryPrintString(ctx, fp, p, pop);
  }
  void BinaryCopyString(int p1, int pop1, int p2, int pop2)
  {
    PGABinaryCopyString(ctx, p1, pop1, p2, pop2);
  }
  int BinaryDuplicate(int p1, int pop1, int p2, int pop2)
  {
    return PGABinaryDuplicate(ctx, p1, pop1, p2, pop2);
  }
  void BinaryInitString(int p, int pop)
  {
    PGABinaryInitString(ctx, p, pop);
  }
  MPI_Datatype BinaryBuildDatatype(int p, int pop)
  {
    return PGABinaryBuildDatatype(ctx, p, pop);
  }
  int BinaryHammingDistance(PGABinary* s1, PGABinary* s2)
  {
    return PGABinaryHammingDistance(ctx, s1, s2);
  }
  void BinaryPrint(FILE* fp, PGABinary* chrom, int nb)
  {
    PGABinaryPrint(ctx, fp, chrom, nb);
  }
  //
  // char.c
  //
  void SetCharacterAllele(int p, int pop, int i, char value)
  {
    PGASetCharacterAllele(ctx, p, pop, i, value);
  }
  char GetCharacterAllele(int p, int pop, int i)
  {
    return PGAGetCharacterAllele(ctx, p, pop, i);
  }
  void SetCharacterInitType(int value)
  {
    PGASetCharacterInitType(ctx, value);
  }
  void CharacterCreateString(int p, int pop, int InitFlag)
  {
    PGACharacterCreateString(ctx, p, pop, InitFlag);
  }
  int CharacterMutation(int p, int pop, double mr)
  {
    return PGACharacterMutation(ctx, p, pop, mr);
  }
  void CharacterOneptCrossover(int p1,
                               int p2,
                               int pop1,
                               int c1,
                               int c2,
                               int pop2)
  {
    PGACharacterOneptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void CharacterTwoptCrossover(int p1,
                               int p2,
                               int pop1,
                               int c1,
                               int c2,
                               int pop2)
  {
    PGACharacterTwoptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void CharacterUniformCrossover(int p1,
                                 int p2,
                                 int pop1,
                                 int c1,
                                 int c2,
                                 int pop2)
  {
    PGACharacterUniformCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void CharacterPrintString(FILE* fp, int p, int pop)
  {
    PGACharacterPrintString(ctx, fp, p, pop);
  }
  void CharacterCopyString(int p1, int pop1, int p2, int pop2)
  {
    PGACharacterCopyString(ctx, p1, pop1, p2, pop2);
  }
  int CharacterDuplicate(int p1, int pop1, int p2, int pop2)
  {
    return PGACharacterDuplicate(ctx, p1, pop1, p2, pop2);
  }
  void CharacterInitString(int p, int pop)
  {
    PGACharacterInitString(ctx, p, pop);
  }
  MPI_Datatype CharacterBuildDatatype(int p, int pop)
  {
    return PGACharacterBuildDatatype(ctx, p, pop);
  }
  //
  // create.c
  //
  void SetRandomInitFlag(int RandomBoolean)
  {
    PGASetRandomInitFlag(ctx, RandomBoolean);
  }
  int GetRandomInitFlag()
  {
    return PGAGetRandomInitFlag(ctx);
  }
  void CreatePop(int pop)
  {
    PGACreatePop(ctx, pop);
  }
  void CreateIndividual(int p, int pop, int initflag)
  {
    PGACreateIndividual(ctx, p, pop, initflag);
  }
  //
  // cross.c
  //
  void Crossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGACrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  int GetCrossoverType()
  {
    return PGAGetCrossoverType(ctx);
  }
  double GetCrossoverProb()
  {
    return PGAGetCrossoverProb(ctx);
  }
  double GetUniformCrossoverProb()
  {
    return PGAGetUniformCrossoverProb(ctx);
  }
  void SetCrossoverType(int crossover_type)
  {
    PGASetCrossoverType(ctx, crossover_type);
  }
  void SetCrossoverProb(double crossover_prob)
  {
    PGASetCrossoverProb(ctx, crossover_prob);
  }
  void SetUniformCrossoverProb(double uniform_cross_prob)
  {
    PGASetUniformCrossoverProb(ctx, uniform_cross_prob);
  }
  //
  // debug.c
  //
//  void SortFuncNameIndex()
//  {
//    PGASortFuncNameIndex(ctx);
//  }
//  void DebugPrint(int level,
//                  char* funcname,
//                  char* msg,
//                  int datatype,
//                  void* data)
//  {
//    PGADebugPrint(ctx, level, funcname, msg, datatype, data);
//  }
//  void SetDebugLevel(int level)
//  {
//    PGASetDebugLevel(ctx, level);
//  }
//  void ClearDebugLevel(int level)
//  {
//    PGAClearDebugLevel(ctx, level);
//  }
//  void SetDebugLevelByName(char* funcname)
//  {
//    PGASetDebugLevelByName(ctx, funcname);
//  }
//  void ClearDebugLevelByName(char* funcname)
//  {
//    PGAClearDebugLevelByName(ctx, funcname);
//  }
//  int GetDebugLevelOfName(char* funcname)
//  {
//    return PGAGetDebugLevelOfName(ctx, funcname);
//  }
//  int GetDebugFlag(char* funcname)
//  {
//    return PGAGetDebugFlag(ctx, funcname);
//  }
//  void SetDebugFlag11(int Flag)
//  {
//    PGASetDebugFlag11(ctx, Flag);
//  }
//  void SetDebugFlag20(int Flag)
//  {
//    PGASetDebugFlag20(ctx, Flag);
//  }
//  void SetDebugFlag21(int Flag)
//  {
//    PGASetDebugFlag21(ctx, Flag);
//  }
//  void SetDebugFlag30(int Flag)
//  {
//    PGASetDebugFlag30(ctx, Flag);
//  }
//  void SetDebugFlag32(int Flag)
//  {
//    PGASetDebugFlag32(ctx, Flag);
//  }
//  void SetDebugFlag34(int Flag)
//  {
//    PGASetDebugFlag34(ctx, Flag);
//  }
//  void SetDebugFlag36(int Flag)
//  {
//    PGASetDebugFlag36(ctx, Flag);
//  }
//  void SetDebugFlag40(int Flag)
//  {
//    PGASetDebugFlag40(ctx, Flag);
//  }
//  void SetDebugFlag42(int Flag)
//  {
//    PGASetDebugFlag42(ctx, Flag);
//  }
//  void SetDebugFlag44(int Flag)
//  {
//    PGASetDebugFlag44(ctx, Flag);
//  }
//  void SetDebugFlag46(int Flag)
//  {
//    PGASetDebugFlag46(ctx, Flag);
//  }
//  void SetDebugFlag48(int Flag)
//  {
//    PGASetDebugFlag48(ctx, Flag);
//  }
//  void SetDebugFlag50(int Flag)
//  {
//    PGASetDebugFlag50(ctx, Flag);
//  }
//  void SetDebugFlag52(int Flag)
//  {
//    PGASetDebugFlag52(ctx, Flag);
//  }
//  void SetDebugFlag54(int Flag)
//  {
//    PGASetDebugFlag54(ctx, Flag);
//  }
//  void SetDebugFlag56(int Flag)
//  {
//    PGASetDebugFlag56(ctx, Flag);
//  }
//  void SetDebugFlag58(int Flag)
//  {
//    PGASetDebugFlag58(ctx, Flag);
//  }
//  void SetDebugFlag60(int Flag)
//  {
//    PGASetDebugFlag60(ctx, Flag);
//  }
//  void SetDebugFlag62(int Flag)
//  {
//    PGASetDebugFlag62(ctx, Flag);
//  }
//  void SetDebugFlag64(int Flag)
//  {
//    PGASetDebugFlag64(ctx, Flag);
//  }
//  void SetDebugFlag66(int Flag)
//  {
//    PGASetDebugFlag66(ctx, Flag);
//  }
//  void PrintDebugOptions()
//  {
//    PGAPrintDebugOptions(ctx);
//  }
  //
  // duplcate.c
  //
  int Duplicate(int p, int pop1, int pop2, int n)
  {
    return PGADuplicate(ctx, p, pop1, pop2, n);
  }
  void Change(int p, int pop)
  {
    PGAChange(ctx, p, pop);
  }
  void SetNoDuplicatesFlag(int no_dup)
  {
    PGASetNoDuplicatesFlag(ctx, no_dup);
  }
  int GetNoDuplicatesFlag()
  {
    return PGAGetNoDuplicatesFlag(ctx);
  }
  //
  // evaluate.c
  //
  void SetEvaluation(int p, int pop, double val)
  {
    PGASetEvaluation(ctx, p, pop, val);
  }
  double GetEvaluation(int p, int pop)
  {
    return PGAGetEvaluation(ctx, p, pop);
  }
  void SetEvaluationUpToDateFlag(int p, int pop, int status)
  {
    PGASetEvaluationUpToDateFlag(ctx, p, pop, status);
  }
  int GetEvaluationUpToDateFlag(int p, int pop)
  {
    return PGAGetEvaluationUpToDateFlag(ctx, p, pop);
  }
  double GetRealFromBinary(int p,
                           int pop,
                           int start,
                           int end,
                           double lower,
                           double upper)
  {
    return PGAGetRealFromBinary(ctx, p, pop, start, end, lower, upper);
  }
  double GetRealFromGrayCode(int p,
                             int pop,
                             int start,
                             int end,
                             double lower,
                             double upper)
  {
    return PGAGetRealFromGrayCode(ctx, p, pop, start, end, lower, upper);
  }
  void EncodeRealAsBinary(int p,
                          int pop,
                          int start,
                          int end,
                          double low,
                          double high,
                          double val)
  {
    PGAEncodeRealAsBinary(ctx, p, pop, start, end, low, high, val);
  }
  void EncodeRealAsGrayCode(int p,
                            int pop,
                            int start,
                            int end,
                            double low,
                            double high,
                            double val)
  {
    PGAEncodeRealAsGrayCode(ctx, p, pop, start, end, low, high, val);
  }
  int GetIntegerFromBinary(int p, int pop, int start, int end)
  {
    return PGAGetIntegerFromBinary(ctx, p, pop, start, end);
  }
  int GetIntegerFromGrayCode(int p, int pop, int start, int end)
  {
    return PGAGetIntegerFromGrayCode(ctx, p, pop, start, end);
  }
  void EncodeIntegerAsBinary(int p, int pop, int start, int end, int val)
  {
    PGAEncodeIntegerAsBinary(ctx, p, pop, start, end, val);
  }
  void EncodeIntegerAsGrayCode(int p, int pop, int start, int end, int val)
  {
    PGAEncodeIntegerAsGrayCode(ctx, p, pop, start, end, val);
  }
  double MapIntegerToReal(int v, int a, int b, double l, double u)
  {
    return PGAMapIntegerToReal(ctx, v, a, b, l, u);
  }
  int MapRealToInteger(double r, double l, double u, int a, int b)
  {
    return PGAMapRealToInteger(ctx, r, l, u, a, b);
  }
  //
  // fitness.c
  //
  void Fitness(int popindex)
  {
    PGAFitness(ctx, popindex);
  }
  int Rank(int p, int* order, int n)
  {
    return PGARank(ctx, p, order, n);
  }
  double GetFitness(int p, int pop)
  {
    return PGAGetFitness(ctx, p, pop);
  }
  int GetFitnessType()
  {
    return PGAGetFitnessType(ctx);
  }
  int GetFitnessMinType()
  {
    return PGAGetFitnessMinType(ctx);
  }
  double GetMaxFitnessRank()
  {
    return PGAGetMaxFitnessRank(ctx);
  }
  void SetFitnessType(int fitness_type)
  {
    PGASetFitnessType(ctx, fitness_type);
  }
  void SetFitnessMinType(int fitness_type)
  {
    PGASetFitnessMinType(ctx, fitness_type);
  }
  void SetMaxFitnessRank(double fitness_rank_max)
  {
    PGASetMaxFitnessRank(ctx, fitness_rank_max);
  }
  void FitnessLinearNormal(PGAIndividual* pop)
  {
    PGAFitnessLinearNormal(ctx, pop);
  }
  void FitnessLinearRank(PGAIndividual* pop)
  {
    PGAFitnessLinearRank(ctx, pop);
  }
  void FitnessMinReciprocal(PGAIndividual* pop)
  {
    PGAFitnessMinReciprocal(ctx, pop);
  }
  void FitnessMinCmax(PGAIndividual* pop)
  {
    PGAFitnessMinCmax(ctx, pop);
  }
  void SetFitnessCmaxValue(double val)
  {
    PGASetFitnessCmaxValue(ctx, val);
  }
  double GetFitnessCmaxValue()
  {
    return PGAGetFitnessCmaxValue(ctx);
  }
  //
  // hamming.c
  //
  double HammingDistance(int popindex)
  {
    return PGAHammingDistance(ctx, popindex);
  }
  //
  // heap.c
  //
  void DblHeapSort(double* a, int* idx, int n)
  {
    PGADblHeapSort(ctx, a, idx, n);
  }
  void IntHeapSort(int* a, int* idx, int n)
  {
    PGAIntHeapSort(ctx, a, idx, n);
  }
  //
  // integer.c
  //
  void SetIntegerAllele(int p, int pop, int i, int value)
  {
    PGASetIntegerAllele(ctx, p, pop, i, value);
  }
  int GetIntegerAllele(int p, int pop, int i)
  {
    return PGAGetIntegerAllele(ctx, p, pop, i);
  }
  void SetIntegerInitPermute(int min, int max)
  {
    PGASetIntegerInitPermute(ctx, min, max);
  }
  void SetIntegerInitRange(int* min, int* max)
  {
    PGASetIntegerInitRange(ctx, min, max);
  }
  int GetIntegerInitType()
  {
    return PGAGetIntegerInitType(ctx);
  }
  int GetMinIntegerInitValue(int i)
  {
    return PGAGetMinIntegerInitValue(ctx, i);
  }
  int GetMaxIntegerInitValue(int i)
  {
    return PGAGetMaxIntegerInitValue(ctx, i);
  }
  void IntegerCreateString(int p, int pop, int InitFlag)
  {
    PGAIntegerCreateString(ctx, p, pop, InitFlag);
  }
  int IntegerMutation(int p, int pop, double mr)
  {
    return PGAIntegerMutation(ctx, p, pop, mr);
  }
  void IntegerOneptCrossover(int p1,
                             int p2,
                             int pop1,
                             int c1,
                             int c2,
                             int pop2)
  {
    PGAIntegerOneptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void IntegerTwoptCrossover(int p1,
                             int p2,
                             int pop1,
                             int c1,
                             int c2,
                             int pop2)
  {
    PGAIntegerTwoptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void IntegerUniformCrossover(int p1,
                               int p2,
                               int pop1,
                               int c1,
                               int c2,
                               int pop2)
  {
    PGAIntegerUniformCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void IntegerPrintString(FILE* fp, int p, int pop)
  {
    PGAIntegerPrintString(ctx, fp, p, pop);
  }
  void IntegerCopyString(int p1, int pop1, int p2, int pop2)
  {
    PGAIntegerCopyString(ctx, p1, pop1, p2, pop2);
  }
  int IntegerDuplicate(int p1, int pop1, int p2, int pop2)
  {
    return PGAIntegerDuplicate(ctx, p1, pop1, p2, pop2);
  }
  void IntegerInitString(int p, int pop)
  {
    PGAIntegerInitString(ctx, p, pop);
  }
  MPI_Datatype IntegerBuildDatatype(int p, int pop)
  {
    return PGAIntegerBuildDatatype(ctx, p, pop);
  }
  //
  // mutation.c
  //
  int Mutate(int p, int pop)
  {
    return PGAMutate(ctx, p, pop);
  }
  void SetMutationType(int mutation_type)
  {
    PGASetMutationType(ctx, mutation_type);
  }
  int GetMutationType()
  {
    return PGAGetMutationType(ctx);
  }
  void SetMutationRealValue(double val)
  {
    PGASetMutationRealValue(ctx, val);
  }
  double GetMutationRealValue()
  {
    return PGAGetMutationRealValue(ctx);
  }
  void SetMutationIntegerValue(int val)
  {
    PGASetMutationIntegerValue(ctx, val);
  }
  int GetMutationIntegerValue()
  {
    return PGAGetMutationIntegerValue(ctx);
  }
  void SetMutationBoundedFlag(int val)
  {
    PGASetMutationBoundedFlag(ctx, val);
  }
  int GetMutationBoundedFlag()
  {
    return PGAGetMutationBoundedFlag(ctx);
  }
  void SetMutationProb(double mutation_prob)
  {
    PGASetMutationProb(ctx, mutation_prob);
  }
  double GetMutationProb()
  {
    return PGAGetMutationProb(ctx);
  }
  //
  // pga.c
  //
  void RunMutationAndCrossover(int oldpop, int newpop)
  {
    PGARunMutationAndCrossover(ctx, oldpop, newpop);
  }
  void RunMutationOrCrossover(int oldpop, int newpop)
  {
    PGARunMutationOrCrossover(ctx, oldpop, newpop);
  }
  void UpdateGeneration(MPI_Comm comm)
  {
    PGAUpdateGeneration(ctx, comm);
  }
  int GetDataType()
  {
    return PGAGetDataType(ctx);
  }
  int GetOptDirFlag()
  {
    return PGAGetOptDirFlag(ctx);
  }
  int GetStringLength()
  {
    return PGAGetStringLength(ctx);
  }
  int GetVariableStringLength(int p, int pop)
  {
    return PGAGetVariableStringLength(ctx, p, pop);
  }
  int GetGAIterValue()
  {
    return PGAGetGAIterValue(ctx);
  }
  void SetMutationOrCrossoverFlag(int flag)
  {
    PGASetMutationOrCrossoverFlag(ctx, flag);
  }
  void SetMutationAndCrossoverFlag(int flag)
  {
    PGASetMutationAndCrossoverFlag(ctx, flag);
  }
  int GetMutationOrCrossoverFlag()
  {
    return PGAGetMutationOrCrossoverFlag(ctx);
  }
  int GetMutationAndCrossoverFlag()
  {
    return PGAGetMutationAndCrossoverFlag(ctx);
  }
  //
  // pop.c
  //
  void SortPop(int pop)
  {
    PGASortPop(ctx, pop);
  }
  int GetPopSize()
  {
    return PGAGetPopSize(ctx);
  }
  int GetNumReplaceValue()
  {
    return PGAGetNumReplaceValue(ctx);
  }
  int GetPopReplaceType()
  {
    return PGAGetPopReplaceType(ctx);
  }
  int GetSortedPopIndex(int n)
  {
    return PGAGetSortedPopIndex(ctx, n);
  }
  void SetPopSize(int popsize)
  {
    PGASetPopSize(ctx, popsize);
  }
  void SetNumReplaceValue(int pop_replace)
  {
    PGASetNumReplaceValue(ctx, pop_replace);
  }
  void SetPopReplaceType(int pop_replace)
  {
    PGASetPopReplaceType(ctx, pop_replace);
  }
  //
  // random.c
  //
  int RandomFlip(double p)
  {
    return PGARandomFlip(ctx, p);
  }
  int RandomInterval(int start, int end)
  {
    return PGARandomInterval(ctx, start, end);
  }
  double Random01(int newseed)
  {
    return PGARandom01(ctx, newseed);
  }
  double RandomUniform(double start, double end)
  {
    return PGARandomUniform(ctx, start, end);
  }
  double RandomGaussian(double mean, double sigma)
  {
    return PGARandomGaussian(ctx, mean, sigma);
  }
  int GetRandomSeed()
  {
    return PGAGetRandomSeed(ctx);
  }
  void SetRandomSeed(int seed)
  {
    PGASetRandomSeed(ctx, seed);
  }
  //
  // real.c
  //
  void SetRealAllele(int p, int pop, int i, double value)
  {
    PGASetRealAllele(ctx, p, pop, i, value);
  }
  double GetRealAllele(int p, int pop, int i)
  {
    return PGAGetRealAllele(ctx, p, pop, i);
  }
  void SetRealInitPercent(double* median, double* percent)
  {
    PGASetRealInitPercent(ctx, median, percent);
  }
  void SetRealInitRange(double* min, double* max)
  {
    PGASetRealInitRange(ctx, min, max);
  }
  double GetMinRealInitValue(int i)
  {
    return PGAGetMinRealInitValue(ctx, i);
  }
  double GetMaxRealInitValue(int i)
  {
    return PGAGetMaxRealInitValue(ctx, i);
  }
  int GetRealInitType()
  {
    return PGAGetRealInitType(ctx);
  }
  void RealCreateString(int p, int pop, int initflag)
  {
    PGARealCreateString(ctx, p, pop, initflag);
  }
  int RealMutation(int p, int pop, double mr)
  {
    return PGARealMutation(ctx, p, pop, mr);
  }
  void RealOneptCrossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGARealOneptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void RealTwoptCrossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGARealTwoptCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void RealUniformCrossover(int p1, int p2, int pop1, int c1, int c2, int pop2)
  {
    PGARealUniformCrossover(ctx, p1, p2, pop1, c1, c2, pop2);
  }
  void RealPrintString(FILE* fp, int p, int pop)
  {
    PGARealPrintString(ctx, fp, p, pop);
  }
  void RealCopyString(int p1, int pop1, int p2, int pop2)
  {
    PGARealCopyString(ctx, p1, pop1, p2, pop2);
  }
  int RealDuplicate(int p1, int pop1, int p2, int pop2)
  {
    return PGARealDuplicate(ctx, p1, pop1, p2, pop2);
  }
  void RealInitString(int p, int pop)
  {
    PGARealInitString(ctx, p, pop);
  }
  MPI_Datatype RealBuildDatatype(int p, int pop)
  {
    return PGARealBuildDatatype(ctx, p, pop);
  }
  //
  // report.c
  //
  void PrintReport(FILE* fp, int pop)
  {
    PGAPrintReport(ctx, fp, pop);
  }
  void SetPrintOptions(int option)
  {
    PGASetPrintOptions(ctx, option);
  }
  void SetPrintFrequencyValue(int print_freq)
  {
    PGASetPrintFrequencyValue(ctx, print_freq);
  }
  int GetPrintFrequencyValue()
  {
    return PGAGetPrintFrequencyValue(ctx);
  }
  void PrintPopulation(FILE* fp, int pop)
  {
    PGAPrintPopulation(ctx, fp, pop);
  }
  void PrintIndividual(FILE* fp, int p, int pop)
  {
    PGAPrintIndividual(ctx, fp, p, pop);
  }
  void PrintString(FILE* file, int p, int pop)
  {
    PGAPrintString(ctx, file, p, pop);
  }
  void PrintContextVariable(FILE* fp)
  {
    PGAPrintContextVariable(ctx, fp);
  }
  //
  // restart.c
  //
  void Restart(int source_pop, int dest_pop)
  {
    PGARestart(ctx, source_pop, dest_pop);
  }
  void SetRestartFlag(int val)
  {
    PGASetRestartFlag(ctx, val);
  }
  int GetRestartFlag()
  {
    return PGAGetRestartFlag(ctx);
  }
  void SetRestartFrequencyValue(int numiter)
  {
    PGASetRestartFrequencyValue(ctx, numiter);
  }
  int GetRestartFrequencyValue()
  {
    return PGAGetRestartFrequencyValue(ctx);
  }
  void SetRestartAlleleChangeProb(double prob)
  {
    PGASetRestartAlleleChangeProb(ctx, prob);
  }
  double GetRestartAlleleChangeProb()
  {
    return PGAGetRestartAlleleChangeProb(ctx);
  }
  //
  // select.c
  //
  void Select(int popix)
  {
    PGASelect(ctx, popix);
  }
  int SelectNextIndex()
  {
    return PGASelectNextIndex(ctx);
  }
  void SetSelectType(int select_type)
  {
    PGASetSelectType(ctx, select_type);
  }
  int GetSelectType()
  {
    return PGAGetSelectType(ctx);
  }
  void SetPTournamentProb(double ptournament_prob)
  {
    PGASetPTournamentProb(ctx, ptournament_prob);
  }
  double GetPTournamentProb()
  {
    return PGAGetPTournamentProb(ctx);
  }
  int SelectProportional(PGAIndividual* pop)
  {
    return PGASelectProportional(ctx, pop);
  }
  void SelectSUS(PGAIndividual* pop)
  {
    PGASelectSUS(ctx, pop);
  }
  int SelectTournament(PGAIndividual* pop)
  {
    return PGASelectTournament(ctx, pop);
  }
  int SelectPTournament(PGAIndividual* pop)
  {
    return PGASelectPTournament(ctx, pop);
  }
  //
  // stop.c
  //
  int Done(MPI_Comm comm)
  {
    return PGADone(ctx, comm);
  }
  int CheckStoppingConditions()
  {
    return PGACheckStoppingConditions(ctx);
  }
  void SetStoppingRuleType(int stoprule)
  {
    PGASetStoppingRuleType(ctx, stoprule);
  }
  int GetStoppingRuleType()
  {
    return PGAGetStoppingRuleType(ctx);
  }
  void SetMaxGAIterValue(int maxiter)
  {
    PGASetMaxGAIterValue(ctx, maxiter);
  }
  int GetMaxGAIterValue()
  {
    return PGAGetMaxGAIterValue(ctx);
  }
  void SetMaxNoChangeValue(int max_no_change)
  {
    PGASetMaxNoChangeValue(ctx, max_no_change);
  }
  void SetMaxSimilarityValue(int max_similarity)
  {
    PGASetMaxSimilarityValue(ctx, max_similarity);
  }
  //
  // system.c
  //
  void Error(char* msg, int level, int datatype, void* data)
  {
    PGAError(ctx, msg, level, datatype, data);
  }
  int GetMaxMachineIntValue()
  {
    return PGAGetMaxMachineIntValue(ctx);
  }
  int GetMinMachineIntValue()
  {
    return PGAGetMinMachineIntValue(ctx);
  }
  double GetMaxMachineDoubleValue()
  {
    return PGAGetMaxMachineDoubleValue(ctx);
  }
  double GetMinMachineDoubleValue()
  {
    return PGAGetMinMachineDoubleValue(ctx);
  }
  void Usage()
  {
    PGAUsage(ctx);
  }
  void PrintVersionNumber()
  {
    PGAPrintVersionNumber(ctx);
  }
  //
  // user.c
  //
  void SetUserFunction(int constant, void* f)
  {
    PGASetUserFunction(ctx, constant, f);
  }

  //
  // utility.c
  //
  double Mean(double *a, int n)
  {
    return PGAMean (ctx, a, n);
  }
  double Stddev ( double *a, int n, double mean)
  {
    return PGAStddev (ctx,  a, n,  mean);
  }
  int Round(double x)
  {
    return PGARound(ctx, x);
  }
  void CopyIndividual( int p1, int pop1, int p2, int pop2)
  {
    PGACopyIndividual(ctx, p1, pop1, p2, pop2);
  }
  int CheckSum(int p, int pop)
  {
    return PGACheckSum(ctx, p, pop);
  }
  int GetWorstIndex(int pop)
  {
    return PGAGetWorstIndex(ctx, pop);
  }
  int GetBestIndex(int pop)
  {
    return PGAGetBestIndex(ctx, pop);
  }
  PGAIndividual *GetIndividual ( int p, int pop)
  {
    return PGAGetIndividual (ctx, p, pop);
  }
  void UpdateAverage(int pop)
  {
    PGAUpdateAverage(ctx, pop);
  }
  void UpdateOnline(int pop)
  {
    PGAUpdateOnline(ctx, pop);
  }
  void UpdateOffline(int pop)
  {
    PGAUpdateOffline(ctx, pop);
  }
  int ComputeSimilarity(PGAIndividual *pop)
  {
    return PGAComputeSimilarity(ctx, pop);
  }
  //
  // parallel.c
  //
  void Evaluate(int pop, func f, MPI_Comm comm)
  {
    PGAEvaluate(ctx, pop, f, comm);
  }
  PGAContext* GetContext()
  {
    return ctx;
  }
};

#endif // PGAPACK_HH
