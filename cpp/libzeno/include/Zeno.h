// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Created: 2019-05-08
//
// ================================================================

#ifndef ZENO_H
#define ZENO_H

#include <string>
#include <vector>

#include "ParametersWalkOnSpheres.h"
#include "ParametersInteriorSampling.h"
#include "ParametersResults.h"

#include "Units.h"
#include "Uncertain.h"

#include "ResultsZeno.h"
#include "ResultsInterior.h"
#include "ResultsCompiler.h"
#include "Results.h"

#include "Geometry/Sphere.h"
#include "Geometry/Cuboid.h"
#include "Geometry/Triangle.h"
#include "Geometry/MixedModel.h"
#include "Geometry/MixedModelProcessed.h"

#include "Virials/ParametersVirial.h"
#include "Virials/ResultsVirial.h"
#include "Virials/IntegratorMSMC.h"

#include "Timer.h"

// ================================================================

#include "RandomNumber/Rand.h"
#include "RandomNumber/SPRNG.h"
#include "RandomNumber/LeapFrog.h"

#include "SpherePoint/RandomSpherePointMarsaglia.h"
#include "SpherePoint/RandomSpherePointPolar.h"

#include "SpherePoint/BiasedSpherePointRejection.h"
#include "SpherePoint/BiasedSpherePointDirect.h"

#include "SpherePoint/RandomBallPointRejection.h"

// ================================================================

namespace zeno {

class Zeno {
  using Model = MixedModelProcessed<double>;

  using BoundingSphere = Sphere<double>;

#if defined USE_RAND_RNG
  using RandomNumberGenerator = Rand;
#elif defined USE_SPRNG_RNG
  using RandomNumberGenerator = SPRNG;
#elif defined USE_LEAP_FROG_RNG
  using RandomNumberGenerator = LeapFrog;
#else
#error "No RandomNumber class selected"
#endif
  
  using RandomSpherePointGenerator = 
    RandomSpherePointMarsaglia<double, 
    RandomNumberGenerator>;
  
  using BiasedSpherePointGenerator = 
    BiasedSpherePointDirect<double, 
    RandomNumberGenerator, 
    RandomSpherePointGenerator>;
  
  using RandomBallPointGenerator =
    RandomBallPointRejection<double,
    RandomNumberGenerator>;
  
 public:
  enum class Status {Success = 0, Error, EmptyModel};

  /// Build the data structure used for the Walk-on-Spheres and 
  /// Interior Sampling algorithms from geometric data.
  ///
  /// The geometry in the MixedModel must be in the unlocked state when the
  /// function is called, and will be locked by the function.  If the geometry
  /// is in the locked state when the function is called, it will not be added.
  ///
  Zeno(MixedModel<double> * modelToProcess);
  
  ~Zeno();

  /// Perform the Walk-on-Spheres calculation with the provided parameters and
  /// store the results internally.
  ///
  /// Parameters that have not been set may have default values computed.  If
  /// so, these values will be written into the parameters object.
  ///
  /// Returns EmptyModel if the model is empty, and Success otherwise.
  ///
  Status doWalkOnSpheres
    (ParametersWalkOnSpheres * parametersWalkOnSpheres,
     ParametersResults * parametersResults);

  /// Perform the Interior Sampling calculation with the provided parameters and
  /// store the results internally.
  ///
  /// Parameters that have not been set may have default values computed.  If
  /// so, these values will be written into the parameters object.
  ///
  /// Returns EmptyModel if the model is empty, and Success otherwise.
  ///
  Status doInteriorSampling
    (ParametersInteriorSampling * parametersInteriorSampling,
     ParametersResults * parametersResults);

  /// Perform the Virial calculation with the provided parameters and
  /// store the results internally.
  ///
  /// Parameters that have not been set may have default values computed.  If
  /// so, these values will be written into the parameters object.
  ///
  /// Returns EmptyModel if the model is empty, and Success otherwise.
  ///
  Status doVirialSampling
    (ParametersVirial * parametersVirial,
     ParametersResults * parametersResults);

  /// Computes final results based on the provided parameters and the results
  /// from the Walk-on-Spheres and Interior Sampling computations.  Different
  /// results are computed depending on which of the computations have been run.
  ///
  /// Parameters that have not been set may have default values computed.  If
  /// so, these values will be written into the parameters object.
  ///
  /// Final results are written into the results object.
  ///
  void getResults(ParametersResults * parametersResults,
		  Results * results) const;

  /// Sets the pointers pointed to by "points" and "charges" to point to vectors
  /// of the hit points and charges from the Walk-on-Spheres algorithm.  If the
  /// Walk-on-Spheres algorithm has not been run with hit point saving 
  /// enabled, these are set to null.
  ///
  void getWalkOnSpheresHitPoints
    (std::vector<Vector3<double> > const * * points,
     std::vector<Vector3<char> > const * * charges) const;

  /// Sets the pointer pointed to by "points" and to point to a vector
  /// of the hit points from the Interior Sampling algorithm.  If the
  /// Interior Sampling algorithm has not been run with hit point saving 
  /// enabled, this is set to null.
  ///
  void getInteriorSamplingHitPoints
    (std::vector<Vector3<double> > const * * points) const;

  double getInitializeTime() const;
  double getPreprocessTime() const;
  double getWalkOnSpheresTime() const;
  double getWalkOnSpheresReductionTime() const;
  double getInteriorSamplingTime() const;
  double getInteriorSamplingReductionTime() const;
  double getVirialTime() const;
  double getVirialReductionTime() const;
  double getTotalTime() const;

 private:
  /// Divides a total number of samples as evenly as possible between a set of
  /// MPI processes of a certain size, and returns how many samples the MPI
  /// process of the given rank is responsible for.
  /// 
  long long computeNumInProcess(int mpiSize, int mpiRank,
				long long totalNumSamples) const;

  /// Sets the default values for parameters that have them if the parameters
  /// have not already been set.
  ///
  void computeDefaultParameters(ParametersWalkOnSpheres * parameters) const;
  void computeDefaultParameters(ParametersInteriorSampling * parameters) const;
  void computeDefaultParameters(ParametersVirial * parameters) const;
  void computeDefaultParameters(ParametersResults * parameters) const;

  /// Allocate a random number generator for each thread, ensuring that each has
  /// a unique stream ID across MPI processes.
  ///
  void setupRNGs(int numThreads,
		 int seed,
		 std::vector<RandomNumberGenerator> * threadRNGs) const;

  /// Estimates the total number of samples that will be required to acheive the
  /// requested error, given the number of samples taken so far and the current
  /// value and error.
  ///
  long long estimateTotalNum(double requestedError,
			     long long numSoFar,
			     Uncertain<double> const & currentValue);
  
  /// Perform Walk-on-Spheres walks until the stopping condition is achieved
  /// (either number of walks or error) and perform a parallel reduction on
  /// the results.
  ///
 void getWalkOnSpheresResults
    (long long numWalksInProcess,			
     ParametersWalkOnSpheres const & parametersWalkOnSpheres,
     ParametersResults const & parametersResults,
     BoundingSphere const & boundingSphere,
     Model const & model,
     std::vector<RandomNumberGenerator> * threadRNGs,
     ResultsZeno * * resultsZeno);

  /// Launches a set of Walk-on-Spheres walks in each of a set of parallel 
  /// threads.
  ///
  void doWalkOnSpheres(ParametersWalkOnSpheres const & parameters,
		       long long numWalksInProcess,
		       BoundingSphere const & boundingSphere, 
		       Model const & nearestSurfacePointFinder,
		       std::vector<RandomNumberGenerator> * threadRNGs,
		       ResultsZeno * resultsZeno);

  /// Launches a given number of Walk-on-Spheres walks and records the results.
  /// Runs in a single thread.
  ///
  static
    void doWalkOnSpheresThread(ParametersWalkOnSpheres const * parameters,
			       BoundingSphere const * boundingSphere, 
			       Model const * nearestSurfacePointFinder,
			       int threadNum,
			       long long numWalks,
			       Timer const * totalTimer,
			       RandomNumberGenerator * randomNumberGenerator,
			       ResultsZeno * resultsZeno);

  /// Perform Interior samples until the stopping condition is achieved
  /// (either number of walks or error) and perform a parallel reduction on
  /// the results.
  ///
  void getInteriorResults
    (long long numSamplesInProcess,			
     ParametersInteriorSampling const & parametersInteriorSampling,
     ParametersResults const & parametersResults,
     BoundingSphere const & boundingSphere,
     Model const & model,
     std::vector<RandomNumberGenerator> * threadRNGs,
     ResultsInterior * * resultsInterior);

  /// Launches a set of Interior samples in each of a set of parallel 
  /// threads.
  ///
  void doInteriorSampling(ParametersInteriorSampling const & parameters,
			  long long numSamplesInProcess,
			  BoundingSphere const & boundingSphere, 
			  Model const & insideOutsideTester,
			  std::vector<RandomNumberGenerator> * threadRNGs,
			  ResultsInterior * resultsInterior);

  /// Performs a given number of Interior samples and records the results.
  /// Runs in a single thread.
  ///
  static 
    void doInteriorSamplingThread(ParametersInteriorSampling const * parameters,
				  BoundingSphere const * boundingSphere, 
				  Model const * insideOutsideTester,
				  int threadNum,
				  long long numSamples,
				  Timer const * totalTimer,
				  RandomNumberGenerator * randomNumberGenerator,
				  ResultsInterior * resultsInterior);

  /// Perform the given number of virial steps
  /// (either number of walks or error) and perform a parallel reduction on
  /// the results.
  ///
  void getVirialResults
    (long long numStepsInProcess,
     ParametersVirial const & parametersVirial,
     ParametersResults const & parametersResults,
     BoundingSphere const & boundingSphere,
     Model const & model,
     std::vector<RandomNumberGenerator> * threadRNGs,
     ResultsVirial * * resultsVirial);

  void doVirialSampling(ParametersVirial const & parameters,
                        long long numStepsInProcess,
                        BoundingSphere const & boundingSphere,
                        Model const & model,
                        std::vector<RandomNumberGenerator> * threadRNGs,
                        ResultsVirial * resultsVirial,
                        double refDiameter);

  static
    void doVirialSamplingThread(ParametersVirial const * parameters,
			        BoundingSphere const & boundingSphere, 
			        Model const & model,
			        int threadNum,
			        long long stepsInThread,
			        Timer const * totalTimer,
			        RandomNumberGenerator * randomNumberGenerator,
			        ResultsVirial * resultsVirial,
			        double refDiameter);

  int mpiSize;
  int mpiRank;

  Model model;

  BoundingSphere modelBoundingSphere;

  ResultsZeno * resultsZeno;
  ResultsInterior * resultsInterior;
  ResultsVirial * resultsVirial;

  Timer initializeTimer;
  Timer preprocessTimer;
  Timer walkOnSpheresTimer;
  Timer walkOnSpheresReductionTimer;
  Timer interiorSamplingTimer;
  Timer interiorSamplingReductionTimer;
  Timer virialTimer;
  Timer virialReductionTimer;
  Timer totalTimer;
};

}

#endif
