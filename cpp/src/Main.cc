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
// Created: Mon Feb 24 11:25:18 2014 EDT
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <thread>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "Parameters.h"
#include "Units.h"
#include "Uncertain.h"

#include "BoundingSphereGenerator.h"

#include "Parser.h"

#include "ResultsZeno.h"
#include "ResultsInterior.h"
#include "ResultsCompiler.h"

#include "Geometry/Sphere.h"
#include "Geometry/Vector3.h"
#include "Geometry/MixedModel.h"

#include "Walker/WalkerExterior.h"
#include "Walker/SamplerInterior.h"

#include "Timer.h"

// ================================================================

using Model = MixedModel<double>;

using BoundingSphere = Sphere<double>;

#if defined USE_RAND_RNG
#include "RandomNumber/Rand.h"
using RandomNumberGenerator = Rand;
#elif defined USE_SPRNG_RNG
#include "RandomNumber/SPRNG.h"
using RandomNumberGenerator = SPRNG;
#else
#error "No RandomNumber class selected"
#endif

#include "SpherePoint/RandomSpherePointMarsaglia.h"
#include "SpherePoint/RandomSpherePointPolar.h"
using RandomSpherePointGenerator = 
  RandomSpherePointMarsaglia<double, 
			     RandomNumberGenerator>;

#include "SpherePoint/BiasedSpherePointRejection.h"
#include "SpherePoint/BiasedSpherePointDirect.h"
using BiasedSpherePointGenerator = 
  BiasedSpherePointDirect<double, 
			  RandomNumberGenerator, 
			  RandomSpherePointGenerator>;

#include "SpherePoint/RandomBallPointRejection.h"
using RandomBallPointGenerator =
  RandomBallPointRejection<double,
			   RandomNumberGenerator>;

// ================================================================

int
getInput(int argc, char **argv,
	 Parameters * parameters,
	 long long * numWalksInProcess,
	 long long * numSamplesInProcess,
	 Model * model,
	 double * initializeTime,
	 double * readTime,
	 double * broadcastTime,
	 std::ofstream * csvOutputFile);

int
preprocess(Parameters const & parameters,
	   BoundingSphere * boundingSphere,
	   Model * model,
	   double * preprocessTime,
	   std::ofstream * csvOutputFile);

void
setupRNGs(Parameters const & parameters,
	  std::vector<RandomNumberGenerator> * threadRNGs);

int
getWalkOnSpheresResults(long long numWalksInProcess,			
			Parameters const & parameters,
			BoundingSphere const & boundingSphere,
			Model const & model,
			Timer const & totalTimer,
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * * resultsZeno,
			double * walkTime,
			double * reduceTime,
			std::ofstream * csvOutputFile);

int
getInteriorResults(long long numSamplesInProcess,			
		   Parameters const & parameters,
		   BoundingSphere const & boundingSphere,
		   Model const & model,
		   Timer const & totalTimer,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * * resultsInterior,
		   double * sampleTime,
		   double * volumeReduceTime,
		   std::ofstream * csvOutputFile);

long long 
computeNumInProcess(int mpiSize, int mpiRank,
		    long long totalNumSamples);

void
parseBodFile(Parameters * parameters,
	     Model * model);

void
getBodData(Parameters * parameters,
	   Model * model,
	   double * readTime,
	   double * broadcastTime);

void
computeDefaultParameters(Parameters * parameters,
			 BoundingSphere * boundingSphere);

long long 
estimateTotalNum(double requestedError,
		 long long numSoFar,
		 Uncertain<double> const & currentValue);

void
doWalkOnSpheres(Parameters const & parameters,
		long long numWalksInProcess,
		BoundingSphere const & boundingSphere, 
		Model const & nearestSurfacePointFinder,
		Timer const & totalTimer,
		std::vector<RandomNumberGenerator> * threadRNGs,
		ResultsZeno * resultsZeno,
		double * walkTime);

void
doWalkOnSpheresThread(Parameters const * parameters,
		      BoundingSphere const * boundingSphere, 
		      Model const * nearestSurfacePointFinder,
		      int threadNum,
		      long long numWalks,
		      Timer const * totalTimer,
		      RandomNumberGenerator * randomNumberGenerator,
		      ResultsZeno * resultsZeno);

void
doInteriorSampling(Parameters const & parameters,
		   long long numSamplesInProcess,
		   BoundingSphere const & boundingSphere, 
		   Model const & insideOutsideTester,
		   Timer const & totalTimer,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * resultsInterior,
		   double * sampleTime);

void
doInteriorSamplingThread(Parameters const * parameters,
			 BoundingSphere const * boundingSphere, 
			 Model const * insideOutsideTester,
			 int threadNum,
			 long long numSamples,
			 Timer const * totalTimer,
			 RandomNumberGenerator * randomNumberGenerators,
			 ResultsInterior * resultsInterior);

void
printOutput(BoundingSphere const & boundingSphere,
	    ResultsInterior const * resultsInterior,
	    ResultsZeno const * resultsZeno,
	    Parameters const & parameters, 
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double sampleTime,
	    double volumeReduceTime,
	    std::ofstream * csvOutputFile);

template <typename T>
void
printExactScalar(std::string const & prettyName,
	         std::string const & csvName,
	         std::string const & units,
	         T property,
	         std::ofstream * csvOutputFile);

void
savePointFiles(ResultsInterior & resultsInterior,
	       ResultsZeno & resultsZeno,
	       Parameters const & parameters);

void
printTime(std::string const & label);

void
printRAM(std::string const & label,
	 std::ofstream * csvOutputFile);

void
writePoints(std::string const & fileName, 
	    std::vector<Vector3<double> > const * points,
	    std::vector<Vector3<char> > const * charges);

// ================================================================

int main(int argc, char **argv) {

  Timer totalTimer;

  totalTimer.start();

  Timer initializeTimer;

  initializeTimer.start();

  int mpiSize = 1, mpiRank = 0;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

  if (mpiRank == 0) {
    printTime("Start time: ");
  }

  Parameters parameters;

  parameters.setMpiSize(mpiSize);
  parameters.setMpiRank(mpiRank);

  long long numWalksInProcess   = 0;
  long long numSamplesInProcess = 0;

  double initializeTime = 0;
  double readTime       = 0;
  double broadcastTime  = 0;

  Model model;

  std::ofstream csvOutputFile;

  initializeTimer.stop();

  int getInputSuccess = getInput(argc, argv,
				 &parameters,
				 &numWalksInProcess,
				 &numSamplesInProcess,
				 &model,
				 &initializeTime,
				 &readTime,
				 &broadcastTime,
				 &csvOutputFile);

  if (getInputSuccess != 0) {
    return getInputSuccess;
  }

  initializeTimer.start();

  BoundingSphere boundingSphere;

  double preprocessTime = 0;

  initializeTimer.stop();

  int preprocessSuccess = preprocess(parameters,
				     &boundingSphere,
				     &model,
				     &preprocessTime,
				     &csvOutputFile);

  if (preprocessSuccess != 0) {
    return preprocessSuccess;
  }

  initializeTimer.start();

  computeDefaultParameters(&parameters, &boundingSphere);

  std::vector<RandomNumberGenerator> threadRNGs;

  setupRNGs(parameters,
	    &threadRNGs);

  ResultsZeno * resultsZeno = NULL;

  double walkTime   = 0;
  double reduceTime = 0;

  initializeTimer.stop();

  int getWalkOnSpheresResultsSuccess = 
    getWalkOnSpheresResults(numWalksInProcess,			
			    parameters,
			    boundingSphere,
			    model,
			    totalTimer,
			    &threadRNGs,
			    &resultsZeno,
			    &walkTime,
			    &reduceTime,
			    &csvOutputFile);

  if (getWalkOnSpheresResultsSuccess != 0) {
    return getWalkOnSpheresResultsSuccess;
  }

  ResultsInterior * resultsInterior = NULL;

  double sampleTime       = 0;
  double volumeReduceTime = 0;

  int getInteriorResultsSuccess = 
    getInteriorResults(numSamplesInProcess,			
		       parameters,
		       boundingSphere,
		       model,
		       totalTimer,
		       &threadRNGs,
		       &resultsInterior,
		       &sampleTime,
		       &volumeReduceTime,
		       &csvOutputFile);

  if (getInteriorResultsSuccess != 0) {
    return getInteriorResultsSuccess;
  }

  initializeTime += initializeTimer.getTime();

  if (parameters.getMaxRunTimeWasSet() &&
      (totalTimer.getTime() > parameters.getMaxRunTime())) {

    std::cerr << std::endl
              << "*** Warning ***" << std::endl
	      << "Max run time was reached.  Not all requested computations "
	      << "may have been performed." << std::endl;
  }

  printOutput(boundingSphere,
	      resultsInterior,
	      resultsZeno,
	      parameters,
	      initializeTime,
	      readTime,
	      broadcastTime,
	      preprocessTime,
	      walkTime,
	      reduceTime,
	      sampleTime,
	      volumeReduceTime,
	      &csvOutputFile);

  savePointFiles(*resultsInterior,
		 *resultsZeno,
		 parameters);

  delete resultsZeno;
  delete resultsInterior;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  totalTimer.stop();

  if (mpiRank == 0) {
    printExactScalar("Total Time", "total_time", "s",
		     totalTimer.getTime(), &csvOutputFile);
    
    std::cout << std::endl;

    printTime("End time: ");
  }

  csvOutputFile.close();

  return 0;
}

/// Get all input to the program, including both parameters and data, from the
/// command line and input files.
///
int
getInput(int argc, char **argv,
	 Parameters * parameters,
	 long long * numWalksInProcess,
	 long long * numSamplesInProcess,
	 Model * model,
	 double * initializeTime,
	 double * readTime,
	 double * broadcastTime,
	 std::ofstream * csvOutputFile) {

  Timer initializeTimer;
  initializeTimer.start();

  parameters->parseCommandLine(argc, argv);

  *numWalksInProcess = 
    computeNumInProcess(parameters->getMpiSize(), 
			parameters->getMpiRank(), 
			parameters->getTotalNumWalks());

  *numSamplesInProcess = 
    computeNumInProcess(parameters->getMpiSize(), 
			parameters->getMpiRank(), 
			parameters->getTotalNumSamples());

  if (parameters->getCsvOutputFileNameWasSet() &&
      parameters->getMpiRank() == 0) {
    
    csvOutputFile->open(parameters->getCsvOutputFileName(),
		        std::ios::out | std::ios::trunc);

    if (!csvOutputFile->is_open()) {
      std::cerr << "*** Warning ***" << std::endl
		<< "Could not open CSV output file "
		<< parameters->getCsvOutputFileName() << std::endl
		<< std::endl;
    }
  }

  if (parameters->getPrintBenchmarks() && 
      parameters->getMpiRank() == 0) {

    printRAM("initialization",
	     csvOutputFile);
  }

  initializeTimer.stop();

  getBodData(parameters,
	     model,
	     readTime,
	     broadcastTime);

  initializeTimer.start();

  if (parameters->getPrintBenchmarks() && 
      parameters->getMpiRank() == 0) {

    printRAM("loading input data",
	     csvOutputFile);
  }

  initializeTimer.stop();
  *initializeTime += initializeTimer.getTime();

  return 0;
}

/// Build the data structure used for the Walk-on-Spheres and Interior Sampling
/// algorithms from geometric data.
///
int
preprocess(Parameters const & parameters,
	   BoundingSphere * boundingSphere,
	   Model * model,
	   double * preprocessTime,
	   std::ofstream * csvOutputFile) {

  Timer preprocessTimer;
  preprocessTimer.start();

  model->preprocess();

  if (model->isEmpty()) {
    std::cerr << "Error: No geometry loaded" << std::endl;

    return 1;
  }
  
  *boundingSphere = BoundingSphereGenerator<double>::generate(*model);

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("building spatial data structure",
	     csvOutputFile);
  }

  preprocessTimer.stop();
  *preprocessTime = preprocessTimer.getTime();

  return 0;
}

/// Allocate a random number generator for each thread, ensuring that each has
/// a unique stream ID across MPI processes.
///
void
setupRNGs(Parameters const & parameters,
	  std::vector<RandomNumberGenerator> * threadRNGs) {

  int numStreams = 
    parameters.getNumThreads() * 
    parameters.getMpiSize();

  threadRNGs->reserve(parameters.getNumThreads());

  for (int threadNum = 0; threadNum < parameters.getNumThreads(); threadNum++) {
    int streamNum = 
      parameters.getMpiRank() * 
      parameters.getNumThreads() + 
      threadNum;

    threadRNGs->emplace_back(streamNum, numStreams, parameters.getSeed());
  }
}

/// Perform Walk-on-Spheres walks until the stopping condition is achieved
/// (either number of walks or error) and perform a parallel reduction on
/// the results.
///
int
getWalkOnSpheresResults(long long numWalksInProcess,			
			Parameters const & parameters,
			BoundingSphere const & boundingSphere,
			Model const & model,
			Timer const & totalTimer,
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * * resultsZeno,
			double * walkTime,
			double * reduceTime,
			std::ofstream * csvOutputFile) {

  Timer reduceTimer;

  bool saveHitPoints = !parameters.getSurfacePointsFileName().empty();

  if (parameters.getTotalNumWalksWasSet() &&
      !parameters.getMaxErrorCapacitanceWasSet() &&
      !parameters.getMaxErrorPolarizabilityWasSet()) {

    *resultsZeno = new ResultsZeno(boundingSphere,
				   parameters.getNumThreads(),
				   saveHitPoints);

    doWalkOnSpheres(parameters,
		    numWalksInProcess,
		    boundingSphere, 
		    model,
		    totalTimer,
		    threadRNGs,
		    *resultsZeno,
		    walkTime);

    reduceTimer.start();
    (*resultsZeno)->reduce();
    reduceTimer.stop();
  }
  else if (parameters.getMaxErrorCapacitanceWasSet() ||
	   parameters.getMaxErrorPolarizabilityWasSet()) {

    *resultsZeno = new ResultsZeno(boundingSphere,
				   parameters.getNumThreads(),
				   saveHitPoints);

    ResultsCompiler resultsCompiler(parameters);

    long long estimatedNumWalksRemaining = parameters.getMinTotalNumWalks();

    while (estimatedNumWalksRemaining > 0) {

      long long estimatedNumWalksRemainingInProcess = 
	computeNumInProcess(parameters.getMpiSize(), 
			    parameters.getMpiRank(), 
			    estimatedNumWalksRemaining);

      doWalkOnSpheres(parameters,
		      estimatedNumWalksRemainingInProcess,
		      boundingSphere, 
		      model,
		      totalTimer,
		      threadRNGs,
		      *resultsZeno,
		      walkTime);

      reduceTimer.start();
      (*resultsZeno)->reduce();
      reduceTimer.stop();

      resultsCompiler.compile(*resultsZeno,
			      NULL,
			      boundingSphere);

      long long estimatedTotalNumWalks = 0;

      long long capacitanceEstimatedTotalNumWalks = 
  	estimateTotalNum(parameters.getMaxErrorCapacitance(),
  			 (*resultsZeno)->getNumWalks(),
  			 resultsCompiler.getCapacitance());

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					capacitanceEstimatedTotalNumWalks);

      long long polarizabilityEstimatedTotalNumWalks = 
  	estimateTotalNum(parameters.getMaxErrorPolarizability(),
  			 (*resultsZeno)->getNumWalks(),
  			 resultsCompiler.getMeanPolarizability());

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					polarizabilityEstimatedTotalNumWalks);

      if (parameters.getTotalNumWalksWasSet()) {
	estimatedTotalNumWalks = std::min(estimatedTotalNumWalks,
					  parameters.getTotalNumWalks());
      }

      estimatedNumWalksRemaining = 
	estimatedTotalNumWalks - (*resultsZeno)->getNumWalks();

      // To avoid repeatedly undershooting, don't take a smaller step than
      // MinTotalNumWalks unless the iterations are about to end
      if (estimatedNumWalksRemaining > 0 &&
	  !(parameters.getTotalNumWalksWasSet() &&
	    estimatedTotalNumWalks == parameters.getTotalNumWalks())) {
	
	estimatedNumWalksRemaining = 
	  std::max(estimatedNumWalksRemaining,
		   parameters.getMinTotalNumWalks());
      }

      if (parameters.getMaxRunTimeWasSet() &&
	  (totalTimer.getTime() > parameters.getMaxRunTime())) {

        break;
      }
    }
  }

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("walk on spheres",
	     csvOutputFile);
  }

  *reduceTime = reduceTimer.getTime();

  return 0;
}

/// Perform Interior samples until the stopping condition is achieved
/// (either number of walks or error) and perform a parallel reduction on
/// the results.
///
int
getInteriorResults(long long numSamplesInProcess,			
		   Parameters const & parameters,
		   BoundingSphere const & boundingSphere,
		   Model const & model,
		   Timer const & totalTimer,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * * resultsInterior,
		   double * sampleTime,
		   double * volumeReduceTime,
		   std::ofstream * csvOutputFile) {

  bool saveInteriorPoints = 
    !parameters.getInteriorPointsFileName().empty();

  Timer volumeReduceTimer;

  if (parameters.getTotalNumSamplesWasSet() &&
      !parameters.getMaxErrorVolumeWasSet()) {

    *resultsInterior = new ResultsInterior(parameters.getNumThreads(),
					   saveInteriorPoints);

    doInteriorSampling(parameters,
		       numSamplesInProcess,
		       boundingSphere, 
		       model,
		       totalTimer,
		       threadRNGs,
		       *resultsInterior,
		       sampleTime);

    volumeReduceTimer.start();
    (*resultsInterior)->reduce();
    volumeReduceTimer.stop();
  }
  else if (parameters.getMaxErrorVolumeWasSet()) {

    *resultsInterior = new ResultsInterior(parameters.getNumThreads(),
					   saveInteriorPoints);

    ResultsCompiler resultsCompiler(parameters);

    long long estimatedNumSamplesRemaining = parameters.getMinTotalNumSamples();

    while (estimatedNumSamplesRemaining > 0) {

      long long estimatedNumSamplesRemainingInProcess = 
	computeNumInProcess(parameters.getMpiSize(), 
			    parameters.getMpiRank(), 
			    estimatedNumSamplesRemaining);

      doInteriorSampling(parameters,
			 estimatedNumSamplesRemainingInProcess,
			 boundingSphere, 
			 model,
			 totalTimer,
			 threadRNGs,
			 *resultsInterior,
			 sampleTime);

      volumeReduceTimer.start();
      (*resultsInterior)->reduce();
      volumeReduceTimer.stop();

      resultsCompiler.compile(NULL,
			      *resultsInterior,
			      boundingSphere);

      long long estimatedTotalNumSamples = 
	estimateTotalNum(parameters.getMaxErrorVolume(),
			 (*resultsInterior)->getNumSamples(),
			 resultsCompiler.getVolume());

      if (parameters.getTotalNumSamplesWasSet()) {
	estimatedTotalNumSamples = std::min(estimatedTotalNumSamples,
					    parameters.getTotalNumSamples());
      }

      estimatedNumSamplesRemaining = 
	estimatedTotalNumSamples - (*resultsInterior)->getNumSamples();

      // To avoid repeatedly undershooting, don't take a smaller step than
      // MinTotalNumSamples unless the iterations are about to end 
      if (estimatedNumSamplesRemaining > 0 &&
	  !(parameters.getTotalNumSamplesWasSet() &&
	    estimatedTotalNumSamples == parameters.getTotalNumSamples())) {
	
	estimatedNumSamplesRemaining = 
	  std::max(estimatedNumSamplesRemaining,
		   parameters.getMinTotalNumSamples());
      }

      if (parameters.getMaxRunTimeWasSet() &&
	  (totalTimer.getTime() > parameters.getMaxRunTime())) {

        break;
      }
    }
  }

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("interior samples",
	     csvOutputFile);
  }

  *volumeReduceTime = volumeReduceTimer.getTime();

  return 0;
}

/// Divides a total number of samples as evenly as possible between a set of
/// MPI processes of a certain size, and returns how many samples the MPI
/// process of the given rank is responsible for.
/// 
long long 
computeNumInProcess(int mpiSize, int mpiRank,
		    long long totalNumSamples) {

  long long numSamplesInProcess = totalNumSamples / mpiSize;

  if (mpiRank < totalNumSamples % mpiSize) {
    numSamplesInProcess ++;
  }

  return numSamplesInProcess;
}

/// Parses the bod file given as input and extracts sphere data, voxel data,
/// and parameters.
///
void
parseBodFile(Parameters * parameters,
	     Model * model) {

  std::string fileName = parameters->getInputFileName();

  std::ifstream inputFile;

  inputFile.open(fileName, std::ifstream::in);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening input file " << fileName << std::endl;
    exit(1);
  }

  Parser parser(inputFile, parameters, model);

  parser.parse();

  inputFile.close();
}

/// Gets the data from the bod file given as input either by parsing the file
/// or by an MPI brodcast, depending on the MPI rank of the process.
///
void
getBodData(Parameters * parameters,
	   Model * model,
	   double * readTime,
	   double * broadcastTime) {

  Timer readTimer;
  Timer broadcastTimer;

  if (parameters->getMpiRank() == 0) {
    readTimer.start();
    parseBodFile(parameters, model);
    readTimer.stop();

    if (parameters->getMpiSize() > 1) {
      broadcastTimer.start();
      model->mpiSend();
      parameters->mpiSend();
      broadcastTimer.stop();
    }
  }
  else {
    broadcastTimer.start();
    model->mpiReceive();
    parameters->mpiReceive();
    broadcastTimer.stop();
  }

  *readTime      = readTimer.getTime();
  *broadcastTime = broadcastTimer.getTime();
}

/// Sets the default values for parameters that have them if the parameters
/// have not already been set.
///
void
computeDefaultParameters(Parameters * parameters,
			 BoundingSphere * boundingSphere) {

  const double defaultSkinThicknessFactor = 0.000001;

  BoundingSphere originalBoundingSphere(*boundingSphere);

  if (parameters->getLaunchCenterWasSet()) {
    boundingSphere->setCenter(parameters->getLaunchCenter());
  }
  else {
    parameters->setLaunchCenter(boundingSphere->getCenter());
  }

  if (parameters->getLaunchRadiusWasSet()) {
    boundingSphere->setRadius(parameters->getLaunchRadius());
  }
  else {
    parameters->setLaunchRadius(boundingSphere->getRadius());
  }

  if (!parameters->getSkinThicknessWasSet()) {
    parameters->setSkinThickness(boundingSphere->getRadius() * 
				 defaultSkinThicknessFactor);
  }

  if (!parameters->getLengthScaleWasSet()) {
    parameters->setLengthScale(1, Units::Length::L);
  }

  if (!boundingSphere->contains(originalBoundingSphere)) {
    std::cerr << std::endl
	      << "*** Warning ***" << std::endl
	      << "User-specified launch sphere may not contain the entire model"
	      << std::endl;
  }
}

/// Estimates the total number of samples that will be required to acheive the
/// requested error, given the number of samples taken so far and the current
/// value and error.
///
long long 
estimateTotalNum(double requestedError,
		 long long numSoFar,
		 Uncertain<double> const & currentValue) {

  double requestedVariance = 
    pow((requestedError / 100) * currentValue.getMean(), 2);

  long long estimatedTotalNum = 
    ceil(currentValue.getVariance() * numSoFar / requestedVariance);

  return estimatedTotalNum;
}

/// Launches a set of Walk-on-Spheres walks in each of a set of parallel 
/// threads.
///
void
doWalkOnSpheres(Parameters const & parameters,
		long long numWalksInProcess,
		BoundingSphere const & boundingSphere, 
		Model const & nearestSurfacePointFinder,
		Timer const & totalTimer,
		std::vector<RandomNumberGenerator> * threadRNGs,
		ResultsZeno * resultsZeno,
		double * walkTime) {

  Timer walkTimer;
  walkTimer.start();

  const int numThreads = parameters.getNumThreads();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numWalksInThread = numWalksInProcess / numThreads;

    if (threadNum < numWalksInProcess % numThreads) {
      numWalksInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doWalkOnSpheresThread,
		      &parameters,
		      &boundingSphere, 
		      &nearestSurfacePointFinder,
		      threadNum,
		      numWalksInThread,
		      &totalTimer,
		      &(threadRNGs->at(threadNum)),
		      resultsZeno);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  walkTimer.stop();
  *walkTime += walkTimer.getTime();
}

/// Launches a given number of Walk-on-Spheres walks and records the results.
/// Runs in a single thread.
///
void
doWalkOnSpheresThread(Parameters const * parameters,
		      BoundingSphere const * boundingSphere, 
		      Model const * nearestSurfacePointFinder,
		      int threadNum,
		      long long numWalks,
		      Timer const * totalTimer,
		      RandomNumberGenerator * randomNumberGenerator,
		      ResultsZeno * resultsZeno) {

  const double shellThickness = parameters->getSkinThickness();
  
  WalkerExterior<double, 
		 RandomNumberGenerator,
		 Model,
		 RandomSpherePointGenerator,
		 BiasedSpherePointGenerator>
    walker(randomNumberGenerator, 
	   *boundingSphere, 
           *nearestSurfacePointFinder,
	   shellThickness);

  for (long long walkNum = 0; walkNum < numWalks; walkNum++) {

    bool hitObject = false;
    int numSteps   = 0;

    Vector3<double> startPoint;
    Vector3<double> endPoint;

    walker.walk(&hitObject, &numSteps,
		&startPoint, &endPoint);

    if (hitObject) {
      resultsZeno->recordHit(threadNum, 
			     startPoint, endPoint,
			     randomNumberGenerator);
    }
    else {
      resultsZeno->recordMiss(threadNum);
    }

    if (parameters->getMaxRunTimeWasSet() &&
	(totalTimer->getTime() > parameters->getMaxRunTime())) {

      break;
    }
  }
}

/// Launches a set of Interior samples in each of a set of parallel 
/// threads.
///
void
doInteriorSampling(Parameters const & parameters,
		   long long numSamplesInProcess,
		   BoundingSphere const & boundingSphere, 
		   Model const & insideOutsideTester,
		   Timer const & totalTimer,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * resultsInterior,
		   double * sampleTime) {

  Timer sampleTimer;
  sampleTimer.start();

  const int numThreads = parameters.getNumThreads();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numSamplesInThread = numSamplesInProcess / numThreads;

    if (threadNum < numSamplesInProcess % numThreads) {
      numSamplesInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doInteriorSamplingThread,
		      &parameters,
		      &boundingSphere, 
		      &insideOutsideTester,
		      threadNum,
		      numSamplesInThread,
		      &totalTimer,
		      &(threadRNGs->at(threadNum)),
		      resultsInterior);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  sampleTimer.stop();
  *sampleTime += sampleTimer.getTime();
}

/// Performs a given number of Interior samples and records the results.
/// Runs in a single thread.
///
void
doInteriorSamplingThread(Parameters const * parameters,
			 BoundingSphere const * boundingSphere, 
			 Model const * insideOutsideTester,
			 int threadNum,
			 long long numSamples,
			 Timer const * totalTimer,
			 RandomNumberGenerator * randomNumberGenerator,
			 ResultsInterior * resultsInterior) {

  SamplerInterior<double, 
		 RandomNumberGenerator,
		 Model,
		 RandomBallPointGenerator>
    sampler(randomNumberGenerator, 
	    *boundingSphere, 
	    *insideOutsideTester);

  for (long long sampleNum = 0; sampleNum < numSamples; sampleNum++) {

    bool hitObject = false;

    Vector3<double> hitPoint;

    sampler.sample(&hitObject,
		   &hitPoint);

    if (hitObject) {
      resultsInterior->recordHit(threadNum,
				 hitPoint);
    }
    else {
      resultsInterior->recordMiss(threadNum);
    }

    if (parameters->getMaxRunTimeWasSet() &&
	(totalTimer->getTime() > parameters->getMaxRunTime())) {

      break;
    }
  }
}

/// Prints parameters, results, and (optionally) detailed running time
/// benchmarks on MPI process 0.
///
void
printOutput(BoundingSphere const & boundingSphere,
	    ResultsInterior const * resultsInterior,
	    ResultsZeno const * resultsZeno,
	    Parameters const & parameters, 
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double sampleTime,
	    double volumeReduceTime,
	    std::ofstream * csvOutputFile) {
  
  if (parameters.getMpiRank() == 0) {
    std::cout << std::endl
	      << "Parameters" << std::endl
	      << "----------" << std::endl
	      << std::endl;

    parameters.print(csvOutputFile);

    std::cout << std::endl
	      << "Results" << std::endl
	      << "-------" << std::endl
	      << std::endl;

    if (resultsZeno != NULL) {
      printExactScalar("Number of walks performed",
		       "num_walks_performed",
		       "1",
		       (long long)resultsZeno->getNumWalks(),
		       csvOutputFile);

      std::cout << std::endl;
    }

    if (resultsInterior != NULL) {
      printExactScalar("Number of interior samples performed",
		       "num_interior_samples_performed",
		       "1",
		       (long long)resultsInterior->getNumSamples(),
		       csvOutputFile);

      std::cout << std::endl;
    }

    ResultsCompiler resultsCompiler(parameters);

    resultsCompiler.compile(resultsZeno,
			    resultsInterior,
			    boundingSphere);

    resultsCompiler.print(parameters.getPrintCounts(),
			  csvOutputFile);

    if (parameters.getPrintBenchmarks()) {
      printExactScalar("Initialize     ", "initialize_time", "s",
		       initializeTime, csvOutputFile);

      printExactScalar("Read           ", "read_time", "s",
		       readTime, csvOutputFile);

      printExactScalar("Broadcast      ", "broadcast_time", "s",
		       broadcastTime, csvOutputFile);
      
      printExactScalar("Preprocess     ", "preprocess_time", "s",
		       preprocessTime, csvOutputFile);

      printExactScalar("Exterior Walk  ", "exterior_walk_time", "s",
		       walkTime, csvOutputFile);

      printExactScalar("Exterior Reduce", "exterior_reduce_time", "s",
		       reduceTime, csvOutputFile);

      printExactScalar("Volume Sample  ", "volume_sample_time", "s",
		       sampleTime, csvOutputFile);

      printExactScalar("Volume Reduce  ", "volume_reduce_time", "s",
		       volumeReduceTime, csvOutputFile);
      
      std::cout << std::endl;
    }
  }
}

template <typename T>
void
printExactScalar(std::string const & prettyName,
	         std::string const & csvName,
	         std::string const & units,
	         T property,
	         std::ofstream * csvOutputFile) {

  std::cout << std::fixed
	    << prettyName << " (" << units << "): " << property << std::endl;
  
  *csvOutputFile << std::scientific
 		 << csvName << ",units," << units
		 << std::endl
		 << csvName << ",mean," << property
		 << std::endl
		 << csvName << ",std_dev," << "0"
		 << std::endl;
}

/// Writes Walk-on-Spheres and Interior Sampling hit points to disk.
/// 
void
savePointFiles(ResultsInterior & resultsInterior,
	       ResultsZeno & resultsZeno,
	       Parameters const & parameters) {

  if (!parameters.getSurfacePointsFileName().empty()) {
    resultsZeno.gatherHitPoints();

    writePoints(parameters.getSurfacePointsFileName(), 
		resultsZeno.getPoints(), 
		resultsZeno.getCharges());
  }

  if (!parameters.getInteriorPointsFileName().empty()) {
    resultsInterior.gatherHitPoints();

    writePoints(parameters.getInteriorPointsFileName(), 
		resultsInterior.getPoints(), 
		NULL);
  }
}

/// Prints the current date and time prefixed by the given label.
///
void
printTime(std::string const & label) {
  const int bufferSize = 256;

  time_t rawtime;
  struct tm * timeinfo;
  char buffer[bufferSize];
    
  time(&rawtime);
  timeinfo = localtime(&rawtime);
    
  strftime(buffer, bufferSize, "%F %T", timeinfo);
  
  std::cout << label << buffer << std::endl;
}

/// Prints the current RAM used by the process prefixed by the given label.
///
void
printRAM(std::string const & label,
	 std::ofstream * csvOutputFile) {
  
  std::ifstream statusFile("/proc/self/status");

  std::string line;

  size_t pos = std::string::npos;

  while (statusFile.good() &&
	 pos == std::string::npos) {

    std::getline(statusFile, line);

    pos = line.find("VmRSS:");
  }

  statusFile.close();

  std::string ram;
  std::string units;
  
  if (pos != std::string::npos) {
    size_t unitsPos = line.find_last_of(" ") + 1;

    size_t ramPos = line.find_last_of("\t :", unitsPos - 2) + 1;
      
    ram = line.substr(ramPos, unitsPos - ramPos);
    
    units = line.substr(unitsPos);
  }
  else {
    ram = "unknown";
    units = "";
  }

  std::string prettyName = "RAM after " + label;
  
  std::string csvName = label + "_ram";
  std::replace(csvName.begin(), csvName.end(), ' ', '_');
  
  printExactScalar(prettyName, csvName, units,
		   ram, csvOutputFile);
}

/// Writes the given set of points and (if not NULL) charges to disk.
///
void
writePoints(std::string const & fileName, 
	    std::vector<Vector3<double> > const * points,
	    std::vector<Vector3<char> > const * charges) {

  std::ofstream outputFile;

  outputFile.open(fileName, std::ofstream::out);

  if (!outputFile.is_open()) {

    std::cerr << "Error opening output file " << fileName << std::endl;
    exit(1);
  }

  for (unsigned int i = 0; i < points->size(); i++) {

    if (charges != NULL) {
      outputFile << charges->at(i).get(0)
		 << charges->at(i).get(1)
		 << charges->at(i).get(2);
    }

    outputFile << std::setw(16) << points->at(i).get(0)
	       << std::setw(16) << points->at(i).get(1)
	       << std::setw(16) << points->at(i).get(2) << std::endl;
  }

  outputFile.close();
}
